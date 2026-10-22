# Drone structure-from-motion (SfM) point cloud to tree-level structure.
#
# Steps
#   1. Read a photogrammetric point cloud (LAS/LAZ) exported from Metashape,
#      Pix4D or WebODM.
#   2. Classify ground with three filters (PMF, CSF, MCC) and build a DTM from
#      each. Score each DTM against surveyed ground control points (GCPs) and
#      keep the filter with the lowest RMSE.
#   3. Normalise heights, build a canopy height model (CHM) and smooth it.
#   4. Detect tree tops with a variable-window local maximum filter.
#   5. Segment crowns (Dalponte 2016) and compute per-tree metrics:
#      height, convex-hull crown area, gap fraction and leaf area density.
#   6. Optionally join trees to field plots for model training.
#
# Distilled from the PhD analysis notebooks (six Victorian eucalypt sites).
# Edit the CONFIG block, then run with: Rscript R/uav_sfm_forest_structure.R
#
# Requires: lidR, terra, sf, dplyr. RCSF is required for CSF.

suppressPackageStartupMessages({
  library(lidR); library(terra); library(sf); library(dplyr)
})

# ---- CONFIG ------------------------------------------------------------------
cfg <- list(
  las_file   = Sys.getenv("UAV_LAS",   unset = "data/uav/site_pointcloud.laz"),
  gcp_file   = Sys.getenv("UAV_GCP",   unset = "data/uav/site_gcps.csv"),  # cols: Easting, Northing, Z
  plots_file = Sys.getenv("UAV_PLOTS", unset = ""),                       # optional polygons with field plot IDs
  plot_id    = "Plot_no",
  out_dir    = Sys.getenv("UAV_OUT",   unset = "outputs/uav_sfm"),
  dtm_res    = 0.5,    # m
  chm_res    = 0.5,    # m
  min_tree_h = 6,      # m, discard tree tops below this height
  max_crown  = 10      # m, Dalponte maximum crown diameter
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Read -----------------------------------------------------------------
las <- readLAS(cfg$las_file, select = "xyzc")
if (is.empty(las)) stop("Point cloud is empty: ", cfg$las_file)
las <- filter_duplicates(las)

# ---- 2. Ground filters and DTM selection by GCP RMSE --------------------------
# Algorithms are built lazily so a missing optional package (RCSF) only drops CSF.
ground_algos <- list(
  PMF = function() pmf(ws = seq(3, 12, 3), th = seq(0.5, 3, length.out = 4)),
  CSF = function() csf(sloop_smooth = TRUE, class_threshold = 1, cloth_resolution = 1, time_step = 1),
  MCC = function() mcc(s = 1.5, t = 0.3)
)
# CSF needs RCSF and MCC needs RMCC; drop a filter if its package is missing.
for (nm in c(CSF = "RCSF", MCC = "RMCC") |> names()) {
  pkg <- c(CSF = "RCSF", MCC = "RMCC")[[nm]]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(pkg, " not installed; skipping ", nm, ".")
    ground_algos[[nm]] <- NULL
  }
}
dtm_score <- function(name) {
  g <- classify_ground(las, ground_algos[[name]](), last_returns = FALSE)
  dtm <- rasterize_terrain(g, res = cfg$dtm_res, algorithm = knnidw(k = 10L, p = 2))
  rmse <- NA_real_
  if (file.exists(cfg$gcp_file)) {
    gcp <- read.csv(cfg$gcp_file)
    z_dtm <- terra::extract(dtm, as.matrix(gcp[, c("Easting", "Northing")]))[, 1]
    rmse <- sqrt(mean((z_dtm - gcp$Z)^2, na.rm = TRUE))
  }
  list(name = name, las = g, dtm = dtm, rmse = rmse)
}
dtms <- lapply(names(ground_algos), dtm_score)
dtm_table <- data.frame(filter = sapply(dtms, `[[`, "name"),
                        gcp_rmse_m = sapply(dtms, `[[`, "rmse"))
print(dtm_table)
write.csv(dtm_table, file.path(cfg$out_dir, "dtm_gcp_rmse.csv"), row.names = FALSE)
best <- dtms[[if (all(is.na(dtm_table$gcp_rmse_m))) length(dtms) else which.min(dtm_table$gcp_rmse_m)]]
message("Using ground filter: ", best$name)
terra::writeRaster(best$dtm, file.path(cfg$out_dir, "dtm.tif"), overwrite = TRUE)

# ---- 3. Normalise and CHM -----------------------------------------------------
nlas <- normalize_height(best$las, knnidw(), dtm = best$dtm)
nlas <- filter_poi(nlas, Z >= -0.5 & Z <= 60)
chm <- rasterize_canopy(nlas, res = cfg$chm_res, algorithm = p2r(subcircle = 0.15))
chm_s <- terra::focal(chm, w = matrix(1, 3, 3), fun = median, na.rm = TRUE)
terra::writeRaster(chm_s, file.path(cfg$out_dir, "chm_smoothed.tif"), overwrite = TRUE)

# ---- 4. Tree tops ------------------------------------------------------------
# Window grows with height. Coefficients were tuned for open eucalypt canopy;
# recalibrate for dense plantations.
ws_fun <- function(h) pmax(3, h * 0.10 + 2)
ttops <- locate_trees(chm_s, lmf(ws_fun))
ttops <- ttops[ttops$Z >= cfg$min_tree_h, ]
message(nrow(ttops), " tree tops detected")

# ---- 5. Crown segmentation and tree metrics ----------------------------------
seg <- segment_trees(nlas, dalponte2016(chm_s, ttops, max_cr = cfg$max_crown))
seg <- filter_poi(seg, !is.na(treeID))

gap_fraction_mean <- function(z, dz = 1, z0 = 2) {
  gf <- lidR::gap_fraction_profile(z, dz = dz, z0 = z0)
  if (nrow(gf) == 0) return(NA_real_)
  mean(gf$gf, na.rm = TRUE)
}
lad_cv <- function(z, dz = 1, k = 0.5, z0 = 2) {
  lad <- lidR::LAD(z, dz = dz, k = k, z0 = z0)$lad
  lad <- lad[is.finite(lad)]
  if (length(lad) < 2) return(NA_real_)
  stats::sd(lad) / mean(lad)
}
tree_fun <- ~list(Zmax = max(Z), Zq95 = stats::quantile(Z, 0.95),
                  Zmean = mean(Z), npts = length(Z),
                  gap_frac = gap_fraction_mean(Z), lad_cv = lad_cv(Z))
crowns <- crown_metrics(seg, func = tree_fun, geom = "convex")
crowns$crown_area_m2 <- as.numeric(sf::st_area(crowns))
sf::st_write(crowns, file.path(cfg$out_dir, "tree_crowns.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---- 6. Optional join to field plots -----------------------------------------
if (nzchar(cfg$plots_file) && file.exists(cfg$plots_file)) {
  plots <- sf::st_transform(sf::st_read(cfg$plots_file, quiet = TRUE), sf::st_crs(crowns))
  trees_in_plots <- sf::st_join(sf::st_centroid(crowns), plots[, cfg$plot_id], left = FALSE)
  write.csv(sf::st_drop_geometry(trees_in_plots),
            file.path(cfg$out_dir, "uav_trees_by_plot.csv"), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "sessionInfo.txt"))
message("Outputs written to ", cfg$out_dir)
