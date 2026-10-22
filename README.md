# uav-sfm-tree-metrics

R workflow that turns a drone structure-from-motion (SfM) or aerial LiDAR point cloud into individual tree crowns and structural metrics ready for biomass or inventory modelling.

**Skills shown:** photogrammetric and LiDAR point cloud processing, ground filtering with accuracy checks, canopy height models, automated tree detection and crown segmentation, geospatial outputs.

## Workflow (`R/uav_sfm_forest_structure.R`)

1. Read a LAS or LAZ point cloud exported from Agisoft Metashape, Pix4D or WebODM.
2. Classify ground with PMF, CSF and MCC and build a DTM from each. Score each DTM against surveyed ground control points by RMSE and keep the best filter.
3. Normalise heights, build a canopy height model and smooth it with a median filter.
4. Detect tree tops with a variable-window local maximum filter.
5. Segment crowns with the Dalponte (2016) algorithm.
6. Calculate per-tree height metrics, convex-hull crown area, mean gap fraction and the coefficient of variation of leaf area density.
7. Optionally join trees to field plot polygons for model training.

Outputs are a DTM and smoothed CHM (GeoTIFF), tree crowns with metrics (GeoPackage), a DTM accuracy table and session information.

## Usage

```r
Sys.setenv(UAV_LAS = "data/site.laz", UAV_GCP = "data/gcps.csv", UAV_PLOTS = "data/plots.gpkg")
source("R/uav_sfm_forest_structure.R")
```
The GCP file needs columns `Easting`, `Northing` and `Z`. Without GCPs the script still runs and uses MCC (or PMF if RMCC is not installed).

The script was tested end to end on the `Topography.laz` sample shipped with lidR (2,264 crowns detected).

## Background

The workflow comes from my PhD analysis of six eucalypt sites in Victoria (plantation and native forest), where UAV tree height and crown area were used to estimate DBH and aboveground biomass.

## Installation

```r
source("install_packages.R")   # lidR, RCSF, RMCC, terra, sf, dplyr
```

## Data

Drone data are not included. `.gitignore` excludes `data/`, `outputs/`, point clouds and rasters.

## Related repositories

- [tls-stem-reconstruction](../tls-stem-reconstruction): TLS stem segmentation and diameter profiles
- [uav-sfm-tree-metrics](../uav-sfm-tree-metrics): drone point cloud tree detection and crown metrics
- [forest-model-transfer-uncertainty](../forest-model-transfer-uncertainty): site-transfer validation and prediction intervals

## Author

Hammed Akano, PhD candidate, Deakin University. Contact hammedaakano@gmail.com.

## Licence

MIT. See `LICENSE`.
