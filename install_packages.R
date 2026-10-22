pkgs <- c("lidR", "RCSF", "RMCC", "terra", "sf", "dplyr")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new)
