# install_R_dependencies.R

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(c("SingleCellExperiment", "zellkonverter"), update = FALSE, ask = FALSE)
install.packages("dplyr", repos = "https://cloud.r-project.org")
