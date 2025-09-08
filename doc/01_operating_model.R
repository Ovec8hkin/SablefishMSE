## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----sable_om, eval=FALSE-----------------------------------------------------
# sable_om <- readRDS(file.path(here::here(), "data", "sablefish_om_big.RDS"))
# 
# sable_om$dem_params
# sable_om$model_options$obs_pars
# sable_om$init_naa
# 

## ----sable_om_dem_params, eval=TRUE, include=FALSE----------------------------
# library(dplyr)
# library(ggplot2)
# library(patchwork)

# sable_om <- readRDS(file.path(here::here(), "data", "sablefish_om_big.RDS"))
# plots <- afscOM::plot_demographic_parameters(sable_om$dem_params, show_plots = TRUE)


