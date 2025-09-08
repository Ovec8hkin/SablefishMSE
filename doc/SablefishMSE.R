## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----setup--------------------------------------------------------------------
# remotes::install_github("ovec8hkin/SablefishMSE")
# library(SablefishMSE)

## ----load, eval=TRUE, echo=FALSE, results="hide", message=FALSE---------------
library(tidyverse)

sablefishMSE_dir <- here::here()
lapply(list.files(file.path(sablefishMSE_dir, "R"), full.names = TRUE), source)

## ----setup-2------------------------------------------------------------------
# remotes::install_github("BenWilliams-NOAA/afscOM")
# remotes::install_github("Craig44/SpatialSablefishAssessment")

## ----mse_func, echo=TRUE------------------------------------------------------
# run_mse_multiple(
#     om,
#     mp,
#     options,
#     seed=seed
# )

## ----echo=TRUE----------------------------------------------------------------
# sable_om <- readRDS("data/sablefish_om_big.RDS") # Read this saved OM from a file
# sable_om$dem_params     # List of demographic parameter arrays
# sable_om$model_options  # List of OM options (e.g. observation process parameters)
# sable_om$init_naa       # Array of numbers-at-age in years 1 of the OM

## ----recruitment, echo=TRUE---------------------------------------------------
# # Define recruitment to occur via historical resampling
# recruitment_obj <- list(
#   # R function for generating future recruitments
#   func=resample_recruits,
#   # List of parameters values to pass to `func`
#   pars <- list(
#       hist_recruits = hist_recruits,
#       nyears = 10*nyears
#   )
# )
# 
# sable_om$recruitment <- recruitment_obj
# 

## ----setup_mp_options, eval=TRUE----------------------------------------------
mp <- setup_mp_options()
mp

## ----tier3, eval=TRUE, include=FALSE------------------------------------------
tier3 <- function(ref_pts, naa, dem_params, avgrec, cutoff_age=1){
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- cutoff_age-1
    ssb <- apply(naa[,a:nages,1,]*dem_params$waa[,a:nages,1,,drop=FALSE]*dem_params$mat[,a:nages,1,,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$Bref, ref_pts$Fref)
    )
}

## ----hcr, eval=TRUE-----------------------------------------------------------
hcr_obj <- list(
    # R function defining HCR 
    func = tier3,
    # List of extra parameters required by `func`
    extra_pars = NA,
    # Specify stability constraints and/or harvest caps
    extra_options = list(
        max_stability = NA,
        harvest_cap = NA
    ),
    # Output units of `func`
    units = "F"
)

hcr_obj

mp$hcr <- hcr_obj

## ----mse_options_list, eval=TRUE----------------------------------------------
mse_options <- setup_mse_options()
mse_options

## ----run_mse------------------------------------------------------------------
# mse <- run_mse(
#   om = sable_om,
#   mp = mp,
#   options = mse_options,
#   seed=seed
# )

## ----bind_mse_outputs, echo=TRUE----------------------------------------------
# model_runs <- list(
#     mse
# )
# extra_columns <- expand.grid(
#     om = c("om1")
#     hcr = c("mp1")
# )
# 
# naa_out <- bind_mse_outputs(model_runs, var=c("naa"), extra_columns)

## ----plot_ssb, eval=FALSE-----------------------------------------------------
# # Plot spawning biomass from OM and EM
# ssb_data <- get_ssb_biomass(
#   model_runs=model_runs,
#   extra_columns=extra_columns,
#   dem_params=sable_om$dem_params,
# )
# 
# plot_ssb(
#   ssb_data,         # SSB Data to plot
#   v1="hcr",         # Variable to associate with color
#   v2="om",          # Variable to facet by
#   v3=NA,
#   show_est = FALSE  # Dont plot SSB estimates
# )
# 

