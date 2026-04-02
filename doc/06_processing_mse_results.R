## ----get_saved_model_runs, eval=FALSE-----------------------------------------
# om_filter <- c("OM1", "OM2")
# hcr_filter <- c("MP1", "MP2")
# 
# model_runs <- get_saved_model_runs(om_filter=om_filter, hcr_filter=hcr_filter)

## ----process_ssb_examps, eval=FALSE-------------------------------------------
# om_filter <- c("OM1", "OM2")
# hcr_filter <- c("MP1", "MP2")
# 
# mse_runs <- get_saved_model_runs(om_filter=om_filter, hcr_filter=hcr_filter)
# model_runs <- mse_runs$model_runs
# extra_columns <- mse_runs$extra_columns
# 
# # Calculate SSB and total biomass across all model runs and simulations
# ssb_biomass <- get_ssb_biomass(model_runs=model_runs, extra_columns=extra_columns, om1$dem_params, om_filter=om_filter, hcr_filter=hcr_filter)
# 
# # Alternatively, if `get_saved_model_runs` is not used
# rm(model_runs, extra_columns)
# mse_runs <- NULL
# 
# # Itertaively read in and process each model run file as appropriate
# ssb_biomass <- get_ssb_biomass(model_runs=NULL, extra_columns=NULL, om1$dem_params, om_filter=om_filter, hcr_filter=hcr_filter)
# 

## ----performance_metric_summary, eval=FALSE-----------------------------------
# performance_metric_summary(
#     model_runs=model_runs,
#     extra_columns=extra_columns,
#     dem_params=om$dem_params,
#     ref_naa=ref_naa,
#     hcr_filter=c("MP1", "MP2"),
#     om_filter=c("OM1", "OM2"),
#     interval_widths=c(0.50, 0.95),
#     # List of metrics to calculate
#     metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "prop_years_lowssb")
# )

