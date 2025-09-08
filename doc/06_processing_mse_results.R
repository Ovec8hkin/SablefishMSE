## ----get_saved_model_runs, eval=FALSE-----------------------------------------
# model_runs <- get_saved_model_runs()

## ----bind_mse_outputs, eval=FALSE---------------------------------------------
# bind_mse_outputs(
#     model_runs=model_runs,
#     var=c("naa"),
#     extra_columns=expand.grid(
#         om=c("OM1", "OM2"),
#         hcr=c("MP1", "MP2")
#     )
# )

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

