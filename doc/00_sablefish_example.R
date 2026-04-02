## ----load_base_om, eval=TRUE--------------------------------------------------
library(SablefishMSE)
data(sable_om)

## ----om1, eval=TRUE-----------------------------------------------------------
# Normal recruitment
data(om_rand_recruit)
om_rand_recruit$recruitment

## ----om2, eval=TRUE-----------------------------------------------------------
data(om_bh_recruit)
om_bh_recruit$recruitment

## ----om3, eval=TRUE-----------------------------------------------------------
data(om_bhcyclic_recruit)
om_bhcyclic_recruit$recruitment

## ----om_list, eval=FALSE------------------------------------------------------
# om_list <- afscOM::listN(om_rand_recruit, om_bh_recruit, om_bhcyclic_recruit, om_immcrash_recruit)

## ----setup_mp_options, eval=TRUE----------------------------------------------
mp_base <- setup_mp_options()
mp_base

## ----mp1, eval=TRUE-----------------------------------------------------------
data(mp_f40)
mp_f40

## ----mp23, eval=TRUE----------------------------------------------------------
data(mp_f50)
mp_f50

## ----mp45, eval=TRUE----------------------------------------------------------
data(mp_5perc)
data(mp_25cap)

mp_5perc
mp_25cap

## ----mp_list, eval=FALSE------------------------------------------------------
# hcr_list <- afscOM::listN(
#     mp_f40, mp_f50,
#     mp_5perc, mp_25cap
# )

## ----mse_options, eval=TRUE---------------------------------------------------
mse_options_base <- setup_mse_options()
mse_options <- mse_options_base
mse_options$n_spinup_years <- 54
mse_options$recruitment_start_year <- 54
mse_options$n_proj_years <- 25

mse_options

mse_options_list <- afscOM::listN(mse_options)

## ----run_mse, eval=FALSE------------------------------------------------------
# nsims <- 20
# seed_list <- sample(1:(1000*nsims), nsims)  # Draw 20 random seeds
# model_runs <- run_mse_multiple(
#     om_list,
#     hcr_list,
#     seed_list,
#     nyears=100,
#     mse_options_list=mse_options_list,
#     diagnostics = TRUE,
#     save=TRUE
# )

## ----extra_columns, eval=FALSE------------------------------------------------
# # Quick Way to get the names of the OMs and HCRs in the same order as they appear in the lists
# om_names <- unlist(lapply(om_list, \(x) x$name))
# hcr_names <- unlist(lapply(hcr_list, \(x) x$name))
# 
# extra_columns <- expand.grid(
#     om = om_names,
#     hcr = hcr_names
# )

## ----data_processing, eval=FALSE----------------------------------------------
# ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=hcr_names, om_filter=om_names)
# 
# f_data <- get_fishing_mortalities(model_runs, extra_columns, hcr_filter=hcr_names, om_filter=om_names)
# 
# abctac <- get_management_quantities(model_runs, extra_columns, spinup_years=mse_options$n_spinup_years, hcr_filter=hcr_names, om_filter=om_names)
# 
# catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=hcr_names, om_filter=om_names)
# 

## ----plotting, eval=FALSE-----------------------------------------------------
# common_trajectory <- mse_options$n_spinup_years # the number of years of historical data before the recruitment and HCR functions begin applying
# 
# plot_ssb(ssb_data, v1="hcr", v2="om", common_trajectory=common_trajectory, show_est = FALSE)
# 
# plot_fishing_mortalities(f_data, v1="hcr", v2="om", common_trajectory = common_trajectory, show_est=FALSE)
# 
# plot_abc_tac(abctac, v1="hcr", v2="om", common_trajectory=common_trajectory)
# 
# plot_landed_catch(catch_data, v1="hcr", v2="om", common_trajectory = common_trajectory)
# 

## ----performance_metrics, eval=FALSE------------------------------------------
# perf_tradeoffs <- compute_performance_metric_summary(
#     model_runs,
#     extra_columns,
#     sable_om$dem_params,
#     ref_naa,
#     hcr_filter=hcr_names,
#     om_filter=om_names,
#     interval_widths=c(0.50, 0.80),
#     time_horizon = c(55, 130),
#     extra_filter = NULL,
#     relative=NULL,
#     summarise_by=c("om", "hcr"),
#     summary_out = FALSE,
#     metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "prop_years_lowssb")
# )
# 
# perf_data <- performance_metrics$perf_data
# plot_performance_metric_summary(perf_data)
# 

