## ----load_base_om, eval=FALSE-------------------------------------------------
# nyears <- 100
# 
# sable_om <- readRDS("data/sablefish_om_big.RDS")
# 
# # Also load the historical recruitment timeseries
# assessment <- dget("data/sablefish_assessment_2023.rdat")
# hist_recruits <- assessment$natage.female[,1]*2

## ----om1, eval=FALSE----------------------------------------------------------
# om1 <- sable_om
# om1$recruitment$func <- resample_recruits
# om1$recruitment$pars <- list(
#     hist_recruits = hist_recruits,
#     nyears = nyears
# )

## ----om2, eval=FALSE----------------------------------------------------------
# om2 <- sable_om
# om2$recruitment$func <- resample_regime_recruits
# om2$recruitment$pars <- list(
#     regime1_recruits = hist_recruits[hist_recruits <= 35],
#     regime2_recruits = hist_recruits[hist_recruits > 35],
#     nyears = 10*nyears,
#     regime_length = c(20, 5),
#     starting_regime = 0
# )

## ----om4, eval=FALSE----------------------------------------------------------
# om4 <- sable_om
# om4$recruitment$func <- regime_recruits
# om4$recruitment$pars <- list(
#     mus = c(mean(hist_recruits[hist_recruits <= 35])),
#     cvs = c(sd(hist_recruits[hist_recruits <= 35])/mean(hist_recruits[hist_recruits <= 35])),
#     nyears = 10*nyears,
#     regime_length = c(nyears*10),
#     starting_regime = 0
# )

## ----om5, eval=FALSE----------------------------------------------------------
# om6 <- sable_om
# om6$recruitment$func <- recruits_crash
# om6$recruitment$pars <- list(
#     crash_start_year = 1,
#     crash_length = 10,
#     crash_value = 0,
#     hist_recruits = hist_recruits,
#     nyears = 10*nyears
# )

## ----om_list, eval=FALSE------------------------------------------------------
# om_list <- listN(om1, om2, om3, om4, om6)

## ----setup_mp_options, eval=FALSE---------------------------------------------
# mp_base <- setup_mp_options()

## ----atainment_function, eval=FALSE-------------------------------------------
# mp_base$management$tac_land_reduction <- list(
#     func = stairstep_attainment,
#     pars = list(
#         breakpoints = c(20, 30),
#         levels = c(0.874, 0.786, 0.647),
#         phase_ins = 0
#     )
# )

## ----mp1, eval=FALSE----------------------------------------------------------
# mp1 <- mp_base
# mp1$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = NA
#     ),
#     units = "F"
# )

## ----mp234, eval=FALSE--------------------------------------------------------
# mp2 <- mp_base
# mp2$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = 0.10,
#         harvest_cap = NA
#     ),
#     units = "F"
# )
# 
# mp3 <- mp_base
# mp3$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = 0.15,
#         harvest_cap = NA
#     ),
#     units = "F"
# )
# 
# mp4 <- mp_base
# mp4$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = 0.25,
#         harvest_cap = NA
#     ),
#     units = "F"
# )

## ----mp5, eval=FALSE----------------------------------------------------------
# mp5 <- mp_base
# mp5$ref_points$spr_target <- 0.50
# mp5$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = NA
#     ),
#     units = "F"
# )

## ----mp67, eval=FALSE---------------------------------------------------------
# mp6 <- mp_base
# mp6$hcr <- list(
#     func = tier3,
#     extra_pars = list(cutoff_age = 7),
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = NA
#     ),
#     units = "F"
# )
# 
# mp7 <- mp_base
# mp7$hcr <- list(
#     func = tier3,
#     extra_pars = list(cutoff_age = 14),
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = NA
#     ),
#     units = "F"
# )

## ----mp8, eval=FALSE----------------------------------------------------------
# mp8 <- mp_base
# mp8$ref_points$spr_target <- 0.45
# mp8$hcr <- list(
#     func = pfmc4010,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = NA
#     ),
#     units = "TAC"
# )

## ----mp91011, eval=FALSE------------------------------------------------------
# mp9 <- mp_base
# mp9$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = 15
#     ),
#     units = "F"
# )
# 
# mp10 <- mp_base
# mp10$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = 20
#     ),
#     units = "F"
# )
# 
# mp11 <- mp_base
# mp11$hcr <- list(
#     func = tier3,
#     extra_pars = NA,
#     extra_options = list(
#         max_stability = NA,
#         harvest_cap = 25
#     ),
#     units = "F"
# )

## ----mp_list, eval=FALSE------------------------------------------------------
# mp_list <- listN(mp1, mp2, mp, mp4, mp5, mp6, mp7, mp8, mp9, mp10, mp11)

## ----mse_options, eval=FALSE--------------------------------------------------
# mse_options <- setup_mse_options()

## ----run_mse, eval=FALSE------------------------------------------------------
# nsims <- 9
# seed_list <- sample(1:(1000*nsims), nsims)
# model_runs <- run_mse_multiple(
#     om_list,
#     hcr_list,
#     seed_list,
#     nyears=nyears,
#     mse_options_list=mse_options_list
# )

## ----extra_columns, eval=FALSE------------------------------------------------
# om_names <- c("Random Recruitment", "Cyclic Recruitment", "Beverton-Holt Recruitment", "Low Recruitment", "Recruitment Crash")
# hcr_names <- c("F40", "F40 +/-10%", "F40 +/-15%", "F40 +/-25%", "F50", "F40 >7yo", "F40 >14yo", "PFMC 40-10", "F40 15k Cap", "F40 20k Cap", "F40 25k Cap")
# 
# extra_columns <- expand.grid(
#     om = om_names,
#     hcr = hcr_names
# )

## ----data_processing, eval=FALSE----------------------------------------------
# ssb_data <- get_ssb_biomass(model_runs, extra_columns2, sable_om$dem_params)
# f_data <- get_fishing_mortalities(model_runs, extra_columns2)
# r_data <- get_recruits(model_runs, extra_columns2)
# abctac <- get_management_quantities(model_runs, extra_columns2)
# catch_data <- get_landed_catch(model_runs, extra_columns2)

## ----plotting, eval=FALSE-----------------------------------------------------
# plot_ssb(ssb_data, v1="hcr", v2="om")
# plot_fishing_mortalities(f_data, v1="hcr", v2="om")
# plot_recruitment(r_data, v1="hcr", v2="om")
# plot_abc_tac(abctac, v1="hcr", v2="om")
# plot_landed_catch(catch_data, v1="hcr", v2="om")

## ----performance_metrics, eval=FALSE------------------------------------------
# performance_metrics <- performance_metric_summary(
#     model_runs,
#     extra_columns2,
#     sable_om$dem_params,
# )
# 
# perf_data <- performance_metrics$perf_data
# plot_performance_metric_summary(perf_data)
# 

