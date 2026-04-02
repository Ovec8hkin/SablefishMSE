## ----source, include=FALSE----------------------------------------------------
source(file.path(here::here(), "R/setup_mse_options.R"))

## ----mse_options2, eval=TRUE, echo=TRUE---------------------------------------
options <- setup_mse_options()
options

## ----run_mse_parallel, eval=FALSE, echo=TRUE----------------------------------
# run_mse_parallel(
#     nsims = 10,
#     seeds = sample(1:1000, 10),
#     om=om,
#     hcr=mp,
#     mse_options=options,
#     nyears = 1
# )

## ----run_mse_multiple, eval=FALSE---------------------------------------------
# run_mse_multiple(
#     # List of OMs to run
#     om_list = list(om1, om2),
#     # List of MPs to run
#     hcr_list = list(mp1, mp2),
#     # List of MSE options to run
#     mse_options_list = list(options)
#     seed_list = sample(1:1000, 10),
#     nyears=1,
#     diagnostics = FALSE,
#     save = FALSE
# )

## ----example, eval=FALSE------------------------------------------------------
# om_list = afscOM::listN(om1, om2)
# mp_list = afscOM::listN(mp_f40, mp_f50)
# 
# mse_options <- setup_mse_options()
# 
# seed_list <- sample(1:1000, 10)
# 
# model_runs <- run_mse_multiple(
#     om_list = om_list,
#     hcr_list = mp_list,
#     seed_list = sample(1:1000, 10),
#     mse_options_list = list(mse_options),
#     nyears=100,
#     diagnostics = FALSE,
#     save = TRUE
# )
# 

