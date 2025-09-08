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
#     nyears=1
# )

