## ----recruitment_obj, eval=FALSE----------------------------------------------
# recruit_func <- function(nyears){
#     return(rep(10, nyears))
# }
# 
# recruitment_obj <- list(
#     # function from which to generate future recruitment events
#     func = recruit_func,
#     # extra parameters to pass to `func`
#     pars = list(
#         nyears = 100
#     )
# )
# 
# om$recruitment <- recruitment_obj

## ----resample_recruitment, eval=FALSE-----------------------------------------
# om$recruitment$func <- resample_recruits
# om$recruitment$pars <- list(
#     # vector of values to resample from
#     hist_recruits = hist_recruits,
#     # number of years to generate future recruitment for
#     nyears = nyears
# )

## ----bevholt_recruitment, eval=FALSE------------------------------------------
# om$recruitment$func <- beverton_holt
# om$recruitment$pars <- list(
#     h = 0.85,   # steepness
#     R0 = 15,    # unfished recruitment
#     S0 = 300,   # unfished spawning biomass
#     sigR = 1.20 # recruitment variability
# )

## ----regime_recruitment1, eval=FALSE------------------------------------------
# om$recruitment$func <- resample_regime_recruits
# om$recruitment$pars <- list(
#     # Vector of recruits to resample from in regime 1
#     regime1_recruits = seq(1, 25, 1)],
#     # Vector of recruits to resample from in regime 2
#     regime2_recruits = seq(60, 100, 10),
#     # Total number of years to generate future recruitment for
#     nyears = nyears,
#     # Lengths of each regime
#     regime_length = c(20, 5),
#     # Start with the first regime
#     starting_regime = 0
# )

## ----regime_recruitment2, eval=FALSE------------------------------------------
# om$recruitment$func <- bevholt_regimes
# om$recruitment$pars <- list(
#     # steepness
#     h = 0.85,
#     # spawning biomass per recruit (to calculate regime specific S0)
#     sbpr = 20,
#     # Regime specfic R0
#     R0 = c(12.5, 50),
#     # Regime specific recruitment variability
#     sigR = c(1.20, 1.20),
#     # Total number of years to generate future recruitment for
#     nyears = nyears,
#     # Lengths of each regime
#     regime_length = c(20, 5),
#     # Start with first regim
#     starting_regime = 0
# )

## ----crash_recruitment, eval=FALSE--------------------------------------------
# om$recruitment$func <- recruits_crash
# om$recruitment$pars <- list(
#     # Start a crash period in simulation years 1
#     crash_start_year = 1,
#     # 20-year crash period
#     crash_length = 20,
#     # Use minimum historical recruitment as average crash recruitment level
#     crash_value = min(hist_recruits),
#     # recruitment vector to resample from outside of the "crash period"
#     hist_recruits = hist_recruits,
#     # Total number of years to generte future recruitment for
#     nyears = nyears
# )

## ----recruit_func2, eval=FALSE------------------------------------------------
# weighted_resample_recruitment <- function(hist_recruits, nyears, weights, seed){
#     set.seed(seed)
#     r <- sample(hist_recruits, nyears, replace=TRUE, prob=weights)
#     return(r)
# }
# 
# 
# recruitment_obj <- list(
#     # function from which to generate future recruitment events
#     func = recruit_func,
#     # extra parameters to pass to `func`
#     pars = list(
#         hist_recruits = seq(1, 100, 5)
#         nyears = 100,
#         weights = c(rep(0.025, 10), rep(0.05, 10))
#     )
# )
# 
# om$recruitment <- recruitment_obj
# 

## ----bevhlt2, eval=FALSE------------------------------------------------------
# beverton_holt_ab <- function(a, b, seed){
#     set.seed(seed)
#     function(ssb, y){
#         bh <- (a*ssb)/(1+b*ssb)
#         return(bh)
#     }
# }
# 
# recruitment_obj <- list(
#     # reference to the beverton_holt_ab function above
#     func = beverton_holt_ab,
#     # parameter values for 'a', and 'b'
#     pars = list(
#         a = 10
#         b = 5
#     )
# )
# 
# om$recruitment <- recruitment_obj
# 

