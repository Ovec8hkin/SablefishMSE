## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(out.height = "\\textheight",  out.width = "\\textwidth")

## ----simple_hcr---------------------------------------------------------------
npfmc_tier3_F <- function(ssb, B_ref, F_ref){
    x <- ssb/B_ref
    f_min <- 0
    f_max = F_ref
    lrp = 0.05
    urp = 1 

    if(x >= urp)  F <- f_max # if stock status >= 1
    if(x > lrp && x < urp) F <- f_max * ((x-lrp)/(urp-lrp)) # if stock status > alpha & stock status < 1
    if(x <= lrp) F <- f_min
    
    return(F)
}

tier3 <- function(ref_pts, naa, dem_params){
    # Calculate ssb as naa*waa*mat (female only)
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    return(
        npfmc_tier3_F(ssb, ref_pts$B40, ref_pts$F40)
    )
}

## ----visualize_tier3_hcr, echo=FALSE, fig.width=7, fig.height=4---------------
library(ggplot2)
ssb <- 0:400
Fs <- sapply(ssb, \(x) npfmc_tier3_F(x, B_ref=119.960, F_ref=0.086))
data <- data.frame(SSB=ssb, F=Fs)

ggplot(data)+
    geom_line(aes(x=SSB, y=F), linewidth=1)+
    geom_vline(xintercept = 119.960, linetype="dashed")+
    geom_hline(yintercept = 0.086, linetype="dashed")+
    scale_y_continuous(limits=c(0, 0.12))+
    ggtitle("NPFMC Tier 3a Harvest Control Rule")+
    coord_cartesian(expand=0)+
    theme_bw()


## ----eval=FALSE---------------------------------------------------------------
#  mp <- list(
#     hcr = list(
#         func = hcr_func, # the HCR function to evaluate
#         extra_pars = NA, # any additional parameter values that are needed by the HCR function
#         extra_options = list(
#             max_stability = NA, # a maximum permissable annual change in TACs (as a proportion)
#             harvest_cap = NA    # a maximum permissable annual TAC
#         ),
#         units = "F" # the output units of the HCR function (should usually be "F" or "TAC")
#     )
#     ref_points = list(
#         spr_target = 0.40   # the target SPR level for SPR-based reference points
#     ),
#     management = list(
#         abc_tac_reduction = 1, # proportion of the ABC to set the TAC too
#         tac_land_reduction = 1 # proprtion of the TAC to set realized landings too
#     )
# )

## ----tac_land_reduction, eval=FALSE-------------------------------------------
# mp$management$tac_land_reduction <- list(
#     func = stairstep_attainment,
#     pars = list(
#         breakpoints = c(20, 30),
#         levels = c(0.9, 0.8, 0.6),
#         phase_ins = 0
#     )
# )

## ----stepwise_hcr_example, eval=FALSE-----------------------------------------
# 
# step_hcr <- function(ref_pts, naa, dem_params){
#     # Calculate ssb as naa*waa*mat (female only) in units of 1,000s tons
#     ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
#     return(ifelse(ssb > 100, 0.10, 0.01))
# }
# 
# hcr <- list(
#     func = step_hcr, # the HCR function to evaluate
#     extra_pars = NA, # no extra parameter needed
#     extra_options = list(
#         max_stability = NA, # no stability constraints
#         harvest_cap = NA # no harvest caps
#     ),
#     units = "F"
# )
# 
# mp <- setup_mp_options()
# mp$hcr <- hcr
# 
# mse <- run_mse(om=om, hcr=hcr, nyears_input=100)
# 

## -----------------------------------------------------------------------------
average_age <- function(naa, ages){
    return(weighted.mean(ages, naa))
}

avgage_threshold_f <- function(ref_pts, naa, dem_params, ref_naa, ages){
    # Calculate ssb as naa*waa*mat (female only) in units of 1,000s tons
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)
    as_stat <- average_age(naa, ages)/average_age(ref_naa, ages)
    as_scalar <- threshold_f(as_stat, f_min=0, f_max=1, lrp=0, urp=1)

    x <- ssb/ref_pts$B40
    f_max <- ref_pts$F40*as_scalar
    f_min <- 0
    lrp <- 0
    urp <- 1

    if(x >= urp)  F <- f_max # if stock status >= 1
    if(x > lrp && x < urp) F <- f_max * ((x-lrp)/(urp-lrp)) # if stock status > alpha & stock status < 1
    if(x <= lrp) F <- f_min

    return(f)
}

## ----eval=FALSE---------------------------------------------------------------
# ages <- 2:31
# ref_naa <- 25*sapply(1:30, \(i) 1*(exp(-0.1))^i)
# 
# hcr <- list(
#     func = avgage_threshold_f, # the HCR function to evaluate
#     extra_pars = list(
#         # extra parameter values required by the "avgage_threshold_f function"
#         ref_naa = ref_naa,
#         ages = age
#     ),
#     extra_options = list(
#         max_stability = NA, # a maximum permissable annual change in TACs (as a proportion)
#         harvest_cap = NA # a maximum permissable annual TAC
#     ),
#     units = "F"
# )
# 
# mp <- setup_mp_options()
# mp$hcr <- hcr
# 
# mse <- run_mse(om=om, hcr=mp, nyears_input=100)
# 

## ----mngmt_rules_list, eval=FALSE---------------------------------------------
# 
# hcr$extra_options = list(
#     max_stability = 0.15,   # maximum allowed percent increase in TAC between successive years
#     harvest_cap = 100000    # maximum allowed TAC
# )
# 

