rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SpatialSablefishAssessment)
library(tictoc)
library(doParallel)

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

devtools::load_all(afscOM_dir)

lapply(list.files("R", full.names = TRUE), source)

#' 1. Set up the OM by defining demographic parameters
#' model options (such as options governing the observation
#' processes), and OM initial conditons
nyears <- 110

sable_om <- readRDS("data/sablefish_om_big.RDS") # Read this saved OM from a file
sable_om$model_options$fleet_apportionment <- matrix(c(0.80, 0.20), nrow=nrow(sable_om$model_options$fleet_apportionment), ncol=2, byrow=TRUE)
sable_om$model_options$obs_pars$catch_cv <- c(1e-5, 1e-5, 1e-5, 1e-5)

# Source all available OM and HCR objects
source("dev/oms.R")
source("dev/hcrs.R")

# om_crash2 <- om_immcrash_recruit
# om_crash2$recruitment$pars$crash_start_year <- 20
# om_crash2$recruitment$pars$crash_length <- 30
# om_crash2$name <- "Crash Recruitment"

om_crash2 <- om_bhcyclic_recruit
om_crash2$recruitment$pars$h <- c(0.85, 0.85)
om_crash2$recruitment$pars$R0 <- c(15, 3.5)
om_crash2$recruitment$pars$regime_length <- c(25, 30)
om_crash2$name <- "Crash Recruitment"

om_cycle_low <- om_bhcyclic_recruit
om_cycle_low$recruitment$pars$h <- c(0.85,0.85)
om_cycle_low$recruitment$pars$R0 <- c(50, 5.5)
om_cycle_low$recruitment$pars$regime_length <- c(5, 20)
om_cycle_low$name <- "Low Regime Recruitment"

mp_f40 # Complete
mp_f50 # Complete
mp_20cap # Complete
mp_10perc # Complete
mp_10perc_up # Complete

mp_nefmc_ramp <- mp_f40
mp_nefmc_ramp$name <- "NEFMC Ramp"
mp_nefmc_ramp$hcr$func <- nefmc_ramp

mp_pfmc4010_pstar <- mp_pfmc4010
mp_pfmc4010_pstar$name <- "PFMC 40-10 P*"
mp_pfmc4010_pstar$hcr$func <- pfmc4010_pstar

#' 3. Run the closed-loop MSE simulation
#' A single MSE simulation can be run using the `run_mse(...)`
#' function, while multiple MSE simulations can be run (serially)
#' using the `run_mse_multiple(...)` function.
#' 
#' It is recommended to always use `run_mse_multiple(...)` even
#' when only a single MSE simulation is required.
#' 

set.seed(895)
nsims <- 200
seed_list <- sample(1:(100*nsims), nsims)  # Draw 10 random seeds

mse_options_base <- setup_mse_options()
mse_options <- mse_options_base
mse_options$n_spinup_years <- 54
mse_options$recruitment_start_year <- 54
mse_options$n_proj_years <- 75

mse_options_list <- listN(mse_options)


om_list <- listN(om_bh_recruit, om_cycle_low, om_crash2)
hcr_list <- listN(
    mp_15cap,
    mp_25cap,
    mp_5perc,
    mp_10perc,
    mp_pfmc4010
)


# tic()
model_runs <- run_mse_multiple(
    om_list, 
    hcr_list, 
    seed_list,
    nyears=100,
    mse_options_list=mse_options_list,
    diagnostics = FALSE,
    save=TRUE
)

filetype <- ".jpeg"
figures_dir <- file.path("~/Desktop/Presentations/jan2026_nsaw/figures")#file.path(here::here(), "figures")
width_small <- 12
height_small <- 8

# all_hcrs <- c("F40", "F40 +/- 10% Up", "20k Harvest Cap", "F40 Hybrid", "F50", "F50 +10% Up", "20k Harvest Cap F50", "F50 Hybrid", "No Fishing")
# all_hcr_names <- c("F40", "F40 Stability Constraint", "F40 Harvest Cap", "F40 Hybrid", "F50", "F50 Stability Constraint", "F50 Harvest Cap", "F50 Hybrid", "No Fishing")
rp_hcrs <- c("F40", "F50")
stab_hcrs <- c("F40", "F40 +/- 5%", "F40 +/- 10%")
cap_hcrs <- c("F40", "15k Harvest Cap", "20k Harvest Cap")
hyb_hcrs <- c("F40", "F40 Hybrid")

all_hcrs <- unique(c(rp_hcrs, stab_hcrs, cap_hcrs, hyb_hcrs))


# publication_hcrs <- c("F40", "F40 +/- 10% Up", "20k Harvest Cap", "F40 Hybrid", "F50")
# publication_hcr_names <- all_hcr_names[match(publication_hcrs, all_hcrs)]
# supplementary_hcrs <- c("F40", "F50", "F50 +10% Up", "20k Harvest Cap F50", "F50 Hybrid", "No Fishing")
# supplementary_hcr_names <- all_hcr_names[match(supplementary_hcrs, all_hcrs)]

all_oms <- c("Beverton-Holt Recruitment","Low Regime Recruitment", "Crash Recruitment")
all_om_names <- c("Equilibrium", "Regimes", "Crash")

# publication_oms <- c("Beverton-Holt Recruitment", "Low Regime Recruitment", "Crash Recruitment")
# publication_om_names <- all_om_names[match(publication_oms, all_oms)]
# supplementary_oms <- c("Random Recruitment", "Beverton-Holt Cyclic Recruitment", "Immediate Crash Recruitment")
# supplementary_om_names <- all_om_names[match(supplementary_oms, all_oms)]


mse_runs <- get_saved_model_runs(om_order=all_oms, hcr_order=all_hcrs)
# mse_runs <- readRDS(file.path(here::here(), "data", "zahneretal2025_hybrid_mseruns_FINAL.RDS"))
model_runs <- mse_runs$model_runs
extra_columns <- mse_runs$extra_columns2
extra_columns <- extra_columns %>% 
    mutate(
        om = factor(om, levels = all_oms, labels = all_om_names),
        hcr = factor(hcr, levels = all_hcrs, labels = all_hcrs)
    )


interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 129)

hcr_colors <- set_hcr_colors2(rp_hcrs)
names(hcr_colors) <- rp_hcrs

om_filter = all_om_names

##############################
## Reference Point Plots
##############################
hcr_filter = rp_hcrs
hcr_colors <- set_hcr_colors2(hcr_filter)
names(hcr_colors) <- hcr_filter

ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)

rp_data <- expand.grid(hcr=hcr_filter, om=om_filter, L1=c("caa", "naa")) 
rp_data <- rp_data %>%
    mutate(
        hcr = factor(hcr, levels=hcr_filter, labels=hcr_filter),
        om = factor(om, levels=om_filter, labels=om_filter),
        L1 = factor(L1, labels=c("Landed Catch", "SSB")),
        RP = case_when(
            L1 == "Landed Catch" ~ 18538/1000,
            L1 == "SSB" ~ 105.935
        )
    )

traj <- plot_ssb_catch(
    ssb_data,
    catch_data,
    v1="hcr",
    v2="om",
    common_trajectory = common_trajectory
) + 
    geom_hline(
        data=rp_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter), 
        aes(yintercept=RP), 
        linetype="dashed", 
        color="grey50"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10)),
            scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))
        )
    )+
    guides(color=guide_legend("Harvest\nControl\nRule", nrow = 1))
ggsave(file.path(figures_dir, paste0("refpt_ssb_catch", filetype)), width=16, height=9)

metric_names <- c("Annual Catch", "Total Catch", "Catch AAV", "SSB", "Average Age", "Years on HCR Ramp")
performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "tot_catch", "avg_variation", "avg_ssb", "time_on_ramp", "avg_age")
)

perf_data <- performance_metrics$perf_data %>%
    mutate(
        name = factor(
            name, 
            levels=metric_names, 
            labels=c("Average Annual Catch", "Total Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp")
        )
    ) %>%
    filter(hcr != "No Fishing") %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

formatted_metric_names <- c(
    "Average Annual\nCatch",
    # "Total\nCatch", 
    "Catch\nAAV", 
    "Average Annual\nSSB", 
    "Average\nAge", 
    "Average Years\nSSB < Bref"
)

labels = data.frame(
    x = 1:5,
    y = 1.25,
    labels = formatted_metric_names
)

yaxs <- data.frame(
    x=2.5,
    y=seq(0, 1, 0.25)
)

radar <- perf_data %>% select(om, hcr, name, scaled) %>%
    filter(name != "Total Catch") %>%
    mutate(
        name=factor(
            name, 
            levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp"), 
            labels=formatted_metric_names
        ),
        title="Performance"
    ) %>%
    ggplot()+
        geom_point(aes(x=name, y=scaled, color=hcr, group=hcr), size=3)+
        geom_line(aes(x=name, y=scaled, color=hcr, group=hcr))+
        geom_polygon(aes(x=name, y=scaled, color=hcr, group=hcr), alpha=0)+
        geom_text(data=labels, aes(x=x, y=y, label=labels), size=4.8)+
        geom_text(data=yaxs, aes(x=x, y=y, label=y), size=4)+
        scale_y_continuous("", breaks=seq(0, 1, 0.25))+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(ylim=c(0, 1))+
        labs(x="", y="")+
        ggiraphExtra::coord_radar()+
        facet_grid(title~om)+
        guides(color="none")+#guide_legend(title="Harvest\nControl\nRule", nrow=1))+
        custom_theme+
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )
ggsave(file.path(figures_dir, paste0("refpt_radar", filetype)), width=16, height=6)

traj+guides(color="none") + radar + plot_layout(ncol=1, guides="collect", heights=c(1.5, 1)) & theme(legend.position='bottom')

ggsave(file.path(figures_dir, paste0("refpt_ssb_catch_radar", filetype)), width=12, height=12)

##############################
## Stability Constraint Plots
##############################
hcr_filter = stab_hcrs
hcr_colors <- set_hcr_colors2(hcr_filter)
names(hcr_colors) <- hcr_filter

ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)

rp_data <- expand.grid(hcr=hcr_filter, om=om_filter, L1=c("caa", "naa")) 
rp_data <- rp_data %>%
    mutate(
        hcr = factor(hcr, levels=hcr_filter, labels=hcr_filter),
        om = factor(om, levels=om_filter, labels=om_filter),
        L1 = factor(L1, labels=c("Landed Catch", "SSB")),
        RP = case_when(
            L1 == "Landed Catch" ~ 18538/1000,
            L1 == "SSB" ~ 105.935
        )
    )

traj <- plot_ssb_catch(
    ssb_data,
    catch_data,
    v1="hcr",
    v2="om",
    common_trajectory = common_trajectory
) + 
    geom_hline(
        data=rp_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter), 
        aes(yintercept=RP), 
        linetype="dashed", 
        color="grey50"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10)),
            scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))
        )
    )+
    guides(color=guide_legend("Harvest\nControl\nRule", nrow = 1))
ggsave(file.path(figures_dir, paste0("stab_ssb_catch", filetype)), width=16, height=9)

metric_names <- c("Annual Catch", "Total Catch", "Catch AAV", "SSB", "Average Age", "Years on HCR Ramp")
performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "tot_catch", "avg_variation", "avg_ssb", "time_on_ramp", "avg_age")
)

perf_data <- performance_metrics$perf_data %>%
    mutate(
        name = factor(
            name, 
            levels=metric_names, 
            labels=c("Average Annual Catch", "Total Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp")
        )
    ) %>%
    filter(hcr != "No Fishing") %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

formatted_metric_names <- c(
    "Average Annual\nCatch",
    # "Total\nCatch", 
    "Catch\nAAV", 
    "Average Annual\nSSB", 
    "Average\nAge", 
    "Average Years\nSSB < Bref"
)

labels = data.frame(
    x = 1:5,
    y = 1.25,
    labels = formatted_metric_names
)

yaxs <- data.frame(
    x=2.5,
    y=seq(0, 1, 0.25)
)

radar <- perf_data %>% select(om, hcr, name, scaled) %>%
    filter(name != "Total Catch") %>%
    mutate(
        name=factor(
            name, 
            levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp"), 
            labels=formatted_metric_names
        ),
        title="Performance"
    ) %>%
    ggplot()+
        geom_point(aes(x=name, y=scaled, color=hcr, group=hcr), size=3)+
        geom_line(aes(x=name, y=scaled, color=hcr, group=hcr))+
        geom_polygon(aes(x=name, y=scaled, color=hcr, group=hcr), alpha=0)+
        geom_text(data=labels, aes(x=x, y=y, label=labels), size=4.8)+
        geom_text(data=yaxs, aes(x=x, y=y, label=y), size=4)+
        scale_y_continuous("", breaks=seq(0, 1, 0.25))+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(ylim=c(0, 1))+
        labs(x="", y="")+
        ggiraphExtra::coord_radar()+
        facet_grid(title~om)+
        guides(color="none")+#guide_legend(title="Harvest\nControl\nRule", nrow=1))+
        custom_theme+
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )
ggsave(file.path(figures_dir, paste0("stab_radar", filetype)), width=16, height=6)

traj+guides(color="none") + radar + plot_layout(ncol=1, guides="collect", heights=c(1.5, 1)) & theme(legend.position='bottom')

ggsave(file.path(figures_dir, paste0("stab_ssb_catch_radar", filetype)), width=12, height=12)


##############################
## Harvest Cap Plots
##############################
hcr_filter = cap_hcrs
hcr_colors <- set_hcr_colors2(hcr_filter)
names(hcr_colors) <- hcr_filter

ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)

rp_data <- expand.grid(hcr=hcr_filter, om=om_filter, L1=c("caa", "naa")) 
rp_data <- rp_data %>%
    mutate(
        hcr = factor(hcr, levels=hcr_filter, labels=hcr_filter),
        om = factor(om, levels=om_filter, labels=om_filter),
        L1 = factor(L1, labels=c("Landed Catch", "SSB")),
        RP = case_when(
            L1 == "Landed Catch" ~ 18538/1000,
            L1 == "SSB" ~ 105.935
        )
    )

traj <- plot_ssb_catch(
    ssb_data,
    catch_data,
    v1="hcr",
    v2="om",
    common_trajectory = common_trajectory
) + 
    geom_hline(
        data=rp_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter), 
        aes(yintercept=RP), 
        linetype="dashed", 
        color="grey50"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10)),
            scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))
        )
    )+
    guides(color=guide_legend("Harvest\nControl\nRule", nrow = 1))
ggsave(file.path(figures_dir, paste0("cap_ssb_catch", filetype)), width=16, height=9)

metric_names <- c("Annual Catch", "Total Catch", "Catch AAV", "SSB", "Average Age", "Years on HCR Ramp")
performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "tot_catch", "avg_variation", "avg_ssb", "time_on_ramp", "avg_age")
)

perf_data <- performance_metrics$perf_data %>%
    mutate(
        name = factor(
            name, 
            levels=metric_names, 
            labels=c("Average Annual Catch", "Total Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp")
        )
    ) %>%
    filter(hcr != "No Fishing") %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

formatted_metric_names <- c(
    "Average Annual\nCatch",
    # "Total\nCatch", 
    "Catch\nAAV", 
    "Average Annual\nSSB", 
    "Average\nAge", 
    "Average Years\nSSB < Bref"
)

labels = data.frame(
    x = 1:5,
    y = 1.25,
    labels = formatted_metric_names
)

yaxs <- data.frame(
    x=2.5,
    y=seq(0, 1, 0.25)
)

radar <- perf_data %>% select(om, hcr, name, scaled) %>%
    filter(name != "Total Catch") %>%
    mutate(
        name=factor(
            name, 
            levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp"), 
            labels=formatted_metric_names
        ),
        title="Performance"
    ) %>%
    ggplot()+
        geom_point(aes(x=name, y=scaled, color=hcr, group=hcr), size=3)+
        geom_line(aes(x=name, y=scaled, color=hcr, group=hcr))+
        geom_polygon(aes(x=name, y=scaled, color=hcr, group=hcr), alpha=0)+
        geom_text(data=labels, aes(x=x, y=y, label=labels), size=4.8)+
        geom_text(data=yaxs, aes(x=x, y=y, label=y), size=4)+
        scale_y_continuous("", breaks=seq(0, 1, 0.25))+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(ylim=c(0, 1))+
        labs(x="", y="")+
        ggiraphExtra::coord_radar()+
        facet_grid(title~om)+
        guides(color="none")+#guide_legend(title="Harvest\nControl\nRule", nrow=1))+
        custom_theme+
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )
ggsave(file.path(figures_dir, paste0("cap_radar", filetype)), width=16, height=6)

traj+guides(color="none") + radar + plot_layout(ncol=1, guides="collect", heights=c(1.5, 1)) & theme(legend.position='bottom')

ggsave(file.path(figures_dir, paste0("cap_ssb_catch_radar", filetype)), width=12, height=12)



##############################
## Hybrid Plots
##############################
hcr_filter = hyb_hcrs
hcr_colors <- set_hcr_colors2(hcr_filter)
names(hcr_colors) <- hcr_filter

ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)

rp_data <- expand.grid(hcr=hcr_filter, om=om_filter, L1=c("caa", "naa")) 
rp_data <- rp_data %>%
    mutate(
        hcr = factor(hcr, levels=hcr_filter, labels=hcr_filter),
        om = factor(om, levels=om_filter, labels=om_filter),
        L1 = factor(L1, labels=c("Landed Catch", "SSB")),
        RP = case_when(
            L1 == "Landed Catch" ~ 18538/1000,
            L1 == "SSB" ~ 105.935
        )
    )

traj <- plot_ssb_catch(
    ssb_data,
    catch_data,
    v1="hcr",
    v2="om",
    common_trajectory = common_trajectory
) + 
    geom_hline(
        data=rp_data %>% filter_hcr_om(hcrs=hcr_filter, oms=om_filter), 
        aes(yintercept=RP), 
        linetype="dashed", 
        color="grey50"
    )+
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10)),
            scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))
        )
    )+
    guides(color=guide_legend("Harvest\nControl\nRule", nrow = 1))
ggsave(file.path(figures_dir, paste0("hyb_ssb_catch", filetype)), width=16, height=9)

metric_names <- c("Annual Catch", "Total Catch", "Catch AAV", "SSB", "Average Age", "Years on HCR Ramp")
performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "tot_catch", "avg_variation", "avg_ssb", "time_on_ramp", "avg_age")
)

perf_data <- performance_metrics$perf_data %>%
    mutate(
        name = factor(
            name, 
            levels=metric_names, 
            labels=c("Average Annual Catch", "Total Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp")
        )
    ) %>%
    filter(hcr != "No Fishing") %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

formatted_metric_names <- c(
    "Average Annual\nCatch",
    # "Total\nCatch", 
    "Catch\nAAV", 
    "Average Annual\nSSB", 
    "Average\nAge", 
    "Average Years\nSSB < Bref"
)

labels = data.frame(
    x = 1:5,
    y = 1.25,
    labels = formatted_metric_names
)

yaxs <- data.frame(
    x=2.5,
    y=seq(0, 1, 0.25)
)

radar <- perf_data %>% select(om, hcr, name, scaled) %>%
    filter(name != "Total Catch") %>%
    mutate(
        name=factor(
            name, 
            levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp"), 
            labels=formatted_metric_names
        ),
        title="Performance"
    ) %>%
    ggplot()+
        geom_point(aes(x=name, y=scaled, color=hcr, group=hcr), size=3)+
        geom_line(aes(x=name, y=scaled, color=hcr, group=hcr))+
        geom_polygon(aes(x=name, y=scaled, color=hcr, group=hcr), alpha=0)+
        geom_text(data=labels, aes(x=x, y=y, label=labels), size=4.8)+
        geom_text(data=yaxs, aes(x=x, y=y, label=y), size=4)+
        scale_y_continuous("", breaks=seq(0, 1, 0.25))+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(ylim=c(0, 1))+
        labs(x="", y="")+
        ggiraphExtra::coord_radar()+
        facet_grid(title~om)+
        guides(color="none")+#guide_legend(title="Harvest\nControl\nRule", nrow=1))+
        custom_theme+
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
        )
ggsave(file.path(figures_dir, paste0("hyb_radar", filetype)), width=16, height=6)

traj+guides(color="none") + radar + plot_layout(ncol=1, guides="collect", heights=c(1.5, 1)) & theme(legend.position='bottom')

ggsave(file.path(figures_dir, paste0("hyb_ssb_catch_radar", filetype)), width=12, height=12)



##############################
## Resiliency Plots
##############################
hcr_filter <- all_hcrs

resilience_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=hcr_filter,
    om_filter=c("Crash"),
    interval_widths=interval_widths,
    time_horizon = c(time_horizon[1]+25, time_horizon[2]), 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("crash_time", "recovery_time")
)

resilience_data <- resilience_metrics$perf_data %>%
    mutate(hcr = factor(hcr, levels=rev(all_hcrs))) %>%
    group_by(.width, om, name) %>%
    mutate(
        scaled = ifelse(
                    name %in% c(
                        "Catch AAV", 
                        "Proportion of Years SSB < B35", 
                        "Average Years on HCR Ramp",
                        "Recovery Time"
                    ), 
                    inf_min(eval(rlang::parse_expr(col_name)))/eval(rlang::parse_expr(col_name)), 
                    eval(rlang::parse_expr(col_name))/inf_max(eval(rlang::parse_expr(col_name)))
                ),
    ) %>%
    arrange(desc(scaled), .by_group=TRUE) %>%
    mutate(
        rank = factor(row_number())
    ) %>%
    print(n=100)

plot_performance_metric_summary(resilience_data)+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("Year", limits=c(0, 30)),
            scale_x_continuous("Year", limits=c(0, 20))
        )
    )+
    theme(
        strip.text.y = element_blank(),
        plot.margin = margin(10, 30, 10, 10),
        panel.spacing.x = unit(30, "pt")
    )
ggsave(file.path(figures_dir, paste0("resilience_crash_recovery", filetype)), width=12, height=6)


##############################
## Recruitmen Plots
##############################
hcr_filter <- c("F40")

rec <- get_recruits(model_runs, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)
examp_rec <- rec %>% 
    as_tibble() %>%
    filter(time > 54, sim %in% sample(seed_list, 5)) %>%
    mutate(
        sim = factor(sim)
    )

mean_rec <- examp_rec %>% group_by(om) %>% summarise(r=median(rec))

summ_rec <- rec %>%
    as_tibble() %>%
    filter(time > 54) %>%
    group_by(time, om) %>%
    median_qi(rec, .width=interval_widths)

ggplot(examp_rec) +
    geom_line(aes(x=time, y=rec, color=om, group=sim), size=0.5, alpha=0.4)+
    geom_line(data=summ_rec, aes(x=time, y=rec), color="black", size=1)+
    geom_hline(data=mean_rec, aes(yintercept=r), linetype="dashed")+
    # geom_text(data=mean_rec, aes(x=123, y=115, label=round(r, 2)), size=6)+
    scale_y_continuous(limits=c(0, 120))+
    scale_x_continuous(breaks=c(seq(55, 130, 20), 134), labels=c(seq(0, 75, 20), 75))+
    guides(color="none")+
    labs(x="Time", y="Recruitment (millions)")+
    coord_cartesian(expand=0)+
    facet_wrap(~om, ncol=3)+
    theme_bw()+
    custom_theme
ggsave(file.path(figures_dir, paste0("recruitment_trajectories", filetype)), width=12, height=5)
