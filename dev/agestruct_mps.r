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

age_structure_hcr <- function(ref_pts, naa, dem_params, avgrec, ref_naa, desired_abi, adjustment, cutoff_age){
    # Compute spawning biomass
    nages <- afscOM::get_model_dimensions(dem_params$sel)$nages
    a <- cutoff_age-1
    ssb <- apply(naa[,,1,]*dem_params$waa[,,1,,drop=FALSE]*dem_params$mat[,,1,,drop=FALSE], 1, sum)

    abi_y <- abi(naa[,,1,], ref_naa, threshold = 0.90, start_age = cutoff_age)

    # adj <- ifelse(abi_y < desired_abi, 1/adjustment, adjustment)

    abi_adjustment <- ifelse(
        abi_y < desired_abi,
        (1+(((abi_y-desired_abi)/desired_abi)))/adjustment,
        (1+(((abi_y-desired_abi)/desired_abi)/adjustment))
    )
    # abi_adjustment <- (abi_y/desired_abi)^adjustment
    abi_adjustment <- min(2, abi_adjustment)  # Limit adjustment factor between 0.01 and 2

    abi_adjustment <- ifelse(ssb <= ref_pts$Bref, 1, abi_adjustment)

    # Compute B40%
    Fmax <- ref_pts$Fref*abi_adjustment
    
    f <- npfmc_tier3_F(ssb, ref_pts$Bref, Fmax)
    
    return(f)
}

mp_agestruct_061 <- setup_mp_options()
mp_agestruct_061$name <- "AS | ref=0.6 | adj=1"
mp_agestruct_061$hcr$func <- age_structure_hcr
mp_agestruct_061$hcr$extra_pars <- list(
    ref_naa = ref_naa0,
    desired_abi = 0.6,
    adjustment = 1,
    cutoff_age = 3
)
mp_agestruct_061$hcr$units = "F"
mp_agestruct_061$management$tac_land_reduction <- 1

mp_agestruct_041 <- setup_mp_options()
mp_agestruct_041$name <- "AS | ref=0.4 | adj=1"
mp_agestruct_041$hcr$func <- age_structure_hcr
mp_agestruct_041$hcr$extra_pars <- list(
    ref_naa = ref_naa0,
    desired_abi = 0.4,
    adjustment = 1,
    cutoff_age = 3
)
mp_agestruct_041$hcr$units = "F"
mp_agestruct_041$management$tac_land_reduction <- 1

mp_agestruct_065 <- setup_mp_options()
mp_agestruct_065$name <- "AS | ref=0.6 | adj=5"
mp_agestruct_065$hcr$func <- age_structure_hcr
mp_agestruct_065$hcr$extra_pars <- list(
    ref_naa = ref_naa0,
    desired_abi = 0.6,
    adjustment = 5,
    cutoff_age = 3
)
mp_agestruct_065$hcr$units = "F"
mp_agestruct_065$management$tac_land_reduction <- 1

mp_agestruct_045 <- setup_mp_options()
mp_agestruct_045$name <- "AS | ref=0.4 | adj=5"
mp_agestruct_045$hcr$func <- age_structure_hcr
mp_agestruct_045$hcr$extra_pars <- list(
    ref_naa = ref_naa0,
    desired_abi = 0.4,
    adjustment = 5,
    cutoff_age = 3
)
mp_agestruct_045$hcr$units = "F"
mp_agestruct_045$management$tac_land_reduction <- 1


mp_fullmat <- mp_f40
mp_fullmat$name <- "F40 Fully Mature"
mp_fullmat$ref_points$rp_start_age <- 13
mp_fullmat$hcr$extra_pars <- list(
    cutoff_age = 14
)

mp_sprhyper <- mp_f40
mp_sprhyper$name <- "F40 Hyperallometric SPR"
mp_sprhyper$ref_points$rp_hyperallometry <- 3
mp_sprhyper$hcr$extra_pars <- list(
    hyperallometry = 3
)


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
mse_options$run_estimation <- TRUE
mse_options_list <- listN(mse_options)


om_list <- listN(om_bh_recruit, om_cycle_low, om_crash2)

hcr_list <- listN(mp_f40, mp_fullmat, mp_sprhyper, mp_agestruct_061, mp_agestruct_041, mp_agestruct_045)# mp_agestruct_061, mp_agestruct_065, mp_agestruct_041, mp_agestruct_045)

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

# toc()


# Data Processing
filetype <- ".jpeg"
figures_dir <- file.path(here::here(), "figures/agestruct_hcrs")#file.path(here::here(), "figures")
width_small <- 12
height_small <- 8

hcr_names <- unlist(lapply(hcr_list, function(x) x$name))

publication_hcrs <- hcr_names

om_names <- unlist(lapply(om_list, function(x) x$name))

publication_oms <- om_names

# publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "")
# metric_names <- c("Annual Catch", "Catch AAV", "Average SSB", "Average Age", "Proportion of Years SSB < B35*")

# mse_runs <- get_saved_model_runs(om_order=all_oms, hcr_order=all_hcrs)
# mse_runs <- readRDS(file.path(here::here(), "data", "zahneretal2025_hybrid_mseruns_FINAL.RDS"))
# model_runs <- mse_runs$model_runs
# extra_columns <- mse_runs$extra_columns2
# extra_columns <- extra_columns %>% 
    # mutate(
        # om = factor(om, levels = all_oms, labels = all_om_names),
        # hcr = factor(hcr, levels = all_hcrs, labels = all_hcr_names)
    # )

extra_columns <- expand.grid(om=om_names, hcr=hcr_names)

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 129)

hcr_colors <- c("black", "red", "blue", "darkgreen", "purple", "#888888")
names(hcr_colors) <- publication_hcrs
# Specify which subset of HCR/OM combinations to get results for
hcr_filter <- publication_hcrs#c("F40", "F50", as_hcrs)
om_filter <- publication_oms
# extra_columns <- expand.grid(om=om_filter, hcr=hcr_filter)

### Spawning Biomass and Catch Plots
ssb_data <- get_ssb_biomass(NULL, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
catch_data <- get_landed_catch(NULL, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)


rp_data <- expand.grid(hcr=hcr_filter, om=om_filter, L1=c("caa", "naa")) 
rp_data <- rp_data %>%
    mutate(
        hcr = factor(hcr),#, levels=all_hcr_names, labels=all_hcr_names),
        om = factor(om),#, levels=all_om_names, labels=all_om_names),
        L1 = factor(L1, labels=c("Landed Catch", "SSB")),
        RP = case_when(
            L1 == "Landed Catch" ~ 18538/1000, #mean(assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"]),
            L1 == "SSB" ~ 105.935
        )
    )

plot_ssb_catch(
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
    ggh4x:::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10)),
            scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))
        )
    )+
    guides(color=guide_legend("Harvest\nControl\nRule", ncol = 3))
    

ggsave(filename=file.path(figures_dir, "ssb_catch.jpeg"), width=16*(3/3), height=10, units="in")


econ_data <- get_dynamic_economic_value(NULL, extra_columns, hcr_filter, om_filter)
plot_econ_value(econ_data, v1="hcr", v2="om", time_horizon=time_horizon)+
    scale_y_continuous("Relative Economic Value", limits=c(0, 12), breaks=seq(0, 12, 2))+
    guides(color=guide_legend("Harvest\nControl\nRule"))

ggsave(filename=file.path(figures_dir, "econ_value.jpeg"), width=16*(3/3), height=8, units="in")

abi_data <- get_abi(NULL, extra_columns, ref_naa0, hcr_filter, om_filter)
plot_average_age(abi_data, v1="hcr", v2="om", time_horizon=time_horizon)+
    scale_y_continuous("ABI0", limits=c(0, 2), breaks=seq(0, 2, 0.25))+
    geom_hline(yintercept=0.40, linetype="dashed")+
    geom_hline(yintercept=0.60, linetype="dashed")+
    geom_hline(yintercept=1.0, linetype="dashed")+
    guides(color=guide_legend("Harvest\nControl\nRule"))

ggsave(filename=file.path(figures_dir, "abi.jpeg"), width=16*(3/3), height=8, units="in")


### Performance Metrics -----------------
#########################################
metric_names <- c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Dynamic Annual Value")

performance_metrics <- performance_metric_summary(
    NULL, 
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
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "dynamic_value")
)
perf_data <- performance_metrics$perf_data %>%
    mutate(
        name = factor(
            name, 
            levels=metric_names, 
            labels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Annual Economic Value")
        ),
        om = factor(om, levels=publication_oms)
    ) %>%
    filter(hcr != "No Fishing") %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")


formatted_metric_names <- c(
    "Average Annual\nCatch", 
    "Catch\nAAV", 
    "Average Annual\nSSB", 
    "Average\nAge",
    "Average Annual\nEconomic Value"
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
    mutate(name=factor(
        name, 
        levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Annual Economic Value"), 
        labels=formatted_metric_names
    )) %>%
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
        facet_wrap(~om)+
        guides(color=guide_legend(title="Harvest\nControl\nRule", ncol=3))+
        custom_theme+
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
ggsave(filename=file.path(figures_dir, paste0("radar_plots", filetype)), width=16, height=7, units="in")


table <- perf_data %>% filter(.width==0.50) %>%
    mutate(scaled = round(scaled, 2)) %>%
    mutate(
        name=factor(
            name, 
            levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Annual Economic Value"), 
            labels=c(
                "Average Annual\nCatch", 
                "Catch AAV", 
                "Average Annual\nSSB", 
                "Average Age", 
                "Average Annual\nEconomic Value"
            )
        ),
        hcr = factor(
            hcr,
            labels = c("F40", "F40\nFully\nMature", "F40\nHyper\nSPR", "AS_0.6_\n0.01", "AS_0.4_\n0.01", "AS_0.4_\n0.05")
        ) 
    ) %>%
    ggplot()+
        geom_tile(aes(x=hcr, y=name, fill=scaled))+
        geom_text(aes(x=hcr, y=name, label=round(median, 2)), size=6, color="black")+
        scale_fill_gradient(oob=scales::squish, limits=c(0.6, 1), low="white", high="#368536")+
        scale_y_discrete(name="", position="right")+
        scale_x_discrete(name="", position="bottom")+
        facet_wrap(~om, nrow=1, strip.position="left")+
        coord_cartesian(expand = 0)+
        custom_theme+
        guides(fill="none")+
        theme(
            strip.text.y = element_blank(),
            legend.postion = "top"
        )

library(patchwork)
radar/table + plot_layout(guides="collect") & theme(legend.position="top")
ggsave(filename=file.path(figures_dir, "performance.jpeg"), width=18, height=14, units="in")

rel_performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    dem_params = sable_om$dem_params, 
    ref_naa = ref_naa,
    hcr_filter=hcr_filter,
    om_filter=om_filter,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative="F40", 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "avg_ssb")
)
rel_perf_data <- rel_performance_metrics$perf_data %>%
    mutate(
        name = factor(
            name, 
            levels=c("Annual Catch", "Average Annual SSB"), 
            labels=c("Average Annual Catch", "Average Annual SSB")
        ),
        hcr = factor(
            hcr,
            levels=hcr_filter,
            labels=hcr_filter
        )
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

catchssb_perf <- ggplot(
        rel_perf_data
    )+
    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=hcr, color=hcr, shape=om), point_size=3, position="dodge")+
    scale_shape_discrete()+
    scale_color_manual(values=hcr_colors)+
    labs(y="", x="", shape="OM", color="Relative MS Order")+
    coord_cartesian(expand=0)+
    guides(shape="none", color=guide_legend(nrow=1))+
    facet_grid(rows=vars(om), cols=vars(name), scales="free_x")+
    custom_theme+
    theme(
        plot.margin = margin(0.25, 1, 0.25, 0.25, "cm"),
        panel.spacing.x = unit(5, "cm"),
        plot.title = element_text(size=18),
        legend.spacing.x = unit(1.5, "cm")
    )+
    theme(
        panel.spacing.x = unit(0.75, "cm"),
    )+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous(limits=c(0, 1.5)),
            scale_x_continuous(limits=c(0, 2), breaks=c(0, 0.5, 1.0, 1.5, 2.0))
        )
    )


library(grid)
library(gtable)
library(ggplotify)
library(ggpmisc)

g1 <- ggplotGrob(catchssb_perf)

insets <- expand.grid(om=om_filter, L1=c("Landed Catch", "SSB"))
insets$x <- 125
insets$y <- c(rep(58, 3), rep(360, 3))
insets$plot <- list(
    as.ggplot(g1[c(10, 15), c(7)]),
    as.ggplot(g1[c(12, 15), c(7)]),
    as.ggplot(g1[c(14, 15), c(7)]),
    as.ggplot(g1[c(10, 15), c(9)]),
    as.ggplot(g1[c(12, 15), c(9)]),
    as.ggplot(g1[c(14, 15), c(9)])
)


plot_ssb_catch(
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
    ggpp::geom_plot(data=insets, aes(x=x, y=y, label=plot))+
    guides(color=guide_legend("Harvest\nControl\nRule", ncol = 3))+
    ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 370))
            )
        )

ggsave(filename=file.path(figures_dir, paste0("ssb_catch_insets", filetype)), width=16, height=10, units="in")


dem_params = sable_om$dem_params 
ref_naa = ref_naa0
extra_filter = NULL
relative=NULL
summarise_by=c("om", "hcr")
summary_out = TRUE

group_columns <- c("sim", summarise_by)


naas <- bind_mse_outputs(model_runs, "naa", extra_columns) %>%
        as_tibble() %>%
        ungroup() %>%
        filter_hcr_om(hcr_filter, om_filter) %>%
        filter_times(time_horizon=time_horizon) %>%
        group_by(time, age, sim, om, hcr) %>%
        mutate(value = sum(value)) %>%
        # round_to_zero("value") %>%
        filter(sex == "F") %>%
        ungroup()

avg_abi <- naas %>%  
        group_by(time, sim, hcr, om) %>%
        summarise(
            avg_abi = abi(value, ref_naa0)
        ) %>%
        relativize_performance(
            rel_column = "hcr",
            value_column = "avg_abi",
            rel_value = NULL,
            grouping = c(group_columns, "time")
        )
        
if(!is.null(extra_filter)){
    avg_abi <- avg_abi %>% filter(eval(extra_filter))
}

out <- avg_abi
if(summary_out){
    avg_abi <- avg_abi %>%
        group_by(across(all_of(c(summarise_by, "time")))) %>%
        median_qi(avg_abi, .width=interval_widths, .simple_names=FALSE)
}

ggplot(avg_abi)+
    geom_line(aes(x=time, y=avg_abi, color=hcr, group=hcr))+
    scale_color_manual(values=hcr_colors)+
    facet_wrap(~om)+
    custom_theme



avg_abi %>% 
    filter_hcr_om(hcrs=as_hcrs, oms=publication_om_names) %>% 
    filter(.width==0.5) %>% 
    select(om, hcr, time, avg_abi) %>%
    separate(hcr, into=c("type", "ref", "exp"), sep=" \\| ", remove=FALSE) %>%
    separate(ref, into=c("ref", "ref_value"), sep="=", remove=TRUE) %>%
    separate(exp, into=c("exp", "exp_value"), sep="=", remove=TRUE) %>%
    select(-type, -ref, -exp) %>%
    mutate(
        ref_value = as.numeric(ref_value),
        exp_value = as.numeric(exp_value),
        adjustment = (avg_abi/ref_value)^exp_value
    ) %>%

    ggplot()+
        geom_point(aes(x=avg_abi, y=adjustment, color=hcr))+
        scale_color_manual(values=c("red", "orange", "blue", "green", "purple", "black"))+
        facet_wrap(~om, scales="free")+
        custom_theme
    













full_performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=publication_oms,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = FALSE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "crash_time", "recovery_time") 
)

agg_perf_data <- full_performance_metrics$perf_data %>%
    ungroup() %>%
    pivot_longer(4:9, names_to="name", values_to="value") %>%
    mutate(
        name = factor(name, levels=c("annual_catch", "aav", "spbio", "avg_age", "crash_time", "recovery_time"), labels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Crash Time", "Recovery Time"))
    ) %>%
    filter(!(om %in% c("Random Recruitment", "Beverton-Holt Regime Recruitment") & name %in% c("Crash Time", "Recovery Time"))) %>%
    group_by(hcr, name) %>%
    summarise(median=median(value, na.rm=TRUE)) %>%
    group_by(name) %>%
    mutate(om="Aggregated") %>%
    select(om, hcr, name, median) 


### Aggregate Performance --------------
#########################################
full_performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=publication_oms,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = FALSE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "crash_time", "recovery_time") 
)

om_aggregated_performance <- full_performance_metrics$perf_data %>%
    group_by(sim, hcr) %>%
    summarise(
        across(2:7, \(x) mean(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(3:8, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "crash_time", "recovery_time"), 
            labels=c("Average Annual Catch", "Catch AAV", "Average SSB", "Average Age", "Crash Time", "Recovery Time")
        ),
        base = ifelse(grepl("F40", hcr), "F40", "F50")
    ) 

### Pairs tradeoff plots -----------------------
################################################
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(size=2, alpha=0.65) +
    geom_smooth(method = method, color = "black", fullrange=TRUE, ...)
  p
}

diagonalFn <- function(data, mapping, ...){
    p <- ggplot(data=data, mapping=mapping)+
        geom_histogram(aes_string(fill="hcr"), linetype="solid", alpha=0.25)
    p
}

upperFn <- function(data, mapping, method = "lm", ...){
    p <- ggplot(data=data, mapping=mapping)+
        geom_smooth(aes_string(color="hcr", linetype="base"), method=method, fullrange=FALSE, ...)
    p
}

GGally::ggpairs(
    om_aggregated_performance %>% 
        pivot_wider(names_from="name", values_from="value") %>%
        filter(hcr != "No Fishing"), 
    mapping=aes(color=hcr, shape=base, linetype=base),
    columns=4:9,
    lower = list(continuous = GGally::wrap(lowerFn, method = "lm", se=FALSE, linewidth=0.5)),
    upper="blank",
    # diag = list(continuous="barDiag"),
    diag = list(continuous = GGally::wrap(diagonalFn)),
    legend=7
)+
scale_color_manual(values=hcr_colors)+
scale_fill_manual(values=hcr_colors)+
scale_shape_manual(values=c(20, 5))+
theme(legend.position = "bottom")+
custom_theme+
guides(
    color=guide_legend("Harvest\nControl\nRule", ncol=4),
    shape=guide_legend("Base\nHCR", nrow=2),
    linetype=guide_legend("Base\nHCR", nrow=2)
)

ggsave(filename=file.path(figures_dir, paste0("performance_frontier_age", filetype)), width=15, height=15, units="in")

################################
### TOPSIS ---------------------
################################
library(MCDA)

compute_topsis <- function(perf_data, topsis_splits, topsis_weights, topsis_minmax){

    new_names <- paste0("Var",1:length(topsis_splits))
    names(new_names) <- topsis_splits

    perf_data %>%
        ungroup() %>%
        # select(c(topsis_splits, hcr, name, value)) %>%
        group_by(across(all_of(topsis_splits))) %>%
        group_split() %>%
        purrr::map(function(df){

            om_name <- df %>% distinct(om)

            table <- df %>% select(!topsis_splits) %>% 
                        # pivot_wider(names_from="name", values_from="value") %>%
                        ungroup() %>% column_to_rownames("hcr")

            if(om_name == "Crash Recruitment"){
                table <- table %>% mutate(across(everything(), ~ ifelse(is.na(.x), 0, .x))) %>% as.matrix
                weights <- topsis_weights
                minmax <- topsis_minmax
            }else{
                table <- table %>% select(-crash_time, -recovery_time) %>% as.matrix
                weights <- topsis_weights[1:4]
                minmax <- topsis_minmax[1:4]
            }

            topsis <- MCDA::TOPSIS(table, weights, minmax)

            return(topsis)
        }) %>% 
        bind_rows() %>%
        bind_cols(
            perf_data[topsis_splits] %>% distinct()
        ) %>%
        select(c(topsis_splits, 1:(ncol(.)-length(topsis_splits))))
        # arrange(region)

}

weights <- c(0.25, 0.125, 0.25, 0.125, 0.065, 0.185)
minmax <- c("max", "min", "max", "min", "max", "min")

topsis <- compute_topsis(
    perf_data = full_performance_metrics$perf_data %>% select(-dyn_annual_value) %>% ungroup(), 
    # full_performance_metrics$perf_data %>% 
    #     select(-crash_time, -recovery_time, -dyn_annual_value) %>% 
    #     group_by(sim, hcr) %>%
    #     summarise(across(annual_catch:prop_years, ~ median(.x))),
    topsis_splits = c("sim", "om"),
    topsis_weights = weights,
    topsis_minmax = minmax
)

topsis %>%
    # mutate(row_sum = rowSums(select(., 3:ncol(.)))) %>% 
    # mutate_at(3:ncol(.), ~ ./row_sum) %>% 
    # select(-row_sum) %>% 
    pivot_longer(3:ncol(.), names_to="hcr", values_to="value") %>%
    mutate(
        hcr = factor(hcr, levels=publication_hcrs, labels=publication_hcrs)
    ) %>%
    group_by(om, hcr) %>%
    median_qi(value, .width=interval_widths) %>%

    ggplot(aes(x=value, xmin=.lower, xmax=.upper, color=hcr, y=hcr))+
        geom_pointinterval()+
        # scale_x_continuous(limits=c(0.09, 0.15))+
        facet_wrap(~om)+
        custom_theme

topsis %>%
    # mutate(row_sum = rowSums(select(., 3:ncol(.)))) %>% 
    # mutate_at(3:ncol(.), ~ ./row_sum) %>% 
    # select(-row_sum) %>% 
    pivot_longer(3:ncol(.), names_to="hcr", values_to="value") %>%
    mutate(
        hcr = factor(hcr, levels=publication_hcrs, labels=publication_hcrs)
    ) %>%
    group_by(om, hcr) %>%
    median_qi(value, .width=interval_widths) %>%
    filter(.width == 0.50) %>%
    select(om, hcr, value) %>%
    pivot_wider(names_from="om", values_from="value") %>%
    mutate(
        agg = rowMeans(select(., where(is.numeric)))
    ) %>%
    arrange(desc(agg))


################
# Kobe Plots 
################

phaseplane_data <- get_phaseplane_data(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    hcr_filter=publication_hcrs,
    om_filter=publication_oms
)
rps <- calculate_ref_points(
    nages = 30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = 0.35
)

plot_phase_diagram(phaseplane_data, rps, common_trajectory = 54)


kobe_tseries_data <- phaseplane_data %>% select(-biomass) %>%
    filter_times(time_horizon) %>%
    mutate(
        Frat = total_F/rps$Fref,
        Brat = spbio/rps$Bref
    ) %>%
    mutate(
        quadrant = case_when(
            Frat < 1 & Brat > 1 ~ 3,
            Frat < 1 & Brat < 1 ~ 4,
            Frat > 1 & Brat < 1 ~ 1,
            Frat > 1 & Brat > 1 ~ 2,
            TRUE ~ 0
        ),
    ) %>%
    group_by(time, om, hcr) %>%
    summarise(
        red = sum(quadrant == 1)/n(),
        yellow = sum(quadrant == 4)/n(),
        green = sum(quadrant == 3)/n(),
        orange = sum(quadrant == 2)/n(),
    ) %>%
    pivot_longer(red:orange, names_to="quadrant", values_to="proportion") %>%
    mutate(
        quadrant = factor(quadrant, levels=c("green", "yellow", "orange", "red")),
    )
    # median_qi(quadrant, .width=interval_widths) 


ggplot(kobe_tseries_data)+
    geom_col(aes(x=time, y=proportion, fill=quadrant), color="black", position="stack")+
    scale_fill_manual(values=c("red"="#ff0000", "orange"="#ffa500", "yellow"="#ffff00", "green"="#00ff00"))+
    facet_grid(hcr~om)+
    coord_cartesian(expand=FALSE)+
    custom_theme

ggsave("~/Desktop/kobe_plot.jpeg", width=16, height=12, units="in")


d <- phaseplane_data %>% select(-biomass) %>%
    filter_times(time_horizon) %>%
    mutate(
        Frat = total_F/rps$Fref,
        Brat = spbio/rps$Bref
    ) %>%
    group_by(time, om, hcr) %>%
    median_qi(Frat, Brat, .width=interval_widths)

ggplot(d)+
    geom_point(aes(x=Brat, y=Frat, color=time))+
    geom_line(aes(x=Brat, y=Frat))+
    facet_grid(hcr~om)+
    geom_hline(yintercept=1, linetype="dashed")+
    geom_vline(xintercept=1, linetype="dashed")+
    scale_y_continuous(limits=c(0, 2))+
    scale_x_continuous(limits=c(0, 2))+
    coord_cartesian(expand=FALSE)+
    custom_theme
ggsave("~/Desktop/kobe_scatter.jpeg", width=16, height=12, units="in")


recruit_data <- get_recruits(model_runs, extra_columns, publication_hcrs, publication_oms)
plot_recruitment(recruit_data, v1="hcr", v2="om")
