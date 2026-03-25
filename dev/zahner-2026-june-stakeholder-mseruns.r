rm(list=ls())

library(tidyverse)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SpatialSablefishAssessment)
library(tictoc)
library(doParallel)
library(afscOM)
# Change to wherever your local copy of afscOM is
library(devtools)
# afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

# devtools::load_all(afscOM_dir)

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

mp_f40_hybrid <- mp_f40
mp_f40_hybrid$hcr$extra_options$max_stability <- mp_10perc_up$hcr$extra_options$max_stability
mp_f40_hybrid$hcr$extra_options$harvest_cap <- mp_20cap$hcr$extra_options$harvest_cap
mp_f40_hybrid$name <- "F40 Hybrid"

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


om_list <- listN(om_bh_recruit, om_crash2)
hcr_list <- listN(
    mp_f40,
    mp_f50,
    mp_5perc
    mp_10perc_up,
    mp_15cap,
    mp_25cap,
    mp_f40_hybrid,
    mp_fullmat,
    mp_sprhyper,
    mp_agestruct_041,
    mp_agestruct_061,
    mp_agestruct_065
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

# toc()


# Data Processing
filetype <- ".jpeg"
figures_dir <- file.path(here::here(), "figures/hybrid_hcrs")#file.path(here::here(), "figures")
width_small <- 12
height_small <- 8

base_hcrs <- c("F40", "F50", "F40 +/- 5%", "F40 +/- 10% Up", "15k Harvest Cap", "25k Harvest Cap", "F40 Hybrid")
base_hcr_names <- c("F40", "F50", "F40 5% Stability Constraint", "F40 15k Harvest Cap", "F40 Hybrid")

as_hcrs <- c("F40", "F40 Fully Mature", "F40 Hyperallometric SPR", "AS | ref=0.4 | adj=1", "AS | ref=0.6 | adj=1")
as_hcr_names <- c("F40", "F40 Fully Mature", "F40 Hyper SPR", "AS_0.4_0.01", "AS_0.6_0.01")

oms <- c("Beverton-Holt Recruitment", "Low Regime Recruitment", "Crash Recruitment")
om_names <- c("Beverton-Holt", "Regimes", "Crash")

mse_runs <- get_saved_model_runs(om_order=oms, hcr_order=unique(c(base_hcrs, as_hcrs)))
# mse_runs <- readRDS(file.path(here::here(), "data", "zahneretal2025_hybrid_mseruns_FINAL2.RDS"))
model_runs <- mse_runs$model_runs
extra_columns <- mse_runs$extra_columns2
extra_columns <- extra_columns %>% 
    mutate(
        om = factor(om, levels = oms, labels = om_names),
        hcr = factor(hcr, levels = c(base_hcrs, as_hcrs), labels = c(base_hcr_names, as_hcr_names))
    )


interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 129)

hcr_filter <- as_hcr_names
om_filter <- om_names

hcr_colors <- set_hcr_colors2(hcr_filter)
names(hcr_colors) <- hcr_filter


### Spawning Biomass and Catch Plots
ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=hcr_filter, om_filter=om_filter)
catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=hcr_filter, om_filter=om_filter)

rp_data <- expand.grid(hcr=hcr_filter, om=om_filter, L1=c("caa", "naa")) 
rp_data <- rp_data %>%
    mutate(
        hcr = factor(hcr)#, levels=stakeholder_hcrs, labels=stakeholder_hcrs),
        om = factor(om)#, levels=stakeholder_oms, labels=stakeholder_oms),
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
    ggh4x::facetted_pos_scales(
        y = list(
            scale_y_continuous(limits=c(0, 60), breaks=seq(0, 60, 10)),
            scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 50))
        )
    )+
    guides(color=guide_legend("Harvest\nControl\nRule", nrow = 1))

ggsave(file.path("~/Desktop", "ssb_catch.jpeg"))



metric_names <- c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Years on HCR Ramp", "Dynamic Annual Value")

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
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "time_on_ramp", "avg_age", "dynamic_value")
)

perf_data <- performance_metrics$perf_data %>%
    mutate(
        name = factor(
            name, 
            levels=metric_names, 
            labels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp", "Average Annual Economic Value")
        )
    ) %>%
    filter(hcr != "No Fishing") %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")


formatted_metric_names <- c(
    "Average Annual\nCatch", 
    "Catch\nAAV", 
    "Average Annual\nSSB", 
    "Average Age", 
    "Average\nYears\nSSB <\n Bref", 
    "Average Annual\nEconomic Value"
)

labels = data.frame(
    x = 1:6,
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
        levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp", "Average Annual Economic Value"), 
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
        guides(color=guide_legend(title="Harvest\nControl\nRule", nrow=1))+
        custom_theme+
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )

table <- perf_data %>% filter(.width==0.50) %>%
    mutate(scaled = round(scaled, 2)) %>%
    mutate(
        name=factor(
            name, 
            levels=c("Average Annual Catch", "Catch AAV", "Average Annual SSB", "Average Age", "Average Years on HCR Ramp", "Average Annual Economic Value"), 
            labels=c(
                "Average Annual\nCatch", 
                "Catch AAV", 
                "Average Annual\nSSB", 
                "Average Age", 
                "Average Years\nSSB < Bref", 
                "Average Annual\nEconomic Value"
            )

        ),
        hcr = factor(
            hcr,
            labels = c("F40", "F50", "F40 +10%\nStability\nConstraint", "F40\n15k\nHarvest\nCap", "F40\n25k\nHarvest\nCap", "F40\nHybrid")
        ) 
    ) %>%
    ggplot()+
        geom_tile(aes(x=hcr, y=name, fill=scaled))+
        geom_text(aes(x=hcr, y=name, label=round(median, 2)), size=5, color="black")+
        scale_fill_gradient(oob=scales::squish, limits=c(0.6, 1), low="white", high="#368536")+
        scale_y_discrete(name="", position="right")+
        scale_x_discrete(name="", position="bottom")+
        facet_wrap(~om, nrow=1, strip.position="left")+
        coord_cartesian(expand = 0)+
        custom_theme+
        guides(fill="none")+
        theme(
            strip.text.y = element_blank(),
            legend.postion = "top",
            axis.text = element_text(size=12)
        )

library(patchwork)
radar/table + plot_layout(guides="collect") & theme(legend.position="top")

ggsave(file.path("~/Desktop", "performance.jpeg"), width=12, units="in")


rank_colors <- rank_colors_small
plot_performance_metric_summary(perf_data, v1="hcr", v2="om")+
    theme(
        panel.spacing.x = unit(0.75, "cm"),
    )+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous(limits=c(0, 55)),
            scale_x_continuous(limits=c(0, 0.08), breaks=c(0, 0.025, 0.05, 0.075)),
            scale_x_continuous(limits=c(0, 550), breaks=c(0, 150, 300, 450)),
            scale_x_continuous(limits=c(0, 15)),
            scale_x_continuous(limits=c(0, 75), breaks=seq(0, 75, 15)),
            scale_x_continuous(limits=c(0, 15), breaks=seq(0, 15, 5))
        )
    )


rec_data <- get_recruits(
    model_runs,
    extra_columns,
    c("F40"),
    om_filter
)



examp_rec <- rec_data %>% 
    as_tibble() %>%
    filter(time > 54, sim %in% sample(seed_list, 5)) %>%
    mutate(
        sim = factor(sim)
    )

mean_rec <- examp_rec %>% group_by(om) %>% summarise(r=median(rec))

summ_rec <- rec_data %>%
    as_tibble() %>%
    filter(time > 54) %>%
    group_by(time, om) %>%
    median_qi(rec, .width=interval_widths)

ggplot(examp_rec) +
    geom_lineribbon(data=summ_rec, aes(x=time, y=rec, ymin=.lower, ymax=.upper), color="black", size=1, alpha=0.75)+
    # geom_line(aes(x=time, y=rec, color=om, group=sim), size=0.5, alpha=0.6)+
    geom_hline(data=mean_rec, aes(yintercept=r), linetype="dashed")+
    geom_text(data=mean_rec, aes(x=120, y=90, label=round(r, 3)), size=6)+
    scale_y_continuous(limits=c(0, 100))+
    scale_x_continuous(breaks=c(seq(55, 130, 20), 134), labels=c(seq(0, 75, 20), 75))+
    scale_fill_brewer(palette = "Blues")+
    guides(color="none", fill="none")+
    labs(x="Time", y="Recruitment (millions)")+
    coord_cartesian(expand=0, ylim=c(0, 100))+
    facet_wrap(~om, ncol=3)+
    theme_bw()+
    custom_theme

ggsave(file.path("~/Desktop/", "recruitment.jpeg"))
