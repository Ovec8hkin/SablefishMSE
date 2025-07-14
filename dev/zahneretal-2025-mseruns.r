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


om_list <- listN(om_rand_recruit, om_bh_recruit, om_bhcyclic_recruit, om_immcrash_recruit)
hcr_list <- listN(
    mp_pfmc4010
)


tic()
model_runs <- run_mse_multiple(
    om_list, 
    hcr_list, 
    seed_list,
    nyears=100,
    mse_options_list=mse_options_list,
    diagnostics = FALSE,
    save=TRUE
)
toc()


 # Data Processing
filetype <- ".jpeg"
figures_dir <- file.path(here::here(), "figures")
width_small <- 12
height_small <- 8

publication_hcrs <- c("F40", "F50", "PFMC 40-10", "British Columbia", "Constant F50", "F40 +/- 5%", "F40 +/- 10%", "15k Harvest Cap", "25k Harvest Cap", "No Fishing")
publication_oms2 <- c("Random Recruitment", "Beverton-Holt Cyclic Recruitment", "Immediate Crash Recruitment")
publication_oms <- c("Random Recruitment", "Beverton-Holt Regime Recruitment", "Crash Recruitment")
publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion of Years with Low SSB")
metric_names <- c("Annual Catch", "Catch AAV", "Average SSB", "Average Age", "Proportion of Years SSB < B35*")

mse_runs <- get_saved_model_runs(om_order=publication_oms2, hcr_order=publication_hcrs)
mse_runs <- readRDS(file.path(here::here(), "data", "zahneretal2025_mseruns_FINAL.RDS"))
model_runs <- mse_runs$model_runs
extra_columns <- mse_runs$extra_columns2
extra_columns <- extra_columns %>% 
    mutate(
        om = factor(om, levels = publication_oms2, labels = publication_oms),
        hcr = factor(hcr, levels = publication_hcrs, labels = publication_hcrs)
    )

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 130)

hcr_colors <- set_hcr_colors2(publication_hcrs)


### Spawning Biomass and Catch Plots
ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=publication_hcrs, om_filter=publication_oms)
plot_ssb(ssb_data, v1="hcr", v2="om", v3=NA, common_trajectory=common_trajectory, show_est = FALSE)
ggsave(filename=file.path(figures_dir, paste0("ssb", filetype)), width=width_small, height=height_small, units=c("in"))
 
# plot_relative_ssb(ssb_data, v1="hcr", v2="om", common_trajectory = common_trajectory, base_hcr = "No Fishing")
# ggsave(filename=file.path(figures_dir, paste0("rel_ssb", filetype)), width=width_small, height=height_small, units=c("in"))

catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms)
plot_landed_catch(catch_data, v1="hcr", v2="om", common_trajectory = common_trajectory)
ggsave(filename=file.path(figures_dir, paste0("catch", filetype)), width=width_small, height=height_small, units=c("in"))


# hist_OFLs <- c(16128, 13397, 15931, 30211, 51726, 61319, 40839, 47857, 55385, 58731)

rp_data <- expand.grid(hcr=publication_hcrs, om=publication_oms, L1=c("caa", "naa")) 
rp_data <- rp_data %>%
    mutate(
        hcr = factor(hcr, levels=publication_hcrs, labels=publication_hcrs),
        om = factor(om, levels=publication_oms, labels=publication_oms),
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
    geom_hline(data=rp_data, aes(yintercept=RP), linetype="dashed", color="grey50")

ggsave(filename=file.path(here::here(), "figures", "ssb_catch.jpeg"), width=16, height=10, units="in")


### Performance Metrics
performance_metrics <- performance_metric_summary(
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
    summary_out = TRUE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "prop_years_lowssb")
)
perf_data <- performance_metrics$perf_data %>% 
    mutate(
        name = factor(name, levels=publication_metrics, labels=metric_names)
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

rank_colors <- rank_colors_small
plot_performance_metric_summary(perf_data, v1="hcr", v2="om")+
    theme(
        panel.spacing.x = unit(0.75, "cm"),
    )
ggsave(filename=file.path(here::here(), "figures", paste0("performance_revised2", filetype)), width=18, height=12, units="in")

### Resilience Metrics
resilience_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=c("Crash Recruitment"),
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("crash_time", "recovery_time")
)

resilience_data <- resilience_metrics$perf_data %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

plot_performance_metric_summary(resilience_data)+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("Year", limits=c(0, 15)),
            scale_x_continuous("Year", limits=c(0, 35))
        )
    )+
    theme(
        strip.text.y = element_blank(),
        plot.margin = margin(10, 30, 10, 10),
        panel.spacing.x = unit(30, "pt")
    )
ggsave(filename=file.path(here::here(), "figures", paste0("resilience", filetype)), width=width_small, height=height_small, units="in")


### Aggregate Performance
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
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "prop_years_lowssb") 
)

om_aggregated_performance <- full_performance_metrics$perf_data %>%
    group_by(sim, hcr) %>%
    summarise(
        across(2:6, \(x) mean(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(3:7, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_years"), 
            labels=c("Average Annual Catch", "Catch AAV", "Average SSB", "Average Age", "Proportion of Years SSB < B35")
        )
    )


hist_abcs <- c(44200, 37100, 33400, 28800, 25200, 25000, 28800, 25300, 19600, 17200, 16800, 15900, 17200, 16900, 17300, 20900, 23000, 21000, 21000, 20100, 18000, 16100, 15200, 16000, 17200, 16200, 13700, 13700, 11800, 13100, 15000, 15100, 22000, 29600, 34500, 40500)

om_means <- om_aggregated_performance %>% group_by(name) %>% summarise(v = mean(value))
om_means$v <- rep(
    c(
            18538/1000,#median(assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"]),
            perf_data %>% ungroup() %>% filter(!is.infinite(median), hcr != "No Fishing", name=="Catch AAV") %>% summarise(median=mean(median)) %>% as.numeric,
            # median(assessment$t.series[,"spbiom"]),
            105.935,
            median(apply(assessment$natage.female, 1, \(x) compute_average_age(x, 2:31))),
            sum(assessment$t.series[,"spbiom"] < 105)/length(assessment$t.series[,"spbiom"])
        ),
    1
)

ggplot(om_aggregated_performance)+
    ggridges::stat_binline(aes(x=value, y=hcr, fill=hcr))+
    geom_vline(data=om_means, aes(xintercept=v), linetype="dashed", color="grey50")+
    scale_fill_manual(values=hcr_colors)+
    facet_wrap(~name, scales="free_x")+
    labs(y="")+
    custom_theme+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("", limits=c(0, 50)),
            scale_x_continuous("", limits=c(0, 0.06)),
            scale_x_continuous("", limits=c(0, 300)),
            scale_x_continuous("", limits=c(5, 10)),
            scale_x_continuous("", limits=c(0, 1), breaks=c(0, 0.25, 0.50, 0.75, 1.0), labels=seq(0, 1, 0.25), expand=c(0, NA))
        )
    )+
    guides(fill=guide_legend("Management \n Strategy", nrow=2))+
    theme(
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.x = unit(1, "cm")
    )
ggsave(file.path(figures_dir, paste0("performance_metrics_agg_hist", filetype)), width=12, height=9, units="in")

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(size=3) +
    geom_smooth(method = method, color = "black", fullrange=TRUE, ...)
  p
}

GGally::ggpairs(
    # om_aggregated_performance %>% 
    #     pivot_wider(names_from="name", values_from="value") %>%
    #     filter(hcr != "No Fishing"), 
    perf_data %>% filter(hcr != "No Fishing", .width==0.50) %>%
        select(om, hcr, name, median) %>%
        group_by(hcr, name) %>%
        summarise(value=mean(median)) %>%
        pivot_wider(names_from="name", values_from="value") %>%
        ungroup() %>%
        mutate(
            color = c("black", "#EA8115", "#83c738", "#FA5CCC", "#f3ac0c", "#30AF6C", "#8115EA", "#29C1D6", "#AC7A53")
        ),
    mapping=aes(color=hcr), 
    columns=2:6,
    lower = list(continuous = GGally::wrap(lowerFn, method = "lm", linewidth=0.5)),
    upper="blank",
    diag = "blank",
    legend=1
)+
scale_color_manual(values=hcr_colors)+
custom_theme

ggsave(filename=file.path(here::here(), "figures", paste0("performance_frontier_1", filetype)), width=12, height=12, units="in")

### Extra Analyses

# Get average performance metric value by HCR across the OM scenarios
om_aggregated_performance %>% group_by(hcr, name) %>%
    summarise(
        value = mean(value)
    ) %>%
    pivot_wider(names_from="name", values_from="value")

# Get number of simulations in which the population went extinct prematurely
ssb_data %>% filter(biomass < 1e-2) %>% pull(sim) %>% unique

# Get number of years in which stability and max harvest constraints were
# applied to their respective HCRs
catch_data %>%
    filter_times(time_horizon) %>%
    filter(fleet == "Fixed") %>%
    filter_hcr_om(hcrs = c("F40 +/- 5%", "F40 +/- 10%", "15k Harvest Cap", "25k Harvest Cap"), oms=publication_oms) %>%
    group_by(sim, om, hcr) %>%
    mutate(
        rel_perc_change = abs(total_catch/lag(total_catch)-1),
        constrained = case_when(
            hcr == "F40 +/- 5%" ~ ifelse(rel_perc_change >= 0.0499, 1, 0),
            hcr == "F40 +/- 10%" ~ ifelse(rel_perc_change >= 0.0999, 1, 0),
            hcr == "15k Harvest Cap" ~ ifelse(total_catch >= 14.99, 1, 0),
            hcr == "25k Harvest Cap" ~ ifelse(total_catch >= 24.99, 1, 0)
        )
    ) %>%
    # arrange(sim, om, time) %>%
    group_by(sim, om, hcr) %>%
    summarise(
        years_constrained = sum(constrained, na.rm=TRUE)/n()
    ) %>%
    group_by(om, hcr) %>%
    median_qi(years_constrained, .width = interval_widths) %>%
    print(n=100)



full_performance_metrics$perf_data %>%
    # filter(hcr != "15k Harvest Cap") %>%
    group_by(sim, om, hcr) %>%
    summarise(
        across(1:5, \(x) median(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(4:8, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_years"), 
            labels=publication_metrics
        )
    ) %>%
    group_by(om, hcr, name) %>%
    summarise(
        value=mean(value)
    ) %>%
    pivot_wider(names_from="hcr", values_from="value") %>%
    mutate(
        across(2:11, \(x) x/`F40`)
    ) %>%
    mutate(
        across(2:11, \(x) case_when(
            abs(x-1) <= 0.05 ~ 1,
            abs(x-1) <= 0.10 ~ 2,
            abs(x-1) <= 0.20 ~ 3,
            TRUE ~ 4
        ))
    ) %>%
    filter(name == "SSB")


recruit_data <- get_recruits(model_runs, extra_columns, hcr_filter=c("F40"), om_filter=publication_oms)

examp_rec <- recruit_data %>% 
    as_tibble() %>%
    filter(time > 54, sim %in% sample(seed_list, 5)) %>%
    mutate(
        sim = factor(sim)
    )

mean_rec <- examp_rec %>% group_by(om) %>% summarise(r=median(rec))

summ_rec <- recruit_data %>%
    as_tibble() %>%
    filter(time > 54) %>%
    group_by(time, om) %>%
    median_qi(rec, .width=interval_widths)

ggplot(examp_rec) +
    geom_line(aes(x=time, y=rec, color=om, group=sim), size=0.5, alpha=0.6)+
    geom_line(data=summ_rec, aes(x=time, y=rec), color="black", size=1)+
    geom_hline(data=mean_rec, aes(yintercept=r), linetype="dashed")+
    geom_text(data=mean_rec, aes(x=123, y=115, label=round(r, 2)), size=6)+
    scale_y_continuous(limits=c(0, 120))+
    scale_x_continuous(breaks=c(seq(55, 130, 20), 134), labels=c(seq(0, 75, 20), 75))+
    guides(color="none")+
    labs(x="Time", y="Recruitment (millions)")+
    coord_cartesian(expand=0)+
    facet_wrap(~om, ncol=3)+
    theme_bw()+
    custom_theme

ggsave(filename=file.path(here::here(), "figures", "example_recruitment.png"), width=width_small, height=6, units = "in")


dp_plot <- afscOM::subset_dem_params(sable_om$dem_params, 64:100, d=1, drop=FALSE)
dimnames(dp_plot$mort)$sex <- c("Female", "Male")
dimnames(dp_plot$waa)$sex <- c("Female", "Male")
dimnames(dp_plot$mat)$sex <- c("Female", "Male")
dimnames(dp_plot$sel)$sex <- c("Female", "Male")

plots <- afscOM::plot_demographic_parameters(
    dp_plot, 
    show_plots = TRUE
)

library(patchwork)

(
    plots$mort_plot+guides(color="none", linetype="none")+guides(linetype="none") +
    plots$waa_plot+guides(color="none", linetype="none")+labs(y="Weight (kg)")+guides(linetype="none")
) / 
(
    plots$mat_plot+guides(color="none", linetype="none")+guides(linetype="none") + 
    plots$sel_plot+ggtitle("Fishery Selectivity")+guides(linetype="none")
) + plot_layout(guides="collect") & theme(legend.position = "bottom")+custom_theme

ggsave(filename=file.path(figures_dir, paste0("demographics", filetype)), width = width_small)

lapply(
    afscOM::subset_dem_params(sable_om$dem_params, 64:100, d=1, drop=FALSE), 
    \(x) 
    dimnames(x)$sex <- c("Female", "Male")
)



mean_rec <- recruit_data %>% 
    filter(L1 == "naa") %>%
    group_by(om) %>% 
    summarise(r=mean(rec))


calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean_rec[3,2]/2,
    spr_target = 0.50
)

b40s <- get_b40_timeseries(
    model_runs, 
    extra_columns,
    c("F40"),
    publication_oms
)

f40_ssbdata<- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=c("F40"), om_filter=publication_oms)
f40_ssbdata <- f40_ssbdata %>%
        select(-c("biomass")) %>%
        # Compute quantiles of SSB distribution
        group_by(time, om, hcr) %>%
        summarise(spbio=median(spbio))

ggplot(b40s%>% filter(time >= 55))+
    geom_lineribbon(aes(x=time, y=B40, ymin=.lower, ymax=.upper, color=om), size=0.85)+
    geom_line(data=b40s %>% filter(time < 56), aes(x=time, y=B40), color="black", size=0.85)+
    geom_line(data=f40_ssbdata, aes(x=time, y=spbio), color="black", alpha=0.5, linetype="longdash", size=0.85)+
    scale_fill_brewer(palette="Blues")+
    scale_y_continuous(limits=c(0, 320))+
    coord_cartesian(expand=0)+
    labs(x="Year")+
    facet_wrap(~om)+
    guides(color="none", fill="none")+
    custom_theme
ggsave(filename=file.path(figures_dir, paste0("b40_timeseries", filetype)), width=width_small, height=6, units="in")



performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=publication_oms,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative="F40", 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "avg_variation", "avg_ssb", "avg_age", "prop_years_lowssb")
)

perf_data <- performance_metrics$perf_data %>% 
    mutate(
        hcr = factor(hcr, levels=publication_hcrs),
        name = factor(name, levels=publication_metrics, labels=metric_names)
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

rank_colors <- rank_colors_small
plot_performance_metric_summary(perf_data, v1="hcr", v2="om", is_relative = TRUE)+
    theme(
        panel.spacing.x = unit(0.75, "cm"),
    )
ggsave(file.path(figures_dir, paste0("performance_relative", filetype)), width=18, height=12)
