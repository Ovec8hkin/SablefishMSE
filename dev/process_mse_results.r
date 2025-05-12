library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggdist)
library(broom)
library(afscOM)

lapply(list.files("R", full.names = TRUE), source)

filetype <- ".jpeg"
figures_dir <- file.path("/Users/jzahner/Desktop/Presentations/may2025_lab_meeting", "figures")

nyears <- 110
sable_om <- readRDS("data/sablefish_om_big.RDS")

# source(file.path(here::here(), "dev", "oms.R"))
# source(file.path(here::here(), "dev", "hcrs.R"))
# om_list <- afscOM::listN(om_rand_recruit, om_bh_recruit, om_bhcyclic_recruit, om_immcrash_recruit)
# hcr_list <- afscOM::listN(
#     mp_f40, mp_f50, mp_b30f40, mp_b40f50,
#     mp_5perc, mp_10perc, mp_10perc_up, mp_15perc,
#     mp_15cap, mp_25cap,
#     mp_f50chr,
#     mp_pfmc4010, mp_bcsable,
#     mp_f00chr
# )

# om_names <- unlist(lapply(om_list, \(x) x$name))
# hcr_names <- unlist(lapply(hcr_list, \(x) x$name))

publication_hcrs <- c("F40", "F50", "F40 +/- 5%", "F40 +/- 10%", "15k Harvest Cap", "25k Harvest Cap", "Constant F50", "PFMC 40-10", "British Columbia", "No Fishing")
# publication_hcrs <- c("F40", "F40 +/- 5%", "20k Harvest Cap", "Complex")
# stakeholder_hcrs <- c("F40", "F50", "F40 +/- 5%", "20k Harvest Cap", "F40 Fully Mature", "No Fishing")
publication_oms <- c("Random Recruitment", "Beverton-Holt Cyclic Recruitment", "Immediate Crash Recruitment")
stakeholder_oms <- c("Historical", "Cyclic", "Crash")
publication_metrics = c("Annual Catch", "Catch AAV", "SSB", "Average Age", "Proportion of Years SSB < B35")
stakeholder_metrics = c("Annual Catch", "Catch AAV", "SSB", "Proportion Large Catch", "Average Age", "Proportion of Years SSB < B35", "Annual Value", "Dynamic Annual Value")

mse_runs <- get_saved_model_runs(om_order=publication_oms, hcr_order=publication_hcrs)
model_runs <- mse_runs$model_runs
extra_columns <- mse_runs$extra_columns2 #%>%
    # mutate(
    #     om = factor(om, levels = publication_oms, labels=stakeholder_oms)
    # )

interval_widths <- c(0.50, 0.80)
common_trajectory <- 54
time_horizon <- c(55, 130)

width_small <- 12
height_small <- 8

# hcr_colors <- set_hcr_colors(publication_hcrs)
hcr_colors <- c(
    "F40" = "#E31C39",
    "F50" = "#EA8115",
    "20k Harvest Cap" = "#1C39E3",
    "No Fishing" = "black",
    "F40 +/- 5%" = "#30AF6C",
    "F40 +/- 10%" = "#8115EA",
    "15k Harvest Cap" = "#29C1D6",
    "25k Harvest Cap" = "#AC7A53",
    "Constant F50" = "#f3ac0c",
    "PFMC 40-10" = "#83c738",
    "British Columbia" = "#FA5CCC"
)

### Spawning Biomass and Catch Plots
ssb_data <- get_ssb_biomass(model_runs, extra_columns, sable_om$dem_params, hcr_filter=publication_hcrs, om_filter=publication_oms)
plot_ssb(ssb_data, v1="hcr", v2="om", v3=NA, common_trajectory=common_trajectory, show_est = FALSE, highlight = NULL)
ggsave(filename=file.path(figures_dir, paste0("ssb", filetype)), width=width_small, height=6, units=c("in"))
 
plot_relative_ssb(ssb_data, v1="hcr", v2="om", common_trajectory = common_trajectory, base_hcr = "No Fishing")
ggsave(filename=file.path(figures_dir, paste0("rel_ssb", filetype)), width=width_small, height=height_small, units=c("in"))

catch_data <- get_landed_catch(model_runs, extra_columns, hcr_filter=publication_hcrs, om_filter=publication_oms)
plot_landed_catch(catch_data, v1="hcr", v2="om", common_trajectory = common_trajectory, highlight = c("F40", "F50", "F40 +/- 5%"))
ggsave(filename=file.path(figures_dir, paste0("catch", filetype)), width=width_small, height=6, units=c("in"))

plot_ssb_catch(
    ssb_data %>% filter(time > 0),
    catch_data %>% filter(time > 0),
    v1="hcr",
    v2="om",
    common_trajectory = common_trajectory,
    highlight = c("F40", "15k Harvest Cap", "25k Harvest Cap")
)#+scale_x_continuous(labels = seq(1960, 2090, 20))
ggsave(filename=file.path(figures_dir, "ssb_catch_small.jpeg"), width=16, height=9, units="in")

### Performance Metrics
fishery_performance_metrics <- performance_metric_summary(
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
    metric_list = c("avg_catch", "avg_variation", "avg_catch_lg")#,  "avg_ssb", "avg_age", "prop_years_lowssb", "annual_value", "dynamic_value")
)
# write_csv(performance_metrics$perf_data, file=file.path(here::here(), "data", "perfs.csv"))
rank_colors <- rank_colors_small#[c(1, length(rank_colors_small))]

perf_data <- fishery_performance_metrics$perf_data %>% 
    filter(!is.na(om), !is.na(hcr)) %>%
    mutate(
        name = factor(name, levels=stakeholder_metrics)
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median") %>%
    filter(!is.na(name))

plot_performance_metric_summary(perf_data, v1="hcr", v2="om")+ 
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous(limits=c(0, 55)),
            scale_x_continuous(limits=c(0, 0.10), breaks=c(0, 0.025, 0.05, 0.075, 0.10), labels=scales::percent_format(suffix="")),
            scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.50, 1.0))
            # scale_x_continuous(limits=c(0, 550), breaks=c(0, 150, 300, 450)),
            # scale_x_continuous(limits=c(0, 1)),
            # scale_x_continuous(limits=c(0, 15), breaks=seq(0, 15, 5)),
            # scale_x_continuous(limits=c(0, 30)),
            # scale_x_continuous(limits=c(0, 12))
        )
    )+
    theme(
        panel.spacing.x = unit(40, "pt")
    )
ggsave(filename=file.path(figures_dir, paste0("fishery_performance", filetype)), width=10, height=8, units="in")


conservation_performance_metrics <- performance_metric_summary(
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
    metric_list = c("avg_ssb", "avg_age", "prop_years_lowssb")#, "annual_value", "dynamic_value")
)
perf_data <- conservation_performance_metrics$perf_data %>% 
    filter(!is.na(om), !is.na(hcr), hcr %in% publication_hcrs) %>%
    mutate(
        name = factor(name, levels=c("SSB", "Average Age", "Proportion of Years with Low SSB"), labels=c("SSB", "Average Age", "Proportion of Years\nSSB < B35"))
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median") 

plot_performance_metric_summary(perf_data, v1="hcr", v2="om", highlight = c("F40", "15k Harvest Cap", "25k Harvest Cap"))+ 
    ggh4x::facetted_pos_scales(
        x = list(
            # scale_x_continuous(limits=c(0, 55)),
            # scale_x_continuous(limits=c(0, 0.10), breaks=c(0, 0.025, 0.05, 0.075, 0.10), labels=scales::percent_format(suffix="")),
            # scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.50, 1.0)),
            scale_x_continuous(limits=c(0, 550), breaks=c(0, 150, 300, 450)),
            scale_x_continuous(limits=c(0, 15)),
            scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25))
            # scale_x_continuous(limits=c(0, 30)),
            # scale_x_continuous(limits=c(5, 12))
        )
    )+
    theme(
        panel.spacing.x = unit(40, "pt")
    )
ggsave(filename=file.path(figures_dir, paste0("conservation_performance", filetype)), width=10, height=8, units="in")


economic_performance_metrics <- performance_metric_summary(
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
    metric_list = c("annual_value", "dynamic_value")
)
perf_data <- economic_performance_metrics$perf_data %>% 
    filter(!is.na(om), !is.na(hcr), hcr %in% publication_hcrs) %>%
    mutate(
        name = factor(name, levels=c("Annual Value", "Dynamic Annual Value"), labels=c("Annual Value", "Dynamic Annual Value"))
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median") 

plot_performance_metric_summary(perf_data, v1="hcr", v2="om")+ 
    ggh4x::facetted_pos_scales(
        x = list(
            # scale_x_continuous(limits=c(0, 55)),
            # scale_x_continuous(limits=c(0, 0.10), breaks=c(0, 0.025, 0.05, 0.075, 0.10), labels=scales::percent_format(suffix="")),
            # scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.50, 1.0)),
            # scale_x_continuous(limits=c(0, 550), breaks=c(0, 150, 300, 450)),
            # scale_x_continuous(limits=c(0, 15)),
            # scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.25))
            scale_x_continuous(limits=c(0, 30)),
            scale_x_continuous(limits=c(5, 12))
        )
    )+
    theme(
        panel.spacing.x = unit(40, "pt")
    )
ggsave(filename=file.path(figures_dir, paste0("economic_performance", filetype)), width=10, height=8, units="in")


perf_data %>% filter(.width == 0.50) %>% 
    select(om, hcr, name, scaled) %>%
    pivot_wider(names_from="name", values_from="scaled") %>%
    ungroup() %>%
    mutate(across(4:8, ~ as.numeric(.))) %>%
    rowwise() %>%
    mutate(
        avg_rank = mean(c_across(4:8), na.rm=TRUE)
    ) %>%
    arrange(om, avg_rank)


### Resilience Metrics
resilience_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=publication_hcrs,
    om_filter=c("Immediate Crash Recruitment"),
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("crash_time", "recovery_time")
)

resilience_data <- resilience_metrics$perf_data %>%
    add_row(
        om = "Immediate Crash Recruitment",
        hcr = "No Fishing",
        .width = c(0.5, 0.8),
        .point = "median",
        .interval = "qi",
        name = "Crash Time",
        median = 20, lower = 20, upper = 20
    ) %>%
    mutate(
        om = factor(om),
        hcr = factor(hcr, levels=publication_hcrs),
        name = factor(name, levels=c("Crash Time", "Recovery Time"))
    ) %>%
    group_by(.width, om, name) %>%
    scale_and_rank("median")

plot_performance_metric_summary(resilience_data)+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("Year", limits=c(0, 21)),
            scale_x_continuous("Year", limits=c(0, 40))
        )
    )+
    theme(
        strip.text.y = element_blank()
    )
ggsave(filename=file.path(figures_dir, paste0("resilience", filetype)), width=width_small, height=height_small, units="in")

### Aggregate Performance Distributions

perf_tradeoffs <- performance_metric_summary(
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

tradeoff_summ <- perf_tradeoffs$perf_data %>%
    # filter(hcr != "15k Harvest Cap") %>%
    group_by(sim, hcr) %>%
    summarise(
        across(2:6, \(x) median(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(3:7, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_years"), 
            labels=publication_metrics
        )
    )

# omavg_performance <- 

perf_tradeoffs$perf_data %>% 
    group_by(sim, hcr) %>%
    summarise(
        across(2:6, \(x) median(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(3:7, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "spbio", "avg_age", "prop_years"), 
            labels=publication_metrics
        )
    ) %>%
    group_by(hcr, name) %>%
    summarise(value=mean(value))

ggplot(tradeoff_summ)+
    ggridges::stat_binline(aes(x=value, y=hcr, fill=hcr))+
    geom_vline(data=omavg_performance, aes(xintercept=value), linetype="dashed")+
    scale_fill_manual(values=hcr_colors)+
    facet_wrap(~name, scales="free_x")+
    labs(y="")+
    custom_theme+
    ggh4x::facetted_pos_scales(
        x = list(
            scale_x_continuous("", limits=c(0, 35)),
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

ggsave(file.path(figures_dir, paste0("performance_metrics_agg_hist", filetype)))

tradeoff_summ %>% filter(hcr == "No Fishing", name == "SSB") %>% pull(value) %>% median


perf_tradeoffs$perf_data %>%
    # filter(hcr != "15k Harvest Cap") %>%
    group_by(sim, om, hcr) %>%
    summarise(
        across(1:5, \(x) median(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(4:8, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "biomass", "avg_age", "prop_years"), 
            labels=publication_metrics
        )
    ) %>%
    group_by(om,hcr, name) %>%
    summarise(
        value=mean(value)
    ) %>%
    pivot_wider(names_from="hcr", values_from="value") %>%
    mutate(
        across(2:12, \(x) x/`F40`)
    ) %>% arrange(name, om)


perf_tradeoffs$perf_data %>%
    # filter(hcr != "15k Harvest Cap") %>%
    group_by(sim, om, hcr) %>%
    summarise(
        across(1:5, \(x) median(x, na.rm=TRUE))
    ) %>% 
    pivot_longer(4:8, names_to="name", values_to="value") %>%
    mutate(
        name = factor(
            name, 
            levels=c("annual_catch", "aav", "biomass", "avg_age", "prop_years"), 
            labels=publication_metrics
        )
    ) %>%
    group_by(om, hcr, name) %>%
    summarise(
        value=mean(value)
    ) %>%
    pivot_wider(names_from="hcr", values_from="value") %>%
    mutate(
        across(2:12, \(x) x/`F40`)
    ) %>%
    mutate(
        across(2:12, \(x) case_when(
            abs(x-1) <= 0.05 ~ 1,
            abs(x-1) <= 0.10 ~ 2,
            abs(x-1) <= 0.20 ~ 3,
            TRUE ~ 4
        ))
    ) %>%
    pivot_longer(3:13, names_to="hcr", values_to="value") %>%
    count(om, name, value) %>% print(n=100)







plots <- afscOM::plot_demographic_parameters(afscOM::subset_dem_params(sable_om$dem_params, 64:100, d=1, drop=FALSE), show_plots = TRUE)

library(patchwork)

(
    plots$mort_plot+guides(color="none", linetype="none") + 
    plots$waa_plot+guides(color="none", linetype="none")+labs(y="Weight (kg)")
) / 
(
    plots$mat_plot+guides(color="none", linetype="none") + 
    plots$sel_plot
) + plot_layout(guides="collect")
ggsave(filename=file.path(figures_dir, paste0("demographics", filetype)))



performance_metrics <- performance_metric_summary(
    model_runs, 
    extra_columns, 
    sable_om$dem_params, 
    ref_naa,
    hcr_filter=stakeholder_hcrs,
    om_filter=stakeholder_oms,
    interval_widths=interval_widths,
    time_horizon = time_horizon, 
    extra_filter = NULL,
    relative=NULL, 
    summarise_by=c("om", "hcr"),
    summary_out = TRUE,
    metric_list = c("avg_catch", "avg_variation", "avg_catch_lg", "avg_ssb", "avg_age", "prop_years_lowssb", "annual_value", "dynamic_value")
)
fishery_performance_metrics$perf_data %>% 
    filter(hcr %in% stakeholder_hcrs, name %in% publication_metrics, om %in% publication_oms) %>%
    mutate(
            hcr = factor(hcr, levels=publication_hcrs),
            om = factor(om, labels=c("Historical", "Cyclic", "Crash")),
            name = factor(name, levels=publication_metrics)
    ) %>% 
    filter(hcr != "No Fishing") %>%
    group_by(om, name, .width) %>%
    mutate(
        scaled = median/max(median)
    ) %>%
    arrange(.width, desc(scaled), .by_group=TRUE) %>%
    mutate(
        rank = ifelse(
                    name %in% c("Catch AAV", "Proportion of Years SSB < B35"), 
                    factor(desc(row_number())), 
                    factor(row_number())
                ),
        # rank = factor(rank)
    ) %>% 
    ungroup() %>% 
    filter(.width == 0.50) %>% 
    select(om, hcr, name, scaled) %>%
    pivot_wider(names_from="name", values_from="scaled") %>%
    # mutate_at(3:ncol(.), as.integer) %>%
    select(-c(`Catch AAV`, `Proportion of Years with Low SSB`)) %>%
    mutate(
        total_rank = rowSums(across(3:ncol(.))),
        average_rank = rowMeans(across(3:ncol(.)))
    ) %>%
    group_by(om) %>%
    arrange(desc(average_rank), .by_group = TRUE) %>%
    print(n=100)



recruit_data <- get_recruits(model_runs, extra_columns, hcr_filter=c("F40"), om_filter=stakeholder_oms)

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
    facet_wrap(~om, nrow=2)+
    theme_bw()+
    custom_theme
ggsave(file.path("~/Desktop/Presentations/apr2025_sable_stakeholders", "figures", "recruit_scenarios.jpeg"), width=10, height=8, unit="in")
