#' Plot Trajectory of Spawning Biomass
#' 
#' Plot median trajectory of SSB across all simulations as a line plot.
#' 
#' @param data tibble output from `get_ssb_biomass()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param show_est whether to show estimated SSB from EM as point ranges
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param relative whether to relativize SSB to a specific HCR (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_ssb
#' 
plot_ssb <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64, base_hcr="F40", relative=NA, highlight=NULL){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "biomass")]
    
    # Plot spawning biomass from OM and EM
    d <- data %>%
        select(-c("biomass")) %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns)) %>%
        mutate(
            color_group = as.character(hcr) 
        )

    if(!is.null(highlight)){
        d <- d %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }

    plot <- plot_timeseries(d  %>% filter(L1 == "naa"), v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="SSB (1000s mt)", highlight=highlight)
    if(show_est){
        plot <- plot + geom_pointrange(data = d %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    return(plot)
}

plot_relative_ssb <- function(data, v1="hcr", v2=NA, common_trajectory=64, base_hcr="No Fishing"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "biomass")]
    
    base_ssb_data <- data %>% filter(hcr == base_hcr, L1 == "spbio")
    rel_ssb <- data %>% left_join(base_ssb_data, by=c("time", "sim", "L1", "om"), suffix=c("", ".nofish")) %>%
        filter(time > common_trajectory) %>%
        mutate(
            rel_ssb = spbio/spbio.nofish
        ) %>%
        filter(L1 == "naa") %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(rel_ssb, .width=interval_widths)

    ylim <- c(0, 1.1)

    plot <- ggplot(rel_ssb) +
        geom_line(aes(x=time, y=rel_ssb, color=.data[[v1]], group=.data[[v1]]), size=0.85)+
        scale_y_continuous(limits=ylim)+
        scale_color_manual(values=hcr_colors)+
        coord_cartesian(expand=0)+
        labs(x="Year", y="Relative SSB")+
        guides(color=guide_legend(title="Management \n Strategy", nrow=2))+
        facet_wrap(~.data[[v2]])

    return(plot+custom_theme)

}

#' Plot Trajectory of Fishing Mortality
#' 
#' Plot median trajectory of total fishing mortality across all simulations as a line plot.
#' 
#' @param data tibble output from `get_fishing_mortalities()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param show_est whether to show estimated SSB from EM as point ranges
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_fishing_mortalities
#' 
plot_fishing_mortalities <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64, interval_widths=c(0.50, 0.80), base_hcr="F40", highlight=NULL){
    # Plot fishing mortality rates from OM and EM
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "F", "total_F")]

    f <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(F, total_F, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns)) %>%
        filter(name == "total_F")

    if(!is.null(highlight)){
        f <- f %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }

    plot <- plot_timeseries(f  %>% filter(L1 == "faa"), v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Fishing Mortality", highlight=highlight)
    if(show_est){
        plot <- plot + geom_pointrange(data = f %>% filter(L1 == "faa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    return(plot)
}

#' Plot Trajectory of Recruitment
#' 
#' Plot median trajectory of recruitment across all simulations as a line plot.
#' 
#' @param data tibble output from `get_recruits()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param show_est whether to show estimated SSB from EM as point ranges
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_recruitment
#' 
plot_recruitment <- function(data, v1="hcr", v2=NA,v3=NA,  show_est=FALSE, common_trajectory=64, interval_widths=c(0.50, 0.80), base_hcr="F40", highlight=NULL){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "rec")]

    r <- data %>%
        # summarise SSB across year and sim 
        group_by(across(all_of(group_columns))) %>%
        mean_qi(rec, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))

    mean_rec <- r %>% pull(median) %>% mean

    if(!is.null(highlight)){
        r <- r %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }

    plot <- plot_timeseries(r %>% filter(L1 == "naa"), v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Recruits (millions)", highlight=highlight)+
        geom_hline(yintercept = mean_rec, linetype="dashed")
    if(show_est){
        plot <- plot + geom_pointrange(data = f %>% filter(L1 == "naa_est"), aes(x=time, y=median, ymin=lower, ymax=upper, color=hcr), alpha=0.35)
    }

    return(plot)
}

#' Plot Trajectory of Landed Catch
#' 
#' Plot median trajectory of total landed catch across all simulations as a line plot.
#' 
#' @param data tibble output from `get_landed_catch()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_landed_catch
#' 
plot_landed_catch <- function(data, v1="hcr", v2=NA, v3=NA, by_fleet=FALSE, common_trajectory=64, base_hcr="F40", highlight=NULL){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    c <- data %>%
        select(-c("catch")) %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns)) %>%
        mutate(
            color_group = as.character(hcr) 
        )

    if(!is.null(highlight)){
        c <- c %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }
    
    plot <- plot_timeseries(c, v1, v2, v3, common_trajectory, interval_widths, base_hcr, ylab="Landings (1000s mt)", highlight=highlight)

    return(plot)

}

#' Plot Trajectory of Dynamic Economic Value
#' 
#' Plot median trajectory of economic value for the fixed gear fleet across 
#' all simulations as a line plot. See `get_dynamic_economic_value()` for 
#' details on how economic value is calculated.
#' 
#' @param data tibble output from `get_dynamic_economic_value()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param relative whether to relativize economic value to a specific HCR (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_dynamic_economic_value
#' 
plot_dynamic_economic_value <- function(data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40", highlight=NULL, relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "total_value")]
    d <- data %>% distinct(.keep_all=TRUE)

    # Relativize to specific HCR
    hcrs <- d %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        d <- d %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column="total_value", 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    # Plot spawning biomass from OM and EM
    d <- d %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(total_value, .width=interval_widths, .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    if(!is.null(highlight)){
        d <- d %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }

    plot <- plot_timeseries(d, v1, v2, v3, common_trajectory, interval_widths, base_hcr, highlight, ylab="Economic Value")
    if(!is.na(relative)){
        plot <- plot+geom_hline(yintercept=1, linetype="dashed")
    }

    return(plot)
}

#' Plot Trajectory of Average Population Age
#' 
#' Plot median trajectory of average population age across all simulations as a line plot.
#' 
#' @param data tibble output from `get_average_age()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param relative whether to relativize average age to a specific HCR (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_average_age
#' 
plot_average_age <- function(data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, time_horizon=NULL, interval_widths=c(0.50, 0.80), base_hcr="F40", highlight=NULL, relative=NA){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "avg_age")]
    d <- data %>% distinct(.keep_all=TRUE)

    if(!is.null(time_horizon)){
        d <- d %>% filter_times(time_horizon)
    }

    # Relativize to specific HCR
    hcrs <- d %>% distinct(hcr) %>% pull
    if(!is.na(relative) && relative %in% hcrs){
        d <- d %>% 
            relativize_performance(
                rel_column="hcr", 
                value_column="avg_age", 
                rel_value=relative, 
                grouping=c("sim", group_columns)
            )
    }

    # Plot spawning biomass from OM and EM
    d <- d %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(avg_age, .width=interval_widths, .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    if(!is.null(highlight)){
        d <- d %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }

    plot <- plot_timeseries(d, v1, v2, v3, common_trajectory, interval_widths, base_hcr, highlight, ylab="Average Age (Years)")
    if(!is.na(relative)){
        plot <- plot+geom_hline(yintercept=1, linetype="dashed")
    }

    return(plot)
}

#' Plot Trajectories of Spawning Biomass and Landed Catch
#' 
#' Plot median trajectory of total landed catch and SSB across all simulations as line plots.
#' 
#' @param ssb_data tibble output from `get_ssb_biomass()`
#' @param catch_data tibble output from `get_landed_catch()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_ssb_catch
#' 
plot_ssb_catch <- function(ssb_data, catch_data, v1="hcr", v2=NA, v3=NA, common_trajectory=64, base_hcr="F40", flip_facet=FALSE, highlight=NULL){

    group_columns <- colnames(ssb_data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "biomass")]
    # Plot spawning biomass from OM and EM
    ssb_d <- ssb_data %>%
        select(-c("biomass")) %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns)) %>%
        mutate(
            color_group = as.character(hcr) 
        )

    if(!is.null(highlight)){
        ssb_d <- ssb_d %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }

    hcr1 <- as.character((ssb_d %>% pull(hcr) %>% unique)[1])

    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- ssb_d %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    ssb_common <- ssb_d %>% left_join(traj, by=traj_column) %>% filter(L1=="naa", hcr==hcr1) %>% group_by(om) %>% filter(time <= common)


    group_columns <- colnames(catch_data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    catch_d <- catch_data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns)) %>%
        mutate(
            color_group = as.character(hcr) 
        )

    if(!is.null(highlight)){
        catch_d <- catch_d %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }
    
    hcr1 <- as.character((catch_d %>% pull(hcr) %>% unique)[1])
    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- catch_d %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    catch_common <- catch_d %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(om) %>% filter(time <= common)

    d <- bind_rows(ssb_d, catch_d) %>% filter(L1 != "naa_est") %>% 
            mutate(L1 = factor(L1, labels=c("Landed Catch", "SSB")))
    common <- bind_rows(ssb_common, catch_common) %>% filter(L1 != "naa_est") %>% 
            mutate(L1 = factor(L1, labels=c("Landed Catch", "SSB")))

    base_hcr <- d %>% filter(hcr == base_hcr)

    colors <- hcr_colors
    sizes <- rep(0.85, length(colors))
    names(sizes) <- names(colors)
    if(!is.null(highlight)){
        colors <- hcr_colors[highlight]
        colors <- c(colors, "Other" = "grey70")

        sizes <- c(rep(1.2, length(highlight)), 0.85)
        names(sizes) <- c(highlight, "Other")
    }

    plot <- ggplot(d) + 
        geom_line(data = base_hcr, aes(x=time, y=median, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, group=.data[[v1]], color=hcr), size=0.85)+
        # geom_line(
        #     data = d %>% filter(color_group != "Other", time > common_trajectory-1), 
        #     aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=color_group, size=color_group)
        # )+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        # geom_hline(yintercept=121.4611, linetype="dashed")+
        scale_fill_brewer(palette="Blues")+
        scale_color_manual(values=colors)+
        # scale_size_manual(values=sizes)+
        # scale_y_continuous(limits=c(0, 320))+
        labs(x="Year", y="1000s Metric Tons")+
        coord_cartesian(expand=0)+
        guides(color=guide_legend("Management \n Strategy", nrow=2), fill="none", size="none")+
        facet_grid(rows=vars(L1), cols=vars(.data[[v2]]), scales="free_y")+
        ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 500))
            )
        )

    if(flip_facet){
        plot <- plot + facet_grid(rows=vars(.data[[v2]]), cols=vars(L1), scales="free_y")
    }
    
    return(plot+custom_theme)
}

plot_atage_trajectory_ternary <- function(data, segments, col_names){
    axis_names = names(data)[6:8]
    return(
        ggplot(data, aes(x=.data[[col_names[1]]], y=.data[[col_names[2]]], z=.data[[col_names[3]]], color=hcr))+
            coord_tern(Tlim=c(0, 1), Llim=c(0, 1), Rlim=c(0, 1))+
            geom_point()+
            geom_segment(
                data = segments, 
                aes(x=x, y=y, z=z, xend=xend, yend=yend, zend=zend, group=hcr),
                arrow=arrow(length = unit(3, "mm"))
            )+
            scale_T_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_L_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_R_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            facet_grid(rows=vars(om), cols=vars(hcr))+
            labs(x=axis_names[1], y=axis_names[2], z=axis_names[3])+
            theme_bw()+
            theme(
                legend.position="bottom",
                tern.axis.arrow.show = TRUE,
                panel.spacing.x = unit(1, "cm"),
                panel.grid.minor = element_line(color="white")
            )
    
    )

}

plot_atage_density_ternary <- function(data, col_names){
    axis_names <- names(data)[6:8]
    return(
        ggplot(data, aes(x=.data[[col_names[1]]], y=.data[[col_names[2]]], z=.data[[col_names[3]]]))+
            coord_tern()+
            stat_density_tern(
                geom='polygon',
                aes(fill=..level..),
                bins=100,
            )+
            geom_mean_ellipse(color="white")+
            scale_fill_viridis(limits=c(0, 65))+
            scale_T_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_L_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            scale_R_continuous(breaks = seq(0, 1, 0.5), labels=seq(0, 100, 50))+
            facet_grid(rows=vars(om), cols=vars(hcr))+
            labs(fill="Simulation Years", x=axis_names[1], y=axis_names[2], z=axis_names[3])+
            theme_bw()+
            theme(
                legend.position="bottom",
                tern.axis.arrow.show = TRUE,
                panel.spacing.x = unit(1, "cm"),
                panel.grid.minor = element_line(color="white")
            )
    )
}

# plot_atage <- function(data, v1="hcr", v2=NA){
#     group_columns <- colnames(data)
#     group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]
    
#     d <- data %>%
#         group_by(time, class, hcr, om, L1) %>%
#         median_qi(value, .width=c(0.50))

#     plot <- ggplot()+
#         geom_bar(aes(x=time, y=value, fill=class), position="fill", stat="identity")+
#         geom_vline(xintercept = 2022-1960+1, color="white", size=1.25, linetype="dashed")+
#         scale_fill_viridis(direction=-1, discrete=TRUE, option="magma")+
#         scale_x_continuous(breaks=seq(1, nyears+1, 20), labels=seq(1960, 1960+nyears+1, 20))+
#         coord_cartesian(expand=0)+
#         labs(x="Year", fill="Age Group")+
#         guides(fill=guide_legend(reverse=TRUE))+
#         # facet_wrap(~L1, ncol=1, scales="free_x")+
#         facet_grid(rows=vars(hcr), cols=vars(L1))+
#         theme_bw()+
#         theme(
#             axis.text = element_text(size=12),
#             axis.title.y=element_blank(), 
#             strip.background = element_blank(),
#             strip.text.x = element_text(size=16, hjust=0),
#             panel.spacing.y = unit(0.4, "in"),
#             legend.position = "bottom"
#         )

#     if(!is.na(v2)){
#         plot <- plot + facet_wrap(~.data[[v2]])
#     }else{
    
#     }
    
#     # if(by_fleet){
#     #     plot <- plot + facet_wrap(~fleet)
#     # }

#     return(plot)

# }

#' Plot Trajectories of ABC, TAC, Expected Landings, and Attainment
#' 
#' Plot median trajectory of management quantitoes across all simulations as line plots.
#' 
#' @param data tibble output from `get_management_quantities()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_ssb_catch
#' 
plot_abc_tac <- function(data, v1="hcr", v2=NA, common_trajectory=64, interval_widths=c(0.50, 0.80), base_hcr="F40", highlight=NULL){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "value")]

    q <- data %>%
        mutate(
            L1 = factor(L1, levels=c("abc", "tac", "exp_land", "attainment"), labels=c("ABC", "TAC", "Expected Landings", "Attainment"))
        ) %>%
        # filter(L1 != "Expected Landings") %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(value, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    if(!is.null(highlight)){
        q <- q %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }

    plot <- plot_timeseries(q, v1, v2="L1", v3="om", common_trajectory, interval_widths, base_hcr, ylab="Quantity", highlight=highlight)

    plot <- plot + ggh4x::facetted_pos_scales(
            y = list(
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 60)),
                scale_y_continuous(limits=c(0, 50)),
                scale_y_continuous(limits=c(0.5, 1.5))
            )
        )

    return(plot)
}

plot_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=64){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "total_F")]

    d <- data %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, total_F, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
            filter(.width == 0.50)

    segments <- d %>% as_tibble() %>% 
        select(spbio, total_F, hcr, om) %>% 
        rename(x=spbio, y=total_F) %>%
        group_by(hcr, om) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    plot <- ggplot(d, aes(x=spbio, y=total_F, color=hcr, group=hcr))+
        geom_point(size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, group=hcr),
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_hline(data=ref_pts, aes(yintercept=Fref), linetype="dashed")+
        geom_vline(data=ref_pts, aes(xintercept=Bref), linetype="dashed")+
        scale_x_continuous(limits=c(0, 200))+
        scale_y_continuous(limits=c(0, 0.125))+
        coord_cartesian(expand=0)+
        facet_grid(cols=vars(hcr), rows=vars(om))

    return(plot)
}

plot_catch_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=64){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "total_catch")]

    d <- data %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, total_catch, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
            filter(.width == 0.50)

    segments <- d %>% as_tibble() %>% 
        select(spbio, total_catch, hcr, om) %>% 
        rename(x=spbio, y=total_catch) %>%
        group_by(hcr, om) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    plot <- ggplot(d, aes(x=spbio, y=total_catch, color=hcr, group=hcr))+
        geom_point(size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, group=hcr),
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_hline(data=ref_pts, aes(yintercept=Fref), linetype="dashed")+
        geom_vline(data=ref_pts, aes(xintercept=Bref), linetype="dashed")+
        scale_x_continuous(limits=c(0, 200))+
        scale_y_continuous(limits=c(0, 35))+
        coord_cartesian(expand=0)+
        facet_grid(cols=vars(hcr), rows=vars(om))

    return(plot)
}

plot_hcr_phase_diagram <- function(data, ref_pts, v1="hcr", v2=NA, common_trajectory=64){

    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio", "value")]

    d <- data %>%
            rename(out_F=value) %>%
            group_by(across(all_of(group_columns))) %>%
            filter(time > common_trajectory) %>%
            median_qi(spbio, out_F, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
            filter(.width == 0.50)

    segments <- d %>% as_tibble() %>% 
        select(spbio, out_F, hcr, om) %>% 
        rename(x=spbio, y=out_F) %>%
        group_by(hcr, om) %>%
        mutate(
            xend = lead(x, 1),
            yend = lead(y, 1)
        ) %>%
        ungroup() %>%
        arrange(hcr) %>%
        drop_na()

    plot <- ggplot(d, aes(x=spbio, y=out_F, color=hcr, group=hcr))+
        geom_point(size=1.5)+
        geom_segment(
            data = segments, 
            aes(x=x, y=y, xend=xend, yend=yend, group=hcr),
            arrow=arrow(length = unit(3, "mm"))
        )+
        geom_hline(data=ref_pts, aes(yintercept=Fref), linetype="dashed")+
        geom_vline(data=ref_pts, aes(xintercept=Bref), linetype="dashed")+
        scale_x_continuous(limits=c(0, 200))+
        scale_y_continuous(limits=c(0, 0.125))+
        coord_cartesian(expand=0)+
        facet_grid(cols=vars(hcr), rows=vars(om))

    return(plot)

}


#' Plot Summary of Performance Metrics
#' 
#' Plot summary of performance metrics across all simulation as point intervals.
#' 
#' @param perf_data tibble output from `performance_metric_summary()`
#' @param v1 variable to map to y-axis (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param is_relative whether the performance metrics have been relativized to a specific HCR
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' 
#' @export plot_performance_metric_summary
#' 
plot_performance_metric_summary <- function(perf_data, v1="hcr", v2="om", is_relative=FALSE, highlight=NULL){

    metric_minmax = perf_data %>% group_by(name) %>% summarise(min=min(lower), max=max(upper))
    axis_scalar <- c(0.9, 1.1)

    # hist_abcs <- c(44200, 37100, 33400, 28800, 25200, 25000, 28800, 25300, 19600, 17200, 16800, 15900, 17200, 16900, 17300, 20900, 23000, 21000, 21000, 20100, 18000, 16100, 15200, 16000, 17200, 16200, 13700, 13700, 11800, 13100, 15000, 15100, 22000, 29600, 34500, 40500)

    summary <- perf_data %>% filter(!is.infinite(median), hcr != "No Fishing") %>% summarise(median=mean(median))
    # summary$median <- rep(
    #     c(
    #         18538/1000,#median(assessment$t.series[,"Catch_HAL"]+assessment$t.series[,"Catch_TWL"]),
    #         perf_data %>% ungroup() %>% filter(!is.infinite(median), hcr != "No Fishing", name=="Catch AAV") %>% summarise(median=mean(median)) %>% as.numeric,
    #         # median(assessment$t.series[,"spbiom"]),
    #         105.935,
    #         median(apply(assessment$natage.female, 1, \(x) compute_average_age(x, 2:31))),
    #         sum(assessment$t.series[,"spbiom"] < 105)/length(assessment$t.series[,"spbiom"])
    #     ),
    #     length(perf_data$om %>% unique)*2
    # )
    if(is_relative){
        summary$median <- rep(1, nrow(summary))
    }

    perf_data <- perf_data %>% mutate(color_group = as.character(hcr))
    if(!is.null(highlight)){
        perf_data <- perf_data %>% mutate(
            color_group = ifelse(hcr %in% highlight, color_group, "Other")
        )
    }
    
    if(!is.null(highlight)){
        colors <- c(hcr_colors, "Other" = "grey70")
        colors <- colors[perf_data$color_group %>% unique]
    }

    color_var <- ifelse(is.null(highlight), "rank", "color_group")

    plot <- ggplot(perf_data)+
                geom_vline(data=summary, aes(xintercept = median), color="grey50", linetype="dashed")+
                scale_shape_discrete()+
                labs(y="", x="", shape="OM", color="Relative MS Order")+
                coord_cartesian(expand=0)+
                guides(shape="none", color=guide_legend(nrow=1))+
                theme(
                    plot.margin = margin(0.25, 1, 0.25, 0.25, "cm"),
                    panel.spacing.x = unit(5, "cm"),
                    plot.title = element_text(size=18),
                    legend.spacing.x = unit(1.5, "cm")
                )

    if(is.character(v2)){
        plot <- plot + 
                    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=.data[[v1]], color=.data[[color_var]], shape=.data[[v2]]), point_size=3, position="dodge")+
                    facet_grid(rows=vars(.data[[v2]]), cols=vars(name), scales="free_x")
    }else{
        plot <- plot + 
                    geom_pointinterval(aes(x=median, xmin=lower, xmax=upper, y=.data[[v1]], color=.data[[color_var]]), point_size=3, position="dodge")+
                    facet_wrap(~name, scales="free_x")
    }

    return(plot+custom_theme)
}

plot_ssb_paginate <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64, base_hcr="F40"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "spbio")]
    # Plot spawning biomass from OM and EM
    d <- data %>%
        # Compute quantiles of SSB distribution
        group_by(across(all_of(group_columns))) %>%
        median_qi(spbio, .width=c(0.50, 0.80), .simple_names=FALSE) %>%
        # Reformat ggdist tibble into long format
        reformat_ggdist_long(n=length(group_columns))

    hcr1 <- as.character((d %>% pull(hcr) %>% unique)[1])

    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- d %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- d %>% left_join(traj, by=traj_column) %>% filter(L1=="naa", hcr==hcr1) %>% group_by(om) %>% filter(time <= common, time >= 40) %>% select(-hcr)

    base_hcr_d <- d %>% filter(L1 == "naa", hcr == base_hcr) %>% select(-hcr)

    ps <- lapply(unlist(unique(c(om_names))), function(o){

        d2 <- d %>% filter(om == o)
        base_hcr_d2 <- base_hcr_d %>% filter(om == o)

        ggplot(d2 %>% filter(L1 == "naa")) + 
            geom_lineribbon(data = base_hcr_d2, aes(x=time, y=median, ymin=lower, ymax=upper, group=1), color="black", size=0.85)+
            geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
            geom_line(data = common, aes(x=time, y=median), size=0.85)+
            geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
            # geom_hline(yintercept=121.4611, linetype="dashed")+
            scale_fill_brewer(palette="Blues")+
            scale_y_continuous(limits=c(0, max(d2 %>% pull(median))*1.2))+
            scale_x_continuous(limits=c(40, max(base_hcr_d2 %>% pull(time))))+
            facet_wrap(~ hcr, ncol=4, nrow=7)+
            labs(x="Year", y="SSB", title=o)+
            coord_cartesian(expand=0)+
            guides(fill="none", color="none")
            custom_theme+
            theme(
                plot.title = element_text(size=20)
            )
    })

    ggsave(
        filename = "~/Desktop/ssb_paginated.pdf", 
        plot = marrangeGrob(ps, nrow=1, ncol=1), 
        width = 8.5, height = 11
    )

    return(ps)
}

plot_catch_paginate <- function(data, v1="hcr", v2=NA, v3=NA, show_est=FALSE, common_trajectory=64, base_hcr="F40"){
    group_columns <- colnames(data)
    group_columns <- group_columns[! group_columns %in% c("sim", "catch", "total_catch")]

    c2 <- data %>%
        group_by(across(all_of(group_columns))) %>%
        median_qi(catch, total_catch, .width=c(0.50, 0.80), .simple_names=TRUE) %>%
        reformat_ggdist_long(n=length(group_columns))
    
    hcr1 <- as.character((c2 %>% pull(hcr) %>% unique)[1])
    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- c2 %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    c2 <- c2 %>% filter(name == "total_catch")

    common <- c2 %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(om) %>% filter(time <= common, time >= 40) %>% select(-hcr)

    base_hcr_c <- c2 %>% filter(hcr == base_hcr) %>% select(-hcr)

    ps <- lapply(unlist(unique(c(om_names))), function(o){

        c3 <- c2 %>% filter(om == o)
        base_hcr_c2 <- base_hcr_c %>% filter(om == o)

        ggplot(c3 %>% filter(L1 == "land_caa")) + 
            geom_lineribbon(data = base_hcr_c2, aes(x=time, y=median, ymin=lower, ymax=upper, group=1), color="black", size=0.85)+
            geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
            geom_line(data = common, aes(x=time, y=median), size=0.85)+
            geom_vline(data=common, aes(xintercept=common_trajectory), linetype="dashed")+
            # geom_hline(yintercept=121.4611, linetype="dashed")+
            scale_fill_brewer(palette="Blues")+
            # scale_y_continuous(limits=c(0, max(c3 %>% pull(median))*1.2))+
            scale_x_continuous(limits=c(40, max(base_hcr_c2 %>% pull(time))))+
            facet_wrap(~ hcr, ncol=4, nrow=7)+
            labs(x="Year", y="Catch (mt)", title=o)+
            coord_cartesian(expand=0, ylim=c(0, 60))+
            guides(fill="none", color="none")
            custom_theme+
            theme(
                plot.title = element_text(size=20)
            )
    })

    ggsave(
        filename = "~/Desktop/catch_paginated.pdf", 
        plot = marrangeGrob(ps, nrow=1, ncol=1), 
        width = 8.5, height = 11
    )

    return(ps)
}

set_hcr_colors <- function(hcrs){
    hcr_colors <- scales::hue_pal()(length(hcrs))
    hcr_colors[which(hcrs == "No Fishing")] <- "#000000"
    names(hcr_colors) <- hcrs
    return(hcr_colors)
}

set_hcr_colors2 <- function(hcrs){
    hcr_colors <- c(
        "F40" = "black",
        "F50" = "#EA8115",
        "B40/F50" = "#1C39E3",
        "No Fishing" = "#E31C39",
        "F40 +/- 5%" = "#30AF6C",
        "F40 +/- 10%" = "#8115EA",
        "15k Harvest Cap" = "#29C1D6",
        "25k Harvest Cap" = "#AC7A53",
        "Constant F50" = "#f3ac0c",
        "PFMC 40-10" = "#83c738",
        "British Columbia" = "#FA5CCC"
    )
    return(hcr_colors)
}


rank_colors_small <- c(
    "#D55E00",
    "#FF740A",
    "#FF8B33",
    "#FFA35C",
    "#FFBA85",
    "#AAAAAA",
    "#5CC3FF",
    "#33B4FF",
    "#0AA5FF",
    "#008EE0",
    "#0072B2"
)

rank_colors_large <- c(
    "#8F3E00",
    "#B85000",
    "#D55E00",
    "#FF740A",
    "#FF8B33",
    "#FFA35C",
    "#FFBA85",
    "#AAAAAA",
    "#5CC3FF",
    "#33B4FF",
    "#0AA5FF",
    "#008EE0",
    "#0072B2",
    "#005A8F",
    "#004166"
)

#' Custom GGPlot2 theme
#'
#' @export custom_theme
custom_theme <- ggplot2::theme_bw()+ggplot2::theme(
    panel.spacing.y = ggplot2::unit(0.5, "cm"),
    panel.grid.minor = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size=14),
    axis.text = ggplot2::element_text(size=14),
    strip.text = ggplot2::element_text(size=14),
    legend.text = ggplot2::element_text(size=14),
    legend.position = "bottom"
)

#' Generialize Timeseries Plotting Function
#' 
#' General function for plotting timeseries of management quantities (e.g. SSB, catch, ABC, TAC) 
#' with options for faceting and highlighting specific HCRs.
#' 
#' @param data tibble summarised with `median_qi()` and reformatted with `reformat_ggdist_long()`
#' @param v1 variable to map to color (e.g. "hcr")
#' @param v2 variable to facet by (e.g. "om")
#' @param v3 variable to facet by in addition to v2 (e.g. "om")
#' @param common_trajectory year to plot vertical dashed line to indicate where MSE projection period begins
#' @param interval_widths vector of credible interval widths to plot (e.g. c(0.50, 0.80))
#' @param base_hcr HCR to plot as a thicker line for reference (must match names in `extra_columns`)
#' @param highlight vector of HCRs to highlight with color (must match names in `extra_columns`)
#' @param ylab label for y-axis
#' 
#' @export plot_timeseries
#' 
plot_timeseries <- function(data, v1="hcr", v2=NA, v3=NA, common_trajectory=54, interval_widths=c(0.50, 0.80), base_hcr="F40", ylab="", highlight=NULL){

    hcr1 <- as.character((data %>% pull(hcr) %>% unique)[1])

    traj_column <- ifelse(is.na(v3), v2, v3)
    traj <- data %>% distinct(eval(rlang::parse_expr(traj_column))) %>% mutate(common=common_trajectory) %>% rename(!!traj_column := 1)

    common <- data %>% left_join(traj, by=traj_column) %>% filter(hcr==hcr1) %>% group_by(eval(rlang::parse_expr(v2))) %>% filter(time <= common)

    base_hcr_d <- data %>% filter(hcr == base_hcr)

    colors <- hcr_colors
    sizes <- rep(0.85, length(colors))
    names(sizes) <- names(colors)
    if(!is.null(highlight)){
        colors <- hcr_colors[highlight]
        colors <- c(colors, "Other" = "grey70")

        sizes <- c(rep(1.2, length(highlight)), 0.85)
        names(sizes) <- c(highlight, "Other")
    }

    plot <- ggplot(data) + 
        geom_line(data = base_hcr_d, aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(aes(x=time, y=median, ymin=lower, ymax=upper, group=.data[[v1]], color=.data[[v1]]), size=0.85)+
        geom_line(data = common, aes(x=time, y=median), size=0.85)+
        geom_vline(data=common, aes(xintercept=common), linetype="dashed")+
        scale_fill_brewer(palette="Greys")+
        scale_color_manual(values=colors)+
        scale_size_manual(values=sizes)+
        labs(x="Year", y=ylab)+
        coord_cartesian(expand=0)+
        guides(color=guide_legend(title="Management \n Strategy", nrow=2), fill="none")

    if(!is.na(v2) && is.na(v3)){
        plot <- plot + facet_wrap(~.data[[v2]])+guides(fill="none")
        plot <- plot + scale_y_continuous(limits=c(0, max(data %>% pull(median))*1.2))
    }else if(!is.na(v2) && !is.na(v3)){
        plot <- plot + facet_grid(rows=vars(.data[[v2]]), cols=vars(.data[[v3]]), scales="free_y")+guides(fill="none")
    }
    
    return(plot+custom_theme)
}
