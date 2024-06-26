#' Get Spawning Biomass and Total Biomass
#' 
#' Process MSE simulations for spawning biomass,
#' and total stock biomass.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_ssb_biomass
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_ssb_biomass(model_runs, extra_columns)
#' }
#'
get_ssb_biomass <- function(model_runs, extra_columns, dem_params){
    group_columns <- c("time", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            # join WAA and maturity-at-age for computing SSB
            left_join(
                melt(sable_om$dem_params$waa, value.name="weight"), 
                by=c("time", "age", "sex")
            ) %>%
            left_join(
                melt(sable_om$dem_params$mat, value.name="maturity"), 
                by=c("time", "age", "sex")
            ) %>%
            drop_na() %>%
            # compute derived quantities
            mutate(
                biomass = value*weight,
                spbio = value*weight*maturity
            ) %>%
            # SSB is females only
            filter(sex == "F") %>%
            # summarise SSB across year and sim 
            group_by(across(all_of(group_columns))) %>%
            summarise(spbio=sum(spbio))
    )
}

#' Get Annual Fishing Mortality
#' 
#' Process MSE simulations for fishing mortality by fleet.
#' Total fishing mortality across fleets is alsoc computed.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_fishing_mortalities
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_fishing_mortalities(model_runs, extra_columns)
#' }
#'
get_fishing_mortalities <- function(model_runs, extra_columns){
    group_columns <- c("time", "fleet", "sim", "L1", names(extra_columns))
    
    return(
        bind_mse_outputs(model_runs, c("faa", "faa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            group_by(across(all_of(group_columns))) %>%
            # compute fleet-based F as the maximum F across age classes
            summarise(
                F = max(value)
            ) %>%
            ungroup() %>%
            group_by(across(all_of(group_columns[-2]))) %>%
            # total F is the sum of fleet-based Fs
            mutate(
                total_F = sum(F)
            ) %>%
            ungroup()
    )
}

#' Get Annual Recruits
#' 
#' Process MSE simulations for annual recruits.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_recruits
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_recruits(model_runs, extra_columns)
#' }
#'
get_recruits <- function(model_runs, extra_columns){
    group_columns <- c("time", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("naa", "naa_est"), extra_columns) %>% 
            as_tibble() %>%
            drop_na() %>%
            filter(age == 2) %>%
            group_by(across(all_of(group_columns))) %>%
            summarise(rec=sum(value))
    )
}

#' Get Landed Catches
#' 
#' Process MSE simulations for landed catches by fleet.
#' Total landed catch across fleets is alsoc computed.
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_landed_catch
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_landed_catch(model_runs, extra_columns)
#' }
#'
get_landed_catch <- function(model_runs, extra_columns){
    group_columns <- c("time", "fleet", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("land_caa"), extra_columns) %>%
            as_tibble() %>%
            drop_na() %>%
            group_by(across(all_of(group_columns))) %>%
            # compute fleet-based F as the maximum F across age classes
            summarise(
                catch = sum(value)
            ) %>%
            ungroup() %>%
            group_by(across(all_of(group_columns[-2]))) %>%
            # total F is the sum of fleet-based Fs
            mutate(
                total_catch = sum(catch)
            ) %>%
            ungroup()
    )
}

#' Get ABC, TAC, and Expected Landings
#' 
#' Process MSE simulations for ABC, TAC, and expected 
#' landings quantities
#'
#' @param model_runs list of completed MSE simualtion objects
#' @param extra_columns additional columns to append to output
#'
#' @export get_management_quantities
#'
#' @example \dontrun{
#'      mse1 <- run_mse(om, hcr1, ...)
#'      mse2 <- run_mse(om, hcr2, ...)
#' 
#'      model_runs <- list(mse1, mse2)
#'      extra_columns <- list(hcr=c("hcr1", "hcr2"))
#'      get_management_quantities(model_runs, extra_columns)
#' }
#'
get_management_quantities <- function(model_runs, extra_columns){
    cols <- c("time", "sim", "value", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs, c("abc", "tac", "exp_land"), extra_columns) %>%
            as_tibble() %>%
            drop_na() %>%
            select(cols)
    )
}

get_numbers_at_age <- function(model_runs, extra_columns){
    group_columns <- c("time", "class", "sim", "L1", names(extra_columns))
    return(
        bind_mse_outputs(model_runs2, c("naa"), extra_columns2) %>%
            as_tibble() %>%
            mutate(
                class = factor(
                    case_when(age < 3 ~ "1/2", age < 5 ~ "2/3", age < 7 ~ "3/4", age < 9 ~ "4/5", age < 15 ~ "5/7", age > 14 ~ "7+"), 
                    levels=c("1/2", "2/3", "3/4", "4/5", "5/7", "7+"), 
                    labels=c("Grade 1/2 (1-2yo)", "Grade 2/3 (3-4yo)", "Grade 3/4 (5-6yo)", "Grade 4/5 (7-8yo)", "Grade 5/7 (9-14yo)", "Grade 7+ (15+yo)")
                ),
                L1 = factor(L1, levels=c("caa", "naa"), labels=c("Catch-at-Age", "Numbers-at-Age"))
            ) %>%
            group_by(across(all_of(group_columns))) %>%
            summarise(value=sum(value))
    )
}

#' Get Reference Points from MSE Simulations
#' 
#' Derive fishing mortality and biomass reference points
#' from completed MSE simulations.
#' 
#' Note that this function hasn't been tested when multiple
#' MSE simulations are present in the `model_runs` list.
#'
#' @param model_runs list of completed MSE simulations
#' @param extra_columns additional columns that should be
#' appended to the resultant data frame
#' @param dem_params list of demographic parameter matrices
#' @param year the simualtion year to calculate reference 
#' points at
#'
#' @export get_reference_points
#'
#' @example
#' \dontrun{
#'      model_runs <- list(mse1)
#'      extra_columns <- list(hcr="hcr1")
#'      get_reference_points(model_runs, extra_columns, om$dem_params, nyears)
#' }
#'
get_reference_points <- function(model_runs, extra_columns){

    get_rps <- function(om_name, recruitment, prop_fs){
        om <- om_list[[which(names(om_list) == om_name)]]
        year <- 64
        joint_selret <- calculate_joint_selret(
            sel=om$dem_params$sel[year,,,,,drop=FALSE],
            ret=om$dem_params$ret[year,,,,,drop=FALSE],
            prop_fs = prop_fs
        )
        ref_pts <- calculate_ref_points(
            30,
            om$dem_params$mort[year,,1,1],
            om$dem_params$mat[year,,1,1],
            om$dem_params$waa[year,,1,1],
            joint_selret$sel[,,1,,drop=FALSE],
            joint_selret$ret[,,1,,drop=FALSE],
            recruitment/2,
            spr_target = 0.40
        )
        return(c(ref_pts$Fref, ref_pts$Bref, ref_pts$B0))
    }

    avg_recruitment <- get_recruits(model_runs, extra_columns) %>%
        group_by(sim, om) %>%
        summarise(rec=median(rec))

    prop_fs_df <- get_fishing_mortalities(model_runs, extra_columns) %>%
        filter(L1 != "faa_est") %>%
        group_by(time, sim, om, hcr, fleet) %>%
        mutate(
            prop_f = F/total_F
        ) %>%
        select(time, sim, fleet, om, hcr, prop_f) %>%
        distinct() %>%
        pivot_wider(names_from = "fleet", values_from="prop_f") %>%
        group_by(sim, om, hcr) %>%
        summarise(Fixed = mean(Fixed), Trawl=mean(Trawl))


    ref_pts_df <- prop_fs_df %>% 
        left_join(avg_recruitment, by=c("sim", "om")) %>%
        group_by(sim, om, hcr) %>%
        reframe(rps = get_rps(om, rec, c(Fixed, Trawl))) %>%
        mutate(rp_name = rep(c("F40", "B40", "B0"), length(hcr_list)*length(om_list)*length(seed_list))) %>%
        pivot_wider(names_from="rp_name", values_from="rps") %>%
        group_by(om) %>%
        median_qi(F40, B40, B0, .width=interval_widths, .simple_names=TRUE) %>%
        # mutate(bad=1) %>%
        # select(bad, 1:12) %>%
        reformat_ggdist_long(n=1)
        # select(-bad)

    return(ref_pts_df)

}
