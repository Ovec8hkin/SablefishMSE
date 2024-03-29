#' Format OM Observation Data for TMB Assessment Model
#' #'
#' Format data and observations from an `afscOM` operating model
#' and update the model data and parameters required by the 
#' `SpatialSablefishAssessment` TMB assessment model. 
#' 
#' This function requires the existence of 'data/sablefish_em_data_2022.RDS',
#' and 'data/sablefish_em_par_2022.RDS' files, which come packaged with the
#' `SablefishMSE` codebase. 
#'
#' @param nyears the number of years of data that will be passed to the model
#' @param dem_params demographic parameter matrices subsetted to 1 year
#' @param land_caa one year worth of landed catch-at-age data (dimensins [1, nages, nsexes, nregions, nfleets])
#' @param survey_indices one year worth of survey indices (LL RPN, LL RPW, and TW RPW)
#' @param fxfish_caa_obs one year worth of catch-at-age observation from the fixed gear fishery
#' @param ll_ac_obs one year worth of age composition observations frmo the longline survey
#' @param model_options list of model options provided to the OM
#' @param added_year the number of new years of data being added (should usually be 1) 
#'
#' @export format_em_data
#'
#' @example
#'
format_em_data <- function(nyears, dem_params, land_caa, survey_indices, fxfish_caa_obs, ll_ac_obs, model_options, added_years=1){

    # Read the most current set of data from an RDS file
    # If those files don't exist, then default back to 
    # using data through 2022.
    if(!file.exists("data/sablefish_em_data_curr.RDS")){
        data_file <- "data/sablefish_em_data_2022.RDS"
        param_file <- "data/sablefish_em_par_2022.RDS"
    }else{
        data_file <- "data/sablefish_em_data_curr.RDS"
        param_file <- "data/sablefish_em_par_curr.RDS"
    }

    new_data <- readRDS(data_file)

    extra_years <- added_years

    # extract ages
    ages <- as.numeric(colnames(dem_params$waa))
    new_data$ages <- ages

    # extract years
    years <- as.double(1960:(1960+nyears-1))
    new_data$years <- years

    # extract maturity
    maturity <- dem_params$mat[,,1,1]
    new_data$maturity <- cbind(new_data$maturity, maturity)

    # extract mortality
    mortality <- dem_params$mort[,,1,1]
    new_data$M <- cbind(new_data$M, mortality)

    # extract WAA
    waa <- dem_params$waa[,,,1, drop=FALSE]
    waa_m <- waa[,,2,]
    waa_f <- waa[,,1,]

    new_data$male_mean_weight_by_age <- cbind(new_data$male_mean_weight_by_age, waa_m)
    new_data$female_mean_weight_by_age <- cbind(new_data$female_mean_weight_by_age, waa_f)

    # age-length transition matrix
    male_al_trans_mat <- SpatialSablefishAssessment::extend_3darray_last_dim(new_data$male_age_length_transition, n=extra_years)
    female_al_trans_mat <- SpatialSablefishAssessment::extend_3darray_last_dim(new_data$female_age_length_transition, n=extra_years)

    new_data$male_age_length_transition <- male_al_trans_mat
    new_data$female_age_length_transition <- female_al_trans_mat

    # proportion male
    prop_male <- rep(0.5, nyears)
    new_data$proportion_male <- prop_male
    new_data$proportion_male2 <- prop_male

    # spawning time proportion
    spawn_time <- rep(0, nyears)
    new_data$spawning_time_proportion <- spawn_time

    # fishery catch vectors
    #landed_caa <- subset_matrix(land_caa[,,,,,drop=FALSE], 1, 6)
    ll_catch <- sum(land_caa[,,,,1], na.rm=TRUE)
    tw_catch <- sum(land_caa[,,,,2], na.rm=TRUE)

    new_data$ll_fishery_catch <- c(new_data$ll_fishery_catch, ll_catch)
    new_data$trwl_fishery_catch <- c(new_data$trwl_fishery_catch, tw_catch)

    # fishery selectivity
    ll_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_sel_by_year_indicator, n=extra_years)
    tw_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$trwl_sel_by_year_indicator, n=extra_years)
    jp_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_sel_by_year_indicator, n=extra_years)

    new_data$ll_sel_by_year_indicator <- ll_sel_yr_indic
    new_data$trwl_sel_by_year_indicator <- tw_sel_yr_indic
    new_data$srv_jap_fishery_ll_sel_by_year_indicator <- jp_sel_yr_indic

    # survey selectivity and catchability
    srv_ll_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_sel_by_year_indicator, n=extra_years)
    srv_tw_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_nmfs_trwl_sel_by_year_indicator, n=extra_years)
    srv_jp_sel_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_sel_by_year_indicator, n=extra_years)

    new_data$srv_dom_ll_sel_by_year_indicator <- srv_ll_sel_yr_indic
    new_data$srv_nmfs_trwl_sel_by_year_indicator <- srv_tw_sel_yr_indic
    new_data$srv_jap_ll_sel_by_year_indicator <- srv_jp_sel_yr_indic

    srv_ll_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_q_by_year_indicator, n=extra_years)
    srv_tw_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_nmfs_trwl_q_by_year_indicator, n=extra_years)
    srv_jp_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_q_by_year_indicator, n=extra_years)
    jp_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_q_by_year_indicator, n=extra_years)

    new_data$srv_dom_ll_q_by_year_indicator <- srv_ll_q_yr_indic
    new_data$srv_nmfs_trwl_q_by_year_indicator <- srv_tw_q_yr_indic
    new_data$srv_jap_ll_q_by_year_indicator <- srv_jp_q_yr_indic
    new_data$srv_jap_fishery_ll_q_by_year_indicator <- jp_q_yr_indic

    ll_cpue_q_yr_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_cpue_q_by_year_indicator, n=extra_years)
    new_data$ll_cpue_q_by_year_indicator <- ll_cpue_q_yr_indic

    # catch at age new_data
    ll_caa_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_catchatage_indicator, n=extra_years)
    ll_caa_indic[(length(ll_caa_indic)-1)] <- 1
    new_data$ll_catchatage_indicator <- ll_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_ll_caa <- apply(fxfish_caa_obs[,,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_ll_catchatage <- cbind(new_data$obs_ll_catchatage, apply(obs_ll_caa, 1, as.double))

    srv_ll_caa_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_age_indicator, n=extra_years)
    srv_ll_caa_indic[(length(years)-1)] <- 1
    new_data$srv_dom_ll_age_indicator <- srv_ll_caa_indic

    # Be really careful here, because we actually need to pass the PREVIOUS
    # years age composition observations, because they are always one year
    # delayed.
    obs_srv_ll_caa <- apply(ll_ac_obs[,,,,drop=FALSE], c(1, 2), sum)
    new_data$obs_srv_dom_ll_age <- cbind(new_data$obs_srv_dom_ll_age, apply(obs_srv_ll_caa, 1, as.double))

    srv_jpll_caa_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_age_indicator, n=extra_years)
    new_data$srv_jap_ll_age_indicator <- srv_jpll_caa_indic

    srv_tw_caa_indic <- c(new_data$srv_nmfs_trwl_age_indicator, rep(0, length.out=extra_years))
    new_data$srv_nmfs_trwl_age_indicator <- srv_tw_caa_indic

    # catch at length
    ll_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_catchatlgth_indicator, n=extra_years)
    tw_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$trwl_catchatlgth_indicator, n=extra_years)

    new_data$ll_catchatlgth_indicator <- ll_cal_indic
    new_data$trwl_catchatlgth_indicator <- tw_cal_indic

    srv_ll_cal_indic <- c(new_data$srv_dom_ll_lgth_indicator, rep(0, length.out=extra_years))
    srv_tw_cal_indic <- c(new_data$srv_nmfs_trwl_lgth_indicator, rep(0, length.out=extra_years))

    new_data$srv_dom_ll_lgth_indicator <- srv_ll_cal_indic
    new_data$srv_nmfs_trwl_lgth_indicator <- srv_tw_cal_indic

    srv_jpll_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_lgth_indicator, n=extra_years)
    jp_ll_cal_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_lgth_indicator, n=extra_years)

    new_data$srv_jap_ll_lgth_indicator <- srv_jpll_cal_indic
    new_data$srv_jap_fishery_ll_lgth_indicator <- jp_ll_cal_indic

    # Survey new_data
    srv_ll_rpn_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_dom_ll_bio_indicator, n=extra_years)
    srv_ll_rpn_obs <- survey_indices$ll_rpn
    srv_ll_rpn_ses <- model_options$obs_pars$surv_ll$rpn_cv*survey_indices$ll_rpn

    new_data$srv_dom_ll_bio_indicator <- srv_ll_rpn_indic
    new_data$obs_dom_ll_bio <- c(new_data$obs_dom_ll_bio, srv_ll_rpn_obs)
    new_data$se_dom_ll_bio <- c(new_data$se_dom_ll_bio, srv_ll_rpn_ses)

    srv_jpll_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_ll_bio_indicator, n=extra_years)
    new_data$srv_jap_ll_bio_indicator <- srv_jpll_indic

    # Be careful here. The trawl survey technically only happens ever other
    # year, but its being simulated occurring in every year.
    srv_tw_rpw_indic <- c(new_data$srv_nmfs_trwl_bio_indicator, rep(1, length.out=extra_years))#c(new_data$srv_nmfs_trwl_bio_indicator, rep(c(0, 1), length.out=extra_years))
    srv_tw_rpw_obs <- survey_indices$tw_rpw
    srv_tw_rpw_ses <- model_options$obs_pars$surv_tw$rpw_cv*survey_indices$tw_rpw

    new_data$srv_nmfs_trwl_bio_indicator <- srv_tw_rpw_indic
    new_data$obs_nmfs_trwl_bio <- c(new_data$obs_nmfs_trwl_bio, srv_tw_rpw_obs)
    new_data$se_nmfs_trwl_bio <- c(new_data$se_nmfs_trwl_bio, srv_tw_rpw_ses)

    ll_cpue_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$ll_cpue_indicator, n=extra_years)
    jp_ll_rpw_indic <- SpatialSablefishAssessment::extend_vec_last_val(new_data$srv_jap_fishery_ll_bio_indicator, n=extra_years)

    new_data$ll_cpue_indicator <- ll_cpue_indic
    new_data$srv_jap_fishery_ll_bio_indicator <- jp_ll_rpw_indic

    # Parameters

    new_parameters <- readRDS(param_file)

    ln_M_year_devs <- SpatialSablefishAssessment::extend_vec_last_val(new_parameters$ln_M_year_devs, n=extra_years)
    ln_M_age_devs  <- SpatialSablefishAssessment::extend_vec_last_val(new_parameters$ln_M_age_devs, n=extra_years)
    ln_rec_dev <- c(new_parameters$ln_rec_dev, rep(0, extra_years))
    ln_ll_F <- c(new_parameters$ln_ll_F_devs, rep(0, extra_years))
    ln_tw_F <- c(new_parameters$ln_trwl_F_devs, rep(0, extra_years))

    new_parameters$ln_M_year_devs <- ln_M_year_devs
    new_parameters$ln_M_age_devs <- ln_M_age_devs
    new_parameters$ln_rec_dev <- ln_rec_dev
    new_parameters$ln_ll_F_devs <- ln_ll_F
    new_parameters$ln_trwl_F_devs <- ln_tw_F

    # my_model = TMB::MakeADFun(data = new_data,
    #                            parameters = new_parameters,
    #                            DLL = "SpatialSablefishAssessment_TMBExports")

    # mle_optim = nlminb(start = my_model$par, objective = my_model$fn, gradient  = my_model$gr, control = list(iter.max = 10000, eval.max = 10000))
    # mle_report = my_model$report(mle_optim$par)  
    # assessment_ssb <- SpatialSablefishAssessment::get_SSB(mle_report) %>% filter(Year == max(Year)) %>% pull(SSB)

    saveRDS(new_data, "data/sablefish_em_data_curr.RDS")
    saveRDS(new_parameters, "data/sablefish_em_par_curr.RDS")

    capture.output(valid <- SpatialSablefishAssessment::validate_input_data_and_parameters(new_data, new_parameters))

    if(valid){
        return(afscOM::listN(new_data, new_parameters))
    }else{
        print("Something was wrong with the EM data.")
        return(FALSE)
    }
}
