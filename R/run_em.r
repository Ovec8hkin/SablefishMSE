#' Run Estimation Method 
#' 
#' Given an OM object and a succesfull MSE run, apply
#' the EM to the observation data generated by the MSE.
#'
#' @param s simulation seed
#' @param sable_om om object
#' @param mse_tier3 mse object
#' @param om string identifed for OM
#'
#' @export
#'
#' @example
#'
run_EM <- function(s, sable_om, mse_tier3, om, fix_pars=NA){
    afscOM_dir <- "~/Desktop/Projects/afscOM"
    
    library(devtools)
    library(TMB) 
    devtools::load_all(afscOM_dir)

    source("R/format_em_data.R")
    source("R/fit_TMB_model.R")

    i=160-63
    y <- 63+i

    assess_inputs <- format_em_data(
                        nyears = y,
                        dem_params = afscOM::subset_dem_params(sable_om$dem_params, 64:y, d=1, drop=FALSE),
                        land_caa = afscOM::subset_matrix(mse_tier3$land_caa[64:y,,,,,s,drop=FALSE], 1, d=6, drop=TRUE),
                        survey_indices = afscOM::subset_dem_params(afscOM::subset_dem_params(mse_tier3$survey_obs, 64:y, d=1, drop=FALSE), s, d=5, drop=TRUE),
                        fxfish_caa_obs = afscOM::subset_matrix(mse_tier3$survey_obs$fxfish_acs[63:(y-1),,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                        ll_ac_obs = afscOM::subset_matrix(mse_tier3$survey_obs$ll_acs[63:(y-1),,,,s,drop=FALSE], 1, d=5, drop=TRUE), # Age comp data is one year delayed
                        model_options = sable_om$model_options,
                        added_years = i,
                        file_suffix = s
                    )
    file.remove(paste0("data/sablefish_em_data_curr_",s,".RDS"))
    file.remove(paste0("data/sablefish_em_par_curr_",s,".RDS"))
    mod_out <- fit_TMB_model(assess_inputs$new_data, assess_inputs$new_parameters, fix_pars = fix_pars)
    #mod_out$opt$par["ln_mean_rec"]

    om_ssb <- apply(mse_tier3$naa[1:y,,1,1,s]*sable_om$dem_params$mat[1:y,,1,1]*sable_om$dem_params$waa[1:y,,1,1], 1, sum)
    om_rec <- apply(mse_tier3$naa[1:y, 1, ,1,s], 1, sum)
    om_f <- mse_tier3$out_f[1:y,1,1,1,s]

    assess_df <- SpatialSablefishAssessment::get_SSB(mod_out$report)
    assess_rec <- mod_out$report$natage_f[1,]#SpatialSablefishAssessment::get_recruitment(mod_out$report)
    assess_f <- SpatialSablefishAssessment::get_fishing_mortalities(mod_out$report) %>% group_by(Year) %>% dplyr::summarise(F=sum(F))

    #assess_df <- data.frame(Year=1960:(1960+160-1))
    #assess_df$ssb <- assess_ssb
    assess_df$om_ssb <- om_ssb
    assess_df$sim <- s
    assess_df$om <- om
    assess_df$rec <- assess_rec[-length(assess_rec)]
    assess_df$om_rec <- om_rec
    assess_df$f <- assess_f$F
    assess_df$om_f <- om_f
    return(assess_df)
}

# s <- 1
# sable_om <- om_sexed_small
# mse_tier3 <- mse_sexed_small
# om <- "smsx"
# fix_pars <- list(
#     "ln_mean_rec" = mean(log(mse_tier3$naa[,1,1,1,s])),
#     "ln_rec_dev" = log(mse_tier3$naa[,1,1,1,s]) - mean(log(mse_tier3$naa[,1,1,1,s])) + 1.04^2/2
# )

# assess_ssb <- run_EM(s, sable_om, mse_tier3, om, fix_pars)

# ggplot(assess_ssb)+
#     geom_line(aes(x=Year, y=om_rec))+
#     geom_line(aes(x=Year, y=rec), color="red")

# log(mean(recruitment))

# log(recruitment)-mean(log(recruitment)) + 

# plot(1:64, log(recruitment), type="l")
# points(1:64, log(recruitment))
# lines(1:66, rec, col="blue")
# abline(h=mean(log(recruitment)), col="red")


# ln_rec_devs <- c(-0.0201727896823, -0.0211653628442, -0.0222762796367, -0.0235167887645, -0.0248376531671, -0.0263140331771, -0.0279775591988, -0.0298335784895, -0.0318628708855, -0.0341045755794, -0.0365927819101, -0.0393071338025, -0.0422739844026, -0.0455070747877, -0.0490012744867, -0.0527273187832, -0.0566564781594, -0.0607215085617, -0.0647839822855, -0.0686725000884, -0.0721473710816, -0.0747356469988, -0.0756711664860, -0.0736932654205, -0.0674020989957, -0.0545598257407, -0.0326432856417, 0.000831967151461, 0.0477478898283, 0.106462820860, 0.170899944240, 0.228214073308, 0.260537688760, 0.272295671697, 0.223888837810, 0.138327569143, 0.0409977654323, -0.0504299349341, -0.139913190910, -0.225833906653, -0.298717396287, -0.351735855526, -0.363169725158, -1.08102049603, -0.951503149105, -0.900354062579, -0.671126924421, 1.14538462882, 0.993817636166, 0.167556465633, 1.45911595657, 0.686962952511, -0.296119130703, -0.240150999009, 0.105223647861, -0.806134517718, -1.20159915722, -1.10370002507, -0.280362194813, 0.385861443905, -0.677219759900, 0.451877334609, -0.816981612496, -0.754631289491, -0.241530322750, 0.269235393826, -0.431497255133, 0.853730546539, 0.0317112429988, 0.0425907663835, 1.01280894948, -0.151876454291, -0.447976422128, -0.272519813754, -0.712160267655, -0.413806228641, -0.468138784792, 0.0159327161347, 0.319484365520, -0.407168000172, -0.321847602343, -1.19277559360, -0.726547291688, -0.0723645963570, 1.16309263388, 0.337204769385, 1.81806716488, 1.85659954644, 1.23682827751, 1.98814135583, 0.468638127271)
# ln_mean_rec <- 3.28521376476

# # bias_ramp <- rep(0, length(1933:2023))
# # for(i in 1:length(bias_ramp)){
# #     if(i < 50){
# #         bias_ramp[i] <- 0
# #     }else if(i >=50 & i < 60){
# #         y <- i + 1930
# #         bias_ramp[i] <- 1*(1-((y-1980)/(1990-1980)))
# #     }else if(i >= 60 & i < 91){
# #         bias_ramp[i] <- 1
# #     }else{
# #         y <- i + 1930
# #         bias_ramp[i] <- 1*(1-((y-2020)/(2023-2020)))
# #     }
# # }

# sig_rs <- c(rep(0.15, length(1930:1975)), rep(exp(0.0435597461633), 5), rep(exp(0.0435597461633), length(1981:2020)))

# years <- 1930:2023
# rec <- log(exp(ln_mean_rec + ln_rec_devs - (sig_rs^2/2))[years >= 1958])


# # log_mean_rec: