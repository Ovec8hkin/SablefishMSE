rm(list=ls())

library(tidyverse)
library(ggplot2)
library(ggdist)
library(ggh4x)
library(reshape2)
library(SpatialSablefishAssessment)
library(tictoc)
library(doParallel)
# library(afscOM)
# library(afscOM) # may work but not certain

# Change to wherever your local copy of afscOM is
library(devtools)
afscOM_dir <- "~/Desktop/Projects/afscOM"
sablefishMSE_dir <- here::here()

devtools::load_all(afscOM_dir)

lapply(list.files("R", full.names = TRUE), source)

sable_om <- readRDS(file.path(here::here(), "data", "sablefish_om_big.RDS"))
assessment <- dget(file.path(here::here(), "data", "sablefish_assessment_2023.rdat"))
hist_recruits <- assessment$natage.female[,1]*2

dp_y <- afscOM::subset_dem_params(sable_om$dem_params, 64, d=1, drop=FALSE)
joint_selret <- calculate_joint_selret(
    sel = dp_y$sel,
    ret = dp_y$ret,
    prop_fs = c(0.80, 0.20)
)

spr40_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.40, 0.40)
)

spr50_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.50, 0.50)
)

spr4050_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.50, 0.40)
)

sprMSY_rp <- calculate_ref_points(
    nages=30,
    mort = dp_y$mort[,,1,],
    mat = dp_y$mat[,,1,],
    waa = dp_y$waa[,,1,],
    sel =  joint_selret$sel[,,1,,drop=FALSE],
    ret = joint_selret$ret[,,1,,drop=FALSE],
    avg_rec = mean(hist_recruits)/2,
    spr_target = c(0.45, 0.45)
)


convert_to_catch <- function(f, ssb){

    if(f == 0 || ssb == 0){
        return(0)
    }

    print(f)
    naa <- array(NA, dim=c(nages, 2))
    naa[1,] <- 1
    zaa <- dp_y$mort + joint_selret$sel*joint_selret$ret*F
    for(a in 2:(nages)){
        naa[a,] <- naa[a-1,]*exp(-zaa[1,a-1,,1])
    }
    equil_naa <- array(naa, dim=c(1, 30, 2, 1))
    equil_ssb <- sum(equil_naa[1,,1,1]*dp_y$mat[1,,1,1]*dp_y$waa[1,,1,1])

    ssb_frac <- ssb/equil_ssb

    equil_naa <- equil_naa*ssb_frac

    # sum(equil_naa*dp_y$waa[1,,,1]*joint_selret$sel[1,,,1])
    
    catch <- afscOM::baranov(f, equil_naa, dp_y$waa, dp_y$mort, joint_selret$sel)
    return(catch)
}


convert_to_f <- function(catch, ssb){
    catch_min <- function(f, catch,ssb){
        c <- convert_to_catch(exp(f),ssb)
        return(abs(c - catch))
    }

    if(catch == 0 || ssb == 0){
        return(0)
    }

    F <- bbmle::mle2(catch_min,
        start = list(f=log(0.05)),
        vecpar = TRUE,
        data=list(
            catch=catch,
            ssb=ssb
        ),
        method="Nelder-Mead",
        optimizer="nlminb",
        control=list(eval_max=5e6, iter_max=5e6)
    )
    return(exp(F@coef))
}

threshold_f_hcr <- function(ssb, FRPs, BRPs){
    return(threshold_f(ssb, f_min=FRPs[1], f_max=FRPs[2], lrp=BRPs[1], urp=BRPs[2]))
}

threshold_f_dep_hcr <- function(ssb, FRPs, BRPs, Bref){
    dep <- ssb/Bref
    return(threshold_f_hcr(dep, FRPs, BRPs))
}

pfmc_catch_special <- function(ssb, B0, Fref){
    dep <- ssb/B0
    # From Maia
    # sel_bio <- ssb_to_selbio(ssb)
    ABC <- convert_to_catch(Fref, ssb)#*sum(exploit_bio)
    # ABC <- OFL*exp(qnorm(pstar, 0, OFLsigma))
    if(dep < 0.1){TAC <- 0}
    if(dep > 0.40) {TAC <- ABC}
    if(dep <= 0.40 & dep >= 0.1){
        TAC <- ((convert_to_catch(Fref, 0.4*B0))/(B0*0.3))*(ssb - B0*0.1)
        # TAC_factor <- ((0.4*B0)/(B0*0.3))*(ssb - B0*0.1)
        # TAC <- Fref*ssb*((ssb-0.1*B0)/(0.4*B0-0.1*B0))
    }
    # TAC <- ABC
    return(TAC)
}

threshold_cap_hcr <- function(ssb, FRPs, BRPs, cap){
    f1 <- threshold_f_hcr(ssb, FRPs, BRPs)
    # std_sel <- joint_selret$sel[,,1,]/sum(joint_selret$sel[,,1,])

    # catch <- afscOM::F_to_mu(f1)*sum(ssb*std_sel)
    catch <- convert_to_catch(f1, ssb)
    if(catch > cap) catch <- cap

    # mu <- catch/sum(ssb*std_sel)
    f2 <- convert_to_f(catch, ssb)
    return(f2)
}

ssb_to_selbio <- function(ssb){
    ssb_to_bio_factor <- mean(assessment$t.series[,"totbiom"]/assessment$t.series[,"spbiom"])
    tot_bio <- ssb*ssb_to_bio_factor
    bio_atage <- tot_bio*(dp_y$waa/sum(dp_y$waa))
    sel_bio_atage <- bio_atage*joint_selret$sel

    return(sum(sel_bio_atage))

}


ssbs <- seq(0, 300, 1)

f40_hcr <- sapply(ssbs, \(x) threshold_f_hcr(x, FRPs=c(0, spr40_rp$Fref), BRPs=c(spr40_rp$Bref*0.05, spr40_rp$Bref)))
f50_hcr <- sapply(ssbs, \(x) threshold_f_hcr(x, FRPs=c(0, spr50_rp$Fref), BRPs=c(spr50_rp$Bref*0.05, spr50_rp$Bref)))
pfmc4010_hcr <- sapply(ssbs, \(x) threshold_f_dep_hcr(x, FRPs=c(0, sprMSY_rp$Fref), BRPs=c(0.1, 0.40), Bref=sprMSY_rp$B0))
pfmc4010_hcr_catch <- sapply(ssbs, \(x) pfmc_catch_special(x, B0 = sprMSY_rp$B0, Fref=sprMSY_rp$Fref))
bcsable_hcr <- sapply(ssbs, \(x) threshold_f_dep_hcr(x, FRPs=c(0, 0.055), BRPs=c(0.40, 0.60), Bref=sprMSY_rp$Bref))
cap15_hcr <- sapply(ssbs, \(x) threshold_cap_hcr(x, FRPs=c(0, spr40_rp$Fref), BRPs=c(spr40_rp$Bref*0.05, spr40_rp$Bref), cap=15))
cap20_hcr <- sapply(ssbs, \(x) threshold_cap_hcr(x, FRPs=c(0, spr40_rp$Fref), BRPs=c(spr40_rp$Bref*0.05, spr40_rp$Bref), cap=20))
cap25_hcr <- sapply(ssbs, \(x) threshold_cap_hcr(x, FRPs=c(0, spr40_rp$Fref), BRPs=c(spr40_rp$Bref*0.05, spr40_rp$Bref), cap=25))
chrf50_hcr <- rep(spr50_rp$Fref, length(ssbs))

hcr_df <- data.frame(
        ssb=ssbs, 
        f40=f40_hcr,
        perc5=f40_hcr,
        perc10=f40_hcr, 
        f50=f50_hcr,
        pfmc4010=pfmc4010_hcr, 
        bcsable=bcsable_hcr, 
        chr50=chrf50_hcr,
        cap15=cap15_hcr,
        cap25=cap25_hcr,
        # cap20=cap20_hcr,
        nofish=0
)

hcr_df_long <- hcr_df %>% pivot_longer(2:ncol(hcr_df), names_to="hcr", values_to="F") %>%
    mutate(
        hcr=factor(
            hcr, 
            levels=c("f40", "f50", "pfmc4010", "bcsable", "chr50", "perc5", "perc10", "cap15", "cap20", "cap25", "nofish"),
            labels=c("F40", "F50", "PFMC 40-10", "British Columbia", "Constant F50", "F40 +/- 5%", "F40 +/- 10%", "15k Harvest Cap", "20k Harvest Cap", "25k Harvest Cap", "No Fishing")
        )
    ) %>%
    rowwise() %>%
    mutate(
        catch = case_when(
            hcr %in% c("PFMC 40-10") ~ pfmc_catch_special(ssb, B0 = sprMSY_rp$B0, Fref=sprMSY_rp$Fref),
            TRUE ~ convert_to_catch(F, ssb)
        ),
        F = case_when(
            hcr %in% c("PFMC 40-10") ~ convert_to_f(catch, ssb),
            TRUE ~ F
        ),
        F = F/spr40_rp$Fref,
        ssb = ssb/spr40_rp$Bref
    )

hcr_colors <- set_hcr_colors2(hcr_df_long %>% pull(hcr) %>% unique)
# hcr_colors["No Fishing"] <- "#F8766D"
# hcr_colors["F40"] <- "black"

ggplot(hcr_df_long)+
    geom_line(aes(x=ssb, y=F, color=hcr, linetype="F"), linewidth=1)+
    geom_line(aes(x=ssb, y=catch/60, color=hcr, linetype="Catch"), linewidth=1)+
    # geom_vline(xintercept = spr40_rp$Bref, linetype="dashed", color="grey")+
    # geom_vline(xintercept = spr50_rp$B0, linetype="dashed", color="grey")+
    # annotate("text", x=spr40_rp$Bref-40, y=0.11, label="B[40]", parse=TRUE, size=6)+
    # annotate("text", x=spr40_rp$B0-40, y=0.11, label="B[100]", parse=TRUE, size=6)+
    scale_x_continuous("Relative SSB (SSB/B40)", expand=c(0, 0), breaks = seq(0, 2.5, 0.50))+
    scale_y_continuous("Relative Fishing Mortality (F/F40)", limits=c(0, 1.10), breaks=seq(0, 1.25, 0.25), sec.axis = sec_axis(transform = ~.*60, breaks=seq(0, 60, 15), name="Catch (1000 mt)"))+
    scale_linetype_manual(values=c("F"="solid", "Catch"="dotdash"))+
    scale_color_manual(values=hcr_colors)+
    coord_cartesian(xlim=c(0, 2.5))+
    # coord_cartesian(expand=0)+
    facet_wrap(~hcr, nrow=2)+
    guides(color="none", linetype=guide_legend("Unit"))+
    custom_theme+
    theme(panel.spacing.x = unit(1, "cm"))
ggsave(file.path(here::here(), "figures", "hcr_diagrams.jpeg"), width=14, height=8, units="in")


