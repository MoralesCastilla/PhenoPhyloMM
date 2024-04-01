#'#############################################################
#' Phylogenetic estimates of species-level phenology improve ecological forecasting 
#' 
#' * Script #1 Phylogenetic Mixed Model fitting
#'
#'  
#'  by Morales-Castilla, I., et al. 
#'  feb 2024
#'#############################################################


## Runs (or reads) the phylogeny models, extracts some output
## Does some basic plotting

rm(list=ls())
options(stringsAsFactors = FALSE)

# Setting working directory. Add in your own path in an if statement for your file structure
setwd("~/your/local/address/")


# Loading packages
library(caper)
library(pez)
library(phytools)
library(rstan)
library(shinystan)
library(plyr)
library(dplyr)

options(mc.cores = parallel::detectCores())



#'######################################
#### load data and phylogeny ####
#'######################################


  d = read.csv("data/ospreebbphyloms_forknb.csv")
  phylo = read.tree("data/phyloforphyloms.tre")
  




#'###################################
# Run  the models      ####
#'###################################

## Fit model here
  fitlambest <- stan("stan_code/PhenoPhyloMM_PMM.stan",
              data=list(N=nrow(d),
                        n_sp=nspecies,
                        sp=d$sppnum,
                        x1=d$force.z,
                        x2 = d$chill.z,
                        x3=d$photo.z,
                        y=d$resp,
                        Vphy=vcv(phylo, corr = TRUE)),
              iter = 4000,
              warmup = 2000, # half the iter as warmp is default, but leaving in case we want to change
              chains = 4,
              seed = 117 
  )
  
  ## Save fitted posterior
  saveRDS(fitlambest, "output/fit_model_PMM.rds")
  
  fitlamb0 <- stan("stan/PhenoPhyloMM_HMM.stan",
                   data=list(N=nrow(d),
                             n_sp=nspecies,
                             sp=d$sppnum,
                             x1=d$force.z,
                             x2 = d$chill.z,
                             x3=d$photo.z,
                             y=d$resp,
                             Vphy=vcv(phylo, corr = TRUE)),
                   iter = 4000,
                   warmup = 2000, 
                   chains = 4,
                   seed = 117 
  )
  saveRDS(fitlamb0, "output/fit_model_HMM.rds")
  


#'###################################
# Explore model fit            ####
#'###################################

## Summarize full fit
# summary(fit)$summary

## Summarize lambdas, b_zf, b_zc, , b_zp, intercept mean, and sigmas
fitsum <- summary(fitlambest, pars = list("a_z", "sigma_interceptsa", 
                                   "b_zf", "sigma_interceptsbf", "lam_interceptsbf", 
                                   "b_zc", "sigma_interceptsbc", "lam_interceptsbc",
                                   "b_zp", "sigma_interceptsbp", "lam_interceptsbp","sigma_y"))$summary

fitsumdf <- as.data.frame(fitsum)

source("source/stan_utility.R")
check_all_diagnostics(fitlambest)





#'###############################################
#### comparing estimates lambda est vs 1 vs 0 ####
#'###############################################


## load models



## Summarize lambdas, b_zf, b_zc, , b_zp, intercept mean, and sigmas
tableresults.0 = summary(fitlam0, pars = list("a_z", "sigma_interceptsa", "b_zf", "sigma_interceptsbf", "b_zc", "sigma_interceptsbc", "b_zp", "sigma_interceptsbp", "sigma_y"))$summary
tableresults.est = summary(fitlambest, pars = list("a_z", "lam_interceptsa", "sigma_interceptsa", "b_zf", "lam_interceptsbf", "sigma_interceptsbf", "b_zc", "lam_interceptsbc", "sigma_interceptsbc", "b_zp", "lam_interceptsbp", "sigma_interceptsbp", "sigma_y"))$summary



## rename model to include species names
names(fitlambest)[grep(pattern = "^a\\[", x = names(fitlambest))] <- phylo$tip.label
names(fitlambest)[grep(pattern = "^b_force", x = names(fitlambest))] <- phylo$tip.label
names(fitlambest)[grep(pattern = "^b_chill", x = names(fitlambest))] <- phylo$tip.label
names(fitlambest)[grep(pattern = "^b_photo", x = names(fitlambest))] <- phylo$tip.label

names(fitlam0)[grep(pattern = "^a\\[", x = names(fitlam0))] <- phylo$tip.label
names(fitlam0)[grep(pattern = "^b_force", x = names(fitlam0))] <- phylo$tip.label
names(fitlam0)[grep(pattern = "^b_chill", x = names(fitlam0))] <- phylo$tip.label
names(fitlam0)[grep(pattern = "^b_photo", x = names(fitlam0))] <- phylo$tip.label




# get model estimates per species ----

## where species are

posspsindata.est <- list(10:200,202:392,394:584)
posspsindata.01 <- list(6:196,198:388,390:580)


## forcing
cueforce = summary(fitlambest)$summary[posspsindata.est[[1]],"mean"]
cueforcesdup = summary(fitlambest)$summary[posspsindata.est[[1]],"75%"]
cueforcesdlow = summary(fitlambest)$summary[posspsindata.est[[1]],"25%"]

cueforce0 = summary(fitlam0)$summary[posspsindata.01[[1]],"mean"]
cueforcesdup0 = summary(fitlam0)$summary[posspsindata.01[[1]],"75%"]
cueforcesdlow0 = summary(fitlam0)$summary[posspsindata.01[[1]],"25%"]


## chill
cuechill = summary(fitlambest)$summary[posspsindata.est[[2]],"mean"]
cuechillsdup = summary(fitlambest)$summary[posspsindata.est[[2]],"75%"]
cuechillsdlow = summary(fitlambest)$summary[posspsindata.est[[2]],"25%"]

cuechill0 = summary(fitlam0)$summary[posspsindata.01[[2]],"mean"]
cuechillsdup0 = summary(fitlam0)$summary[posspsindata.01[[2]],"75%"]
cuechillsdlow0 = summary(fitlam0)$summary[posspsindata.01[[2]],"25%"]


## photo
cuephoto = summary(fitlambest)$summary[posspsindata.est[[3]],"mean"]
cuephotosdup = summary(fitlambest)$summary[posspsindata.est[[3]],"75%"]
cuephotosdlow = summary(fitlambest)$summary[posspsindata.est[[3]],"25%"]

cuephoto0 = summary(fitlam0)$summary[posspsindata.01[[3]],"mean"]
cuephotosdup0 = summary(fitlam0)$summary[posspsindata.01[[3]],"75%"]
cuephotosdlow0 = summary(fitlam0)$summary[posspsindata.01[[3]],"25%"]





### plot correlations angio ----
plotting = F
lambdazero = F

if(plotting){
  
  dev.off()
  par(mfrow=c(1,3))
  
  virid <-  colorRampPalette(c("yellow","darkcyan","purple"))
  
  colschill <- virid(30)[as.numeric(cut(c(cuechill0, cuechill),breaks = 30))]
  colschillhmm <- colschill[1:length(cuechill0)]
  colschillpmm <- colschill[(length(cuechill0)+1):length(colschill)]
  
  
  plot(cuechill0, cuechill, 
       xlab="sensitivity to chilling HMM",
       ylab="sensitivity to chilling PMM", 
       pch=16, col=adjustcolor(colschillpmm,0.4),cex=1.2, cex.lab=1.5,
       xlim=c(-30,5),ylim=c(-30,5))
  abline(v=mean(cuechill0), col='grey', lty=2, lwd=2)  
  
  for(i in 1:length(cueforce0)){
    lines(c(cuechillsdlow0[i],cuechillsdup0[i]),
          rep(cuechill[i],2), col=adjustcolor(colschillpmm[i],0.2))
    lines(rep(cuechill0[i],2),
          c(cuechillsdlow[i],cuechillsdup[i]),
          col=adjustcolor(colschillhmm[i],0.2))
  }
  points(cuechill0, cuechill,pch=16, col=adjustcolor(colschillpmm,0.4),cex=1.2)
  
  abline(a=0,b=1, col='darkgrey', lty=2, lwd=1.5)  
  #abline(lm(cuechill~cuechill0), lwd=1.5)
  mtext("a", side = 3, adj = 0.05,line=-2,cex=1.5)
  
  
  colsforce <- virid(30)[as.numeric(cut(c(cueforce0, cueforce),breaks = 30))]
  colsforcehmm <- colsforce[1:length(cueforce0)]
  colsforcepmm <- colsforce[(length(cueforce0)+1):length(colsforce)]
  
  plot(cueforce0, cueforce, 
       xlab="sensitivity to forcing HMM",
       ylab="sensitivity to forcing PMM", 
       pch=16, col=adjustcolor(colsforcepmm,0.4),cex=1.2, cex.lab=1.5,
       xlim=c(-20,5),ylim=c(-20,5))
  abline(v=mean(cueforce0), col='grey', lty=2, lwd=2)  
  
  for(i in 1:length(cueforce0)){
    lines(c(cueforcesdlow0[i],cueforcesdup0[i]),
          rep(cueforce[i],2), col=adjustcolor(colsforcepmm[i],0.2))
    lines(rep(cueforce0[i],2),
          c(cueforcesdlow[i],cueforcesdup[i]),
          col=adjustcolor(colsforcehmm[i],0.2))
    
  }
  points(cueforce0, cueforce,pch=16, col=adjustcolor(colsforcepmm,0.4),cex=1.2)
  
  abline(a=0,b=1, col='darkgrey', lty=2, lwd=1.5)  
  #abline(lm(cueforce~cueforce0), lwd=1.5)
  mtext("b", side = 3, adj = 0.05,line=-2,cex=1.5)
  
  colsphoto <- virid(30)[as.numeric(cut(c(cuephoto0, cuephoto),breaks = 30))]
  colsphotohmm <- colsphoto[1:length(cuephoto0)]
  colsphotopmm <- colsphoto[(length(cuephoto0)+1):length(colsphoto)]
  
  plot(cuephoto0, cuephoto, 
       xlab="sensitivity to photoperiod HMM",
       ylab="sensitivity to photoperiod PMM", 
       pch=16, col=adjustcolor(colsphotohmm,0.4),cex=1.2, cex.lab=1.5,
       xlim=c(-10,3),ylim=c(-10,3))
  abline(v=mean(cuephoto0), col='grey', lty=2, lwd=2)  
  
  for(i in 1:length(cuephoto0)){
    lines(c(cuephotosdlow0[i],cuephotosdup0[i]),
          rep(cuephoto[i],2), col=adjustcolor(colsphotohmm[i],0.2))
    
    lines(rep(cuephoto0[i],2),
          c(cuephotosdlow[i],cuephotosdup[i]),
          col=adjustcolor(colsphotohmm[i],0.2))
  }
  points(cuephoto0, cuephoto,pch=16, col=adjustcolor(colsphotopmm,0.4),cex=1.2)
  
  abline(a=0,b=1, col='darkgrey', lty=2, lwd=1.5)  
  #abline(lm(cuephoto~cuephoto0), lwd=1.5)
  mtext("c", side = 3, adj = 0.05,line=-2,cex=1.5)
  
  
}



# end ----
