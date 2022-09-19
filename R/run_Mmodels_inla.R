################################################################################
## Title: Bayesian inference in multivariate spatio-temporal areal models     ##
##        using INLA: analysis of gender-based violence in small areas        ##
##                                                                            ##
## Authors: Vicente, G. - Goicoa, T. -  Ugarte, M.D.                          ##
##                                                                            ##
## https://doi.org/10.1007/s00477-020-01808-x                                 ##
##                                                                            ##
################################################################################
##                      Spatio-temporal M-models (INLA)                       ##
################################################################################
rm(list=ls())
library(INLA)

## Some necessary packages ##
packages <- c("dplyr","sf","spdep","spatialreg","fastDummies")

## Install packages not yet installed ##
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## Packages loading ##
invisible(lapply(packages, library, character.only = TRUE))

## Folder to save results ##
if(!file.exists("results")) {dir.create("results")}


########################################
## 1) Load data and cartography files ##
########################################
load("./dataMmodel.RData")

head(carto_UP)
str(data_UP)

plot(carto_india$geometry)
plot(carto_UP$geometry, add=T, col="lightblue")

J <- length(data_UP)
S <- length(unique(data_UP[[1]]$dist))
T <- length(unique(data_UP[[1]]$year))

data <- do.call(rbind,data_UP)
data$Crime <- rep(1:J,each=S*T)


####################################
## 2) Fit spatio-temporal M-model ##
####################################
source("functions/INLA_Mmodel_icar.R")  ## Intrinsic multivariate CAR latent effect
source("functions/INLA_Mmodel_lcar.R")  ## Leroux multivariate CAR latent effect
source("functions/INLA_Mmodel_pcar.R")  ## Proper multivariate CAR latent effect
source("functions/MCAR_INLA_ST.R")      ## Fits INLA models using different priors for the spatial, temporal and spatio/temporal random effects


############################################
## 2.1) Multivariate intrinsic CAR models ##
############################################

## Additive 
icar.ad <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="intrinsic", temporal="rw1", interaction="none",
                        strategy="simplified.laplace")

## Type I
icar.t1 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="intrinsic", temporal="rw1", interaction="TypeI",
                        strategy="simplified.laplace")

## Type II
icar.t2 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="intrinsic", temporal="rw1", interaction="TypeII",
                        strategy="simplified.laplace")

## Type III
icar.t3 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="intrinsic", temporal="rw1", interaction="TypeIII",
                        strategy="simplified.laplace")

## Type IV
icar.t4 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="intrinsic", temporal="rw1", interaction="TypeIV",
                        strategy="simplified.laplace")


## Save the models ##
MODELS.inla.icar <- list(Additive=icar.ad, TypeI=icar.t1, TypeII=icar.t2, TypeIII=icar.t3, TypeIV=icar.t4)

save(MODELS.inla.icar, file=paste0("./results/", gsub("\\.", "_", "MODELS.inla.icar"), ".RData"))


#########################################
## 2.2) Multivariate Leroux CAR models ##
#########################################

## Additive 
lcar.ad <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="Leroux", temporal="rw1", interaction="none",
                        strategy="simplified.laplace")

## Type I
lcar.t1 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="Leroux", temporal="rw1", interaction="TypeI",
                        strategy="simplified.laplace")

## Type II
lcar.t2 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="Leroux", temporal="rw1", interaction="TypeII",
                        strategy="simplified.laplace")

## Type III
lcar.t3 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="Leroux", temporal="rw1", interaction="TypeIII",
                        strategy="simplified.laplace")

## Type IV
lcar.t4 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="Leroux", temporal="rw1", interaction="TypeIV",
                        strategy="simplified.laplace")


## Save the models ##
MODELS.inla.lcar <- list(Additive=lcar.ad, TypeI=lcar.t1, TypeII=lcar.t2, TypeIII=lcar.t3, TypeIV=lcar.t4)

save(MODELS.inla.lcar, file=paste0("./results/", gsub("\\.", "_", "MODELS.inla.lcar"), ".RData"))



#########################################
## 2.2) Multivariate proper CAR models ##
#########################################

## Additive 
pcar.ad <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="proper", temporal="rw1", interaction="none",
                        strategy="simplified.laplace")

## Type I
pcar.t1 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="proper", temporal="rw1", interaction="TypeI",
                        strategy="simplified.laplace")

## Type II
pcar.t2 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="proper", temporal="rw1", interaction="TypeII",
                        strategy="simplified.laplace")

## Type III
pcar.t3 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="proper", temporal="rw1", interaction="TypeIII",
                        strategy="simplified.laplace")

## Type IV
pcar.t4 <- MCAR_INLA_ST(carto=carto_UP, data=data, ID.area="dist", ID.year="year", ID.disease="Crime",
                        O="obs", E="exp", spatial="proper", temporal="rw1", interaction="TypeIV",
                        strategy="simplified.laplace")


## Save the models ##
MODELS.inla.pcar <- list(Additive=pcar.ad, TypeI=pcar.t1, TypeII=pcar.t2, TypeIII=pcar.t3, TypeIV=pcar.t4)

save(MODELS.inla.pcar, file=paste0("./results/", gsub("\\.", "_", "MODELS.inla.pcar"), ".RData"))
