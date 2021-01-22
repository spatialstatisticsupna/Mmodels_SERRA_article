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
########################
## Some necessary packages (This only installs packages not yet installed)
########################
## Package names
packages <- c("spdep", "spData", "rgdal","fastDummies")

## Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## Packages loading
invisible(lapply(packages, library, character.only = TRUE))


## Folder to save results
if(!file.exists("resul")) {dir.create("resul")}

########################################
## Data organization for INLA         ##
########################################
## Load data and Uttar Pradesh sf
load("./dataMmodel.RData")
rm(list = c("Mzz","carto_india"))

crimes<- c("rape", "dowry")

n <- length(unique(data$dist))
t <- length(unique(data$year))
J <- length(crimes)

data <- data[order(data[,"ID_year"], data[,"ID_area"]),]
obs<- reshape2::melt(data, id.vars = c("dist","year"), measure.vars = crimes,
                     variable.name = "crime", value.name = "observed")
exp<- reshape2::melt(data, id.vars = c("dist","year"), measure.vars = paste0("e_",crimes),
                     variable.name = "crime", value.name = "expected")

df.inla <- data.frame(dist= obs[,"dist"], year= obs[,"year"],
                      observed= obs[,"observed"], expected= exp[,"expected"],
                      crime=obs[,"crime"])

Data_UP<- vector("list",2)
names(Data_UP)<- crimes
Data_UP$rape<- df.inla[df.inla$crime==crimes[1], c("dist", "year", "observed", "expected")]
Data_UP$dowry<- df.inla[df.inla$crime==crimes[2], c("dist", "year", "observed", "expected")]

################################################################################
##  1) Fitting models                                                         ##
################################################################################
source("./functions/Mmodel_icar.R")  ## Intrinsic multivariate CAR latent effect
source("./functions/Mmodel_lcar.R")  ## Leroux multivariate CAR latent effect
source("./functions/Mmodel_pcar.R")  ## Proper multivariate CAR latent effect
source("./functions/Mmodel_iid.R")   ## iid

## Fit a spatial-temporal multivariate Poisson mixed model to areal count data, 
## where dependence between spatial/temporal patterns of the diseases is addressed
## through the use of M-models
source("./functions/MCAR_INLA_st.R") 

## The strategy to use for the approximations
strategy <- "gaussian"

########################################
## 1.1) iCAR                          ##
########################################
####################
## 1.1.1) FE      ##
####################
## icar.ad.fe
icar.ad.fe <- MCAR_INLA_st(carto=carto_up, data=Data_UP, ID.area="dist", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic", 
                           prior.disease="FE", prior.interaction=0, strategy=strategy)

## icar.t1.fe
icar.t1.fe <- MCAR_INLA_st(carto=carto_up, data=Data_UP, ID.area="dist", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="FE", prior.interaction=1, strategy=strategy)

## icar.t2.fe
icar.t2.fe <- MCAR_INLA_st(carto=carto_up, data=Data_UP, ID.area="dist", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="FE", prior.interaction=2, strategy=strategy)

## icar.t3.fe
icar.t3.fe <- MCAR_INLA_st(carto=carto_up, data=Data_UP, ID.area="dist", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="FE", prior.interaction=3, strategy=strategy)

## icar.t4.fe
icar.t4.fe <- MCAR_INLA_st(carto=carto_up, data=Data_UP, ID.area="dist", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="FE", prior.interaction=4, strategy=strategy)

## save
resulta.inla.icar.fe <- list(icar.ad.fe=icar.ad.fe, icar.t1.fe=icar.t1.fe,
                             icar.t2.fe=icar.t2.fe, icar.t3.fe=icar.t3.fe,
                             icar.t4.fe=icar.t4.fe)

save(resulta.inla.icar.fe, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.icar.fe"), ".RData" ) )

####################
## 1.1.2) RE      ##
####################
## icar.ad.re
icar.ad.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="RE", prior.interaction=0, strategy=strategy)

## icar.t1.re
icar.t1.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="RE", prior.interaction=1, strategy=strategy)

## icar.t2.re
icar.t2.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="RE", prior.interaction=2, strategy=strategy)

## icar.t3.re
icar.t3.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic",
                           prior.disease="RE", prior.interaction=3, strategy=strategy)

## icar.t4.re
icar.t4.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="intrinsic", 
                           prior.disease="RE", prior.interaction=4, strategy=strategy)

## save
resulta.inla.icar.re <- list(icar.ad.re=icar.ad.re, icar.t1.re=icar.t1.re,
                             icar.t2.re=icar.t2.re, icar.t3.re=icar.t3.re,
                             icar.t4.re=icar.t4.re)

save(resulta.inla.icar.re, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.icar.re"), ".RData" ) )

########################################
## 1.2) LCAR                          ##
########################################
####################
## 1.2.1) FE      ##
####################
## lcar.ad.fe
lcar.ad.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="FE", prior.interaction=0, strategy=strategy)
## lcar.t1.fe
lcar.t1.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="FE", prior.interaction=1, strategy=strategy)
## lcar.t2.fe
lcar.t2.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="FE", prior.interaction=2, strategy=strategy)
## lcar.t3.fe
lcar.t3.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="FE", prior.interaction=3, strategy=strategy)
## lcar.t4.fe
lcar.t4.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="FE", prior.interaction=4, strategy=strategy)

## save
resulta.inla.lcar.fe <- list(lcar.ad.fe=lcar.ad.fe, lcar.t1.fe=lcar.t1.fe,
                             lcar.t2.fe=lcar.t2.fe, lcar.t3.fe=lcar.t3.fe,
                             lcar.t4.fe=lcar.t4.fe)

save(resulta.inla.lcar.fe, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.lcar.fe"), ".RData" ) )

####################
## 1.2.2) RE      ##
####################
## lcar.ad.re
lcar.ad.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="RE", prior.interaction=0,strategy=strategy)
## lcar.t1.re
lcar.t1.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="RE", prior.interaction=1, strategy=strategy)
## lcar.t2.re
lcar.t2.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="RE", prior.interaction=2, strategy=strategy)
## lcar.t3.re
lcar.t3.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="RE", prior.interaction=3, strategy=strategy)
## lcar.t4.re
lcar.t4.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="Leroux",
                           prior.disease="RE", prior.interaction=4, strategy=strategy)

## save
resulta.inla.lcar.re <- list(lcar.ad.re=lcar.ad.re, lcar.t1.re=lcar.t1.re,
                             lcar.t2.re=lcar.t2.re, lcar.t3.re=lcar.t3.re,
                             lcar.t4.re=lcar.t4.re)

save(resulta.inla.lcar.re, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.lcar.re"), ".RData" ) )

########################################
## 1.3) pCAR                          ##
########################################
####################
## 1.3.1) FE      ##
####################
## pcar.ad.fe
pcar.ad.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="FE", prior.interaction=0, strategy=strategy)
## pcar.t1.fe
pcar.t1.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="FE", prior.interaction=1, strategy=strategy)
## pcar.t2.fe
pcar.t2.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="FE", prior.interaction=2, strategy=strategy)
## pcar.t3.fe
pcar.t3.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="FE", prior.interaction=3, strategy=strategy)
## pcar.t4.fe
pcar.t4.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="FE", prior.interaction=4, strategy=strategy)

## save
resulta.inla.pcar.fe <- list(pcar.ad.fe=pcar.ad.fe, pcar.t1.fe=pcar.t1.fe,
                             pcar.t2.fe=pcar.t2.fe, pcar.t3.fe=pcar.t3.fe,
                             pcar.t4.fe=pcar.t4.fe)

save(resulta.inla.pcar.fe, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.pcar.fe"), ".RData" ) )

####################
## 1.3.2) RE      ##
####################
## pcar.ad.re
pcar.ad.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="RE", prior.interaction=0, strategy=strategy)
## pcar.t1.re
pcar.t1.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="RE", prior.interaction=1, strategy=strategy)
## pcar.t2.re
pcar.t2.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="RE", prior.interaction=2, strategy=strategy)
## pcar.t3.re
pcar.t3.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="RE", prior.interaction=3, strategy=strategy)
## pcar.t4.re
pcar.t4.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                           O="observed", E="expected", prior.spatial="proper",
                           prior.disease="RE", prior.interaction=4, strategy=strategy)

## save
resulta.inla.pcar.re <- list(pcar.ad.re=pcar.ad.re, pcar.t1.re=pcar.t1.re,
                             pcar.t2.re=pcar.t2.re, pcar.t3.re=pcar.t3.re,
                             pcar.t4.re=pcar.t4.re)

save(resulta.inla.pcar.re, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.pcar.re"), ".RData" ) )


########################################
## 1.4) BYM                           ##
########################################
####################
## 1.4.1) FE      ##
####################
## bym.ad.fe
bym.ad.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="FE", prior.interaction=0, strategy=strategy)
## bym.t1.fe
bym.t1.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="FE", prior.interaction=1, strategy=strategy)
## bym.t2.fe
bym.t2.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="FE", prior.interaction=2, strategy=strategy)
## bym.t3.fe
bym.t3.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="FE", prior.interaction=3, strategy=strategy)
## bym.t4.fe
bym.t4.fe <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="FE", prior.interaction=4, strategy=strategy)
## save
resulta.inla.bym.fe <- list(bym.ad.fe=bym.ad.fe, bym.t1.fe=bym.t1.fe,
                            bym.t2.fe=bym.t2.fe, bym.t3.fe=bym.t3.fe,
                            bym.t4.fe=bym.t4.fe)

save(resulta.inla.bym.fe, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.bym.fe"), ".RData" ) )

####################
## 1.4.2) RE      ##
####################
## bym.ad.re
bym.ad.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="RE", prior.interaction=0, strategy=strategy)
## bym.t1.re
bym.t1.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="RE", prior.interaction=1, strategy=strategy)
## bym.t2.re
bym.t2.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="RE", prior.interaction=2, strategy=strategy)
## bym.t3.re
bym.t3.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="RE", prior.interaction=3, strategy=strategy)
## bym.t4.re
bym.t4.re <- MCAR_INLA_st(carto=Carto_UP, data=Data_UP, ID.area="district", ID.year="year",
                          O="observed", E="expected", prior.spatial="BYM",
                          prior.disease="RE", prior.interaction=4, strategy=strategy)

## save
resulta.inla.bym.re <- list(bym.ad.re=bym.ad.re, bym.t1.re=bym.t1.re,
                            bym.t2.re=bym.t2.re, bym.t3.re=bym.t3.re,
                            bym.t4.re=bym.t4.re)

save(resulta.inla.bym.re, file=paste0("./resul/", gsub("\\.", "_", "resulta.inla.bym.re"), ".RData" ) )

################################################################################
################################################################################