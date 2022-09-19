################################################################################
## Title: Bayesian inference in multivariate spatio-temporal areal models     ##
##        using INLA: analysis of gender-based violence in small areas        ##
##                                                                            ##
## Authors: Vicente, G. - Goicoa, T. -  Ugarte, M.D.                          ##
##                                                                            ##
## https://doi.org/10.1007/s00477-020-01808-x                                 ##
##                                                                            ##
################################################################################
##                    Spatio-temporal M-models (WinBUGS)                      ##
################################################################################
rm(list=ls())

## libraries
library(spdep); library(INLA); library(abind)

# For running the models in parallel calls to WinBUGS
# devtools::install_github("fisabio/pbugs")
library(pbugs) 

## Folder to save results
if(!file.exists("results")) {dir.create("results")}

## Load data and Uttar Pradesh SpatialPolygonsDataFrame
load("./dataMmodel.RData")

################################################################################
## Data organization for WinBUGS                                              ##
################################################################################

## Define number of: areas, crimes and time periods
Ndiseases <- length(data_UP)
Nareas <- length(unique(data_UP$dowry$dist))
Nyears <- length(unique(data_UP$dowry$year))

## Array with observed cases
O <- array(NA, dim=c(Nareas,Ndiseases,Nyears),
           dimnames=list(unique(data_UP$dowry$dist),names(data_UP),unique(data_UP$dowry$year)))

O[,"rape",] <- matrix(data_UP$rape$obs,Nareas,Nyears)
O[,"dowry",] <- matrix(data_UP$dowry$obs,Nareas,Nyears)


## Array with observed cases
E <- array(NA, dim=c(Nareas,Ndiseases,Nyears),
           dimnames=list(unique(data_UP$dowry$dist),names(data_UP),unique(data_UP$dowry$year)))

E[,"rape",] <- matrix(data_UP$rape$exp,Nareas,Nyears)
E[,"dowry",] <- matrix(data_UP$dowry$exp,Nareas,Nyears)


## Spatial neighbourhood structure
nb <- spdep::poly2nb(carto_UP)

## for LCAR model
off <- c(0,sapply(nb,length))
for(i in 2:length(off)){
  off[i]<- off[i]+off[i-1]
}

## Temporal neighbourhood structure
nbt <- list()
nbt[[1]] <- c(2)
nbt[[Nyears]] <- c(Nyears-1)
for(l in 2:(Nyears-1)){
  nbt[[l]]<-c(l-1,l+1)
}

################################################################################
## load functions                                                             ##
################################################################################
source("functions/WINBUGS_nva_nva_bym_ST.R")

## bugs.directory
bugs.dir <- c("C:/WinBUGS14") # Set an appropriate directory

## Set 'n.chains', 'n.iter' and 'n.burnin' parameters
num.chains <- 3
num.iter <- 30000
num.burnin <- 5000

################################################################################
## Run models                                                                 ##
## RE M-models - BYM                                                          ##
################################################################################

########################################
## Additive                           ##
########################################
data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E,
             adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length),
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length))

initials <- function(){list(mu=rnorm(Ndiseases,0,0.1),
                            sdstruct.het=runif(1,0,1),
                            Het=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstruct.sp=runif(1,0,1),  
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstructg=runif(1,0,1), 
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears) ) }

param <- c("smr.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", "c.Het",
           "Spatial", "M", "tPhi","Theta","sdstruct.sp","sdstruct.het", "Sigma.s",
           "Corre.s", "Gam", "Mg", "sdstructg", "Sigma.t", "Corre.t", "Espat", "Etemp")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials,
                                        parameters.to.save = param, model.file = NVA.NVA.ad.bym,
                                        n.chains = num.chains, n.iter = num.iter,
                                        n.burnin = num.burnin, bugs.directory = bugs.dir,
                                        DIC = FALSE))

results.nva.nva.bym.ad <- result 
t.results.nva.nva.bym.ad <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Type I                             ##
########################################
data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E,
             adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length),
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length))

initials <- function(){list(mu=rnorm(Ndiseases,0,0.1),
                            sdstruct.het=runif(1,0,1),
                            Het=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstruct.sp=runif(1,0,1), 
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstructg=runif(1,0,1),
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears),
                            sdZet=runif(Ndiseases,0,1),
                            Zet.aux=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) ) }

param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda",
           "mu", "c.Het", "Spatial", "M", "tPhi","Theta","sdstruct.sp","sdstruct.het",
           "Sigma.s", "Corre.s", "Gam", "Mg", "sdstructg", "Sigma.t", "Corre.t", "Espat",
           "Etemp", "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials,
                                        parameters.to.save = param, model.file = NVA.NVA.t1.bym,
                                        n.chains = num.chains, n.iter = num.iter,
                                        n.burnin = num.burnin, bugs.directory = bugs.dir,
                                        DIC = FALSE))

results.nva.nva.bym.t1 <- result 
t.results.nva.nva.bym.t1 <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Type II                            ##
########################################
data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E,
             adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length),
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length),
             adj.zt=unlist(nbt), weights.zt=rep(1,length(unlist(nbt))), num.zt=sapply(nbt,length))

initials <- function(){list(mu=rnorm(Ndiseases,0,0.1),
                            sdstruct.het =runif(1,0,1),
                            Het=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstruct.sp  =runif(1,0,1),
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstructg=runif(1,0,1),
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears),
                            sdZet=runif(Ndiseases,0,1),
                            Temporal.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) ) }

param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda",
           "mu", "c.Het", "Spatial", "M", "tPhi","Theta","sdstruct.sp","sdstruct.het",
           "Sigma.s", "Corre.s", "Gam", "Mg", "sdstructg", "Sigma.t", "Corre.t", "Espat",
           "Etemp", "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials,
                                        parameters.to.save = param, model.file = NVA.NVA.t2.bym,
                                        n.chains = num.chains, n.iter = num.iter,
                                        n.burnin = num.burnin, bugs.directory = bugs.dir,
                                        DIC = FALSE))
results.nva.nva.bym.t2 <- result 
t.results.nva.nva.bym.t2 <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Type III                           ##
########################################
data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E,
             adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length),
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length),
             adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length))

initials <- function(){list(mu=rnorm(Ndiseases,0,0.1),
                            sdstruct.het =runif(1,0,1),
                            Het=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstruct.sp =runif(1,0,1),
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstructg=runif(1,0,1),
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears),
                            sdZet=runif(Ndiseases,0,1),
                            Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) )}

param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR","lambda", "mu",
           "c.Het","Spatial","M", "tPhi","Theta","sdstruct.sp", "sdstruct.het","Sigma.s",
           "Corre.s", "Gam","Mg","sdstructg", "Sigma.t", "Corre.t", "Espat","Etemp",
           "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = NVA.NVA.t3.bym,
                                        n.chains = num.chains, n.iter = num.iter,
                                        n.burnin = num.burnin, bugs.directory = bugs.dir,
                                        DIC = FALSE))

results.nva.nva.bym.t3 <- result 
t.results.nva.nva.bym.t3 <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Type IV                            ##
########################################
D <- diff(diag(Nyears), differences=1)
Rt <- t(D)%*%D
Rt.ginv <- MASS::ginv(Rt)
Mzz <- chol(Rt.ginv)

data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E,
             adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length),
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length),
             adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length),
             Mzz=Mzz)

initials <- function(){list(mu=rnorm(Ndiseases,0,0.1),
                            sdstruct.het =runif(1,0,1),
                            Het=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstruct.sp =runif(1,0,1),
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            sdstructg=runif(1,0,1),
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears),
                            sdZet=runif(Ndiseases,0,1),
                            Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) )}

param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR","lambda","mu",
           "c.Het","Spatial","M", "tPhi","Theta","sdstruct.sp","sdstruct.het", "Sigma.s",
           "Corre.s", "Gam","Mg","sdstructg", "Sigma.t", "Corre.t", "Espat","Etemp",
           "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials,
                                        parameters.to.save = param, model.file = NVA.NVA.t4.bym,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir,
                                        DIC = FALSE))

results.nva.nva.bym.t4 <- result 
t.results.nva.nva.bym.t4 <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Save results                       ##
########################################
results.winbugs.bym.re <- list(bym.ad.re=results.nva.nva.bym.ad,
                               bym.t1.re=results.nva.nva.bym.t1,
                               bym.t2.re=results.nva.nva.bym.t2,
                               bym.t3.re=results.nva.nva.bym.t3,
                               bym.t4.re=results.nva.nva.bym.t4)

t.results.winbugs.bym.re <- list(bym.ad.re=t.results.nva.nva.bym.ad,
                                 bym.t1.re=t.results.nva.nva.bym.t1,
                                 bym.t2.re=t.results.nva.nva.bym.t2,
                                 bym.t3.re=t.results.nva.nva.bym.t3,
                                 bym.t4.re=t.results.nva.nva.bym.t4)

save(results.winbugs.bym.re, t.results.winbugs.bym.re,
     file=paste0("results/",gsub("\\.", "_", "results.winbugs.bym.re"),".RData"))
