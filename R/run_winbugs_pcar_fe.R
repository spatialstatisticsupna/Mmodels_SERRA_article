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
source("functions/WINBUGS_fe_fe_pcar_ST.R")

## bugs.directory
bugs.dir <- c("C:/WinBUGS14") # Set an appropriate directory

## Set 'n.chains', 'n.iter' and 'n.burnin' parameters
num.chains <- 3
num.iter <- 30000
num.burnin <- 5000

################################################################################
## Run models                                                                 ##
## FE M-models - pCAR                                                         ##
################################################################################

########################################
## Additive                           ##
########################################
data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, 
             adj=unlist(nb), num=sapply(nb,length), M.car.s= 1/sapply(nb,length), 
             C.car.s= rep(1/sapply(nb,length), sapply(nb,length)), 
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length))

## FE.FE.ad.pcar
initials <- function(){list(mu=rnorm(Ndiseases,0,0.1),
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears) )}

param <- c("smr.prob", "espat.prob", "etemp.prob", "gamma","SMR", "lambda", "mu",
           "Theta", "M", "Sigma.s", "Corre.s", "Gam", "Mg", "Sigma.t", "Corre.t",
           "Espat", "Etemp")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = FE.FE.ad.pcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

results.fe.fe.pcar.ad <- result 
t.results.fe.fe.pcar.ad <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type I                             ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, 
           adj=unlist(nb), num=sapply(nb,length), M.car.s= 1/sapply(nb,length), 
           C.car.s= rep(1/sapply(nb,length), sapply(nb,length)), 
           adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length))

## FE.FE.t1.pcar
initials <- function(){list(mu=rnorm(Ndiseases,0,0.1), 
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), 
                            M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                            sdZet=runif(Ndiseases,0,1), 
                            Zet.aux=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) )}

param <- c("smr.prob", "espat.prob", "etemp.prob", "eint.prob","gamma","SMR", "lambda", 
           "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam", "Mg", "Sigma.t", "Corre.t", 
           "Espat", "Etemp", "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = FE.FE.t1.pcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

results.fe.fe.pcar.t1 <- result 
t.results.fe.fe.pcar.t1 <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Type II                            ##
########################################
data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, 
             adj=unlist(nb), num=sapply(nb,length), M.car.s= 1/sapply(nb,length),
             C.car.s= rep(1/sapply(nb,length), sapply(nb,length)), 
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length),
             adj.zt=unlist(nbt), weights.zt=rep(1,length(unlist(nbt))), num.zt=sapply(nbt,length))

initials <- function(){list(mu=rnorm(Ndiseases,0,0.1), 
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                            M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                            sdZet=runif(Ndiseases,0,1), 
                            Temporal.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) )}

param <- c("smr.prob", "espat.prob", "etemp.prob", "eint.prob","gamma","SMR", "lambda", 
           "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam", "Mg", "Sigma.t", "Corre.t",
           "Espat", "Etemp", "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = FE.FE.t2.pcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, DIC = FALSE))

results.fe.fe.pcar.t2 <- result 
t.results.fe.fe.pcar.t2 <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Type III                           ##
########################################
data <- list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, 
             adj=unlist(nb), num=sapply(nb,length), M.car.s= 1/sapply(nb,length),
             C.car.s= rep(1/sapply(nb,length), sapply(nb,length)), 
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), 
             adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length))

initials <- function(){list(mu=rnorm(Ndiseases,0,0.1), 
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), 
                            M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases),
                            Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                            sdZet=runif(Ndiseases,0,1), 
                            Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) )}

param <- c("smr.prob", "espat.prob", "etemp.prob", "eint.prob","gamma","SMR", "lambda", 
           "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam","Mg", "Sigma.t", "Corre.t", 
           "Espat","Etemp", "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = FE.FE.t3.pcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

results.fe.fe.pcar.t3 <- result 
t.results.fe.fe.pcar.t3 <- t.result

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
             adj=unlist(nb), num=sapply(nb,length), M.car.s= 1/sapply(nb,length), 
             C.car.s= rep(1/sapply(nb,length), sapply(nb,length)), 
             adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), 
             adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length), 
             Mzz=Mzz)
  
initials <- function(){list(mu=rnorm(Ndiseases,0,0.1), 
                            Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), 
                            M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), 
                            Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                            sdZet=runif(Ndiseases,0,1), 
                            Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) )}

param <- c("smr.prob", "espat.prob", "etemp.prob", "eint.prob","gamma","SMR", "lambda",
           "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam","Mg", "Sigma.t", "Corre.t",
           "Espat","Etemp", "Eint", "Zet", "sdZet")

t.result <- system.time(result <- pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = FE.FE.t4.pcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

results.fe.fe.pcar.t4 <- result
t.results.fe.fe.pcar.t4 <- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")


########################################
## Save results                       ##
########################################
results.winbugs.pcar.fe <- list(pcar.ad.fe=results.fe.fe.pcar.ad,
                                pcar.t1.fe=results.fe.fe.pcar.t1,
                                pcar.t2.fe=results.fe.fe.pcar.t2,
                                pcar.t3.fe=results.fe.fe.pcar.t3,
                                pcar.t4.fe=results.fe.fe.pcar.t4)

t.results.winbugs.pcar.fe <- list(pcar.ad.fe=t.results.fe.fe.pcar.ad,
                                  pcar.t1.fe=t.results.fe.fe.pcar.t1,
                                  pcar.t2.fe=t.results.fe.fe.pcar.t2,
                                  pcar.t3.fe=t.results.fe.fe.pcar.t3,
                                  pcar.t4.fe=t.results.fe.fe.pcar.t4)
 
save(results.winbugs.pcar.fe, t.results.winbugs.pcar.fe, 
     file =paste0("results/",gsub("\\.", "_", "results.winbugs.pcar.fe"),".RData"))
