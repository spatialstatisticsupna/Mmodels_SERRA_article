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
library(pbugs) # For running the models in parallel calls to WinBUGS

## Folder to save results
if(!file.exists("resul")) {dir.create("resul")}

## Load data and Uttar Pradesh SpatialPolygonsDataFrame
load("./dataMmodel.RData")

################################################################################
## Data organization for WinBUGS                                              ## 
################################################################################
datos<- data
carto<- carto_up

## crimes
crimes<- c("rape", "dowry")

## Number of areas and number of time periods
n<- length(unique(datos$ID_area))
t<- length(unique(datos$ID_year))

## array with observed cases
Obs<- datos[,c(crimes,"ID_area","ID_year")]   
aux<- Obs[,c(crimes, "ID_year")]
O<- abind::abind(split(aux[,1:length(crimes)], aux[,dim(aux)[2]]),along = 3)
rownames(O)<-NULL
rm(list=c("aux", "Obs")) 

## array with expected cases
Esp<-datos[,c(paste0("e_",crimes),"ID_area","ID_year")]
aux<- Esp[,c(paste0("e_",crimes), "ID_year")]
E<- abind::abind(split(aux[,1:length(crimes)], aux[,dim(aux)[2]]),along = 3)
rownames(E)<-NULL
rm(list=c("aux", "Esp")) 

## Define number of: areas, crimes and time periods
Nareas<-dim(O)[1]
Ndiseases<-dim(O)[2]
Nyears<-dim(O)[3]

## Neighborhood structure (spatial)
nb<-spdep::poly2nb(carto)
## for lcar
off<- c(0,sapply(nb,length))
for(i in 2:length(off)){off[i]<- off[i]+off[i-1]}

## Neighborhood structure (temporal)
nbt<-list()
nbt[[1]]<-c(2)
nbt[[Nyears]]<-c(Nyears-1)
for(l in 2:(Nyears-1)){
  nbt[[l]]<-c(l-1,l+1)
}
rm(l) 

################################################################################
## load functions                                                             ##
################################################################################
source("functions/nva_nva_lcar_spaciotemporal_functions.R") # functions (models)

## bugs.directory
bugs.dir<- c("") # Set an appropiate directory

## n.chains, n.iter, n.burnin
num.chains<- 3
num.iter<- 30000
num.burnin<- 5000

################################################################################
## Run models                                                                 ##
## RE M-models - LCAR                                                         ##
################################################################################

########################################
## Additive                           ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, 
           nnv=length(unlist(nb)), adj=unlist(nb), num=sapply(nb,length), off=off, 
           adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length) )

## NVA.NVA.ad.lcar
initials<-function(){list(mu=rnorm(Ndiseases,0,0.1), 
                          sdstruct=runif(1,0,1), 
                          Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), 
                          gamma1=runif(Ndiseases,0.1,0.9),
                          sdstructg=runif(1,0,1), 
                          Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears) ) }

param<-c("smr.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", "Theta", 
         "M", "sdstruct", "prec",  "Sigma.s", "Corre.s", "Gam", "Mg", "sdstructg", 
         "precg", "Sigma.t", "Corre.t", "Espat", "Etemp", "gamma1")

t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = NVA.NVA.ad.lcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

resulta.nva.nva.lcar.ad<- result 
t.resulta.nva.nva.lcar.ad<- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type I                             ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, 
           nnv=length(unlist(nb)), adj=unlist(nb), num=sapply(nb,length), off=off, 
           adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length) )

## NVA.NVA.t1.lcar
initials<-function(){list(mu=rnorm(Ndiseases,0,0.1), 
                          sdstruct=runif(1,0,1), 
                          Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), 
                          gamma1=runif(Ndiseases,0.1,0.9),
                          sdstructg=runif(1,0,1), 
                          Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                          sdZet=runif(Ndiseases,0,1), 
                          Zet.aux=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) ) }

param<-c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda", 
         "mu", "Theta", "M", "sdstruct", "prec",  "Sigma.s", "Corre.s", "Gam", 
         "Mg", "sdstructg", "precg", "Sigma.t", "Corre.t", "Espat", "Etemp", 
         "Eint", "Zet", "sdZet", "gamma1")

t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = NVA.NVA.t1.lcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

resulta.nva.nva.lcar.t1<- result 
t.resulta.nva.nva.lcar.t1<- t.result
## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type II                            ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E,
           nnv=length(unlist(nb)), adj=unlist(nb), num=sapply(nb,length), off=off, 
           adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), 
           adj.zt=unlist(nbt), weights.zt=rep(1,length(unlist(nbt))), num.zt=sapply(nbt,length))

## NVA.NVA.t2.lcar
initials<-function(){list(mu=rnorm(Ndiseases,0,0.1), 
                          sdstruct=runif(1,0,1), 
                          Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), 
                          gamma1=runif(Ndiseases,0.1,0.9),
                          sdstructg=runif(1,0,1), 
                          Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                          sdZet=runif(Ndiseases,0,1), 
                          Temporal.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) ) }

param<-c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", 
         "Theta", "M", "sdstruct", "prec",  "Sigma.s", "Corre.s", "Gam", "Mg", 
         "sdstructg", "precg", "Sigma.t", "Corre.t", "Espat", "Etemp", "Eint", 
         "Zet", "sdZet", "gamma1")

t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = NVA.NVA.t2.lcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

resulta.nva.nva.lcar.t2<- result 
t.resulta.nva.nva.lcar.t2<- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type III                           ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E,
           nnv=length(unlist(nb)), adj=unlist(nb), num=sapply(nb,length), off=off, 
           adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), 
           adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length))

## NVA.NVA.t3.lcar
initials<-function(){list(mu=rnorm(Ndiseases,0,0.1), 
                          sdstruct=runif(1,0,1), 
                          Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas),
                          gamma1=runif(Ndiseases,0.1,0.9),
                          sdstructg=runif(1,0,1), 
                          Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                          sdZet=runif(Ndiseases,0,1), 
                          Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) ) }

param<-c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", 
         "Theta", "M", "sdstruct", "prec",  "Sigma.s", "Corre.s", "Gam", "Mg", 
         "sdstructg", "precg", "Sigma.t", "Corre.t", "Espat", "Etemp", "Eint", "Zet",
         "sdZet", "gamma1")

t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = NVA.NVA.t3.lcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

resulta.nva.nva.lcar.t3<- result 
t.resulta.nva.nva.lcar.t3<- t.result

## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type IV                            ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, 
           nnv=length(unlist(nb)), adj=unlist(nb), num=sapply(nb,length), off=off, 
           adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), 
           adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length), 
           Mzz=Mzz)

## NVA.NVA.t4.lcar
initials<-function(){list(mu=rnorm(Ndiseases,0,0.1), 
                          sdstruct=runif(1,0,1), 
                          Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), 
                          gamma1=runif(Ndiseases,0.1,0.9),
                          sdstructg=runif(1,0,1), 
                          Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), 
                          sdZet=runif(Ndiseases,0,1),
                          Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) ) }

param<-c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR","lambda","mu", 
         "Theta","M","sdstruct", "prec",  "Sigma.s", "Corre.s", "Gam","Mg", 
         "sdstructg", "precg", "Sigma.t", "Corre.t", "Espat","Etemp", "Eint", 
         "Zet", "sdZet", "gamma1")

t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, 
                                        parameters.to.save = param, model.file = NVA.NVA.t4.lcar,
                                        n.chains = num.chains, n.iter = num.iter, 
                                        n.burnin = num.burnin, bugs.directory = bugs.dir, 
                                        DIC = FALSE))

resulta.nva.nva.lcar.t4<- result 
t.resulta.nva.nva.lcar.t4<- t.result
## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## save                               ##
########################################
resulta.winbugs.lcar.re<- list(lcar.ad.re=resulta.nva.nva.lcar.ad,
                               lcar.t1.re=resulta.nva.nva.lcar.t1,
                               lcar.t2.re=resulta.nva.nva.lcar.t2,
                               lcar.t3.re=resulta.nva.nva.lcar.t3,
                               lcar.t4.re=resulta.nva.nva.lcar.t4)

t.resulta.winbugs.lcar.re<- list(lcar.ad.re=t.resulta.nva.nva.lcar.ad,
                                 lcar.t1.re=t.resulta.nva.nva.lcar.t1,
                                 lcar.t2.re=t.resulta.nva.nva.lcar.t2,
                                 lcar.t3.re=t.resulta.nva.nva.lcar.t3,
                                 lcar.t4.re=t.resulta.nva.nva.lcar.t4)

save(resulta.winbugs.lcar.re, t.resulta.winbugs.lcar.re, 
     file =paste0("resul/",gsub("\\.", "_", "resulta.winbugs.lcar.re"),".RData"))

################################################################################
################################################################################