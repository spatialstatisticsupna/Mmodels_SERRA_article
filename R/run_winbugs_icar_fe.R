################################################################################
###############        Spatio-temporal M-models (WinBUGS)        ###############
################################################################################
rm(list=ls())

## Working directory
DirMain<-""     # Set an appropiate directory
setwd(DirMain)

## libraries
library(spdep); library(INLA)
library(pbugs) # For running the models in parallel calls to WinBUGS

### save results
tdir <- paste(getwd(), "/resul", sep = "", collapse = "")
if(!file.exists(tdir)) {dir.create(tdir)}
rm(tdir)

## Load data and Uttar Pradesh SpatialPolygonsDataFrame
load("dataMmodel.RData")

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
source("functions_WinBUGS/fe_fe_icar_spaciotemporal_functions.R") # functions (models)

## bugs.directory
bugs.dir<- c("") # Set an appropiate directory

## n.chains, n.iter, n.burnin
num.chains<- 3
num.iter<- 30000
num.burnin<- 5000

################################################################################
## Run models                                                                 ##
## FE M-models - iCAR                                                         ##
################################################################################

########################################
## Additive                           ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length), adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length))
initials <- function(){ list(mu=rnorm(Ndiseases,0,0.1), Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears) ) }
param <- c("smr.prob","espat.prob", "etemp.prob", "SMR", "lambda", "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam", "Mg", "Sigma.t", "Corre.t", "Espat", "Etemp")
t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, parameters.to.save = param, model.file = FE.FE.ad.icar,
                                        n.chains = num.chains, n.iter = num.iter, n.burnin = num.burnin, bugs.directory = bugs.dir, DIC = FALSE))
resulta.fe.fe.icar.ad<- result 
t.resulta.fe.fe.icar.ad<- t.result
## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type I                             ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length), adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length))
initials <- function(){ list(mu=rnorm(Ndiseases,0,0.1), Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), sdZet=runif(Ndiseases,0,1), Zet.aux=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) ) }
param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam", "Mg", "Sigma.t", "Corre.t", "Espat", "Etemp", "Eint", "Zet", "sdZet")
t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, parameters.to.save = param, model.file = FE.FE.t1.icar,
                                        n.chains = num.chains, n.iter = num.iter, n.burnin = num.burnin, bugs.directory = bugs.dir, DIC = FALSE))
resulta.fe.fe.icar.t1<- result 
t.resulta.fe.fe.icar.t1<- t.result
## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type II                            ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length), adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), adj.zt=unlist(nbt), weights.zt=rep(1,length(unlist(nbt))), num.zt=sapply(nbt,length))
initials <- function(){ list(mu=rnorm(Ndiseases,0,0.1), Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), sdZet=runif(Ndiseases,0,1), Temporal.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nareas,Ndiseases,Nyears)) )}
param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam", "Mg", "Sigma.t", "Corre.t", "Espat", "Etemp", "Eint", "Zet", "sdZet")
t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, parameters.to.save = param, model.file = FE.FE.t2.icar,
                                        n.chains = num.chains, n.iter = num.iter, n.burnin = num.burnin, bugs.directory = bugs.dir, DIC = FALSE))
resulta.fe.fe.icar.t2<- result 
t.resulta.fe.fe.icar.t2<- t.result
## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type III                           ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length), adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length))
initials <- function(){list(mu=rnorm(Ndiseases,0,0.1), Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), sdZet=runif(Ndiseases,0,1), Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) ) }
param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam","Mg", "Sigma.t", "Corre.t", "Espat","Etemp", "Eint", "Zet", "sdZet")
t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, parameters.to.save = param, model.file = FE.FE.t3.icar,
                                        n.chains = num.chains, n.iter = num.iter, n.burnin = num.burnin, bugs.directory = bugs.dir, DIC = FALSE))
resulta.fe.fe.icar.t3<- result 
t.resulta.fe.fe.icar.t3<- t.result
## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## Type IV                            ##
########################################
data<-list(Nyears=Nyears, Ndiseases=Ndiseases, Nareas=Nareas, O=O, E=E, adj=unlist(nb), weights=rep(1,length(unlist(nb))), num=sapply(nb,length), adjt=unlist(nbt), weightst=rep(1,length(unlist(nbt))), numt=sapply(nbt,length), adj.zs=unlist(nb), weights.zs=rep(1,length(unlist(nb))), num.zs=sapply(nb,length), Mzz=Mzz)
initials <- function(){list(mu=rnorm(Ndiseases,0,0.1), Spatial=matrix(rnorm(Nareas*Ndiseases), nrow=Ndiseases, ncol=Nareas), M=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Mg=matrix(rnorm(Ndiseases), nrow=Ndiseases, ncol=Ndiseases), Temporal=matrix(rnorm(Ndiseases*Nyears), nrow=Ndiseases, ncol=Nyears), sdZet=runif(Ndiseases,0,1), Spatial.z=array(rnorm(Nareas*Ndiseases*Nyears),dim=c(Nyears,Ndiseases,Nareas)) )}
param <- c("smr.prob", "eint.prob", "espat.prob", "etemp.prob", "SMR", "lambda", "mu", "Theta", "M", "Sigma.s", "Corre.s", "Gam","Mg", "Sigma.t", "Corre.t", "Espat","Etemp", "Eint", "Zet", "sdZet")
t.result <- system.time(result <- Pbugs(program="winbugs", data = data, inits = initials, parameters.to.save = param, model.file = FE.FE.t4.icar,
                                        n.chains = num.chains, n.iter = num.iter, n.burnin = num.burnin, bugs.directory = bugs.dir, DIC = FALSE))
resulta.fe.fe.icar.t4<- result 
t.resulta.fe.fe.icar.t4<- t.result
## rm
rm(list = c("cl", "initials", "param", "result", "t.result"))
rm("data")

########################################
## save                               ##
########################################
resulta.winbugs.icar.fe<- list(icar.ad.fe=resulta.fe.fe.icar.ad,
                               icar.t1.fe=resulta.fe.fe.icar.t1,
                               icar.t2.fe=resulta.fe.fe.icar.t2,
                               icar.t3.fe=resulta.fe.fe.icar.t3,
                               icar.t4.fe=resulta.fe.fe.icar.t4)

t.resulta.winbugs.icar.fe<- list(icar.ad.fe=t.resulta.fe.fe.icar.ad,
                                 icar.t1.fe=t.resulta.fe.fe.icar.t1,
                                 icar.t2.fe=t.resulta.fe.fe.icar.t2,
                                 icar.t3.fe=t.resulta.fe.fe.icar.t3,
                                 icar.t4.fe=t.resulta.fe.fe.icar.t4)

save(resulta.winbugs.icar.fe, t.resulta.winbugs.icar.fe, file= paste0("resul/",gsub("\\.", "_", "resulta.winbugs.icar.fe"),".RData"))
################################################################################
################################################################################