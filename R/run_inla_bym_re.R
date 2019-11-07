################################################################################
###############          Spatio-temporal M-models (INLA)         ###############
################################################################################
rm(list=ls())

## Working directory
DirMain<- ""     # Set an appropiate directory
setwd(DirMain)

## libraries
library(INLA); library(spdep); library(spData); library(rgdal)

### save results
tdir <- paste(getwd(), "/resul", sep = "", collapse = "")
if(!file.exists(tdir)) {dir.create(tdir)}
rm(tdir)

## Load data and Uttar Pradesh SpatialPolygonsDataFrame
load("dataMmodel.RData")

################################################################################
## Data organization for INLA                                                 ##
################################################################################
datos<- data
carto<- carto_up

## crimes
crimes<- c("rape", "dowry")
e_crimes<- paste0("e_", crimes)

## Number of areas and number of time periods
n<- length(unique(datos$ID_area))
t<- length(unique(datos$ID_year))

## data.frame INLA
n.rep <- 1
datos<- datos[order(datos$ID_year, datos$ID_area),]
d1<- as.list(rep(NA,length(crimes)))
for(i in 1:length(crimes)){ d1[[i]]<- list()}
for(i in 1:length(crimes)){
  for(j in 1:t){
    d1[[i]][[j]]<- datos[datos$ID_year==j, c(crimes[i], e_crimes[i], "ID_area", "ID_year" )]
    d1[[i]][[j]]$ID_disease<- i
    colnames(d1[[i]][[j]])<- c("obs", "esp", "ID_area", "ID_year", "ID_disease")
  }
}

d.data_frame<- NULL
for(j in 1:t){d.data_frame<- rbind(d.data_frame, d1[[1]][[j]], d1[[2]][[j]])}   

## Intercept for each crime
intercepts<- paste0("I",1:length(crimes))
d.data_frame[intercepts]<- NA
for(i in 1:length(crimes)){
  d.data_frame[d.data_frame$ID_disease==i, intercepts[i]]<- 1
}

## Parameters of the Mmodel
k <- length(crimes)
alpha.min <- 0
alpha.max <- 1  

## Define idx, idx.u, idx.v  (different ID_area for different crimes)
d.data_frame$idx<- (d.data_frame$ID_disease-1)*n + d.data_frame$ID_area
d.data_frame$idx.u<- d.data_frame$idx
d.data_frame$idx.v<- d.data_frame$idx

## Define idy  (different ID_year for different crimes)
d.data_frame$idy<- (d.data_frame$ID_disease-1) *t +  d.data_frame$ID_year

## Define idxy.j (different ID_area_year for different crimes)
d.data_frame$idxy<- d.data_frame$ID_area + (d.data_frame$ID_year-1)*n
## idxy.
idxy.<- paste0("idxy.", 1:k)
d.data_frame[idxy.]<- NA
for(i in 1:k){
  d.data_frame[d.data_frame$ID_disease==i, idxy.[i]]<- (d.data_frame[d.data_frame$ID_disease==i, c("ID_year")]-1)*n+d.data_frame[d.data_frame$ID_disease==i, c("ID_area")]
}

## data list for INLA
d <- list(OBS = d.data_frame$obs, EXP = d.data_frame$esp, 
          ID_area = d.data_frame$ID_area, ID_year=d.data_frame$ID_year, ID_disease=d.data_frame$ID_disease,
          idx= d.data_frame$idx, idx.u= d.data_frame$idx.u, idx.v=d.data_frame$idx.v,
          idy= d.data_frame$idy,
          idxy=d.data_frame$idxy,
          idxy.1=d.data_frame$idxy.1, idxy.2=d.data_frame$idxy.2,
          I1=d.data_frame$I1, I2=d.data_frame$I2)


## W: adjacency matrix (spatial)
adj<-spdep::poly2nb(carto)
W <- as(nb2mat(adj, style = "B"), "Matrix")

## Spatial neighborhood matrix (Q_{xi})
spdep::nb2INLA("carto_nb.graph", spdep::poly2nb(carto))
g <- INLA::inla.read.graph("carto_nb.graph")
Q_xi <- matrix(0, g$n, g$n)
for (i in 1:g$n){
  Q_xi[i,i]=g$nnbs[[i]]
  Q_xi[i,g$nbs[[i]]]=-1
}

## W.t: adjacency matrix (temporal)
adj.t<-list(); adj.t[[1]]<-c(2); adj.t[[t]]<-c(t-1)
for(l in 2:(t-1)){adj.t[[l]]<-c(l-1,l+1)}
for(x in 1:length(adj.t)){adj.t[[x]]<-as.integer(adj.t[[x]])}
names(adj.t)<- 1:length(adj.t)
b=list(class="nb", call="nuestro_nb", cell=TRUE, rook=TRUE, sym=TRUE)
attributes(adj.t)<-b
W.t <- as(nb2mat(adj.t, style = "B"), "Matrix")

## Temporal structure matrix for a RW1 prior
D1 <- diff(diag(t), differences=1)
Q_gammaRW1 <- t(D1)%*%D1

## Define appropriate hyperprior distributions
sdunif="expression:
logdens=-log_precision/2;
return(logdens)"

lunif = "expression:
a = 1;
b = 1;
beta = exp(theta)/(1+exp(theta));
logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
log_jacobian = log(beta*(1-beta));
return(logdens+log_jacobian)"

################################################################################
## load functions                                                             ##
################################################################################
source("functions_INLA/inla_rgeneric_Mmodel_model_icar_normal.R") # icar
model.s.u <- inla.rgeneric.define(inla.rgeneric.Mmodel.model.icar.n, debug = TRUE, k = k, W = W, alpha.min = alpha.min, alpha.max = alpha.max)
model.t <- inla.rgeneric.define(inla.rgeneric.Mmodel.model.icar.n, debug = TRUE, k = k, W = W.t, alpha.min = alpha.min, alpha.max = alpha.max)

source("functions_INLA/inla_rgeneric_Mmodel_model_iid_normal.R") # iid
model.s.v <- inla.rgeneric.define(inla.rgeneric.Mmodel.model.iid.n, debug = TRUE, k = k, W = W, alpha.min = alpha.min, alpha.max = alpha.max)

################################################################################
## Run models                                                                 ##
## RE M-models - BYM                                                          ##
################################################################################
## constraints
A_constr.s<- kronecker(diag(k), matrix(1,1,n))
A_constr.t<- kronecker(diag(k), matrix(1,1,t))

R_1_2 <- kronecker(Q_gammaRW1,diag(n))
r_def_1_2 <- n
A_constr_1_2<- kronecker(matrix(1,1,t), diag(n))

R_1_3 <- kronecker(diag(t),Q_xi)
r_def_1_3 <- t
A_constr_1_3<- kronecker(diag(t),matrix(1,1,n))

R_1_4 <- kronecker(Q_gammaRW1,Q_xi)
r_def_1_4 <- n+t-1
A.1.1 <- kronecker(matrix(1,1,t),diag(n))
A.1.2 <- kronecker(diag(t),matrix(1,1,n))
A_constr_1_4 <- rbind(A.1.1,A.1.2)

########################################
## Additive                           ##
########################################
bym.ad.re.re<- inla(OBS ~ -1 + I1 + I2 +
                      f(idx.u, model = model.s.u, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) +
                      f(idx.v, model = model.s.v, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) +
                      f(idy, model = model.t, constr=FALSE, extraconstr=list(A=A_constr.t, e=rep(0,k))), 
                    data = d, E = EXP,
                    family = "poisson", control.predictor = list(compute = TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.fixed(list(mean=list(I1=0, I2=0), prec=list(I1=0.01, I2=0.01))),
                    verbose = TRUE)

########################################
## Type I                             ##
########################################
bym.t1.re.re<- inla(OBS ~ -1 + I1 + I2 +
                      f(idx.u, model = model.s.u, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) +
                      f(idx.v, model = model.s.v, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) + 
                      f(idy, model = model.t, constr=FALSE, extraconstr=list(A=A_constr.t, e=rep(0,k))) +
                      f(idxy.1, model="iid", constr=TRUE, hyper=list(prec=list(prior=sdunif)))+
                      f(idxy.2, model="iid", constr=TRUE, hyper=list(prec=list(prior=sdunif))),
                    data = d, E = EXP,
                    family = "poisson", control.predictor = list(compute = TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.fixed(list(mean=list(I1=0, I2=0), prec=list(I1=0.01, I2=0.01))),
                    verbose = TRUE)

########################################
## Type II                            ##
########################################
bym.t2.re.re<- inla(OBS ~ -1 + I1 + I2 +
                      f(idx.u, model = model.s.u, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) +
                      f(idx.v, model = model.s.v, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) +
                      f(idy, model = model.t, constr=FALSE, extraconstr=list(A=A_constr.t, e=rep(0,k)) ) +
                      f(idxy.1, model="generic0", Cmatrix=R_1_2, rankdef=r_def_1_2, hyper=list(prec=list(prior=sdunif)), constr=FALSE, extraconstr=list(A=A_constr_1_2, e=rep(0,n)))+
                      f(idxy.2, model="generic0", Cmatrix=R_1_2, rankdef=r_def_1_2, hyper=list(prec=list(prior=sdunif)), constr=FALSE, extraconstr=list(A=A_constr_1_2, e=rep(0,n))),
                    data = d, E = EXP,
                    family = "poisson", control.predictor = list(compute = TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.fixed(list(mean=list(I1=0, I2=0), prec=list(I1=0.01, I2=0.01))),
                    verbose = TRUE)

########################################
## Type III                           ##
########################################
bym.t3.re.re<- inla(OBS ~ -1 + I1 + I2 +
                      f(idx.u, model = model.s.u, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) +
                      f(idx.v, model = model.s.v, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) + 
                      f(idy, model = model.t, constr=FALSE, extraconstr=list(A=A_constr.t, e=rep(0,k)) ) +
                      f(idxy.1, model="generic0", Cmatrix=R_1_3, rankdef=r_def_1_3, hyper=list(prec=list(prior=sdunif)), constr=FALSE, extraconstr=list(A=A_constr_1_3, e=rep(0,t))) +
                      f(idxy.2, model="generic0", Cmatrix=R_1_3, rankdef=r_def_1_3, hyper=list(prec=list(prior=sdunif)), constr=FALSE, extraconstr=list(A=A_constr_1_3, e=rep(0,t))),
                    data = d, E = EXP,
                    family = "poisson", control.predictor = list(compute = TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.fixed(list(mean=list(I1=0, I2=0), prec=list(I1=0.01, I2=0.01))),
                    verbose = TRUE)

########################################
## Type IV                            ##
########################################
bym.t4.re.re<- inla(OBS ~ -1 + I1 + I2 +
                      f(idx.u, model = model.s.u, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) +
                      f(idx.v, model = model.s.v, constr=FALSE, extraconstr=list(A=A_constr.s, e=rep(0,k))) + 
                      f(idy, model = model.t, constr=FALSE, extraconstr=list(A=A_constr.t, e=rep(0,k)) ) +
                      f(idxy.1, model="generic0", Cmatrix=R_1_4, rankdef=r_def_1_4, hyper=list(prec=list(prior=sdunif)), constr=FALSE, extraconstr=list(A=A_constr_1_4, e=rep(0,n+t)))+
                      f(idxy.2, model="generic0", Cmatrix=R_1_4, rankdef=r_def_1_4, hyper=list(prec=list(prior=sdunif)), constr=FALSE, extraconstr=list(A=A_constr_1_4, e=rep(0,n+t))),
                    data = d, E = EXP,
                    family = "poisson", control.predictor = list(compute = TRUE, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                    control.fixed(list(mean=list(I1=0, I2=0), prec=list(I1=0.01, I2=0.01))),
                    verbose = TRUE)

########################################
## save                               ##
########################################
resulta.inla.bym.re<- list(bym.ad.re=bym.ad.re.re,
                           bym.t1.re=bym.t1.re.re,
                           bym.t2.re=bym.t2.re.re,
                           bym.t3.re=bym.t3.re.re,
                           bym.t4.re=bym.t4.re.re)
save(resulta.inla.bym.re, file=paste0(gsub("\\.", "_", "resulta.inla.bym.re"), ".RData" ) )

################################################################################
################################################################################