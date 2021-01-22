################################################################################
## Title: Bayesian inference in multivariate spatio-temporal areal models     ##
##        using INLA: analysis of gender-based violence in small areas        ##
##                                                                            ##
## Authors: Vicente, G. - Goicoa, T. -  Ugarte, M.D.                          ##
##                                                                            ##
## https://doi.org/10.1007/s00477-020-01808-x                                 ##
##                                                                            ##
################################################################################
##                          Spatio-temporal M-models                          ##
################################################################################
rm(list=ls())

## libraries
library(sp); library(sf); library(tmap); library(tmaptools); library(INLA)
library(RColorBrewer); library(grid)

## Load data and Uttar Pradesh SpatialPolygonsDataFrame
load("./dataMmodel.RData")

####################################
##  Data organization
####################################
## Crimes
crime <- c("rape", "dowry")
e_crime <- paste0("e_", crime)
crime_name <- c("Rapes", "Dowry deaths")

## Number of crimes
k <- length(crime)

## Number of areas and number of time periods
n <- length(unique(data$ID_area))
t <- length(unique(data$ID_year))
x <- 1:t

## Initial and final time periods
t.from <- min(data$year)
t.to <- max(data$year)

## data.frame INLA
data<- data[order(data$ID_year, data$ID_area),]
d1<- as.list(rep(NA,length(crime)))
for(i in 1:length(crime)){ d1[[i]]<- list()}
for(i in 1:length(crime)){
  for(j in 1:t){
    d1[[i]][[j]]<- data[data$ID_year==j, c(crime[i], e_crime[i], "ID_area", "ID_year" )]
    d1[[i]][[j]]$ID_disease<- i
    colnames(d1[[i]][[j]])<- c("obs", "esp", "ID_area", "ID_year", "ID_disease")
  }
}
d.data_frame<- NULL
for(j in 1:t){d.data_frame<- rbind(d.data_frame, d1[[1]][[j]], d1[[2]][[j]])}   

## Intercept for each crime
intercepts<- paste0("I",1:length(crime))
d.data_frame[intercepts]<- NA
for(i in 1:length(crime)){
  d.data_frame[d.data_frame$ID_disease==i, intercepts[i]]<- 1
}

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
  d.data_frame[d.data_frame$ID_disease==i, idxy.[i]]<- 
    (d.data_frame[d.data_frame$ID_disease==i, c("ID_year")]-1)*n +
    d.data_frame[d.data_frame$ID_disease==i, c("ID_area")]
}

## rm
rm(list = c("d1", "intercepts", "idxy.", "i", "j"))

####################################
##  Load results
####################################
## INLA
load("resul/resulta_inla_bym_fe.RData")   # FE (bym)
load("resul/resulta_inla_bym_re.RData")   # RE (bym)
load("resul/resulta_inla_icar_fe.RData")  # FE (icar)
load("resul/resulta_inla_icar_re.RData")  # RE (icar)
load("resul/resulta_inla_lcar_fe.RData")  # FE (lcar)
load("resul/resulta_inla_lcar_re.RData")  # RE (lcar)
load("resul/resulta_inla_pcar_fe.RData")  # FE (pcar)
load("resul/resulta_inla_pcar_re.RData")  # RE (pcar)

resulta.inla<- list(
  ## icar-fe
  icar.ad.fe=resulta.inla.icar.fe$icar.ad.fe, icar.t1.fe=resulta.inla.icar.fe$icar.t1.fe,
  icar.t2.fe=resulta.inla.icar.fe$icar.t2.fe, icar.t3.fe=resulta.inla.icar.fe$icar.t3.fe,
  icar.t4.fe=resulta.inla.icar.fe$icar.t4.fe,
  ## icar-re
  icar.ad.re=resulta.inla.icar.re$icar.ad.re, icar.t1.re=resulta.inla.icar.re$icar.t1.re,
  icar.t2.re=resulta.inla.icar.re$icar.t2.re, icar.t3.re=resulta.inla.icar.re$icar.t3.re,
  icar.t4.re=resulta.inla.icar.re$icar.t4.re,
  ## pcar-fe
  pcar.ad.fe=resulta.inla.pcar.fe$pcar.ad.fe, pcar.t1.fe=resulta.inla.pcar.fe$pcar.t1.fe, 
  pcar.t2.fe=resulta.inla.pcar.fe$pcar.t2.fe,  pcar.t3.fe=resulta.inla.pcar.fe$pcar.t3.fe, 
  pcar.t4.fe=resulta.inla.pcar.fe$pcar.t4.fe,
  ## pcar-re
  pcar.ad.re=resulta.inla.pcar.re$pcar.ad.re, pcar.t1.re=resulta.inla.pcar.re$pcar.t1.re, 
  pcar.t2.re=resulta.inla.pcar.re$pcar.t2.re, pcar.t3.re=resulta.inla.pcar.re$pcar.t3.re, 
  pcar.t4.re=resulta.inla.pcar.re$pcar.t4.re,
  ## lcar-fe
  lcar.ad.fe=resulta.inla.lcar.fe$lcar.ad.fe, lcar.t1.fe=resulta.inla.lcar.fe$lcar.t1.fe, 
  lcar.t2.fe=resulta.inla.lcar.fe$lcar.t2.fe, lcar.t3.fe=resulta.inla.lcar.fe$lcar.t3.fe, 
  lcar.t4.fe=resulta.inla.lcar.fe$lcar.t4.fe,
  ## lcar-re
  lcar.ad.re=resulta.inla.lcar.re$lcar.ad.re, lcar.t1.re=resulta.inla.lcar.re$lcar.t1.re, 
  lcar.t2.re=resulta.inla.lcar.re$lcar.t2.re, lcar.t3.re=resulta.inla.lcar.re$lcar.t3.re, 
  lcar.t4.re=resulta.inla.lcar.re$lcar.t4.re,
  ## bym-fe
  bym.ad.fe=resulta.inla.bym.fe$bym.ad.fe, bym.t1.fe=resulta.inla.bym.fe$bym.t1.fe, 
  bym.t2.fe=resulta.inla.bym.fe$bym.t2.fe, bym.t3.fe=resulta.inla.bym.fe$bym.t3.fe,
  bym.t4.fe=resulta.inla.bym.fe$bym.t4.fe,
  ## bym-re
  bym.ad.re=resulta.inla.bym.re$bym.ad.re, bym.t1.re=resulta.inla.bym.re$bym.t1.re, 
  bym.t2.re=resulta.inla.bym.re$bym.t2.re, bym.t3.re=resulta.inla.bym.re$bym.t3.re, 
  bym.t4.re=resulta.inla.bym.re$bym.t4.re)

## rm
rm(list = c(paste0("resulta.inla.", c("icar", "lcar", "pcar", "bym"),
                   rep(c(".fe", ".re"),each=4))))

## WinBUGS
load("resul/resulta_winbugs_bym_fe.RData")   # FE (bym)
load("resul/resulta_winbugs_bym_re.RData")   # RE (bym)
load("resul/resulta_winbugs_icar_fe.RData")   # FE (icar)
load("resul/resulta_winbugs_icar_re.RData")   # RE (icar)
load("resul/resulta_winbugs_pcar_fe.RData")   # FE (pcar)
load("resul/resulta_winbugs_pcar_re.RData")   # RE (pcar)
load("resul/resulta_winbugs_lcar_fe.RData")   # FE (lcar)
load("resul/resulta_winbugs_lcar_re.RData")   # RE (lcar)


resulta.winbugs<- list(
  ## icar-fe
  icar.ad.fe=resulta.winbugs.icar.fe$icar.ad.fe, icar.t1.fe=resulta.winbugs.icar.fe$icar.t1.fe, 
  icar.t2.fe=resulta.winbugs.icar.fe$icar.t2.fe, icar.t3.fe=resulta.winbugs.icar.fe$icar.t3.fe, 
  icar.t4.fe=resulta.winbugs.icar.fe$icar.t4.fe,
  ## icar-re
  icar.ad.re=resulta.winbugs.icar.re$icar.ad.re, icar.t1.re=resulta.winbugs.icar.re$icar.t1.re, 
  icar.t2.re=resulta.winbugs.icar.re$icar.t2.re, icar.t3.re=resulta.winbugs.icar.re$icar.t3.re, 
  icar.t4.re=resulta.winbugs.icar.re$icar.t4.re,
  ## pcar-fe
  pcar.ad.fe=resulta.winbugs.pcar.fe$pcar.ad.fe, pcar.t1.fe=resulta.winbugs.pcar.fe$pcar.t1.fe,
  pcar.t2.fe=resulta.winbugs.pcar.fe$pcar.t2.fe, pcar.t3.fe=resulta.winbugs.pcar.fe$pcar.t3.fe, 
  pcar.t4.fe=resulta.winbugs.pcar.fe$pcar.t4.fe,
  ## pcar-re
  pcar.ad.re=resulta.winbugs.pcar.re$pcar.ad.re, pcar.t1.re=resulta.winbugs.pcar.re$pcar.t1.re, 
  pcar.t2.re=resulta.winbugs.pcar.re$pcar.t2.re, pcar.t3.re=resulta.winbugs.pcar.re$pcar.t3.re,
  pcar.t4.re=resulta.winbugs.pcar.re$pcar.t4.re,
  ## lcar-fe
  lcar.ad.fe=resulta.winbugs.lcar.fe$lcar.ad.fe, lcar.t1.fe=resulta.winbugs.lcar.fe$lcar.t1.fe, 
  lcar.t2.fe=resulta.winbugs.lcar.fe$lcar.t2.fe, lcar.t3.fe=resulta.winbugs.lcar.fe$lcar.t3.fe,
  lcar.t4.fe=resulta.winbugs.lcar.fe$lcar.t4.fe,
  ## lcar-re
  lcar.ad.re=resulta.winbugs.lcar.re$lcar.ad.re, lcar.t1.re=resulta.winbugs.lcar.re$lcar.t1.re, 
  lcar.t2.re=resulta.winbugs.lcar.re$lcar.t2.re, lcar.t3.re=resulta.winbugs.lcar.re$lcar.t3.re, 
  lcar.t4.re=resulta.winbugs.lcar.re$lcar.t4.re,
  ## bym-fe
  bym.ad.fe=resulta.winbugs.bym.fe$bym.ad.fe, bym.t1.fe=resulta.winbugs.bym.fe$bym.t1.fe, 
  bym.t2.fe=resulta.winbugs.bym.fe$bym.t2.fe, bym.t3.fe=resulta.winbugs.bym.fe$bym.t3.fe,
  bym.t4.fe=resulta.winbugs.bym.fe$bym.t4.fe,
  ## bym-re
  bym.ad.re=resulta.winbugs.bym.re$bym.ad.re, bym.t1.re=resulta.winbugs.bym.re$bym.t1.re, 
  bym.t2.re=resulta.winbugs.bym.re$bym.t2.re, bym.t3.re=resulta.winbugs.bym.re$bym.t3.re, 
  bym.t4.re=resulta.winbugs.bym.re$bym.t4.re)


t.resulta.winbugs<- list(
  ## icar-fe
  icar.ad.fe=t.resulta.winbugs.icar.fe$icar.ad.fe, icar.t1.fe=t.resulta.winbugs.icar.fe$icar.t1.fe,
  icar.t2.fe=t.resulta.winbugs.icar.fe$icar.t2.fe, icar.t3.fe=t.resulta.winbugs.icar.fe$icar.t3.fe, 
  icar.t4.fe=t.resulta.winbugs.icar.fe$icar.t4.fe,
  ## icar-re
  icar.ad.re=t.resulta.winbugs.icar.re$icar.ad.re, icar.t1.re=t.resulta.winbugs.icar.re$icar.t1.re, 
  icar.t2.re=t.resulta.winbugs.icar.re$icar.t2.re, icar.t3.re=t.resulta.winbugs.icar.re$icar.t3.re, 
  icar.t4.re=t.resulta.winbugs.icar.re$icar.t4.re,
  ## pcar-fe
  pcar.ad.fe=t.resulta.winbugs.pcar.fe$pcar.ad.fe, pcar.t1.fe=t.resulta.winbugs.pcar.fe$pcar.t1.fe, 
  pcar.t2.fe=t.resulta.winbugs.pcar.fe$pcar.t2.fe, pcar.t3.fe=t.resulta.winbugs.pcar.fe$pcar.t3.fe, 
  pcar.t4.fe=t.resulta.winbugs.pcar.fe$pcar.t4.fe,
  ## pcar-re
  pcar.ad.re=t.resulta.winbugs.pcar.re$pcar.ad.re, pcar.t1.re=t.resulta.winbugs.pcar.re$pcar.t1.re, 
  pcar.t2.re=t.resulta.winbugs.pcar.re$pcar.t2.re, pcar.t3.re=t.resulta.winbugs.pcar.re$pcar.t3.re, 
  pcar.t4.re=t.resulta.winbugs.pcar.re$pcar.t4.re,
  ## lcar-fe
  lcar.ad.fe=t.resulta.winbugs.lcar.fe$lcar.ad.fe, lcar.t1.fe=t.resulta.winbugs.lcar.fe$lcar.t1.fe, 
  lcar.t2.fe=t.resulta.winbugs.lcar.fe$lcar.t2.fe, lcar.t3.fe=t.resulta.winbugs.lcar.fe$lcar.t3.fe, 
  lcar.t4.fe=t.resulta.winbugs.lcar.fe$lcar.t4.fe,
  ## lcar-re
  lcar.ad.re=t.resulta.winbugs.lcar.re$lcar.ad.re, lcar.t1.re=t.resulta.winbugs.lcar.re$lcar.t1.re, 
  lcar.t2.re=t.resulta.winbugs.lcar.re$lcar.t2.re, lcar.t3.re=t.resulta.winbugs.lcar.re$lcar.t3.re, 
  lcar.t4.re=t.resulta.winbugs.lcar.re$lcar.t4.re,
  ## bym-fe
  bym.ad.fe=t.resulta.winbugs.bym.fe$bym.ad.fe, bym.t1.fe=t.resulta.winbugs.bym.fe$bym.t1.fe, 
  bym.t2.fe=t.resulta.winbugs.bym.fe$bym.t2.fe, bym.t3.fe=t.resulta.winbugs.bym.fe$bym.t3.fe,
  bym.t4.fe=t.resulta.winbugs.bym.fe$bym.t4.fe,
  ## bym-re
  bym.ad.re=t.resulta.winbugs.bym.re$bym.ad.re, bym.t1.re=t.resulta.winbugs.bym.re$bym.t1.re, 
  bym.t2.re=t.resulta.winbugs.bym.re$bym.t2.re, bym.t3.re=t.resulta.winbugs.bym.re$bym.t3.re, 
  bym.t4.re=t.resulta.winbugs.bym.re$bym.t4.re)

## rm
rm(list = c(paste0("resulta.winbugs.", c("icar", "lcar", "pcar", "bym"), 
                   rep(c(".fe", ".re"),each=4))))
rm(list = c(paste0("t.resulta.winbugs.", c("icar", "lcar", "pcar", "bym"), 
                   rep(c(".fe", ".re"),each=4))))

####################################
##  Selected model
####################################
resulta<- resulta.inla$icar.t2.fe 

##################################################################################
################                    Figures                       ################
##################################################################################
####################################
##  Figure 1. Map of the administrative division of Uttar Pradesh into districts
##            and its location in India (top right corner)                      
####################################
## india map
carto_india$color<- rep(c(0,2,0),c(30,1,2))
map_india<- tm_shape(carto_india) + 
  tm_polygons(col="color", legend.show=FALSE) +
  tm_layout(frame.lwd = 1.5)+
  tm_style("classic")+
  tm_shape(carto_up) + tm_polygons(col="palegreen3")

## uttar pradesh map 
map_up<- tm_shape(carto_up, bbox=c(77,23.5,85,30.5)) + tm_polygons(col="palegreen3") +
  tm_text(text="ID_area", size=1, fontface=5) + 
  tm_grid(n.x=4, n.y=4, alpha=0.5) +
  tm_style("classic")+ tm_format("World")


## pdf figure 1
pdf("figure_1.pdf", height=7.5, width=7.5, onefile=FALSE)
map_up
print(map_india, vp = grid::viewport(0.825, 0.8, width = 0.32, height = 0.32))
dev.off()

## rm
rm(list = c("map_india", "map_up"))

####################################
##  Figure 2. Evolution of the crude rates (per 100000 women) of rapes and
##            dowry deaths in Uttar Pradesh in the period 2001-2014
####################################
## crude rates (per 100000 women)
crude_rates <- 100000*
  aggregate(data[,c("rape", "dowry")], by=list(data$ID_year), FUN=sum)[,c("rape", "dowry")]/ 
  aggregate(data[,c("pop")], by=list(data$ID_year), FUN=sum)[,2]

## minimum and maximum values
inf <- round(min(crude_rates))
top <- round(max(crude_rates))

## color
selected_colors<- c("chocolate1", "tomato4")

## pdf figure 2
pdf("figure_2.pdf", height=5, width=7.5, onefile=FALSE)
plot(t.from:t.to, crude_rates[,"rape"], type="l", xlab ="",ylab ="", 
     ylim=c(inf, top), col=selected_colors[1], cex.axis=1, lwd=4)
lines(t.from:t.to, crude_rates[,"dowry"], col=selected_colors[2],lwd=4)
legend("topleft",  crime_name, ncol=1, pch=c("-", "-"), 
       col=selected_colors, bty="n",lwd=c(4,4), cex=1.2)
dev.off()

## rm
rm(list = c("crude_rates", "inf", "top", "selected_colors"))

####################################
##  Figure 3. Dispersion plots of the final relative risks for rapes and dowry
##            deaths obtained with the Type II interaction RE M-model with in
##            INLA (y-axis) vs. WinBUGS (x-axis), using the iCAR (first row),
##            pCAR (second row), LCAR (third row) and the BYM (last row)     
##            spatial priors
####################################
## pdf figure 3
pdf("figure_3.pdf", height=10, width=5, onefile=FALSE)
par(mfrow=c(4,2), oma = c(0, 2, 2, 0))
spatial <-c("icar", "pcar", "lcar", "bym")
selected_colors<- c("darkorchid4", "olivedrab4") # color
for(i in 1:length(spatial)){
  eval(parse(text= paste0("d.data_frame$risks<-resulta.inla$", 
                          spatial[i],".t2.re$summary.fitted.values$mean")))
  aux<-c()
  eval(parse(text = paste0("aux<-resulta.winbugs$",spatial[i], ".t2.re$mean$SMR")))
  for(j in 1:k){
    plot(as.vector(aux[,j,]),
         as.vector(d.data_frame[d.data_frame$ID_disease==j,c("risks")]),
         xlab=expression(paste(hat(R)[ijt], " WinBUGS")), ylab="",
         col=selected_colors[j], cex.lab=1, main=crime_name[j])
    title(ylab=expression(paste(hat(R)[ijt], " INLA")), line=2.5, cex.lab=1)
    abline(a=0,b=1,col =selected_colors[j])
  }
  rm(aux)
  d.data_frame<- d.data_frame[-which(names(d.data_frame)%in%"risks")]
}
mtext("iCAR spatial priors", outer=TRUE,  cex=0.9, line=-1.5, font=2)
mtext("pCAR spatial priors", outer=TRUE,  cex=0.9, line=-19.5, font=2)
mtext("lCAR spatial priors", outer=TRUE,  cex=0.9, line=-38.0, font=2)
mtext("BYM spatial priors", outer=TRUE,  cex=0.9, line=-56.5, font=2)
dev.off()

## rm
rm(list = c("spatial", "selected_colors", "i", "j"))

####################################
##  Figure 4. Posterior mean of the district-specific spatial risk, exp(theta_ij),
##            and the exceedence probabilities,  i.e., P(exp(theta_ij) > 1|O), for
##            rapeand dowry deaths
####################################
carto_use<- carto_up

## spatial effects 
spatial<- unlist(lapply(resulta$marginals.random$idx, function(x) inla.emarginal(exp,x)))
paleta <- brewer.pal(8,"YlOrRd")
inf<- min(spatial)
top<- max(spatial)
values<- c(round(seq(inf, 1, length.out = 5),2), round(seq(1, top, length.out = 5),2)[-1] )
for(j in 1:length(crime)){
  eval(parse(text=paste0("carto_use$spatial_",crime[j], "<- spatial[(n*(",j-1,")+1):(n*",j,")]")))
}

Rates.spat<- list()
for(j in 1:length(crime)){
  Rates.spat[[j]]<- tm_shape(carto_use) +
    tm_polygons(col=paste("spatial_",crime[j],sep=""), 
                palette=paleta, title="", legend.show=T, legend.reverse=T, 
                style="fixed", breaks=values, interval.closure="left") +
    tm_layout(main.title=crime_name[j], main.title.position="center", legend.text.size=1)
}

## probs spatial effects > 1
prob_spatial<- unlist(lapply(resulta$marginals.random$idx, function(x){1-inla.pmarginal(0, x)}))
paleta.p <- brewer.pal(5,"PuBu")
values.p <- c(0,0.1,0.2,0.8,0.9,1)
for(j in 1:length(crime)){
  eval(parse(text=paste0("carto_use$prob_spatial_",crime[j], "<- prob_spatial[(n*(",j-1,")+1):(n*",j,")]")))
}

Rates.prob.spat<- list()
for(j in 1:length(crime)){
  Rates.prob.spat[[j]]<- tm_shape(carto_use) +
    tm_polygons(col=paste("prob_spatial_",crime[j],sep=""), palette=paleta.p,
                title="", legend.show=T, legend.reverse=T, style="fixed", 
                breaks=values.p, interval.closure="left",
                labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) +
    tm_layout(main.title=crime_name[j], main.title.position="center", legend.text.size=1) 
}

## pdf figure 4
pdf("figure_4.pdf", height=10, width=10, onefile=FALSE)
pushViewport(viewport(layout=grid.layout(2,2)))
print(Rates.spat[[1]], vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Rates.prob.spat[[1]], vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(Rates.spat[[2]], vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(Rates.prob.spat[[2]], vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

## rm
rm(list=c("carto_use","spatial","paleta","top","inf","values","Rates.spat",
          "prob_spatial","paleta.p","values.p","Rates.prob.spat","j"))

####################################
##  Figure 5. Temporal pattern of incidence risks (posterior means of exp(gamma_tj))
##            for rape and dowry deaths in Uttar Pradesh
####################################
temp <- matrix(unlist(lapply(resulta$marginals.random$idy,
                             function(x) inla.emarginal(exp,x))), ncol=length(crime), byrow = FALSE)
aux <- lapply(resulta$marginals.random$idy, function(x) inla.tmarginal(exp,x))
q1 <- matrix(unlist(lapply(aux, function(x) inla.qmarginal(0.025,x))), 
             ncol=length(crime), byrow = FALSE)
q2 <- matrix(unlist(lapply(aux, function(x) inla.qmarginal(0.975,x))), 
             ncol=length(crime), byrow = FALSE)

## minimum and maximum values
inf <-min(q1)-0.05
top <-max(q2)+0.05

## color
selected_colors <- c(rgb(154,192,205,alpha=150, maxColorValue=255),
                     rgb(69,139,116,alpha=150, maxColorValue=255))

## pdf figure 5
pdf("figure_5.pdf", height=5, width=7.5, onefile=FALSE)
plot(range(x), c(inf, top), type="n", xlab="Year", ylab="",
     xaxt="n", cex.lab=1, cex.axis=1, cex.main=1, main=NULL)
title(ylab=expression(exp(gamma[t])), line=2.5, cex.lab=1.1)
axis(1, at=seq(1,t), labels=seq(t.from,t.to), las=0, cex.axis=1)
for(i in 1:length(crime)){
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(q1[,i], tail(q2[,i], 1), rev(q2[,i]), q1[1,i])
  polygon(X.Vec, Y.Vec, col = selected_colors[i], border = NA)
  lines(temp[,i])
  abline(h=1,lty=2)
  ## rm
  rm(list = c("X.Vec", "Y.Vec"))
}
legend("topleft", inset=.02, crime_name, fill=selected_colors,
       horiz=FALSE, cex=1,box.lty=0)
dev.off()

# rm
rm(list=c("temp", "aux", "q1", "q2","inf", "top", "selected_colors","i"))

####################################
##  Figure 6. Map of estimated incidence risks for rape (top) and posterior probabilities
##            that the relative risk is greater than 1 (P(Ritj>1|O)) (bottom) in Uttar 
##            Pradesh
##  Figure 7. Map of estimated incidence risks for dowry deaths (top) and posterior 
##            probabilities that the relative risk isgreater than one (P(Ritj > 1|O)) 
##            (bottom) in Uttar Pradesh
####################################
carto_use<- carto_up

## risks
d.data_frame$risks<- resulta$summary.fitted.values$mean
risks<- list()
for(i in 1:length(crime)){
  risks[[i]]<- data.frame(matrix(d.data_frame[d.data_frame$ID_disease==i,c("risks")], 
                                 nrow=n, ncol=t, byrow=FALSE))
  colnames(risks[[i]])<- paste0("risk_", crime[i],"_",t.from:t.to)
  for(l in t.from:t.to){
    eval(parse(text=  paste0("carto_use$risks_",crime[i], "_",
                             l,"<- risks[[",i,"]]$risk_",crime[i],"_", l) ))
  }
}

paleta <- brewer.pal(8,"YlOrRd")
minimo<- min(d.data_frame$risks)
maximo<- max(d.data_frame$risks)
values<- c(round(seq(minimo, 1, length.out = 5),2),
           round(seq(1, maximo, length.out = 5),2)[-1] )
d.data_frame<- d.data_frame[-which(names(d.data_frame)%in% c("risks"))]


Rates.Risks<- list()
for(j in 1:length(crime)){
  Rates.Risks[[j]]<- tm_shape(carto_use) + 
    tm_polygons(col=paste("risks_",crime[j],"_",seq(t.from,t.to), sep=""),
                palette=paleta, title="", legend.show=T, legend.reverse=T,
                style="fixed", breaks=values, interval.closure="left") +
    tm_layout(main.title=crime_name[j] , main.title.position="center",
              legend.text.size=1,
              panel.labels=as.character(round(seq(t.from,t.to,length.out=t))),
              panel.show=T, panel.label.size=2, panel.label.bg.color="lightskyblue",
              legend.outside=T, legend.outside.position="right", legend.frame=F,
              legend.outside.size=0.25, outer.margins=c(0.02,0.01,0.02,0.01)) +
    tm_facets(ncol=5, nrow=3) 
}

## prob risks > 1
d.data_frame$prp<- 1-resulta$summary.linear.predictor[,"0 cdf"]
p.risks<- list()
for(i in 1:length(crime)){
  p.risks[[i]]<- data.frame(matrix(d.data_frame[d.data_frame$ID_disease==i,c("prp")], 
                                   nrow=n, ncol=t, byrow=FALSE))
  colnames(p.risks[[i]])<- paste0("p.risk_", crime[i],"_",t.from:t.to)
  for(l in t.from:t.to){
    eval(parse(text=  paste0("carto_use$p.risks_",crime[i], "_", 
                             l,"<- p.risks[[",i,"]]$p.risk_",crime[i],"_", l) ))
  }
}

paleta.p <- brewer.pal(5,"PuBu")
values.p <- c(0,0.1,0.2,0.8,0.9,1)
d.data_frame<- d.data_frame[-which(names(d.data_frame)%in% c("prp"))]

Prob.Risks<- list()
for(j in 1:length(crime)){
  Prob.Risks[[j]]<- tm_shape(carto_use) + 
    tm_polygons(col=paste("p.risks_",crime[j],"_",seq(t.from,t.to), sep=""), 
                palette=paleta.p, title="", legend.show=T, legend.reverse=T, 
                style="fixed", breaks=values.p, interval.closure="left") +
    tm_layout(main.title=crime_name[j], main.title.position="center", 
              legend.text.size=1,
              panel.labels=as.character(round(seq(t.from,t.to,length.out=t))),
              panel.show=T, panel.label.size=2, panel.label.bg.color="lightskyblue",
              legend.outside=T, legend.outside.position="right", legend.frame=F, 
              legend.outside.size=0.25, outer.margins=c(0.02,0.01,0.02,0.01)) +
    tm_facets(ncol=5, nrow=3) 
}

## pdf figure 6
pdf("figure_6.pdf", height=10, width=7.5, onefile=FALSE)
pushViewport(viewport(layout=grid.layout(2,1)))
print(Rates.Risks[[1]], vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Prob.Risks[[1]], vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

## pdf figure 7
pdf("figure_7.pdf", height=10, width=7.5, onefile=FALSE)
pushViewport(viewport(layout=grid.layout(2,1)))
print(Rates.Risks[[2]], vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(Prob.Risks[[2]], vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

## rm
rm(list = c("carto_use", "risks", "paleta", "minimo", "maximo", "values", 
            "Rates.Risks", "p.risks", "paleta.p", "values.p", "Prob.Risks", 
            "i","l","j"))

####################################
##  Figure 8. Temporal evolution of final risk estimates for rape and dowry deaths
##            in some districts in Uttar Pradesh:  Ghazlabad, Kheri, Mainpuri, 
##            Sant Kabir Nagar, and Varanasi
####################################
## selected districts
id.area<- c(2,28,43,49,61,70)


d.data_frame$risks.mean <- resulta$summary.fitted.values[,c("mean")]
d.data_frame$risks.q025 <- resulta$summary.fitted.values[,c("0.025quant")]
d.data_frame$risks.q975 <- resulta$summary.fitted.values[,c("0.975quant")]

## minimum and maximum values
inf <- round(min(d.data_frame[d.data_frame$ID_area%in%id.area, c("risks.q025")]),2)-0.05
top <- round(max(d.data_frame[d.data_frame$ID_area%in%id.area, c("risks.q975")]),2)+1.5

## color
selected_colors<-c(rgb(154,192,205,alpha=150, maxColorValue=255),
                   rgb(69,139,116,alpha=150, maxColorValue=255))

## pdf figure 8
pdf("figure_8.pdf", height=12, width=12, onefile=FALSE)
par(mfrow=c(3,2))

for (i in id.area){
  plot(range(x),c(inf, top), type="n",xlab="Year", ylab="", xaxt="n", 
       cex.main=2.5, cex.lab=1.5, cex.axis=1.5, 
       main=paste0(unique(data[data$ID_area==i, c("dist")])," (ID area ", i,")")) 
  title(ylab=expression(R[ijt]), line=2.0, cex.lab=1.7)
  axis(1, at=seq(1,t), labels=seq(t.from,t.to), las=0, cex.axis=1.5)
  for(j in 1:length(crime)){
    X.Vec <- c(x, tail(x, 1), rev(x), x[1])
    Y.Vec <- c(d.data_frame[d.data_frame$ID_area==i & d.data_frame$ID_disease==j, c("risks.q025")],
               tail(d.data_frame[d.data_frame$ID_area==i & d.data_frame$ID_disease==j, c("risks.q975")], 1), 
               rev(d.data_frame[d.data_frame$ID_area==i & d.data_frame$ID_disease==j, c("risks.q975")]),
               d.data_frame[d.data_frame$ID_area==i & d.data_frame$ID_disease==j, c("risks.q025")][1])
    polygon(X.Vec, Y.Vec, col = selected_colors[j], border = NA)
    lines(d.data_frame[d.data_frame$ID_area==i & d.data_frame$ID_disease==j, c("risks.mean")], lwd=1.5)
    abline(h=1,lty=2)
    legend("topleft", inset=.02, crime_name, fill=selected_colors, horiz=FALSE, cex=2.2,box.lty=0)
    rm(list = c("X.Vec", "Y.Vec"))
  }
}
dev.off()

## rm
d.data_frame<- d.data_frame[-which(names(d.data_frame)%in% c("risks.mean", "risks.q025", "risks.q975"))]
rm(list=c("id.area", "inf", "top", "selected_colors", "i", "j"))

##################################################################################
################                    Tables                        ################
##################################################################################
####################################
##  Table 1.  Descriptive statistics. Minimun, first quartile (q1), mean, third
##            quartile (q3), maximun, standard desviation (sd), and coefficient of
##            variation of the number of rapes and dowry deaths in the districts 
##            of Uttar Pradesh per year
####################################

table.1<- matrix(NA, nrow=t, ncol=15)
colnames(table.1)<- c("Year", "min", "q1", "mean", "q3", "max", "sd", "|cv|",
                      "min", "q1", "mean", "q3", "max", "sd", "|cv|")
for(i in 1:t){
  table.1[i,]<- c(2000+i, summary(data[data$ID_year==i, crime[1]])[-3],
                  sd(data[data$ID_year==i, crime[1]]), NA,
                  summary(data[data$ID_year==i, crime[2]])[-3],
                  sd(data[data$ID_year==i, crime[2]]), NA)
  table.1[i,8]<- abs(table.1[i,7]/ table.1[i,4])
  table.1[i,15]<- abs(table.1[i,14]/ table.1[i,11])
}
table.1<- round(table.1, 1)

## latex
latex_table<-xtable::xtable(table.1,
                            caption="Descriptive statistics. Minimun, first quartile ($q_1$), mean, third quartile ($q_3$), maximun, standard desviation (sd), and coefficient of variation of the number of rapes and dowry deaths in the districts of Uttar Pradesh per year.",
                            label="t_crime", digits=c(1),
                            display=c("d", "d", "d","f","f","f","d","f","f", "d","f","f","f","d","f","f"))
xtable::print.xtable(latex_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))

## rm
rm(list = c("table.1", "i", "latex_table"))

####################################
##  Table 2.  Correlations between spatial (by year) and temporal patterns
##            (by district) of rape and dowry deaths 
####################################
cor.spat<- c()
for(i in 1:t){ cor.spat<- c(cor.spat, cor(data[data$ID_year==i,c("smr_rape","smr_dowry")])[1,2]) }
cor.temp<- c()
for(i in 1:n){ cor.temp<- c(cor.temp, cor(data[data$ID_area==i,c("smr_rape","smr_dowry")])[1,2]) }


table.2<- rbind(c(summary(cor.spat), sd(cor.spat), abs(sd(cor.spat)/mean(cor.spat))), 
                c(summary(cor.temp), sd(cor.temp), abs(sd(cor.temp)/mean(cor.temp))))
colnames(table.2)<- c("min", "q1", "median", "mean", "q3", "max", "sd", "|cv|")
table.2<- as.data.frame(table.2)
table.2$Correlation <- c("spatial patterns", "temporal trends")
table.2<- table.2[,c("Correlation", "min", "q1", "median", "mean", "q3", "max", 
                     "sd", "|cv|")]


## latex
latex_table<-xtable::xtable(table.2,
                            caption="Correlations between spatial (by year) and temporal patterns (by district) of rape and dowry deaths.",
                            label="t_corre", digits=c(3))
xtable::print.xtable(latex_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))

## rm
rm(list = c("cor.spat", "cor.temp", "table.2", "i", "latex_table"))

####################################
##  Table 3.  Model selection criteria, DIC, WAIC and LS, for different models
####################################
table.3<- matrix(NA, nrow=length(resulta.inla), ncol=3)
colnames(table.3)<-c("DIC","WAIC","LS")
random.effect<- strsplit(names(resulta.inla),".", fixed = TRUE)

spatemp<-spat<-temp<-c()
for(i in 1:length(resulta.inla)){
  aux<-resulta.inla[[i]]
  table.3[i,]<-round(c(aux$dic$dic, aux$waic$waic, -mean(log(aux$cpo$cpo),na.rm=T)), 3)
  spatemp[i]<- random.effect[[i]][3]
  spat[i]<- random.effect[[i]][1]
  temp[i]<- random.effect[[i]][2]
}
table.3<-as.data.frame(table.3)
table.3<-cbind(spat, spatemp, temp, table.3)


## latex
x_table<-xtable::xtable(table.3, 
                        caption="Model selection criteria, DIC, WAIC and LS, for different models.", 
                        label="t_dic", digits=c(3))
xtable::print.xtable(x_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))

## rm
rm(list=c("table.3", "random.effect","spatemp", "spat", "temp", "i", "aux", "x_table"))

##################################################################################
################                  Appendix                        ################
##################################################################################
####################################
##  Table A.1.  District identifiers (ID) of Uttar Pradesh
####################################
ID<- data.frame(dist=carto_up$dist, ID_area=carto_up$ID_area)
table.a1<- cbind(ID[1:24,c("ID_area", "dist")], ID[25:48,c("ID_area", "dist")],
                 ID[49:72,c("ID_area", "dist")])
colnames(table.a1)<- c("ID", "Dist","ID", "Dist", "ID", "Dist")


## latex
x_table<-xtable::xtable(table.a1, 
                        caption="District identifiers (ID) of Uttar Pradesh.",
                        label="table_a1", digits=c(3))
xtable::print.xtable(x_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))

## rm
rm(list = c("ID", "table.a1", "x_table"))

####################################
##  Table A.2.  Posterior means, standard deviations, and 95% credible intervals 
##              for the crime-specific intercepts (alpha_j , j=1,2) of the models
##              with a spatio-temporal Type II interaction term
####################################
table.r<- table.l<- c()
spatial <-c("icar", "pcar", "lcar", "bym")
for(i in 1:length(spatial)){
  aux<- aux1<- aux2<- c()
  aux1<-eval(parse(text= paste0("resulta.winbugs$", spatial[i], ".t2.fe") ))
  aux2<-eval(parse(text= paste0("resulta.inla$", spatial[i], ".t2.fe")))
  
  aux.r<- aux1.r<- aux2.r<- c()
  aux1.r<- eval(parse(text= paste0("resulta.winbugs$", spatial[i], ".t2.re") ))
  aux2.r<- eval(parse(text= paste0("resulta.inla$", spatial[i], ".t2.re")))
  
  for(j in 1:k){
    aux<- rbind(aux,
                aux1$summary[grepl("mu", rownames(aux1$summary), fixed=TRUE),][j,c(1,2,3,7)],
                aux2$summary.fixed[j,c(1,2,3,5)])
    aux.r<- rbind(aux.r,
                aux1.r$summary[grepl("mu", rownames(aux1.r$summary), fixed=TRUE),][j,c(1,2,3,7)],
                aux2.r$summary.fixed[j,c(1,2,3,5)])
  }
  table.l<- rbind(table.l, aux)
  table.r<- rbind(table.r, aux.r)
  rm(list = c("aux", "aux1", "aux2", "aux.r", "aux1.r", "aux2.r"))
}
table.a2<- as.data.frame(cbind(table.l, table.r))
colnames(table.a2)<- c("mean.fe", "sd.fe", "q.025.fe", "q.975.fe", "mean.re", 
                       "sd.re", "q.025.re", "q.975.re")
rownames(table.a2)<- NULL
table.a2$spatial<- rep(c("iCAR", "pCAR", "LCAR", "BYM"), each=length(spatial))
table.a2$crime<- rep(rep(crime_name, each=2), times= length(spatial))
table.a2$met<- rep(c("MCMC", "INLA"), times= length(spatial)*k)
table.a2<- table.a2[,c("spatial", "crime", "met", "mean.fe", "sd.fe", "q.025.fe",
                       "q.975.fe", "mean.re", "sd.re", "q.025.re", "q.975.re" )]


## latex
latex_table<-xtable::xtable(table.a2,
                            caption="Posterior means, standard deviations, and 95% credible intervals for the crime-specific intercepts (alpha_{j}, j=1,2) of the models with a spatio-temporal Type II interaction term.",
                            label="t_mu", digits=3)
xtable::print.xtable(latex_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))

# rm
rm(list = c("table.a2", "spatial", "table.l", "table.r", "latex_table", "i", "j"))

####################################
##  Table A.3.  Posterior means, standard deviations, and 95% credible intervals 
##              for the hyperparameters of the models with a spatio-temporal Type 
##              II interaction term
####################################
hyperparam<- c()
##############
## icar
##############
## fe
m<- c("icar"); m.e<- c("FE")

model.inla<- resulta.inla$icar.t2.fe
hyper.i<- matrix(NA, nrow = length(crime), ncol=3)
colnames(hyper.i)<- c("mean", "q.025", "q.975" )
rownames(hyper.i)<- c(paste0("sigma.idxy.",1:length(crime)))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))

model.winbugs<- resulta.winbugs$icar.t2.fe
names_summ<-rownames(model.winbugs$summary)
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),] # sigma spatio-temporal
hyper.w<- rbind(sd_st)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )
hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)

hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep(paste0("sigma.idxy",1:2), times=2)
rownames(hyper)<- NULL
hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)
rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ", 
            "sd_st", "hyper.w", "hyper", "m.e"))

m.e<- c("RE")
model.inla<- resulta.inla$icar.t2.re
hyper.i<- matrix(NA, nrow = length(crime)+2, ncol=3)
colnames(hyper.i)<- c("mean", "q.025", "q.975" )
rownames(hyper.i)<- c("sigma.spat", "sigma.temp", paste0("sigma.idxy.",1:length(crime)))
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta5 for idx`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta5 for idy`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[3,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[4,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))

model.winbugs<- resulta.winbugs$icar.t2.re
names_summ<-rownames(model.winbugs$summary)
sd_spat<- model.winbugs$summary[grepl("sdstruct",names_summ, fixed=TRUE),] # sigma spatial
sd_spat<-sd_spat[1,]
sd_temp<- model.winbugs$summary[grepl("sdstructg",names_summ, fixed=TRUE),] # sigma spatial
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),] # sigma spatio-temporal
hyper.w<-rbind(sd_spat,sd_temp,sd_st)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )

hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)
hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep( c("sigma.spat", "sigma.temp", paste0("sigma.idxy.",1:length(crime))), times=2)
rownames(hyper)<- NULL
hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)

rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ", 
            "sd_spat","sd_temp","sd_st", "hyper.w", "hyper", "m.e"))

##############
## pcar
##############
m<- c("pcar")
## FE
m.e<- c("FE")

model.inla<- resulta.inla$pcar.t2.fe
hyper.i<- matrix(NA, nrow = length(crime)+2, ncol=3)
colnames(hyper.i)<- c("mean", "q.025", "q.975" )
rownames(hyper.i)<- c(paste0("sigma.idxy.",1:2), paste0("gamma.",1:length(crime)))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta1 for idx`)
hyper.i[3,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta2 for idx`)
hyper.i[4,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))

model.winbugs<- resulta.winbugs$pcar.t2.fe
names_summ<-rownames(model.winbugs$summary)
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),] # sigma spatio-temporal
gamma<- model.winbugs$summary[grepl("gamma",names_summ, fixed=TRUE),] # sigma spatio-temporal
hyper.w<-rbind(sd_st,gamma)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )

hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)

hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep(c(paste0("sigma.idxy",1:2), paste0("lambda",1:2)), times=2)

rownames(hyper)<- NULL
hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)

rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ",
            "sd_st", "gamma", "hyper.w", "hyper", "m.e"))

## RE
m.e<- c("RE")
model.inla<- resulta.inla$pcar.t2.re
hyper.i<- matrix(NA, nrow = length(crime)+4, ncol=3)
colnames(hyper.i)<- c("mean", "q.025","q.975" )
rownames(hyper.i)<- c("sigma.spat", "sigma.temp", paste0("sigma.idxy.",1:length(crime)), paste0("gamma.",1:length(crime)) )
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta7 for idx`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta5 for idy`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[3,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[4,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta1 for idx`)
hyper.i[5,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta2 for idx`)
hyper.i[6,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))

model.winbugs<- resulta.winbugs$pcar.t2.re
names_summ<-rownames(model.winbugs$summary)
sd_spat<- model.winbugs$summary[grepl("sdstruct",names_summ, fixed=TRUE),]
sd_spat<-sd_spat[1,]
sd_temp<-model.winbugs$summary[grepl("sdstructg",names_summ, fixed=TRUE),]
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),]
gamma<- model.winbugs$summary[grepl("gamma",names_summ, fixed=TRUE),]
hyper.w<- rbind(sd_spat,sd_temp,sd_st,gamma)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )

hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)
hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep( c("sigma.spat", "sigma.temp", paste0("sigma.idxy.",1:length(crime)), paste0("lambda",1:2) ), times=2)
rownames(hyper)<- NULL

hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)

##
rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ", 
            "sd_spat","sd_temp","sd_st","gamma","hyper.w", "hyper", "m.e"))
rm(m)
##############
## lcar
##############
m<- c("lcar")
## FE
m.e<- c("FE")

model.inla<- resulta.inla$lcar.t2.fe
hyper.i<- matrix(NA, nrow = length(crime)+2, ncol=3)
colnames(hyper.i)<- c("mean", "q.025", "q.975" )
rownames(hyper.i)<- c(paste0("sigma.idxy.",1:2), paste0("gamma.",1:length(crime)))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta1 for idx`)
hyper.i[3,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta2 for idx`)
hyper.i[4,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))

model.winbugs<- resulta.winbugs$lcar.t2.fe
names_summ<-rownames(model.winbugs$summary)
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),] # sigma spatio-temporal
gamma<- model.winbugs$summary[grepl("gamma",names_summ, fixed=TRUE),] # sigma spatio-temporal
hyper.w<-rbind(sd_st,gamma)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )

hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)
hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep(c(paste0("sigma.idxy",1:2), paste0("lambda",1:2)), times=2)
rownames(hyper)<- NULL
hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)

rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ",
            "sd_st", "gamma", "hyper.w", "hyper", "m.e"))

## RE
m.e<- c("RE")

model.inla<- resulta.inla$lcar.t2.re
hyper.i<- matrix(NA, nrow = length(crime)+4, ncol=3)
colnames(hyper.i)<- c("mean", "q.025","q.975" )
rownames(hyper.i)<- c("sigma.spat", "sigma.temp", paste0("sigma.idxy.",1:length(crime)), paste0("gamma.",1:length(crime)) )
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta7 for idx`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta5 for idy`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[3,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[4,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta1 for idx`)
hyper.i[5,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/(1+exp(-x))}, model.inla$marginals.hyperpar$`Theta2 for idx`)
hyper.i[6,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))

model.winbugs<- resulta.winbugs$lcar.t2.re
names_summ<-rownames(model.winbugs$summary)
sd_spat<- model.winbugs$summary[grepl("sdstruct",names_summ, fixed=TRUE),] # sigma spatial
sd_spat<-sd_spat[1,]
sd_temp<-model.winbugs$summary[grepl("sdstructg",names_summ, fixed=TRUE),] # sigma temp
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),] # sigma spatio-temporal
gamma<- model.winbugs$summary[grepl("gamma",names_summ, fixed=TRUE),] # sigma spatio-temporal
hyper.w<- rbind(sd_spat,sd_temp,sd_st,gamma)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )

hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)
hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep( c("sigma.spat", "sigma.temp", paste0("sigma.idxy.",1:length(crime)), paste0("lambda",1:2) ), times=2)
rownames(hyper)<- NULL

hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)

rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ",
            "sd_spat","sd_temp","sd_st","gamma","hyper.w", "hyper", "m.e"))
rm(m)
##############
## bym
##############
m<- c("bym")
## FE
m.e<- c("FE")

model.inla<- resulta.inla$bym.t2.fe
hyper.i<- matrix(NA, nrow = length(crime), ncol=3)
colnames(hyper.i)<- c("mean", "q.025","q.975" )
rownames(hyper.i)<- c(paste0("sigma.idxy.",1:length(crime)))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025, 0.975), marg))

model.winbugs<- resulta.winbugs$bym.t2.fe
names_summ<-rownames(model.winbugs$summary)
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),] # sigma spatio-temporal
hyper.w<-rbind(sd_st)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )

hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)
hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep(c(paste0("sigma.idxy",1:2)), times=2)
rownames(hyper)<- NULL
hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)

rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ",
            "sd_st", "hyper.w", "hyper", "m.e"))

## RE
m.e<- c("RE")

model.inla<- resulta.inla$bym.t2.re
hyper.i<- matrix(NA, nrow = length(crime)*2+1, ncol=3)
colnames(hyper.i)<- c("mean","q.025","q.975" )
rownames(hyper.i)<- c("sigma.idx.u", "sigma.idx.v", "sigma.idy", paste0("sigma.idxy.",1:length(crime)))
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta5 for idx.u`)
hyper.i[1,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025,0.975), marg))
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta5 for idx.v`)
hyper.i[2,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025,0.975), marg))
marg<- inla.tmarginal(fun= function(x){exp(-(1/2)*x)}, model.inla$marginals.hyperpar$`Theta5 for idy`)
hyper.i[3,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025,0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.1`)
hyper.i[4,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025,0.975), marg))
marg<- inla.tmarginal(fun= function(x){ 1/sqrt(x)}, model.inla$marginals.hyperpar$`Precision for idxy.2`)
hyper.i[5,]<- c(inla.emarginal(fun= function(x){ x}, marg), inla.qmarginal(c(0.025,0.975), marg))

model.winbugs<- resulta.winbugs$bym.t2.re
names_summ<-rownames(model.winbugs$summary)
sd_spat<- model.winbugs$summary[grepl("sdstruct",names_summ, fixed=TRUE),] # sigma spatial
sd_spat<-sd_spat[1:2,]
sd_temp<- model.winbugs$summary[grepl("sdstructg",names_summ, fixed=TRUE),] # sigma spatial
sd_st<- model.winbugs$summary[grepl("sdZet",names_summ, fixed=TRUE),] # sigma spatio-temporal
hyper.w<-rbind(sd_spat,sd_temp,sd_st)
hyper.w<- hyper.w[,c("mean","2.5%", "97.5%")]
colnames(hyper.w)<- c("mean", "q.025", "q.975" )

hyper<- rbind(hyper.i, hyper.w)
hyper<- as.data.frame(hyper)
hyper$model_names<- m
hyper$Mmodel<- m.e
hyper$metod<- rep(c("INLA", "MCMC"), each=dim(hyper)[1]/2)
hyper$param<- rep(c("sigma.idx.u","sigma.idx.v", "sigma.temp",paste0("sigma.idxy", 1:2)), times=2)
rownames(hyper)<- NULL
hyper<- hyper[,c("model_names", "Mmodel", "metod", "param", "mean", "q.025", "q.975")]

hyperparam<- rbind(hyperparam, hyper)

##
rm(list = c("model.inla", "hyper.i", "marg", "model.winbugs", "names_summ", 
            "sd_spat","sd_temp","sd_st","hyper.w", "hyper", "m.e"))
rm(m)

##############
## presentation table A.3
##############
table.a3<- cbind(hyperparam[hyperparam$metod=="INLA",c("model_names", "Mmodel", "param", "mean", "q.025", "q.975")], 
                hyperparam[hyperparam$metod=="MCMC",c("mean", "q.025", "q.975")])
## latex
latex_table<-xtable::xtable(table.a3, 
                            caption="Posterior means, standard deviations, and 95% credible intervals for the hyperparameters of the models with a spatio-temporal Type II interaction term.",
                            label="t_hyper", digits=c(3))
xtable::print.xtable(latex_table, include.rownames = FALSE, comment=FALSE,
                     caption.placement = getOption("xtable.caption.placement", "top"))

## rm
rm(list = c("table.a3", "latex_table", "hyperparam"))
##################################################################################
##################################################################################