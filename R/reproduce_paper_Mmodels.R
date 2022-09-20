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

library(grid)
library(INLA)
library(RColorBrewer)
library(sf)
library(tmap)
library(tmaptools)
 

########################################
## 1) Load data and cartography files ##
########################################
load("./dataMmodel.RData")

head(carto_UP)
str(data_UP)

J <- length(data_UP)
S <- length(unique(data_UP[[1]]$dist))
T <- length(unique(data_UP[[1]]$year))
t.from <- min(data_UP[[1]]$year)
t.to <- max(data_UP[[1]]$year)

data <- do.call(rbind,data_UP)
data$Crime <- rep(1:J,each=S*T)


######################################
## 2) Load previously fitted models ##
######################################

## INLA ##
load("results/MODELS_inla_icar.RData")
load("results/MODELS_inla_lcar.RData")
load("results/MODELS_inla_pcar.RData")

## WinBUGS ##
load("resul/results_winbugs_bym_fe.RData")
load("resul/results_winbugs_bym_re.RData")
load("resul/results_winbugs_icar_fe.RData")
load("resul/results_winbugs_icar_re.RData")
load("resul/results_winbugs_pcar_fe.RData")
load("resul/results_winbugs_pcar_re.RData")
load("resul/results_winbugs_lcar_fe.RData")
load("resul/results_winbugs_lcar_re.RData")


##################################################################################
################                    Figures                       ################
##################################################################################

## Folder to save the figures ##
if(!file.exists("figures")) {dir.create("figures")}


####################################
##  Figure 1. Map of the administrative division of Uttar Pradesh into districts
##            and its location in India (top right corner)                      
####################################
carto_india$color <- rep(c(0,2,0),c(30,1,2))

map_india<- tm_shape(carto_india) + 
  tm_polygons(col="color", legend.show=FALSE) +
  tm_layout(frame.lwd=1.5) +
  tm_style("classic") +
  tm_shape(carto_UP) + tm_polygons(col="palegreen3")

map_UP <- tm_shape(carto_UP, bbox=c(77,23.5,85,30.5)) + tm_polygons(col="palegreen3") +
  tm_text(text="ID_area", size=1, fontface=5) + 
  tm_grid(n.x=4, n.y=4, alpha=0.5) +
  tm_style("classic") + tm_format("World")


## File: figure_1.pdf
graphics.off()
pdf("figures/figure_1.pdf", height=7.5, width=7.5, onefile=FALSE)
print(map_UP)
print(map_india, vp = grid::viewport(0.825, 0.8, width = 0.32, height = 0.32))
dev.off()



####################################
##  Figure 2. Evolution of the crude rates (per 100000 women) of rapes and
##            dowry deaths in Uttar Pradesh in the period 2001-2014
####################################
crude_rates <- data.frame(dowry=1e+5*aggregate(data_UP$dowry$obs, by=list(data_UP$dowry$year), sum)$x/aggregate(data_UP$dowry$pop, by=list(data_UP$dowry$year), sum)$x,
                          rape=1e+5*aggregate(data_UP$rape$obs, by=list(data_UP$rape$year), sum)$x/aggregate(data_UP$rape$pop, by=list(data_UP$rape$year), sum)$x)

inf <- round(min(crude_rates))
top <- round(max(crude_rates))

selected_colors <- c("chocolate1", "tomato4")


## File: figure_2.pdf
pdf("figures/figure_2.pdf", height=5, width=7.5, onefile=FALSE)
plot(t.from:t.to, crude_rates$rape, type="l", xlab ="",ylab ="", 
     ylim=c(inf, top), col=selected_colors[1], cex.axis=1, lwd=4)
lines(t.from:t.to, crude_rates$dowry, col=selected_colors[2],lwd=4)
legend("topleft",  c("Rapes","Dowry deaths"), ncol=1, pch=c("-", "-"), 
       col=selected_colors, bty="n",lwd=c(4,4), cex=1.2)
dev.off()


####################################
##  Figure 3. Dispersion plots of the final relative risks for rapes and dowry
##            deaths obtained with the Type II interaction RE M-model with in
##            INLA (y-axis) vs. WinBUGS (x-axis), using the iCAR (first row),
##            pCAR (second row), LCAR (third row) and the BYM (last row)     
##            spatial priors
####################################
Models.INLA <- list(iCAR=MODELS.inla.icar$TypeII,LCAR=MODELS.inla.lcar$TypeII,pCAR=MODELS.inla.pcar$TypeII)
Models.WinBUGS <- list(iCAR=results.winbugs.icar.re$icar.t2.re,
                       LCAR=results.winbugs.lcar.re$lcar.t2.re,
                       pCAR=results.winbugs.pcar.re$pcar.t2.re)

risks.INLA <- lapply(Models.INLA, function(x) matrix(Models.INLA$iCAR$summary.fitted.values$mean,S*T,J,byrow=F))
risks.WinBUGS <- lapply(Models.WinBUGS, function(x) cbind(c(x$mean$SMR[,1,]),c(x$mean$SMR[,2,])))


## File: figure_3.pdf
graphics.off()

pdf("figures/figure_3.pdf", height=10, width=7, onefile=FALSE)

par(mfrow=c(3,2), oma = c(0, 2, 2, 0))
selected_colors <- c("darkorchid4", "olivedrab4") # color

for(i in names(risks.INLA)){
  plot(risks.INLA[[i]][,1], risks.WinBUGS[[i]][,1], col=selected_colors[1], cex.lab=1, main="Rapes",
       xlab=expression(paste(hat(R)[ijt], " WinBUGS")), ylab=expression(paste(hat(R)[ijt], " INLA")))
  abline(a=0, b=1, col=selected_colors[1])

  plot(risks.INLA[[i]][,2], risks.WinBUGS[[i]][,2], col=selected_colors[2], cex.lab=1, main="Dowry deaths",
       xlab=expression(paste(hat(R)[ijt], " WinBUGS")), ylab=expression(paste(hat(R)[ijt], " INLA")))
  abline(a=0, b=1, col=selected_colors[2])
}

mtext("iCAR spatial priors", outer=TRUE,  cex=0.9, line=-1.5, font=2)
mtext("pCAR spatial priors", outer=TRUE,  cex=0.9, line=-25, font=2)
mtext("lCAR spatial priors", outer=TRUE,  cex=0.9, line=-50, font=2)

dev.off()


####################################
##  Figure 4. Posterior mean of the district-specific spatial risk, exp(theta_ij),
##            and the exceedence probabilities,  i.e., P(exp(theta_ij) > 1|O), for
##            rapeand dowry deaths
####################################
Model <- MODELS.inla.icar$TypeII
carto <- carto_UP

## Posterior means of district-specific spatial risks ##
spatial <- matrix(unlist(lapply(Model$marginals.random$idx, function(x) inla.emarginal(exp,x))),S,J,byrow=F)

inf <- min(spatial)
top <- max(spatial)
values <- c(round(seq(inf,1,length.out=5),2), round(seq(1,top,length.out=5),2)[-1])

paleta <- brewer.pal(8,"YlOrRd")
carto$spatial.crime1 <- spatial[,1]
carto$spatial.crime2 <- spatial[,2]

Map.spatial <- tm_shape(carto) + 
  tm_polygons(col=c("spatial.crime1","spatial.crime2"), palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(panel.labels=c("Rapes","Dowry deaths"), legend.position=c("right","top"))


## Posterior exceedence probabilities ##
probs <- matrix(unlist(lapply(Model$marginals.random$idx, function(x){1-inla.pmarginal(0, x)})),S,J,byrow=F)

paleta <- brewer.pal(5,"PuBu")
values <- c(0,0.1,0.2,0.8,0.9,1)

carto$prob.crime1 <- probs[,1]
carto$prob.crime2 <- probs[,2]

Map.prob <- tm_shape(carto) + 
  tm_polygons(col=c("prob.crime1","prob.crime2"), palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) + 
  tm_layout(panel.labels=c("Rapes","Dowry deaths"), legend.position=c("right","top"))


## File: figure_4.pdf
tmap_save(tmap_arrange(Map.spatial, Map.prob, ncol=2), file="figures/figure_4.pdf")


####################################
##  Figure 5. Temporal pattern of incidence risks (posterior means of exp(gamma_tj))
##            for rape and dowry deaths in Uttar Pradesh
####################################
Model <- MODELS.inla.icar$TypeII

## Posterior means of year-specific temporal risks ##
temporal <- matrix(unlist(lapply(Model$marginals.random$idy, function(x) inla.emarginal(exp,x))),T,J,byrow=F)

aux <- lapply(Model$marginals.random$idy, function(x) inla.tmarginal(exp,x))
q1 <- matrix(unlist(lapply(aux, function(x) inla.qmarginal(0.025,x))),T,J,byrow=F) 
q2 <- matrix(unlist(lapply(aux, function(x) inla.qmarginal(0.975,x))),T,J,byrow=F) 

inf <- min(q1)-0.05
top <- max(q2)+0.05
selected_colors <- c(rgb(154,192,205,alpha=150, maxColorValue=255),
                     rgb(69,139,116,alpha=150, maxColorValue=255))

## File: figure_5.pdf
graphics.off()

pdf("figures/figure_5.pdf", height=5, width=7.5, onefile=FALSE)
plot(range(x), c(inf, top), type="n", xlab="Year", ylab="",
     xaxt="n", cex.lab=1, cex.axis=1, cex.main=1, main=NULL)
title(ylab=expression(exp(gamma[t])), line=2.5, cex.lab=1.1)
axis(1, at=seq(1,T), labels=seq(t.from,t.to), las=0, cex.axis=1)

x <- seq(1,T)
for(i in seq(J)){
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(q1[,i], tail(q2[,i], 1), rev(q2[,i]), q1[1,i])
  polygon(X.Vec, Y.Vec, col = selected_colors[i], border = NA)
  lines(temporal[,i])
  abline(h=1,lty=2)
}
legend("topleft", inset=.02, legend=c("Rapes","Dowry deaths"), fill=selected_colors, horiz=FALSE, cex=1,box.lty=0)
dev.off()


####################################
##  Figure 6. Map of estimated incidence risks for rape (top) and posterior probabilities
##            that the relative risk is greater than 1 (P(Ritj>1|O)) (bottom) in Uttar 
##            Pradesh
##  Figure 7. Map of estimated incidence risks for dowry deaths (top) and posterior 
##            probabilities that the relative risk isgreater than one (P(Ritj > 1|O)) 
##            (bottom) in Uttar Pradesh
####################################
Model <- MODELS.inla.icar$TypeII

## Maps of posterior mean estimates of relative risks ##
risks <- matrix(Model$summary.fitted.values$mean,S*T,J,byrow=F)

inf <- min(risks)
top <- max(risks)
values <- c(round(seq(inf,1,length.out=5),2), round(seq(1,top,length.out=5),2)[-1])

paleta <- brewer.pal(8,"YlOrRd")

carto.rapes <- carto_UP
carto.dowry <- carto_UP
for(i in seq(T)){
  carto.rapes$var <- matrix(risks[,1],S,T,byrow=F)[,i]
  carto.dowry$var <- matrix(risks[,2],S,T,byrow=F)[,i]
  
  names(carto.rapes)[ncol(carto.rapes)] <- paste("Year",i,sep=".")
  names(carto.dowry)[ncol(carto.dowry)] <- paste("Year",i,sep=".")
}

Map.risk1 <- tm_shape(carto.rapes) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Rapes", main.title.position="center", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=3, ncol=5)

Map.risk2 <- tm_shape(carto.dowry) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Dowry deaths", main.title.position="center", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=3, ncol=5)


## Maps of posterior exceedence probabilities ##
probs <- matrix(1-Model$summary.fitted.values$`1 cdf`,S*T,J,byrow=F)

paleta <- brewer.pal(6,"Blues")[-1]
values <- c(0,0.1,0.2,0.8,0.9,1)

carto.rapes <- carto_UP
carto.dowry <- carto_UP
for(i in seq(T)){
  carto.rapes$var <- matrix(probs[,1],S,T,byrow=F)[,i]
  carto.dowry$var <- matrix(probs[,2],S,T,byrow=F)[,i]
  
  names(carto.rapes)[ncol(carto.rapes)] <- paste("Year",i,sep=".")
  names(carto.dowry)[ncol(carto.dowry)] <- paste("Year",i,sep=".")
}

Map.prob1 <- tm_shape(carto.rapes) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Rapes", main.title.position="center", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=3, ncol=5)

Map.prob2 <- tm_shape(carto.dowry) +
  tm_polygons(col=paste("Year",1:T,sep="."), palette=paleta, title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Dowry deaths", main.title.position="center", panel.label.size=1.5,
            panel.labels=seq(t.from,t.to), legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=3, ncol=5)


## File: figure_6.pdf
tmap_save(tmap_arrange(Map.risk1, Map.prob1, nrow=2), file="figures/figure_6.pdf")

## File: figure_7.pdf
tmap_save(tmap_arrange(Map.risk2, Map.prob2, nrow=2), file="figures/figure_7.pdf")


####################################
##  Figure 8. Temporal evolution of final risk estimates for rape and dowry deaths
##            in some districts in Uttar Pradesh:  Ghazlabad, Kheri, Mainpuri, 
##            Sant Kabir Nagar, and Varanasi
####################################
Model <- MODELS.inla.icar$TypeII

## Selected districts ##
id.area<- c(2,28,43,49,61,70)

RR.rape <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==1,]
RR.dowry <- Model$summary.fitted.values[(Model$.args$data$ID.area %in% id.area) & Model$.args$data$ID.disease==2,]

aux <- lapply(list(rape=RR.rape,dowry=RR.dowry), function(x){
  mean <- data.frame(mean=matrix(x[,"mean"],ncol=T,byrow=F))
  q1 <- data.frame(mean=matrix(x[,"0.025quant"],ncol=T,byrow=F))
  q2 <- data.frame(mean=matrix(x[,"0.975quant"],ncol=T,byrow=F))
  colnames(mean) <- colnames(q1) <- colnames(q1) <- seq(1,T)
  
  list(mean=mean, q1=q1, q2=q2)
})


inf <- min(unlist(aux))-0.05
top <- max(unlist(aux))+1.5
selected_colors <- c(rgb(154,192,205,alpha=150, maxColorValue=255),
                     rgb(69,139,116,alpha=150, maxColorValue=255))


## File: figure_8.pdf
graphics.off()

pdf("figures/figure_8.pdf", height=12, width=12, onefile=FALSE)
par(mfrow=c(3,2))

x <- seq(1,T)
main.title <- unique(paste(Model$.args$data$Area," (ID area ", Model$.args$data$ID.area,")",sep="")[Model$.args$data$ID.area %in% id.area])
  
for(i in seq(length(id.area))){
  
  plot(range(x), c(inf, top), type="n", xlab="Year", ylab="",
       xaxt="n", cex.lab=1.5, cex.axis=1.5, cex.main=2, main=main.title[i])
  title(ylab=expression(R[ijt]), line=2.0, cex.lab=1.7)
  axis(1, at=seq(1,T,2), labels=seq(t.from,t.to,2), las=0, cex.axis=1)
  
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(aux$rape$q1[i,], tail(aux$rape$q2[i,1]), rev(aux$rape$q2[i,]), aux$rape$q1[i,1])
  polygon(X.Vec, Y.Vec, col=selected_colors[1], border = NA)
  lines(1:T, aux$rape$mean[i,], lwd=1.5)
    
  X.Vec <- c(x, tail(x, 1), rev(x), x[1])
  Y.Vec <- c(aux$dowry$q1[i,], tail(aux$dowry$q2[i,1]), rev(aux$dowry$q2[i,]), aux$dowry$q1[i,1])
  polygon(X.Vec, Y.Vec, col=selected_colors[2], border = NA)
  lines(1:T, aux$dowry$mean[i,], lwd=1.5)
    
  abline(h=1,lty=2)
  legend("topleft", inset=.02, c("Rapes","Dowry deaths"), fill=selected_colors, horiz=FALSE, cex=2.2, box.lty=0)
}

dev.off() 
  

##################################################################################
################                    Tables                        ################
##################################################################################

####################################
##  Table 1.  Descriptive statistics. Minimun, first quartile (q1), mean, third
##            quartile (q3), maximun, standard desviation (sd), and coefficient of
##            variation of the number of rapes and dowry deaths in the districts 
##            of Uttar Pradesh per year
####################################

table.1 <- matrix(NA, nrow=T, ncol=15)
colnames(table.1) <- c("Year", "min", "q1", "mean", "q3", "max", "sd", "|cv|",
                       "min", "q1", "mean", "q3", "max", "sd", "|cv|")

k <- 1
for(i in seq(t.from,t.to)){
  table.1[k,]<- c(i, summary(data$obs[data$year==i & data$Crime==1])[-3],
                    sd(data$obs[data$year==i & data$Crime==1]), NA,
                    summary(data$obs[data$year==i & data$Crime==2])[-3],
                    sd(data$obs[data$year==i & data$Crime==2]), NA)
  k <- k+1
}
table.1[,8] <- abs(table.1[,7]/ table.1[,4])
table.1[,15]<- abs(table.1[,14]/ table.1[,11])
table.1 <- round(table.1, 1)
print(table.1)

## latex
latex_table<-xtable::xtable(table.1,
                            caption="Descriptive statistics. Minimun, first quartile ($q_1$), mean, third quartile ($q_3$), maximun, standard desviation (sd), and coefficient of variation of the number of rapes and dowry deaths in the districts of Uttar Pradesh per year.",
                            label="t_crime", digits=c(1),
                            display=c("d", "d", "d","f","f","f","d","f","f", "d","f","f","f","d","f","f"))
xtable::print.xtable(latex_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))


####################################
##  Table 2.  Correlations between spatial (by year) and temporal patterns
##            (by district) of rape and dowry deaths 
####################################
cor.spat <- c()
for(i in unique(data$year)){
  cor.spat <- c(cor.spat, cor(data[data$year==i & data$Crime==1,"SMR"],data[data$year==i & data$Crime==2,"SMR"]))
}

cor.temp <- c()
for(i in unique(data$dist)){
  cor.temp <- c(cor.temp, cor(data[data$dist==i & data$Crime==1,"SMR"],data[data$dist==i & data$Crime==2,"SMR"]))
}

table.2 <- rbind(c(summary(cor.spat), sd(cor.spat), abs(sd(cor.spat)/mean(cor.spat))), 
                c(summary(cor.temp), sd(cor.temp), abs(sd(cor.temp)/mean(cor.temp))))
colnames(table.2) <- c("min", "q1", "median", "mean", "q3", "max", "sd", "|cv|")
table.2 <- as.data.frame(table.2)
table.2$Correlation <- c("spatial patterns", "temporal trends")
table.2 <- table.2[,c("Correlation", "min", "q1", "median", "mean", "q3", "max", "sd", "|cv|")]
print(table.2)

## latex
latex_table<-xtable::xtable(table.2,
                            caption="Correlations between spatial (by year) and temporal patterns (by district) of rape and dowry deaths.",
                            label="t_corre", digits=c(3))
xtable::print.xtable(latex_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))


####################################
##  Table 3.  Model selection criteria, DIC, WAIC and LS for INLA models
####################################
table3.iCAR <- do.call(rbind,lapply(MODELS.inla.icar, function(x) data.frame(DIC=x$dic$dic, WAIC=x$waic$waic, LS=-sum(log(x$cpo$cpo)))))
rownames(table3.iCAR) <- paste("iCAR",rownames(table3.iCAR),sep=".")

table3.LCAR <- do.call(rbind,lapply(MODELS.inla.lcar, function(x) data.frame(DIC=x$dic$dic, WAIC=x$waic$waic, LS=-sum(log(x$cpo$cpo)))))
rownames(table3.LCAR) <- paste("LCAR",rownames(table3.LCAR),sep=".")

table3.pCAR <- do.call(rbind,lapply(MODELS.inla.pcar, function(x) data.frame(DIC=x$dic$dic, WAIC=x$waic$waic, LS=-sum(log(x$cpo$cpo)))))
rownames(table3.pCAR) <- paste("pCAR",rownames(table3.pCAR),sep=".")

table.3 <- rbind(table3.iCAR,table3.LCAR,table3.pCAR)
print(table.3)

## latex
x_table<-xtable::xtable(table.3, 
                        caption="Model selection criteria, DIC, WAIC and LS, for different models.", 
                        label="t_dic", digits=c(3))
xtable::print.xtable(x_table, include.rownames = T, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))


##################################################################################
################                  Appendix                        ################
##################################################################################

####################################
##  Table A.1.  District identifiers (ID) of Uttar Pradesh
####################################
ID <- data.frame(dist=carto_UP$dist, ID_area=carto_UP$ID_area)
table.a1 <- cbind(ID[1:24,c("ID_area", "dist")], ID[25:48,c("ID_area", "dist")],
                 ID[49:72,c("ID_area", "dist")])
colnames(table.a1) <- c("ID", "Dist","ID", "Dist", "ID", "Dist")
print(table.a1)

## latex
x_table <-xtable::xtable(table.a1, 
                        caption="District identifiers (ID) of Uttar Pradesh.",
                        label="table_a1", digits=c(3))
xtable::print.xtable(x_table, include.rownames = FALSE, comment=FALSE, 
                     caption.placement = getOption("xtable.caption.placement", "top"))


####################################
##  Table A.2.  Posterior means, standard deviations, and 95% credible intervals 
##              for the crime-specific intercepts (alpha_j , j=1,2) of the models
##              with a spatio-temporal Type II interaction term
####################################
alpha.INLA <- lapply(Models.INLA, function(x) x$summary.fixed[,c(1,2,3,5)])
alpha.WinBUGS <- lapply(Models.WinBUGS, function(x) x$summary[grepl("mu",rownames(x$summary), fixed=TRUE),c(1,2,3,7)])

table.a2 <- mapply(function(x,y){
  aux <- rbind(x["I1",],y["I1",],x["I2",],y["I2",])
  rownames(aux) <- paste(rep(c("Rapes","Dowry_deaths"),each=2),c("INLA","MCMC"),sep=".")
  aux <- round(aux,3)
  aux
  }, x=alpha.INLA, y=alpha.WinBUGS, SIMPLIFY=F)

print(table.a2)

