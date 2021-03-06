##################################################################################################
## Spatial: FE (lcar); Temporal: FE (rw1); Spatio-temporal: ad
## FE.FE.ad.lcar
##################################################################################################
FE.FE.ad.lcar <- function(){
  ########################################
  for(k in 1:nnv){
    for(j in 1:Ndiseases){
      uneigh[j,k]<- Spatial[j, adj[k] ]
    }
  }
  ########################################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      for(l in 1:Nyears){
        ## Observed of the mean for each area, disease and time
        O[i,j,l] ~ dpois(lambda[i,j,l])
        ## Modeling of the mean for each area, disease and time
        log(lambda[i,j,l]) <- log(E[i,j,l]) + mu[j] + Theta[i,j] + Gam[l,j]
        ## Risk for each area, disease and time
        SMR[i,j,l] <- exp(mu[j] + Theta[i,j] + Gam[l,j])
        smr.prob[i,j,l]<- step(SMR[i,j,l]-1)
      }
      ## Spatial effects
      Espat[i,j] <- exp(Theta[i,j])
      espat.prob[i,j]<- step(Espat[i,j]-1)
    } 
    ## Temporal effects
    for(l in 1:Nyears){
      Etemp[l,j] <- exp(Gam[l,j])
      etemp.prob[l,j] <-  step(Etemp[l,j]-1)
    }
  }
  ####################################
  ## Prior distribution for the mean risk for all areas
  ####################################
  for (j in 1:Ndiseases){ mu[j] ~ dflat() }
  ####################################
  ## spatial 
  ####################################
  ##################
  ## Theta = Phi * M;  (IxJ) = (IxK)*(KxJ)
  ##################
  for (i in 1:Nareas){
    for (j in 1:Ndiseases){
      Theta[i,j] <- inprod2(tPhi[,i], M[,j])
    }
  }
  ##################
  ## Phi (IxK, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      tPhi[j,i] <-  c.Spatial[j,i]
      c.Spatial[j,i]<- Spatial[j,i] - mean(Spatial[j,])
      Spatial[j,i] ~ dnorm(umed[j,i],uvar[j,i])
      umed[j,i]<- gamma1[j] * (sum( uneigh[j, (off[i]+1):off[i+1] ])) / (1- gamma1[j] + gamma1[j] *num[i])
      uvar[j,i]<-  tau1[j]* (1- gamma1[j] +gamma1[j] *num[i])
    }
    gamma1[j] ~ dunif(0,1)
    tau1[j]<- 1  # lo fijamos en 1 para que M sea quien se ocupa de la variabilidad
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (i in 1:(1 * Ndiseases)) {
    for (j in 1:Ndiseases) {
      M[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlations
  ##################
  # Sigma_b =t(M)%*%M
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.s[i,j] <- inprod2(M[,i], M[,j])
      Corre.s[i,j] <- Sigma.s[i,j]/(pow(Sigma.s[i,i],0.5)*pow(Sigma.s[j,j],0.5))
    }
  }
  ####################################
  # temporal
  ####################################
  ##################
  # Gamma = Phig * Mg; (txJ) = (txK)*(KxJ)
  ##################
  for(l in 1:Nyears){
    for(j in 1:Ndiseases){
      Gam[l,j]<-  inprod2(tPhiG[,l], Mg[,j])
    }
  }
  ##################
  ## PhiG (txK, K=2J)
  ##################
  for (j in 1:(Ndiseases)){
    Temporal[j, 1:Nyears] ~ car.normal(adjt[], weightst[], numt[],1)
    for (l in 1:Nyears){
      tPhiG[j,l] <-Temporal[j,l]
    }
  }
  ##################
  ## Mg-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:(Ndiseases)){
      Mg[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlations
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.t[i,j] <- inprod2(Mg[,i], Mg[,j])
      Corre.t[i,j] <- Sigma.t[i,j]/(pow(Sigma.t[i,i],0.5)*pow(Sigma.t[j,j],0.5))
    }
  }
  ######################################
  ######################################
}
##################################################################################################
## Spatial: FE (lcar); Temporal: FE (rw1); Spatio-temporal: Type I
## FE.FE.t1.lcar
############################################
FE.FE.t1.lcar <- function(){
  ########################################
  for(k in 1:nnv){
    for(j in 1:Ndiseases){
      uneigh[j,k]<- Spatial[j, adj[k] ]
    }
  }
  ########################################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      for(l in 1:Nyears){
        ## Observed of the mean for each area, disease and time
        O[i,j,l] ~ dpois(lambda[i,j,l])
        ## Modeling of the mean for each area, disease and time
        log(lambda[i,j,l]) <- log(E[i,j,l]) + mu[j] + Theta[i,j] + Gam[l,j] + sdZet[j] * Zet[i,j,l]
        ## Risk for each area, disease and time
        SMR[i,j,l] <- exp(mu[j] + Theta[i,j] + Gam[l,j] +  sdZet[j] * Zet[i,j,l])
        smr.prob[i,j,l]<- step(SMR[i,j,l]-1)
        ## Spatio-temporal effect
        Eint[i,j,l] <- exp(sdZet[j] * Zet[i,j,l])
        eint.prob[i,j,l]<- step(Eint[i,j,l]-1)
      }
      ## Spatial effects
      Espat[i,j] <- exp(Theta[i,j])
      espat.prob[i,j]<- step(Espat[i,j]-1)
    } 
    ## Temporal effects
    for(l in 1:Nyears){
      Etemp[l,j] <- exp(Gam[l,j])
      etemp.prob[l,j] <-  step(Etemp[l,j]-1)
    }
  }
  ####################################
  ## Prior distribution for the mean risk for all areas
  ####################################
  for (j in 1:Ndiseases){ mu[j] ~ dflat() }
  ####################################
  ## spatial 
  ####################################
  ##################
  ## Theta = Phi * M;  (IxJ) = (IxK)*(KxJ)
  ##################
  for (i in 1:Nareas){
    for (j in 1:Ndiseases){
      Theta[i,j] <- inprod2(tPhi[,i], M[,j])
    }
  }
  ##################
  ## Phi (IxK, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      tPhi[j,i] <-  c.Spatial[j,i]
      c.Spatial[j,i]<- Spatial[j,i] - mean(Spatial[j,])
      Spatial[j,i] ~ dnorm(umed[j,i],uvar[j,i])
      umed[j,i]<- gamma1[j] * (sum( uneigh[j, (off[i]+1):off[i+1] ])) / (1- gamma1[j] + gamma1[j] *num[i])
      uvar[j,i]<-  tau1[j]* (1- gamma1[j] +gamma1[j] *num[i])
    }
    gamma1[j] ~ dunif(0,1)
    tau1[j]<- 1  # lo fijamos en 1 para que M sea quien se ocupa de la variabilidad
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (i in 1:(1 * Ndiseases)) {
    for (j in 1:Ndiseases) {
      M[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlations
  ##################
  # Sigma_b =t(M)%*%M
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.s[i,j] <- inprod2(M[,i], M[,j])
      Corre.s[i,j] <- Sigma.s[i,j]/(pow(Sigma.s[i,i],0.5)*pow(Sigma.s[j,j],0.5))
    }
  }
  ####################################
  # temporal
  ####################################
  ##################
  # Gamma = Phig * Mg; (txJ) = (txK)*(KxJ)
  ##################
  for(l in 1:Nyears){
    for(j in 1:Ndiseases){
      Gam[l,j]<-  inprod2(tPhiG[,l], Mg[,j])
    }
  }
  ##################
  ## PhiG (txK, K=2J)
  ##################
  for (j in 1:(Ndiseases)){
    Temporal[j, 1:Nyears] ~ car.normal(adjt[], weightst[], numt[],1)
    for (l in 1:Nyears){
      tPhiG[j,l] <-Temporal[j,l]
    }
  }
  ##################
  ## Mg-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:(Ndiseases)){
      Mg[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlations
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.t[i,j] <- inprod2(Mg[,i], Mg[,j])
      Corre.t[i,j] <- Sigma.t[i,j]/(pow(Sigma.t[i,i],0.5)*pow(Sigma.t[j,j],0.5))
    }
  }
  ######################################
  # spatio-temporal (M-model)
  ######################################
  ##################
  # Zeta = Phiz * Mz  #(IxT) = (IxK)*(KxT) K=T
  ##################
  for(j in 1:Ndiseases){
    for(i in 1:Nareas){
      for(l in 1:Nyears){
        Zet.aux[i,j,l] ~ dnorm(0,1)
        Zet[i,j,l] <- Zet.aux[i,j,l]
      }
    }
  }
  ##################
  # Prior distribution for the standard deviations
  ##################
  for(j in 1:Ndiseases){
    sdZet[j] ~ dunif(0,100)
  }
  ######################################
  ######################################
}
##################################################################################################
# Spatial: FE (lcar); Temporal: FE (rw1); Spatio-temporal: Type II
# FE.FE.t2.lcar
##################################################################################################
FE.FE.t2.lcar <- function(){
  ########################################
  for(k in 1:nnv){
    for(j in 1:Ndiseases){
      uneigh[j,k]<- Spatial[j, adj[k] ]
    }
  }
  ########################################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      for(l in 1:Nyears){
        ## Observed of the mean for each area, disease and time
        O[i,j,l] ~ dpois(lambda[i,j,l])
        ## Modeling of the mean for each area, disease and time
        log(lambda[i,j,l]) <- log(E[i,j,l]) + mu[j] + Theta[i,j] + Gam[l,j] + sdZet[j] * Zet[i,j,l]
        ## Risk for each area, disease and time
        SMR[i,j,l] <- exp(mu[j] + Theta[i,j] + Gam[l,j] + sdZet[j] * Zet[i,j,l])
        smr.prob[i,j,l]<- step(SMR[i,j,l]-1)
        ## Spatio-temporal effect
        Eint[i,j,l] <- exp(sdZet[j] * Zet[i,j,l])
        eint.prob[i,j,l]<- step(Eint[i,j,l]-1)
      }
      ## Spatial effect
      Espat[i,j] <- exp(Theta[i,j])
      espat.prob[i,j]<- step(Espat[i,j]-1)
    }
    ## Temporal effect
    for(l in 1:Nyears){
      Etemp[l,j] <- exp(Gam[l,j])
      etemp.prob[l,j] <-  step(Etemp[l,j]-1)
    }
  }
  ####################################
  ## Prior distribution for the mean risk for all municipalities
  ####################################
  for (j in 1:Ndiseases){ mu[j] ~ dflat() }
  ####################################
  ## spatial 
  ####################################
  ##################
  ## Theta = Phi * M;  (IxJ) = (IxK)*(KxJ)
  ##################
  for (i in 1:Nareas){
    for (j in 1:Ndiseases){
      Theta[i,j] <- inprod2(tPhi[,i], M[,j])
    }
  }
  ##################
  ## Phi (IxK, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      tPhi[j,i] <-  c.Spatial[j,i]
      c.Spatial[j,i]<- Spatial[j,i] - mean(Spatial[j,])
      Spatial[j,i] ~ dnorm(umed[j,i],uvar[j,i])
      umed[j,i]<- gamma1[j] * (sum( uneigh[j, (off[i]+1):off[i+1] ])) / (1- gamma1[j] + gamma1[j] *num[i])
      uvar[j,i]<-  tau1[j]* (1- gamma1[j] +gamma1[j] *num[i])
    }
    gamma1[j] ~ dunif(0,1)
    tau1[j]<- 1  # lo fijamos en 1 para que M sea quien se ocupa de la variabilidad
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (i in 1:(1 * Ndiseases)) {
    for (j in 1:Ndiseases) {
      M[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlations
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.s[i,j] <- inprod2(M[,i], M[,j])
      Corre.s[i,j] <- Sigma.s[i,j]/(pow(Sigma.s[i,i],0.5)*pow(Sigma.s[j,j],0.5))
    }
  }
  ####################################
  ## temporal
  ####################################
  ##################
  ## Gamma = Phig * Mg;  (txJ) = (txK)*(KxJ)
  ##################
  for(l in 1:Nyears){
    for(j in 1:Ndiseases){
      Gam[l,j]<-  inprod2(tPhiG[,l], Mg[,j])
    }
  }
  ##################
  ## PhiG (txK, K=2J)
  ##################
  for (j in 1:(Ndiseases)){
    Temporal[j, 1:Nyears] ~ car.normal(adjt[], weightst[], numt[],1) # structured
    for (l in 1:Nyears){
      tPhiG[j,l] <-Temporal[j,l]
    }
  }
  ##################
  ## Mg-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:(Ndiseases)){
      Mg[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlations
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.t[i,j] <- inprod2(Mg[,i], Mg[,j])
      Corre.t[i,j] <- Sigma.t[i,j]/(pow(Sigma.t[i,i],0.5)*pow(Sigma.t[j,j],0.5))
    }
  }
  ######################################
  ## spatio-temporal (M-model)
  ######################################
  ##################
  ## Zet = Phi.z * M; (IxT) = (IxK)*(KxT)
  ##################
  for(j in 1:Ndiseases){
    for (i in 1:Nareas){
      for (l in 1:Nyears){
        Zet[i,j,l] <- Mz[i,j,l]
      }
    }
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for(j in 1:Ndiseases){
    for (i in 1:(1*Nareas)){
      for (l in 1:Nyears){
        Mz[i,j,l]<- Mzaux[i,j,l]
      }
    }
    ## Mzaux
    for (i in 1:(Nareas)){
      Temporal.z[i,j,1:Nyears] ~ car.normal(adj.zt[], weights.zt[], num.zt[],1) # structured
      for (l in 1:Nyears){
        Mzaux[i,j,l] <- Temporal.z[i,j,l]
      }
    }
  }
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  for(j in 1:Ndiseases){
    sdZet[j] ~ dunif(0,100)
  }
  ######################################
  ######################################
}
##################################################################################################
## Spatial: FE (lcar); Temporal: FE (rw1); Spatio-temporal: Type III
## FE.FE.t3.lcar
##################################################################################################
FE.FE.t3.lcar <- function(){
  ########################################
  for(k in 1:nnv){
    for(j in 1:Ndiseases){
      uneigh[j,k]<- Spatial[j, adj[k] ]
    }
  }
  ########################################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      for(l in 1:Nyears){
        ## Observed of the mean for each area, disease and time
        O[i,j,l] ~ dpois(lambda[i,j,l])
        ## Modeling of the mean for each area, disease and time
        log(lambda[i,j,l]) <- log(E[i,j,l]) + mu[j] + Theta[i,j] + Gam[l,j] + sdZet[j] * Zet[i,j,l]
        ## Risk for each area, disease and time
        SMR[i,j,l] <- exp(mu[j] + Theta[i,j] + Gam[l,j] + sdZet[j] * Zet[i,j,l])
        smr.prob[i,j,l]<- step(SMR[i,j,l]-1)
        ## Spatio-temporal effect
        Eint[i,j,l] <- exp(sdZet[j] * Zet[i,j,l])
        eint.prob[i,j,l]<- step(Eint[i,j,l]-1)
      }
      ## Spatial effect
      Espat[i,j] <- exp(Theta[i,j])
      espat.prob[i,j]<- step(Espat[i,j]-1)
    }
    ## Temporal effect
    for(l in 1:Nyears){
      Etemp[l,j] <- exp(Gam[l,j])
      etemp.prob[l,j] <-  step(Etemp[l,j]-1)
    }
  }
  ####################################
  ## Prior distribution for the mean risk for all areas
  ####################################
  for (j in 1:Ndiseases){ mu[j] ~ dflat() }
  ####################################
  ## spatial 
  ####################################
  ##################
  ## Theta = Phi * M; (IxJ) = (IxK)*(KxJ)
  ##################
  for (i in 1:Nareas){
    for (j in 1:Ndiseases){
      Theta[i,j] <- inprod2(tPhi[,i], M[,j])
    }
  }
  ##################
  ## Phi (IxK, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      tPhi[j,i] <-  c.Spatial[j,i]
      c.Spatial[j,i]<- Spatial[j,i] - mean(Spatial[j,])
      Spatial[j,i] ~ dnorm(umed[j,i],uvar[j,i])
      umed[j,i]<- gamma1[j] * (sum( uneigh[j, (off[i]+1):off[i+1] ])) / (1- gamma1[j] + gamma1[j] *num[i])
      uvar[j,i]<-  tau1[j]* (1- gamma1[j] +gamma1[j] *num[i])
    }
    gamma1[j] ~ dunif(0,1)
    tau1[j]<- 1  # lo fijamos en 1 para que M sea quien se ocupa de la variabilidad
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (i in 1:(1 * Ndiseases)) {
    for (j in 1:Ndiseases) {
      M[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlation
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.s[i,j] <- inprod2(M[,i], M[,j])
      Corre.s[i,j] <- Sigma.s[i,j]/(pow(Sigma.s[i,i],0.5)*pow(Sigma.s[j,j],0.5))
    }
  }
  ####################################
  ## temporal
  ####################################
  ##################
  ## Gamma = Phig * Mg; (txJ) = (txK)*(KxJ)
  ##################
  for(l in 1:Nyears){
    for(j in 1:Ndiseases){
      Gam[l,j]<-  inprod2(tPhiG[,l], Mg[,j])
    }
  }
  ##################
  ## PhiG (txK, K=2J)
  ##################
  for (j in 1:(Ndiseases)){
    Temporal[j, 1:Nyears] ~ car.normal(adjt[], weightst[], numt[],1) # structured
    for (l in 1:Nyears){
      tPhiG[j,l] <-Temporal[j,l]
    }
  }
  ##################
  ## Mg-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:(Ndiseases)){
      Mg[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlation
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.t[i,j] <- inprod2(Mg[,i], Mg[,j])
      Corre.t[i,j] <- Sigma.t[i,j]/(pow(Sigma.t[i,i],0.5)*pow(Sigma.t[j,j],0.5))
    }
  }
  ######################################
  ## spatio-temporal (M-model)
  ######################################
  ##################
  ## Zet = Phi.z * M; (IxT) = (IxK)*(KxT)
  ##################
  for(j in 1:Ndiseases){
    for (i in 1:Nareas){
      for (l in 1:Nyears){
        Zet[i,j,l] <- tPhi.z[l,j,i]
      }
    }
  }
  ##################
  # Phi.z (IxK, K=T)
  ##################
  for(j in 1:Ndiseases){
    for (l in 1:Nyears){
      Spatial.z[l, j, 1:Nareas] ~ car.normal(adj.zs[], weights.zs[], num.zs[],1) # structured
      for (i in 1:Nareas){
        tPhi.z[l,j,i] <- Spatial.z[l,j,i]
      }
    }
  }
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  for(j in 1:Ndiseases){
    sdZet[j] ~ dunif(0,100)
  }
  ######################################
  ######################################
}
##################################################################################################
## Spatial: FE (lcar); Temporal: FE (rw1); Spatio-temporal: Type IV
## FE.FE.t4.lcar
##################################################################################################
FE.FE.t4.lcar <- function(){
  ########################################
  for(k in 1:nnv){
    for(j in 1:Ndiseases){
      uneigh[j,k]<- Spatial[j, adj[k] ]
    }
  }
  ########################################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      for(l in 1:Nyears){
        ## Observed of the mean for each area, disease and time
        O[i,j,l] ~ dpois(lambda[i,j,l])
        ## Modeling of the mean for each area, disease and time
        log(lambda[i,j,l]) <- log(E[i,j,l]) + mu[j] + Theta[i,j] + Gam[l,j] + sdZet[j] * Zet[i,j,l]
        ## Risk for each area, disease and time
        SMR[i,j,l] <- exp(mu[j] + Theta[i,j] + Gam[l,j] + sdZet[j] * Zet[i,j,l])
        smr.prob[i,j,l]<- step(SMR[i,j,l]-1)
        ## Spatio-temporal effect
        Eint[i,j,l] <- exp(sdZet[j] * Zet[i,j,l])
        eint.prob[i,j,l]<- step(Eint[i,j,l]-1)
      }
      ## Spatial effect
      Espat[i,j] <- exp(Theta[i,j])
      espat.prob[i,j]<- step(Espat[i,j]-1)
    }
    ## Temporal effect
    for(l in 1:Nyears){
      Etemp[l,j] <- exp(Gam[l,j])
      etemp.prob[l,j] <-  step(Etemp[l,j]-1)
    }
  }
  ####################################
  ## Prior distribution for the mean risk for all areas
  ####################################
  for (j in 1:Ndiseases){ mu[j] ~ dflat() }
  ####################################
  ## spatial 
  ####################################
  ##################
  ## Theta = Phi * M; (IxJ) = (IxK)*(KxJ)
  ##################
  for (i in 1:Nareas){
    for (j in 1:Ndiseases){
      Theta[i,j] <- inprod2(tPhi[,i], M[,j])
    }
  }
  ##################
  ## Phi (IxK, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Nareas){
      tPhi[j,i] <-  c.Spatial[j,i]
      c.Spatial[j,i]<- Spatial[j,i] - mean(Spatial[j,])
      Spatial[j,i] ~ dnorm(umed[j,i],uvar[j,i])
      umed[j,i]<- gamma1[j] * (sum( uneigh[j, (off[i]+1):off[i+1] ])) / (1- gamma1[j] + gamma1[j] *num[i])
      uvar[j,i]<-  tau1[j]* (1- gamma1[j] +gamma1[j] *num[i])
    }
    gamma1[j] ~ dunif(0,1)
    tau1[j]<- 1  # lo fijamos en 1 para que M sea quien se ocupa de la variabilidad
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (i in 1:(1 * Ndiseases)) {
    for (j in 1:Ndiseases) {
      M[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlation
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.s[i,j] <- inprod2(M[,i], M[,j])
      Corre.s[i,j] <- Sigma.s[i,j]/(pow(Sigma.s[i,i],0.5)*pow(Sigma.s[j,j],0.5))
    }
  }
  ####################################
  ## temporal
  ####################################
  ##################
  ## Gamma = Phig * Mg; (txJ) = (txK)*(KxJ)
  ##################
  for(l in 1:Nyears){
    for(j in 1:Ndiseases){
      Gam[l,j]<-  inprod2(tPhiG[,l], Mg[,j])
    }
  }
  ##################
  ## PhiG (txK, K=2J)
  ##################
  for (j in 1:(Ndiseases)){
    Temporal[j, 1:Nyears] ~ car.normal(adjt[], weightst[], numt[],1)
    for (l in 1:Nyears){
      tPhiG[j,l] <-Temporal[j,l]
    }
  }
  ##################
  ## Mg-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:(Ndiseases)){
      Mg[i,j] ~ dflat()
    }
  }
  ##################
  ## Sigma=M'M and correlation
  ##################
  for(i in 1:Ndiseases){
    for(j in 1:Ndiseases){
      Sigma.t[i,j] <- inprod2(Mg[,i], Mg[,j])
      Corre.t[i,j] <- Sigma.t[i,j]/(pow(Sigma.t[i,i],0.5)*pow(Sigma.t[j,j],0.5))
    }
  }
  ######################################
  ## spatio-temporal (M-model)
  ######################################
  ##################
  ## Zet = Phi.z * M; (IxT) = (IxK)*(KxT)
  ##################
  for(j in 1:Ndiseases){
    for (i in 1:Nareas){
      for (l in 1:Nyears){
        Zet[i,j,l] <- inprod2(tPhi.z[,j,i], Mz[,l])
      }
    }
  }
  ##################
  # Phi.z (IxK, K=T)
  ##################
  for(j in 1:Ndiseases){
    for (l in 1:Nyears){
      Spatial.z[l,j, 1:Nareas] ~ car.normal(adj.zs[], weights.zs[], num.zs[],1) # structured
      for (i in 1:Nareas){
        tPhi.z[l,j,i] <- Spatial.z[l,j,i]
      }
    }
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (i in 1:(1*Nyears)){
    for (l in 1:Nyears){
      Mz[i,l]<- Mzz[i,l]
    }
  }
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  for(j in 1:Ndiseases){
    sdZet[j] ~ dunif(0,100)
  }
  ######################################
  ######################################
}
##################################################################################################
##################################################################################################