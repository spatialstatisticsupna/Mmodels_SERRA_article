##################################################################################################
## Spatial: NVA (bym); Temporal: NVA (rw1); Spatio-temporal: ad
## NVA.NVA.ad.bym
##################################################################################################
NVA.NVA.ad.bym <- function(){
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
    Spatial[j, 1:Nareas] ~ car.normal(adj[], weights[], num[],1) # structured
    for (i in 1:Nareas){
      Het[j, i] ~ dnorm(0, 1) # unstructured
      c.Het[j,i]<- Het[j,i] - mean(Het[j,])
      tPhi[j,i] <- Spatial[j,i]
    }
  }
  for (j in (Ndiseases + 1):(2 * Ndiseases)){
    for (i in 1:Nareas){
      tPhi[j,i] <- c.Het[(j - Ndiseases),i]
    }
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Ndiseases){
      M[i,j] ~ dnorm(0, prec.sp)
    }
    for (i in (Ndiseases+1):(2*Ndiseases)){
      M[i,j] ~ dnorm(0, prec.het)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  prec.sp <- pow(sdstruct.sp, -2)
  sdstruct.sp ~ dunif(0, 100)
  prec.het <- pow(sdstruct.het, -2)
  sdstruct.het ~ dunif(0, 100)
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
  # PhiG (txK, K=2J)
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
  for (i in 1:(Ndiseases)){
    for (j in 1:Ndiseases){
      Mg[i,j] ~ dnorm(0, precg)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  precg <- pow(sdstructg, -2)
  sdstructg ~ dunif(0, 100)
  ######################################
  ######################################
}
##################################################################################################
## Spatial: NVA (bym); Temporal: NVA (rw1); Spatio-temporal: Type I
## NVA.NVA.t1.bym
##################################################################################################
NVA.NVA.t1.bym <- function(){
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
    Spatial[j, 1:Nareas] ~ car.normal(adj[], weights[], num[],1) # structured
    for (i in 1:Nareas){
      Het[j, i] ~ dnorm(0, 1) # unstructured
      c.Het[j,i]<- Het[j,i] - mean(Het[j,])
      tPhi[j,i] <- Spatial[j,i]
    }
  }
  for (j in (Ndiseases + 1):(2 * Ndiseases)){
    for (i in 1:Nareas){
      tPhi[j,i] <- c.Het[(j - Ndiseases),i]
    }
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Ndiseases){
      M[i,j] ~ dnorm(0, prec.sp)
    }
    for (i in (Ndiseases+1):(2*Ndiseases)){
      M[i,j] ~ dnorm(0, prec.het)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  prec.sp <- pow(sdstruct.sp, -2)
  sdstruct.sp ~ dunif(0, 100)
  prec.het <- pow(sdstruct.het, -2)
  sdstruct.het ~ dunif(0, 100)
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
  # PhiG (txK, K=2J)
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
  for (i in 1:(Ndiseases)){
    for (j in 1:Ndiseases){
      Mg[i,j] ~ dnorm(0, precg)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  precg <- pow(sdstructg, -2)
  sdstructg ~ dunif(0, 100)
  ######################################
  ## spatio-temporal (M-model)
  ######################################
  ##################
  ## Zeta = Phiz * Mz  #(IxT) = (IxK)*(KxT) K=T
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
  ## Prior distribution for the standard deviations
  ##################
  for(j in 1:Ndiseases){
    sdZet[j] ~ dunif(0,100)
  }
  ######################################
  ######################################
}
##################################################################################################
## Spatial: NVA (bym); Temporal: NVA (rw1); Spatio-temporal: Type II
## NVA.NVA.t2.bym
##################################################################################################
NVA.NVA.t2.bym <- function(){
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
    Spatial[j, 1:Nareas] ~ car.normal(adj[], weights[], num[],1) # structured
    for (i in 1:Nareas){
      Het[j, i] ~ dnorm(0, 1) # unstructured
      c.Het[j,i]<- Het[j,i] - mean(Het[j,])
      tPhi[j,i] <- Spatial[j,i]
    }
  }
  for (j in (Ndiseases + 1):(2 * Ndiseases)){
    for (i in 1:Nareas){
      tPhi[j,i] <- c.Het[(j - Ndiseases),i]
    }
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Ndiseases){
      M[i,j] ~ dnorm(0, prec.sp)
    }
    for (i in (Ndiseases+1):(2*Ndiseases)){
      M[i,j] ~ dnorm(0, prec.het)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  prec.sp <- pow(sdstruct.sp, -2)
  sdstruct.sp ~ dunif(0,100)
  prec.het <- pow(sdstruct.het, -2)
  sdstruct.het ~ dunif(0,100)
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
  for (i in 1:(Ndiseases)){
    for (j in 1:Ndiseases){
      Mg[i,j] ~ dnorm(0, precg)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  precg <- pow(sdstructg, -2)
  sdstructg ~ dunif(0,100)
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
## Spatial: NVA (bym); Temporal: NVA (rw1); Spatio-temporal: Type III
## NVA.NVA.t3.bym
##################################################################################################
NVA.NVA.t3.bym <- function(){
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
    Spatial[j, 1:Nareas] ~ car.normal(adj[], weights[], num[],1) # structured
    for (i in 1:Nareas){
      Het[j, i] ~ dnorm(0, 1) # unstructured
      c.Het[j,i]<- Het[j,i] - mean(Het[j,])
      tPhi[j,i] <- Spatial[j,i]
    }
  }
  for (j in (Ndiseases + 1):(2 * Ndiseases)){
    for (i in 1:Nareas){
      tPhi[j,i] <- c.Het[(j - Ndiseases),i]
    }
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Ndiseases){
      M[i,j] ~ dnorm(0, prec.sp)
    }
    for (i in (Ndiseases+1):(2*Ndiseases)){
      M[i,j] ~ dnorm(0, prec.het)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  prec.sp <- pow(sdstruct.sp, -2)
  sdstruct.sp ~ dunif(0,100)
  prec.het <- pow(sdstruct.het, -2)
  sdstruct.het ~ dunif(0,100)
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
  for (i in 1:(Ndiseases)){
    for (j in 1:Ndiseases){
      Mg[i,j] ~ dnorm(0, precg)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  precg <- pow(sdstructg, -2)
  sdstructg ~ dunif(0,100)
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
## Spatial: NVA (bym); Temporal: NVA (rw1); Spatio-temporal: Type IV
## NVA.NVA.t4.bym
##################################################################################################
NVA.NVA.t4.bym <- function(){
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
    Spatial[j, 1:Nareas] ~ car.normal(adj[], weights[], num[],1) # structured
    for (i in 1:Nareas){
      Het[j, i] ~ dnorm(0, 1) # unstructured
      c.Het[j,i]<- Het[j,i] - mean(Het[j,])
      tPhi[j,i] <- Spatial[j,i]
    }
  }
  for (j in (Ndiseases + 1):(2 * Ndiseases)){
    for (i in 1:Nareas){
      tPhi[j,i] <- c.Het[(j - Ndiseases),i]
    }
  }
  ##################
  ## M-matrix (KxJ, K=2J)
  ##################
  for (j in 1:Ndiseases){
    for (i in 1:Ndiseases){
      M[i,j] ~ dnorm(0, prec.sp)
    }
    for (i in (Ndiseases+1):(2*Ndiseases)){
      M[i,j] ~ dnorm(0, prec.het)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  prec.sp <- pow(sdstruct.sp, -2)
  sdstruct.sp ~ dunif(0,100)
  prec.het <- pow(sdstruct.het, -2)
  sdstruct.het ~ dunif(0,100)
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
  for (i in 1:(Ndiseases)){
    for (j in 1:Ndiseases){
      Mg[i,j] ~ dnorm(0, precg)
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
  ##################
  ## Prior distribution for the standard deviations of the random effects
  ##################
  precg <- pow(sdstructg, -2)
  sdstructg ~ dunif(0,100)
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