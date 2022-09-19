options(dplyr.summarise.inform = FALSE)

######################################################################################
## Fit a spatial-temporal multivariate Poisson mixed model to areal count data,     ##
## where dependence between spatial/temporal patterns of the diseases is addressed  ##
## through the use of M-models (Botella-Rocamora et al, 2015)                       ##
######################################################################################
MCAR_INLA_ST <- function(carto=NULL, data=NULL, ID.area=NULL, ID.year=NULL, ID.disease=NULL,
                         O=NULL, E=NULL, W=NULL, spatial="intrinsic", temporal="rw1",
                         interaction="TypeIV", strategy="simplified.laplace"){
  
  ## Check for errors ##
  if(is.null(carto))
    stop("the 'carto' argument is missing")
  if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
    stop("the 'carto' argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
  if(is.null(ID.area))
    stop("the 'ID.area' argument is missing")
  if(is.null(ID.year))
    stop("the 'ID.year' argument is missing")
  if(is.null(ID.disease))
    stop("the 'ID.disease' argument is missing")
  if(is.null(O))
    stop("the 'O' argument is missing")
  if(is.null(E))
    stop("the 'E' argument is missing")
  if(!(spatial %in% c("Leroux","intrinsic","proper","BYM")))
    stop("invalid 'spatial' argument")
  if(!(temporal %in% c("rw1","rw2")))
    stop("invalid 'temporal' argument")
  if(!(interaction %in% c("none","TypeI","TypeII","TypeIII","TypeIV")))
    stop("invalid 'interaction' argument")
  if(!(strategy %in% c("gaussian","simplified.laplace","laplace","adaptative")))
    stop("invalid 'strategy' argument")
  
  
  cat("STEP 1: Pre-processing data\n")
  
  ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
  carto <- sf::st_as_sf(carto)
  
  ## Order the data ##
  if(!ID.area %in% colnames(carto))
    stop(sprintf("'%s' variable not found in carto object",ID.area))
  if(!ID.area %in% colnames(data))
    stop(sprintf("'%s' variable not found in data object",ID.area))
  if(!ID.year %in% colnames(data))
    stop(sprintf("'%s' variable not found in data object",ID.year))
  if(!ID.disease %in% colnames(data))
    stop(sprintf("'%s' variable not found in data object",ID.disease))
  if(!O %in% colnames(data))
    stop(sprintf("'%s' variable not found in carto object",O))
  if(!E %in% colnames(data))
    stop(sprintf("'%s' variable not found in carto object",E))
  
  carto <- carto[order(sf::st_set_geometry(carto, NULL)[,ID.area]),]
  data <- merge(data,carto[,c(ID.area)])
  data$geometry <- NULL
  data[,ID.year] <- paste(sprintf("%02d", as.numeric(as.character(data[,ID.year]))))
  data[,ID.disease] <- paste(sprintf("%02d", as.numeric(as.character(data[,ID.disease]))))
  data <- data[order(data[,ID.disease],data[,ID.year],data[,ID.area]),]
  
  
  ## Define spatial adjacency/structure matrix ##
  if(is.null(W)){
    carto.nb <- poly2nb(carto)
    Ws <- as(nb2listw(carto.nb, style="B"),"CsparseMatrix")
  }else{
    carto.nb <- mat2listw(W, style="B")$neighbours
    Ws <- W
  }
  carto.nc <- n.comp.nb(carto.nb)$nc
  if(carto.nc!=1) stop(sprintf("'carto' object has %d disjoint connected subgraphs",carto.nc))
  
  S <- length(carto.nb)
  Rs <- inla.as.sparse(Diagonal(S,colSums(Ws))-Ws)
  Rs.Leroux <- inla.as.sparse(Diagonal(S)-Rs)
  
  
  ## Define temporal adjacency/structure matrix ##
  T <- length(unique(data[,ID.year]))
  if(temporal=="rw1") dif <- 1
  if(temporal=="rw2") dif <- 2
  D <- diff(diag(T), differences=dif)
  Rt <- inla.as.sparse(t(D)%*%D)
  
  Wt <- -Rt
  diag(Wt) <- 0
  Wt <- drop0(Wt)
  

  ## Define hyperprior distributions ##
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
  
  
  ## Identifiability constraints for the interaction term ##
  constraints <- function(Rs, Rt) {
    
    S <- nrow(Rs)
    T <- nrow(Rt)
    
    if(interaction=="none"){
      R <- NULL
      r.def <- NULL
      A.constr.st <- NULL
    }
    if(interaction=="TypeI"){
      R <- Diagonal(S*T)
      r.def <- 0
      if(temporal=="rw1") A.constr.st <- matrix(1,1,S*T)
      if(temporal=="rw2") A.constr.st <- matrix(rep(1:S,each=T),1,S*T)
    }
    if(interaction=="TypeII"){
      R <- kronecker(Rt,Diagonal(S))
      if(temporal=="rw1") r.def <- S
      if(temporal=="rw2") r.def <- 2*S
      A.constr.st <- kronecker(matrix(1,1,T),Diagonal(S))
      A.constr.st <- as(A.constr.st[-1,],"matrix")
    }
    if(interaction=="TypeIII"){
      R <- kronecker(Diagonal(T),Rs)
      r.def <- T
      A.constr.st <- kronecker(Diagonal(T),matrix(1,1,S))
      A.constr.st <- as(A.constr.st[-1,],"matrix")
    }
    if(interaction=="TypeIV"){
      R <- kronecker(Rt,Rs)
      if(temporal=="rw1") r.def <- S+T-1
      if(temporal=="rw2") r.def <- 2*S+T-2
      A1 <- kronecker(matrix(1,1,T),Diagonal(S))
      A2 <- kronecker(Diagonal(T),matrix(1,1,S))
      A.constr.st <- as(rbind(A1[-1,],A2[-1,]),"matrix")
    }
    
    return(list(R=R,r.def=r.def,A.constr.st=A.constr.st))
  }

  
  # Formula for INLA model ##
  J <- length(unique(data[,ID.disease]))
  
  form <- "O ~ -1+"
  form <- paste(form, paste(paste0("I",1:J),collapse="+"),sep="")
  form <- paste(form, "+ f(idx, model=Mmodel.s, constr=FALSE, extraconstr=list(A=A.constr.s, e=rep(0,J)))")
  form <- paste(form, "+ f(idy, model=Mmodel.t, constr=FALSE, extraconstr=list(A=A.constr.t, e=rep(0,J)))")
  if(interaction!="none"){
    form <- paste(form, "+", paste(paste0("f(idxy.",1:J,", model='generic0', Cmatrix=R, rankdef=r.def, constr=FALSE, extraconstr=list(A=A.constr.st, e=rep(0,nrow(A.constr.st))), hyper=list(prec=list(prior=sdunif)))"), collapse=" + "), sep=" ")
  }
  formula <- stats::as.formula(form)
  
  
  
  cat("STEP 2: Fitting the model with INLA (this may take a while...)\n")
  
  ## Identifiability constraints for each disease ##
  A.constr.s <- kronecker(diag(J), matrix(1,1,S))
  A.constr.t <- kronecker(diag(J), matrix(1,1,T))
  constr <- constraints(Rs,Rt)
  R <- constr$R
  r.def <- constr$r.def
  A.constr.st <- constr$A.constr.st
  
  
  ## Define data frame for INLA model ##
  data.INLA <- data.frame(O=data[,O], E=data[,E], Area=data[,ID.area], Year=data[,ID.year], Disease=data[,ID.disease],
                          ID.area=rep(1:S,T*J), ID.year=rep(rep(1:T,each=S),J), ID.disease=rep(1:J,each=S*T), ID.area.year=rep(1:J,each=S*T))
  
  intercepts <- dummy_cols(data.INLA$ID.disease)[,-1]
  intercepts[intercepts==0] <- NA
  colnames(intercepts) <- paste0("I",1:J)
  data.INLA <- cbind(data.INLA, intercepts)
  
  data.INLA$idx <- (data.INLA$ID.disease-1)*S + data.INLA$ID.area
  data.INLA$idy <- (data.INLA$ID.disease-1)*T +  data.INLA$ID.year
  idxy.<- paste0("idxy.", 1:J)
  for(j in 1:J){
    data.INLA[data.INLA$ID.disease==j, idxy.[j]] <-  (data.INLA[data.INLA$ID.disease==j, "ID.year"]-1)*S + data.INLA[data.INLA$ID.disease==j, "ID.area"]
  }
  
  
  ## Initial values for spatial M-model ##
  aux <- as.data.frame(data.INLA %>% group_by(ID.disease, ID.area) %>% summarise(SMR=sum(O)/sum(E)))
  Sigma <- cov(matrix(aux$SMR,S,J,byrow=F))
  N <- t(chol(Sigma))
  initial.values.s <- as.vector(c(log(diag(N)), N[lower.tri(N,diag=FALSE)]))

  ## Initial values for temporal M-model ##
  aux <- as.data.frame(data.INLA %>% group_by(ID.disease, ID.year) %>% summarise(SMR=sum(O)/sum(E)))
  Sigma <- cov(matrix(aux$SMR,T,J,byrow=F))
  N <- t(chol(Sigma))
  initial.values.t <- as.vector(c(log(diag(N)), N[lower.tri(N,diag=FALSE)]))
  
  
  ## Define selected M-model ##
  if(spatial=="intrinsic"){
    Mmodel.s <- inla.rgeneric.define(Mmodel_icar, debug=FALSE, J=J, W=Ws, initial.values=initial.values.s)
  }
  if(spatial=="Leroux"){
    Mmodel.s <- inla.rgeneric.define(Mmodel_lcar, debug=FALSE, J=J, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
  }
  if(spatial=="proper"){
    Mmodel.s <- inla.rgeneric.define(Mmodel_pcar, debug=FALSE, J=J, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
  }
  
  Mmodel.t <- inla.rgeneric.define(Mmodel_icar, debug=FALSE, J=J, W=Wt, initial.values=initial.values.t)
  
  
  ## Fit the INLA model  ##
  Model <- inla(formula, family="poisson", data=data.INLA, E=E,
                control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                control.inla=list(strategy=strategy))
  
  Model$Mmodel <- list(spatial=spatial, temporal=spatial, interaction=interaction)
  
  
  ## Compute spatial/temporal between-disease correlations and marginal variances ##
  Mmodel.compute <- Mmodel_compute_cor(Model, n.sample=1000)
  Model$summary.cor <- Mmodel.compute$summary.cor
  Model$marginals.cor <- Mmodel.compute$marginals.cor
  Model$summary.var <- Mmodel.compute$summary.var
  Model$marginals.var <- Mmodel.compute$marginals.var
  
  return(Model)
}


#########################
## Auxiliary functions ##
#########################
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

compute.summary <- function(marginal){
  aux <- data.frame(INLA::inla.emarginal(function(x) x, marginal),
                    sqrt(INLA::inla.emarginal(function(x) x^2, marginal)-INLA::inla.emarginal(function(x) x, marginal)^2),
                    INLA::inla.qmarginal(0.025,marginal),
                    INLA::inla.qmarginal(0.5,marginal),
                    INLA::inla.qmarginal(0.975,marginal),
                    INLA::inla.mmarginal(marginal),
                    1-INLA::inla.pmarginal(0,marginal))
  colnames(aux) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode","1 cdf")
  aux
}

Mmodel_compute_cor <- function(model, n.sample=10000){
  
  o <- tryCatch.W.E({
    J <- length(unique(model$.args$data$ID.disease))
    hyperpar.sample <- INLA::inla.hyperpar.sample(n.sample, model, improve.marginals=TRUE)
    hyperpar.sample.s <- hyperpar.sample[,grep("idx",colnames(hyperpar.sample))]
    hyperpar.sample.t <- hyperpar.sample[,grep("idy",colnames(hyperpar.sample))]
    
    
    ## Covariance/correlation matrix for spatial M-model ##
    if(model$Mmodel$spatial=="intrinsic"){
      hyperpar.sample.s[,1:J] <- exp(hyperpar.sample.s[,1:J])
      hyperpar.sample.s <- split(hyperpar.sample.s[,seq(J*(J+1)/2)], seq(nrow(hyperpar.sample.s)))
    }else{
      hyperpar.sample.s[,seq(J+1,2*J)]<- exp(hyperpar.sample.s[,seq(J+1,2*J)])
      hyperpar.sample.s <- split(hyperpar.sample.s[,seq(1+J,J+J*(J+1)/2)], seq(nrow(hyperpar.sample.s)))
    }

    param.sample.s <- lapply(hyperpar.sample.s, function(x){
      N <- diag(x[seq(J)]) 
      N[lower.tri(N, diag=FALSE)] <- x[-seq(J)]
      Sigma <- N %*% t(N)
      Rho <- cov2cor(Sigma)
      Rho.values <- Rho[lower.tri(Rho)]
      
      return(list(sigma=diag(Sigma),rho=Rho.values))
    })
    
    cor.sample.s <- do.call(rbind,lapply(param.sample.s, function(x) x$rho))
    cor.density.s <- apply(cor.sample.s, 2, function(x) density(x, n=75, bw="SJ", from=-1, to=1))
      
    marginals.cor.s <- lapply(cor.density.s, function(xx) cbind(x=xx$x, y=xx$y))
    names(marginals.cor.s) <- paste("rho.s",apply(combn(J,2), 2, function(x) paste0(x, collapse="")),sep="")
    summary.cor.s <- do.call(rbind,lapply(marginals.cor.s, compute.summary))[,1:6]
      
    var.sample.s <- do.call(rbind,lapply(param.sample.s, function(x) x$sigma))
    var.density.s <- apply(var.sample.s, 2, function(x) density(x, n=75, bw="SJ", from=0))
      
    marginals.var.s <- lapply(var.density.s, function(xx) cbind(x=xx$x, y=xx$y))
    names(marginals.var.s) <- paste("var.s",1:J,sep="")
    summary.var.s <- do.call(rbind,lapply(marginals.var.s, compute.summary))[,1:6]
    
    
    ## Covariance/correlation matrix for temporal M-model ##
    hyperpar.sample.t[,1:J] <- exp(hyperpar.sample.t[,1:J])
    hyperpar.sample.t <- split(hyperpar.sample.t[,seq(J*(J+1)/2)], seq(nrow(hyperpar.sample.t)))
    
    param.sample.t <- lapply(hyperpar.sample.t, function(x){
      N <- diag(x[seq(J)]) 
      N[lower.tri(N, diag=FALSE)] <- x[-seq(J)]
      Sigma <- N %*% t(N)
      Rho <- cov2cor(Sigma)
      Rho.values <- Rho[lower.tri(Rho)]
      
      return(list(sigma=diag(Sigma),rho=Rho.values))
    })
    
    cor.sample.t <- do.call(rbind,lapply(param.sample.t, function(x) x$rho))
    cor.density.t <- apply(cor.sample.t, 2, function(x) density(x, n=75, bw="SJ", from=-1, to=1))
    
    marginals.cor.t <- lapply(cor.density.t, function(xx) cbind(x=xx$x, y=xx$y))
    names(marginals.cor.t) <- paste("rho.t",apply(combn(J,2), 2, function(x) paste0(x, collapse="")),sep="")
    summary.cor.t <- do.call(rbind,lapply(marginals.cor.t, compute.summary))[,1:6]
    
    var.sample.t <- do.call(rbind,lapply(param.sample.t, function(x) x$sigma))
    var.density.t <- apply(var.sample.t, 2, function(x) density(x, n=75, bw="SJ", from=0))
    
    marginals.var.t <- lapply(var.density.t, function(xx) cbind(x=xx$x, y=xx$y))
    names(marginals.var.t) <- paste("var.t",1:J,sep="")
    summary.var.t <- do.call(rbind,lapply(marginals.var.t, compute.summary))[,1:6]
  })
    
  if(any(class(o[[1]])=="error")){
    
    summary.cor <- data.frame(rep(NA,2*ncol(combn(J,2))),rep(NA,2*ncol(combn(J,2))),rep(NA,2*ncol(combn(J,2))),rep(NA,2*ncol(combn(J,2))),rep(NA,2*ncol(combn(J,2))),rep(NA,2*ncol(combn(J,2))))
    colnames(summary.cor) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
    rownames(summary.cor) <- c(paste("rho.s",apply(combn(J,2), 2, function(x) paste0(x, collapse="")),sep=""),
                               paste("rho.t",apply(combn(J,2), 2, function(x) paste0(x, collapse="")),sep=""))
      
    marginals.cor <- as.list(rep(NA,2*ncol(combn(J,2))))
    names(marginals.cor) <- rownames(summary.cor)
      
    summary.var <- data.frame(rep(NA,2*J),rep(NA,2*J),rep(NA,2*J),rep(NA,2*J),rep(NA,2*J),rep(NA,2*J))
    colnames(summary.var) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
    rownames(summary.var) <- c(paste("var.s",1:J,sep=""),paste("var.t",1:J,sep=""))
      
    marginals.var <- as.list(rep(NA,2*J))
    names(marginals.var) <- rownames(summary.var)
  }
    
  Mmodel.compute <- list(summary.cor=rbind(summary.cor.s,summary.cor.t),
                         marginals.cor=append(marginals.cor.s,marginals.cor.t),
                         summary.var=rbind(summary.var.s,summary.var.t),
                         marginals.var=append(marginals.var.s,marginals.var.t))
  
  return(Mmodel.compute)
}
