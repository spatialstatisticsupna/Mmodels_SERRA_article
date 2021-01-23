MCAR_INLA_st <- function(carto=NULL,
                         data=NULL, 
                         ID.area=NULL,
                         ID.year=NULL,
                         ID.group=NULL,
                         O=NULL,
                         E=NULL,
                         prior.spatial="Leroux",
                         prior.disease="FE",
                         prior.interaction=0,   
                         model="global",
                         strategy="simplified.laplace"
                         ){
        
        ## Check for errors ##
        if(is.null(carto))
                stop("the carto argument is missing")
        if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
                stop("the carto argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
        if(is.null(ID.area))
                stop("the ID.area argument is missing")
        if(is.null(ID.year))
                stop("the ID.year argument is missing")
        
        if(is.null(O))
                stop("the O argument is missing")
        if(is.null(E))
                stop("the E argument is missing")
        

        if(!(prior.spatial %in% c("Leroux","intrinsic","proper","BYM")))
                stop("invalid prior.spatial argument")
        if(!(prior.disease %in% c("FE","RE")))
                stop("invalid prior.disease argument")
        
        if(!(prior.interaction %in% c(0:4)))
                stop("invalid prior.interaction argument")
        
        
        if(!(model %in% c("global","partition")))
                stop("invalid model argument")
        if(!(strategy %in% c("gaussian","simplified.laplace","laplace","adaptative")))
                stop("invalid strategy argument")
        
        
        ## Pre-processing data                  ##
        ## ---------------------------------------
        cat("STEP 1: Pre-processing data\n")
        
        
        ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class ##
        carto <- sf::st_as_sf(carto)
        
        ## Order the data ##
        j <- length(data)
        if(is.null(names(data))) names(data) <- paste("disease",1:j,sep=".")
        
        for(i in names(data)){
                if(!ID.area %in% colnames(carto))
                        stop(sprintf("'%s' variable not found in carto object",ID.area))
                if(!ID.area %in% colnames(data[[i]]))
                        stop(sprintf("'%s' variable not found in '%s' data.frame object",ID.area,i))
                if(!ID.year %in% colnames(data[[i]]))
                        stop(sprintf("'%s' variable not found in '%s' data.frame object",ID.year,i))
                
                if(!O %in% colnames(data[[i]]))
                        stop(sprintf("'%s' variable not found in '%s' data.frame object",O,i))
                if(!E %in% colnames(data[[i]]))
                        stop(sprintf("'%s' variable not found in '%s' data.frame object",E,i))
                if(length(unique(data[[i]][,ID.area]))!=length(unique(sf::st_set_geometry(carto, NULL)[,ID.area])))
                        stop(sprintf("The values of '%s' variable do not match between data and carto objects", ID.area))
                if(!all(as.character(sort(unique(data[[i]][,ID.area])))==as.character(sort(unique(sf::st_set_geometry(carto, NULL)[,ID.area])))))
                        stop(sprintf("The values of '%s' variable do not match between data and carto objects", ID.area))
        }
        
        ## Define hyperprior distributions      ##
        ## ---------------------------------------
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
        
        
        ## Formula for INLA model               ##
        ## ---------------------------------------
        form <- "O ~ -1+"
        form <- paste(form, paste(paste0("I",1:j),collapse="+"),sep="")
        form <- paste(form, "+ f(idx, model=Mmodel, constr=FALSE, extraconstr=list(A=A.constr.s, e=rep(0,j)))")
        if(prior.spatial %in% c("BYM")){
                form <- paste(form, "+ f(idx.v, model=Mmodel.v, constr=FALSE, extraconstr=list(A=A.constr.s, e=rep(0,j)))")
        }
        form <- paste(form, "+ f(idy, model=Mmodel.t, constr=FALSE, extraconstr=list(A=A.constr.t, e=rep(0,j)))")
        if(prior.interaction %in% c(1:4)){
                form <- paste(form, "+", paste(paste0("f(idxy.",1:j,", model='generic0', Cmatrix=R.st, rankdef=r.def.st, constr=FALSE, extraconstr=list(A=A.constr.st, e=rep(0,dim(A.constr.st)[1])), hyper=list(prec=list(prior=sdunif)) )"), collapse=" + "), sep=" ")
        }
                
        formula <- stats::as.formula(form)
        
        
        
        ## Global model                         ##
        ## ---------------------------------------
        if(model=="global"){
                
                ## Fitting global model with INLA       ##
                ## ---------------------------------------
                cat("STEP 2: Fitting global model with INLA (this may take a while...)\n")
                
                ## order
                data <- lapply(data, function(x) x[order(x[,ID.year],x[,ID.area]),] )
                
                n <- length(unique(data[[1]][,ID.area]))
                t <- length(unique(data[[1]][,ID.year]))
                
                
                ## Adjacency matrices                   ##
                ## ---------------------------------------
                
                ## W: adjacency matrix (spatial)
                adj <- spdep::poly2nb(carto)
                Ws <- as(nb2mat(adj, style = "B"), "Matrix")
                
                ## W.t: adjacency matrix (temporal)
                Pt <- crossprod(diff(diag(t),differences=1))
                Wt <- INLA::inla.as.sparse(Matrix::Diagonal(n=nrow(Pt),diag(Pt))- Pt) ## temp.corre = TRUE
                
                ## Precision matrices                   ##
                ## ---------------------------------------
                Rs <- Matrix::Diagonal(n,colSums(Ws))-Ws
                Rt <- Matrix::Diagonal(t,colSums(Wt))-Wt
                
                # R.st: Precision matrices for spatio-temporal interaction
                if(prior.interaction %in% c(1)){
                        R.st <- diag(n*t)
                        r.def.st <- 0
                }
                if(prior.interaction %in% c(2)){
                        R.st <- kronecker(Rt, Matrix::Diagonal(n,1))
                        r.def.st <- n
                }
                if(prior.interaction %in% c(3)){
                        R.st <- kronecker(Matrix::Diagonal(t,1),Rs)
                        r.def.st <- t
                }
                if(prior.interaction %in% c(4)){
                        R.st <- kronecker(Rt, Rs)
                        r.def.st <- n+t-1
                }
                
                ## Define appropriate constraints matrices##
                ## ---------------------------------------
                A.constr.s<- kronecker(diag(j), matrix(1,1,n))
                A.constr.t<- kronecker(diag(j), matrix(1,1,t))
                
                if(prior.interaction %in% c(1)){
                        A.constr.st <- matrix(1, nrow=1, ncol=n*t)
                }
                if(prior.interaction %in% c(2)){
                        A.constr.st <- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
                }
                if(prior.interaction %in% c(3)){
                        A.constr.st <- kronecker(diag(1,t), matrix(1,nrow=1, ncol=n))
                }
                if(prior.interaction %in% c(4)){
                        A1 <- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
                        A2 <- kronecker(diag(1,t), matrix(1,nrow=1, ncol=n))
                        A.constr.st <- rbind(A1[-1,], A2)
                }
                
                
                ## Data for INLA model                  ##
                ## ---------------------------------------
                data.INLA <- data.frame(O=as.numeric(unlist(lapply(data, function(x) x[,O]))),
                                        E=as.numeric(unlist(lapply(data, function(x) x[,E]))),
                                        Area=as.character(unlist(lapply(data, function(x) x[,ID.area]))),
                                        Year=as.character(unlist(lapply(data, function(x) x[,ID.year]))),
                                        ID.area=1:n,
                                        ID.year=rep(1:t, each=n),
                                        ID.disease=rep(1:j, each=n*t)
                                        )
                intercepts <- fastDummies::dummy_cols(data.INLA$ID.disease)[,-1]
                intercepts[intercepts==0] <- NA
                colnames(intercepts) <- paste0("I",1:j)
                data.INLA <- cbind(data.INLA, intercepts)
                
                data.INLA$idx <- (data.INLA$ID.disease-1)*n + data.INLA$ID.area
                data.INLA$idx.v <- data.INLA$idx
                data.INLA$idy<- (data.INLA$ID.disease-1) *t +  data.INLA$ID.year
                
                idxy.<- paste0("idxy.", 1:j)
                for(i in 1:j){
                        data.INLA[data.INLA$ID.disease==i, idxy.[i] ]<- 
                                (data.INLA[data.INLA$ID.disease==i, "ID.year" ]-1)* n + data.INLA[data.INLA$ID.disease==i, "ID.area"]
                }
                
                ## Initial values for hyperparameters   ##
                ## ---------------------------------------
                library(tidyverse)
                ## Spatial initial values (hyperparameters)
                sir.s <- data.INLA %>% 
                        dplyr::group_by(ID.disease, ID.area) %>%
                        dplyr::summarise(SIR = sum(O)/sum(E))
                Sigma.s <- cov(matrix(sir.s$SIR, nrow=n, ncol=j, byrow = FALSE))
                Rho.s <- cor(matrix(sir.s$SIR, nrow=n, ncol=j, byrow = FALSE))
                # Rho.s<- stats::cov2cor(Sigma.s)
                
                ## Temporal initial values (hyperparameters)
                sir.t <- data.INLA %>%
                        dplyr::group_by(ID.disease, ID.year) %>%
                        dplyr::summarise(SIR = sum(O)/sum(E))
                Sigma.t <- cov(matrix(sir.t$SIR, nrow=t, ncol=j, byrow = FALSE))
                Rho.t <- cor(matrix(sir.t$SIR, nrow=t, ncol=j, byrow = FALSE))
                # Rho.t<- stats::cov2cor(Sigma.t)
                
                
                initial.values.s <- as.vector(
                        c(-log(diag(Sigma.s)),
                          log(1+Rho.s[upper.tri(Rho.s, diag = FALSE)])-log(1-Rho.s[upper.tri(Rho.s, diag = FALSE)])
                          ))
                initial.values.t <- as.vector(
                        c(-log(diag(Sigma.t)),
                          log(1+Rho.t[upper.tri(Rho.t, diag = FALSE)])-log(1-Rho.t[upper.tri(Rho.t, diag = FALSE)])
                          ))
                
                
                ## Define selected Mmodel               ##
                ## ---------------------------------------
                ## Spatial
                if(prior.spatial=="intrinsic" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                if(prior.spatial=="intrinsic" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                if(prior.spatial=="Leroux" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="Leroux" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="proper" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_pcar_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="proper" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_pcar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="BYM" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_fe,  debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                        Mmodel.v <- INLA::inla.rgeneric.define(Mmodel_iid_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                if(prior.spatial=="BYM" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                        Mmodel.v <- INLA::inla.rgeneric.define(Mmodel_iid_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                
                
                ## Temporal
                if(prior.disease=="FE"){
                        Mmodel.t <- INLA::inla.rgeneric.define(Mmodel_icar_fe, debug=FALSE, k=j, W=Wt, initial.values=initial.values.t)
                }
                if(prior.disease=="RE"){
                        Mmodel.t <- INLA::inla.rgeneric.define(Mmodel_icar_re, debug=FALSE, k=j, W=Wt, initial.values=initial.values.t)
                }
                
                
                
                ## Fit the INLA model                   ##
                ## ---------------------------------------
                Model <- INLA::inla(formula, family="poisson", data=data.INLA, E=E,
                                    control.predictor=list(compute=TRUE, cdf=c(log(1))),
                                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                                    control.inla=list(strategy=strategy),
                                    debug = F,
                                    verbose = F)
                
                Model$Mmodel <- list(model=model, prior.spatial=prior.spatial, prior.disease=prior.disease, prior.interaction=prior.interaction)
                
        }

        return(Model)
}
