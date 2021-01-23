########################################################################
## Mmodels - LCAR - FE (Whishart prior)
########################################################################
Mmodel_lcar_fe <-
        function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                         "log.prior", "quit"), theta = NULL) {
                envir = parent.env(environment())
                if (!exists("cache.done", envir = envir)) {
                        D <- as.vector(apply(W, 1, sum))
                        DW<- (Matrix::Diagonal(x = D) -  W)
                        DW<- inla.as.sparse(DW)
                        
                        assign("DW", DW, envir = envir)
                        assign("cache.done", TRUE, envir = envir)
                }
                ########################################################################
                ## theta
                ########################################################################
                interpret.theta = function(){
                        alpha <- alpha.min + (alpha.max - alpha.min) / (1 + exp(-theta[as.integer(1:k)]))
                        
                        sigma.j<- sapply(theta[as.integer(k+1:k)], function(x) { exp(-0.5*x) })
                        corre.j<- sapply(theta[as.integer(2*k + 1:(k*(k-1)/2))], function(x){ (2* exp(x))/(1+exp(x)) -1})
                                                            
                        Rho<- diag(1,k) 
                        Rho[lower.tri(Rho, diag = FALSE)] <- corre.j
                        Rho[upper.tri(Rho, diag = FALSE)] <- t(Rho)[upper.tri(Rho)]
                        sigma.mat<- matrix(sigma.j, ncol = 1) %*% matrix(sigma.j, nrow = 1)
                        Covar<- Rho * sigma.mat
                        
                        e<- eigen(Covar)
                        M<- t(e$vectors %*% diag(sqrt(e$values)))
                        
                        return (list(alpha = alpha, Covar=Covar, M = M))
                }
                ########################################################################
                ## Graph of precision function; i.e., a 0/1 representation of precision matrix
                ########################################################################
                graph = function(){
                        return (Q())
                }
                ########################################################################
                ## Precision matrix
                ########################################################################
                Q = function(){
                        param <- interpret.theta()
                        M.inv <- solve(param$M)
                        MI <- kronecker(M.inv, Matrix::Diagonal(nrow(W), 1))
                        BlockIW <- Matrix::bdiag(lapply(1:k, function(i) {
                                param$alpha[i] * DW + (1-param$alpha[i]) * Matrix::Diagonal(nrow(W), 1)
                        }))
                        Q <- (MI %*% BlockIW) %*% t(MI)
                        Q <- inla.as.sparse(Q)
                        return (Q)
                }
                ########################################################################
                ## Mean of model
                ########################################################################
                mu = function() {
                        return(numeric(0))
                }
                ########################################################################
                ## log.norm.const
                ########################################################################
                log.norm.const = function() {
                        val <- numeric(0)
                        return (val)
                }
                ########################################################################
                ## log.prior
                ########################################################################
                log.prior = function() {
                        ## return the log-prior for the hyperparameters.
                        param = interpret.theta()
                        ##############################
                        ## Uniform prior in (alpha.min, alpha.max) on model scale
                        ##############################
                        val =  sum(-theta[as.integer(1:k)] - 2 * log(1 + exp(-theta[as.integer(1:k)])) )
                        ##############################
                        ## Whishart prior for Covar = t(M)%*%M
                        ##############################
                        sigma2 <- 1000 # Wishart parameter
                        val = val + log(MCMCpack::dwish(W = param$Covar, v = k, S = diag(rep(sigma2, k))))
                        ## sigma (1:J)
                        val = val + sum(-(5/2)* theta[as.integer(k+1:k)])
                        ## rho
                        val = val + sum( log(2) + theta[as.integer( 2*k+1:(k*(k-1)/2) )] - 2*log(1 + exp(theta[as.integer( 2*k+1:(k*(k-1)/2) )])) )
                        return (val)
                }
                ########################################################################
                ## initial
                ########################################################################
                initial = function() {
                        ## return initial values
                        p <- (0.65 - alpha.min) / (alpha.max - alpha.min)
                        return ( c(rep(log(p / (1 - p)), k), as.vector(initial.values) ))
                }
                ########################################################################
                ########################################################################
                quit = function() {
                        return (invisible())
                }
                if (!length(theta)) theta = initial()
                val = do.call(match.arg(cmd), args = list())
                return (val)
        }
########################################################################
########################################################################


########################################################################
## Mmodels - LCAR - RE (Whishart prior)
########################################################################
Mmodel_lcar_re <-
        function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                         "log.prior", "quit"), theta = NULL) {
                envir = parent.env(environment())
                if (!exists("cache.done", envir = envir)) {
                        D <- as.vector(apply(W, 1, sum))
                        DW<- (Matrix::Diagonal(x = D) -  W)
                        DW<- inla.as.sparse(DW)
                        
                        assign("DW", DW, envir = envir)
                        assign("cache.done", TRUE, envir = envir)
                }
                ########################################################################
                ## theta
                ########################################################################
                interpret.theta = function(){
                        alpha <- alpha.min + (alpha.max - alpha.min) / (1 + exp(-theta[as.integer(1:k)]))
                        
                        sigma.j<- sapply(theta[as.integer(k+1:k)], function(x) { exp(-0.5*x) })
                        corre.j<- sapply(theta[as.integer( 2*k+1:(k*(k-1)/2) )], function(x){ (2* exp(x))/(1+exp(x)) -1})
                        
                        Rho<- diag(1,k) 
                        Rho[lower.tri(Rho, diag = FALSE)] <- corre.j
                        Rho[upper.tri(Rho, diag = FALSE)] <- t(Rho)[upper.tri(Rho)]
                        sigma.mat<- matrix(sigma.j, ncol = 1) %*% matrix(sigma.j, nrow = 1)
                        Covar<- Rho * sigma.mat
                        
                        e<- eigen(Covar)
                        M<- t(e$vectors %*% diag(sqrt(e$values)))
                        
                        ## The last parameter is sigma spatial
                        sigma<- exp(-(1/2)*theta[as.integer(k*(k+3)/2+1)])
                        
                        return (list(alpha = alpha, Covar=Covar, M = M, sigma=sigma))
                }
                ########################################################################
                ## Graph of precision function; i.e., a 0/1 representation of precision matrix
                ########################################################################
                graph = function(){
                        return (Q())
                }
                ########################################################################
                ## Precision matrix
                ########################################################################
                Q = function(){
                        param <- interpret.theta()
                        M.inv <- solve(param$M)
                        MI <- kronecker(M.inv, Matrix::Diagonal(nrow(W), 1))
                        BlockIW <- Matrix::bdiag(lapply(1:k, function(i) {
                                param$alpha[i] * DW + (1-param$alpha[i]) * Matrix::Diagonal(nrow(W), 1)
                        }))
                        Q <- (MI %*% BlockIW) %*% t(MI)
                        Q <- inla.as.sparse(Q)
                        return (Q)
                }
                ########################################################################
                ## Mean of model
                ########################################################################
                mu = function() {
                        return(numeric(0))
                }
                ########################################################################
                ## log.norm.const
                ########################################################################
                log.norm.const = function() {
                        val <- numeric(0)
                        return (val)
                }
                ########################################################################
                ## log.prior
                ########################################################################
                log.prior = function() {
                        ## return the log-prior for the hyperparameters
                        param = interpret.theta()
                        ###################################
                        ## log-prior for the spatial smoothing parameter
                        ###################################
                        val =  sum(-theta[as.integer(1:k)] - 2 * log(1 + exp(-theta[as.integer(1:k)])) )
                        ###################################
                        ## sigma spatial
                        ###################################
                        val = val + log(0.005)- (0.5)* theta[as.integer(k*(k+3)/2+1)]
                        ##############################
                        ## Whishart prior for Covar = t(M)%*%M
                        ##############################
                        val = val + log(MCMCpack::dwish(W = param$Covar, v = k, S = diag(rep(param$sigma^2, k))))
                        ## sigma (1:J)
                        val = val + sum(-(5/2)* theta[as.integer(k+1:k)])
                        ## rho
                        val = val + sum( log(2) + theta[as.integer( 2*k+1:(k*(k-1)/2) )] - 2*log(1 + exp(theta[as.integer( 2*k+1:(k*(k-1)/2) )])) )
                        
                        return (val)
                }
                ########################################################################
                ## initial
                ########################################################################
                initial = function() {
                        ## return initial values
                        p <- (0.65 - alpha.min) / (alpha.max - alpha.min)
                        return ( c(rep(log(p / (1 - p)), k), as.vector(initial.values), 0.5 ))
                }
                ########################################################################
                ########################################################################
                quit = function() {
                        return (invisible())
                }
                if (!length(theta)) theta = initial()
                val = do.call(match.arg(cmd), args = list())
                return (val)
        }
########################################################################
########################################################################