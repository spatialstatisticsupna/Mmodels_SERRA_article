#' Leroux et al. (1999) multivariate CAR latent effect
#'
#' @description M-model implementation of the Leroux et al. (1999) multivariate CAR latent effect with different spatial smoothing parameters using the \code{rgeneric} model of INLA.
#' In the \code{Mmodel_lcar_fe()} function the cells of the M-matrix are treated as fixed effects with a normal prior with mean 0 and large fixed variance (denoted as \code{FE} model).
#' In contrast, in the \code{Mmodel_lcar_re()} function the cells of the M-matrix are treated as independent normal random variables with mean 0 and standard deviation \eqn{\sigma} (denoted as \code{RE} model). In this case, a uniform prior distribution between 0 and a large number is considered for \eqn{\sigma}.
#'
#' @details These functions defines a latent effect \eqn{\Theta=\Phi M} using the M-model parameterization proposed by \insertCite{botella2015unifying;textual}{bigDM},
#' where the columns of \eqn{\Phi} follow a Leroux et al. (1999) CAR prior distribution (within-disease correlation) and \eqn{M} is a \eqn{k \times k} matrix which defines the correlation structure between-diseases.
#' Uniform prior distributions on the interval [\code{alpha.min}, \code{alpha.max}] are considered for all the spatial smoothing parameters.
#' \cr\cr
#' The following arguments are required to be defined before calling the functions:
#' \itemize{
#' \item \code{W}: binary adjacency matrix of the spatial areal units
#' \item \code{k}: number of diseases
#' \item \code{initial.values}: initial values defined for the cells of the M-matrix
#' \item \code{alpha.min}: lower limit defined for the uniform prior distribution of the spatial smoothing parameters
#' \item \code{alpha.max}: upper limit defined for the uniform prior distribution of the spatial smoothing parameters
#' }
#'
#' @references
#' \insertRef{botella2015unifying}{bigDM}
#'
#' @param cmd Internall functions used by the \code{rgeneric} model to define the latent effect.
#' @param theta Vector of hyperparameters.
#'
#' @return This is used internally by the \code{INLA::inla.rgeneric.define()} function.
#'
#' @importFrom INLA inla.rgeneric.define
#' @importFrom Matrix bdiag Diagonal Matrix
#' @importFrom MCMCpack dwish
#'
#' @examples
#' ## Load the data and cartography of CAW data in Uttar Pradesh ##
#' data(Carto_UP)
#' head(Carto_UP)
#' plot(Carto_UP$geometry, axes=TRUE)
#'
#' data(Data_UP)
#' str(Data_UP)
#'
#' ## Merge the data.frame and sf objects using the 'ID.area' variable ##
#' ## and compute the binary spatial adjacency matrix                  ##
#' ID.area <- "district"
#' k <- length(Data_UP)
#'
#' aux <- lapply(Data_UP, function(x){
#'   carto.aux <- merge(Carto_UP,x)
#'   data.aux <- sf::st_set_geometry(carto.aux, NULL)
#'   carto.aux <- carto.aux[order(data.aux[,ID.area]),]
#'
#'   invisible(utils::capture.output(aux <- connect_subgraphs(carto.aux, ID.area)))
#'   aux
#' })
#'
#' if(length(unique(aux))==1){
#'   aux <- unique(aux)[[1]]
#' }else{
#'   sprintf("Internal error when creating the neighbour object")
#' }
#'
#' W <- aux$W
#' n <- nrow(W)
#'
#' ## Define the data for INLA model ##
#' data.INLA <- data.frame(
#'   O=as.numeric(unlist(lapply(Data_UP, function(x) x$observed))),
#'   E=as.numeric(unlist(lapply(Data_UP, function(x) x$expected))),
#'   Area=as.character(unlist(lapply(Data_UP, function(x) x[,ID.area]))),
#'   ID.area=rep(1:n,2), ID.disease=rep(1:k, each=n), idx=seq(1,n*k)
#' )
#'
#' intercepts <- fastDummies::dummy_cols(data.INLA$ID.disease)[,-1]
#' intercepts[intercepts==0] <- NA
#' colnames(intercepts) <- paste0("I",1:k)
#' data.INLA <- cbind(data.INLA, intercepts)
#'
#' ## Initial values for the cells of the M-matrix ##
#' e <- eigen(cov(do.call(cbind,lapply(Data_UP, function(x) x$observed/x$expected))))
#' M <- diag(sqrt(e$values)) %*% t(e$vectors)
#' initial.values <- as.vector(M)
#'
#' ## Define the M-model ##
#' Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar_fe, debug=FALSE, k=k, W=W,
#'                                      initial.values=initial.values, alpha.min=0, alpha.max=1)
#' Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar_re, debug=FALSE, k=k, W=W,
#'                                      initial.values=initial.values, alpha.min=0, alpha.max=1)
#'
#' ## Call the inla() function ##
#' A_constr <- kronecker(diag(k), matrix(1,1,n))
#'
#' formula <- O ~ -1 + I1 + I2 +
#'   f(idx, model=Mmodel, constr=FALSE, extraconstr=list(A=A_constr, e=rep(0,k)))
#'
#' Model <- INLA::inla(formula, family="poisson", data=data.INLA, E=E,
#'                     control.predictor=list(compute=TRUE, cdf=c(log(1))),
#'                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE))
#'
#' summary(Model)
#'
#' @name Mmodel_lcar
NULL

#' @rdname Mmodel_lcar
#'
#' @export
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


#' @rdname Mmodel_lcar
#'
#' @export
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