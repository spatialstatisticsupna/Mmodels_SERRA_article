########################################################################
## M-models (INLA)
########################################################################
'inla.rgeneric.Mmodel.model.icar' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL)
  {
    ## M-model implementation using a intrinsic CAR with different parameters
    ## k: number of diseases/blocks
    ## W: adjacency matrix
    ########################################################################
    ## theta
    ########################################################################
    ## k * k entries of the M matrix, by columns.
    ## NOTE: The CAR distributions do NOT have a precision parameter.
    interpret.theta = function()
    {
      # Function for changing from internal scale to external scale
      # also, build the inverse of the matrix used to model in the external scale
      # the between-disease variability and then this matrix is inverted.
      
      # The k*k parameters are the entries in the k matrix, by cols.
      M <- matrix(theta, ncol = k)
      return (list(M = M))
    }
    ########################################################################
    ## Graph of precision function; i.e., a 0/1 representation of precision matrix
    ########################################################################
    graph = function()
    {
      # M \kronecker I
      MI <- kronecker(Matrix(1, ncol = k, nrow = k),
                      Diagonal(nrow(W), 1))
      # I + W
      IW <- Diagonal(nrow(W), 1) + W
      # Block diagonal
      BlockIW <- bdiag(replicate(k, IW, simplify = FALSE))
      G <- (MI %*% BlockIW) %*% MI
      return (G)
    }
    ########################################################################
    ## Precision matrix
    ########################################################################
    Q = function()
    {
      #Parameters in model scale
      param <- interpret.theta()
      # M^{-1} \kronecker I
      M.inv <- solve(param$M)
      MI <- kronecker(M.inv, Diagonal(nrow(W), 1))
      # Number of neighbours
      D <- as.vector(apply(W, 1, sum))
      BlockIW <- bdiag(lapply(1:k, function(i) {
        #        Diagonal(x = D) - param$alpha[i] * W
        Diagonal(x = D) -  W
      }))
      Q <- (MI %*% BlockIW) %*% kronecker(t(M.inv), Diagonal(nrow(W), 1))
      return (Q)
    }
    ########################################################################
    ## Mean of model
    ########################################################################
    mu = function() {
      return(numeric(0))
    }
    log.norm.const = function() {
      ## return the log(normalising constant) for the model
      #param = interpret.theta()
      #
      #val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation)) +
      # 0.5 * log(1.0 - param$alpha^2)
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
      ## Whishart prior for M^T * M
      ##############################
      sigma2 <- 1000 #Wishart parameter
      val = log(MCMCpack::dwish(W = crossprod(param$M), v = k, S = diag(rep(sigma2, k))))
      return (val)
    }
    ########################################################################
    ## initial
    ########################################################################
    initial = function() {
      ## return initial values
      ## The Initial values form a diagonal for M
      return ( c(as.vector(diag(rep(1, k)))))
    }
    ########################################################################
    ########################################################################
    quit = function() {
      return (invisible())
    }
    val = do.call(match.arg(cmd), args = list())
    return (val)
  }
########################################################################
########################################################################