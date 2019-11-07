########################################################################
## M-models (INLA)
########################################################################
'inla.rgeneric.Mmodel.model.iid.n' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL)
  {
    ## M-model implementation using a proper CAR with different parameters
    ## k: number of diseases/blocks
    ## W: adjacency matrix
    ## alpha.min: minimum value for alpha
    ## alpha.max: maximum value for alpha
    ########################################################################
    ## theta
    ########################################################################
    ## k * k entries of the M matrix, by columns.
    ## NOTE: The CAR distributions do NOT have a precision parameter.
    interpret.theta = function()
    {
      #Function for changing from internal scale to external scale
      # also, build the inverse of the matrix used to model in the external scale
      # the between-disease variability and then this matrix is inverted.
      
      # The (k * k) parameters are the entries in the k matrix, by cols.
      M <- matrix(theta[as.integer(1:(k*k))], ncol = k)
      
      # The last parameter is sigma spatial
      sigma<- exp(-(1/2)*theta[as.integer((k*k)+1)])
      return (list(M = M, sigma=sigma))
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
      IW <- Diagonal(nrow(W), 1)
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
      BlockIW <- bdiag(lapply(1:k, function(i) {
        Diagonal(nrow(W), 1)
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
    ########################################################################
    ## log.norm.const
    ########################################################################
    log.norm.const = function() {
      val <- numeric(0)
      return (val)
    }
    ########################################################################
    # log.prior
    ########################################################################
    log.prior = function() {
      ## return the log-prior for the hyperparameters.
      param = interpret.theta()
      ###################################
      ## sigma spatial
      ###################################
      val= log(0.005)-0.5*theta[as.integer((k*k)+1)]
      ###################################
      ## normal for mij
      ###################################
      val = val + sum(dnorm(theta[as.integer(c(1:(k*k)))], mean=0, sd=param$sigma, log=TRUE))
      return (val)
    }
    ########################################################################
    ## initial
    ########################################################################
    initial = function() {
      ## return initial values
      # The Initial values form a diagonal matrix
      # 0.5 values for alphas and diagonal for M
      #      p <- (0.5 - alpha.min) / (alpha.max - alpha.min)
      #      return ( c(rep(log(p / (1 - p)), k),
      #                 as.vector(diag(rep(1, k)))))
      return ( c(as.vector(diag(rep(1, k))) ,0.5) )
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