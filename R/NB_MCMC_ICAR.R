#' @title Function for the NB fused model with ICAR prior
#'
#' @param y1 is the expression vector of gene 1, gene 1 will be regressed on gene 2
#' @param y2 is the expression vector of gene 2
#' @param X is the matrix of additional cell-level covariates
#' @param nIter is the number of MCMC iterations
#' @param beta_thres is a hard threshold on the absolute values of beta's
#' @param nug_sig1 is the additional precision on beta0's
#' @param nug_sig2 is the additional precision on beta1's
#' @param which.r.sampler is the choice of sampling r, either CRT or MH
#' @param verbose is TRUE or FALSE, to show or suppress progression updates
#' @return It returns the estimated beta's, r, shrinkage parameters, and deviance
#'
#' @export


NB_ICAR_model <-function(y1, y2, X = NULL, G, nIter = 5000, beta_thres = 10, 
                    nug_sig1 = 0.05, nug_sig2 = 0.1, which.r.sampler = "CRT", 
                    verbose = "TRUE"){
  
  # variable definitions
  y <- y1
  K <- spam::as.spam(diag(scale(log1p(y2))[, 1]))                      # design matrix with scaled log1p(y2)'s on the diagonal
  p <- ifelse(!is.null(X), ncol(X), 0)                  # number of covariates
  N <- length(y1)                                       # number of samples
  p_mst <- length(which(G == 1))                        # twice the number of connections
  half_p_mst <- p_mst/2
  
  # non-zero indices from the adjacency matrix
  W_sp <- spam::as.spam(G)                              # ICAR adjacency is just G
  G[lower.tri(G)] <- 0
  sorted_indices <- which(G == 1, arr.ind = T)          # indices of non-zero elements of the upper triangle of G matrix
  sorted_indices <- sorted_indices[order(sorted_indices[,1], sorted_indices[,2]), ]  
  
  # MCMC features
  thin <- 1				            # thinning interval
  burn <- nIter/2		          # burn-in
  lastit <- (nIter-burn)/thin	# last stored value
  
  # Storage objects
  Beta <- matrix(0, lastit, (2*N))                   # stores betas
  if(p > 0){Beta <- matrix(0, lastit, (2*N + p))}          
  Tau <- matrix(0, lastit, 2)                    # stores ICAR variances
  Deviance <- matrix(0, lastit, 1)                   # stores deviance
  R <- matrix(0,lastit,1)                            # stores over-dispersion parameter
  
  # Initialize
  beta0 <- beta1 <- rep(1, N)      # beta0's and beta1's
  if(p != 0){beta_cov = rep(0, p)  # slopes for extra covariates if any
  beta1 <- c(beta_cov, beta1)      # putting the slopes with beta1's
  beta <- c(beta0, beta1)
  }else{beta<-c(beta0, beta1)}
  
  r <- 1                           # NB over-dispersion parameter
  
  ## Build adjacency and (unnormalized) Laplacian Q = D - W on N nodes
  deg <- as.numeric(spam::rowSums(W_sp))
  Qicar <- spam::diag.spam(deg, N, N) - W_sp
  
  ## Connected components
  g  <- igraph::graph_from_edgelist(as.matrix(sorted_indices), directed = FALSE)
  comp <- igraph::components(g)$membership
  n_comp <- length(unique(comp))
  
  ## ICAR scales and priors
  tau0 <- 1.0; a_tau0 <- 1.0; b_tau0 <- 1.0
  tau1 <- 1.0; a_tau1 <- 1.0; b_tau1 <- 1.0
  
  
  # Initialize precision matrices using spam
  B0_mat <- B1_mat <- spam::spam(0, nrow = N, ncol = N)  
  
  if(p > 0){
    T0_sp <-  spam::spam(0, p, p)           # or diag.spam(rep(lambda, p))
    B1_star_mat <- spam::bdiag.spam(T0_sp, B1_mat)  # stays in spam-land
  }
  
  T0 <- 0.1                                          # prior precision of the additional covariates
  l <- rep(0, N)                                     # latent vector for CRT-based update of r
  a <- b <- 1                                        # gamma prior parameters on r for CRT                      
  Acc <- 0; s <- 0.01                                # counter for MH-based update of r, and proposal SD 
  del1 <- 1; del2 <- 1                               # extra prior parameters on shrinkage parameters, unused
  nugget <- 1e-8
  
  
  if(verbose == "TRUE"){print(paste0("SpaceBF: NB model fitting has started."))}
  if(which.r.sampler == "CRT" & p == 0){
    # algorithm with CRT-based r sampling and no covariates
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the predictor function and probability vector q
      #############################################################################
      Kbeta1 <- K%*%beta1 
      eta <- beta0 + Kbeta1 
      q <- 1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q[q==0] <- 1e-5
      
      #############################################################################
      ## Algorithm 2nd step: Chinese restaurant table-based r update
      #############################################################################
      l <- sapply(1:N, function(j) sum(rbinom(y[j], 1, round(r / (r + 1:y[j] - 1), 6))))
      r <- rgamma(1, a + sum(l), b - sum(log(q)))
      
      
      #############################################################################
      ## Algorithm 3rd step: Polya-gamma weights and latent response update
      #############################################################################
      w <- BayesLogit::rpg(N,y+r,eta)                   # Polya weights
      z <- (y-r)/(2*w)                                  # latent response
      
      
      #########################################################
      ## Step 4: ICAR precisions for beta0 and beta1  
      #########################################################
      
      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + w +  nug_sig1
      
      # K^T W^(1/2), needed to construct K^TWK for B1_mat                     
      sqrt_w_K <- K
      sqrt_w_K@entries <- sqrt(w) * sqrt_w_K@entries

      B1_mat <- tau1 * Qicar + t(sqrt_w_K) %*% sqrt_w_K
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2
      
      
      #############################################################################
      ## Step 5: sample beta0 and beta1 
      #############################################################################
      beta0 <- as.vector(spam::rmvnorm.canonical(1, b = w * (z - as.vector(K %*% beta1)), Q = B0_mat))
      b1 <- as.vector(t(sqrt_w_K) %*% (sqrt(w) * (z - beta0)))
      beta1 <- as.vector(spam::rmvnorm.canonical(1, b = b1, Q = B1_mat))

      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      
      beta <- c(beta0, beta1)
      
      #############################################################################
      ## Step 6: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      quad1 <- as.numeric(t(beta1) %*% (Qicar %*% beta1))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        R[j] <- r
        Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(which.r.sampler == "CRT" & p != 0){
    # same algorithm as above but with additional covariate adjustment
    K = spam::cbind(X, K)      # redefine K to include the covariates matrix X
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the predictor function and probability vector q
      #############################################################################
      Kbeta1 <- K%*%beta1 
      eta <- beta0 + Kbeta1 
      q <- 1/(1+exp(eta)) 
      q[q==0] <- 1e-5
      
      #############################################################################
      ## Algorithm 2nd step: Chinese restaurant table-based r update
      #############################################################################
      l <- sapply(1:N, function(j) sum(rbinom(y[j], 1, round(r / (r + 1:y[j] - 1), 6))))
      r <- rgamma(1, a + sum(l), b - sum(log(q)))
      
      #############################################################################
      ## Algorithm 3rd step: Polya-gamma weights and latent response update
      #############################################################################
      w <- BayesLogit::rpg(N,y+r,eta)                   
      z <- (y-r)/(2*w)                                  
      
      
      #########################################################
      ## Step 4: ICAR precisions for beta0 and beta1  
      #########################################################
      
      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + w +  nug_sig1
      
      # K^T W^(1/2), needed to construct K^TWK for B1_mat                     
      sqrt_w_K <-spam::diag(sqrt(w)) %*% K
      
      B1_mat <- tau1 * Qicar 
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2
      
      B1_star_mat <- spam::bdiag.spam(spam::spam(T0, p, p, sparse = T, doDiag=FALSE),  B1_mat)
      B1_star_mat <- (B1_star_mat) + t(sqrt_w_K) %*% sqrt_w_K
      
      #############################################################################
      ## Step 5: sample beta0 and beta1 
      #############################################################################
      beta0 <- as.vector(spam::rmvnorm.canonical(1, b = w * (z - as.vector(K %*% beta1)), Q = B0_mat))
      b1 <- as.vector(t(sqrt_w_K) %*% (sqrt(w) * (z - beta0)))
      beta1 <- as.vector(spam::rmvnorm.canonical(1, b = b1, Q = B1_star_mat))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta  <- c(beta0, beta1)
      
      #############################################################################
      ## Step 6: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   # f' Q f
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      beta1_spatial <- beta1[-c(1:p)]                   # extracting spatial slopes: beta1(s)'s and not the covariates
      quad1 <- as.numeric(t(beta1_spatial) %*% (Qicar %*% beta1_spatial))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        R[j] <- r
        Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  } else if(which.r.sampler == "MH" & p == 0){
    # algorithm with MH-based r sampling and no covariates
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the predictor function and probability vector q
      #############################################################################
      Kbeta1 <- K%*%beta1 
      eta <- beta0 + Kbeta1 
      q <- 1/(1+exp(eta)) # dnbinom fn uses q=1-psi
      q[q==0] <- 1e-5
      
      #############################################################################
      ## Algorithm 2nd step: Metropolis Hastings-based r update
      #############################################################################
      rnew <- msm::rtnorm(1, r, sqrt(s), lower=0) # Proposal
      ratio <- sum(dnbinom(y, rnew ,q, log=T), na.rm = T) - 
        sum(dnbinom(y, r, q, log=T), na.rm = T) +      # Likelihood difference, diffuse prior for r
        msm::dtnorm(rnew, r, sqrt(s), 0, log=T) - msm::dtnorm(r, rnew, sqrt(s), 0, log=T)  # Proposal not symmetric
      if (is.na(ratio) == F & log(runif(1))<ratio) {
        r <- rnew                                        # Update r
        Acc <- Acc+1
      }
      
      #############################################################################
      ## Algorithm 3rd step: Polya-gamma weights and latent response update
      #############################################################################
      w <- BayesLogit::rpg(N,y+r,eta)                   # Polya weights
      z <- (y-r)/(2*w)                                  # latent response
      
      
      #########################################################
      ## Step 4: ICAR precisions for beta0 and beta1  
      #########################################################
      
      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + w +  nug_sig1
      
      # K^T W^(1/2), needed to construct K^TWK for B1_mat                     
      sqrt_w_K <- K
      sqrt_w_K@entries <- sqrt(w) * sqrt_w_K@entries
      
      B1_mat <- tau1 * Qicar + t(sqrt_w_K) %*% sqrt_w_K
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2
      
      
      #############################################################################
      ## Step 5: sample beta0 and beta1 
      #############################################################################
      beta0 <- as.vector(spam::rmvnorm.canonical(1, b = w * (z - as.vector(K %*% beta1)), Q = B0_mat))
      b1 <- as.vector(t(sqrt_w_K) %*% (sqrt(w) * (z - beta0)))
      beta1 <- as.vector(spam::rmvnorm.canonical(1, b = b1, Q = B1_mat))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      
      beta <- c(beta0, beta1)
      
      #############################################################################
      ## Step 6: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   # f' Q f
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      quad1 <- as.numeric(t(beta1) %*% (Qicar %*% beta1))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        R[j] <- r
        Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(which.r.sampler == "MH" & p != 0){
    # same algorithm as above but with additional covariate adjustment
    K = spam::cbind(X, K)      # redefine K to include the covariates matrix X
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the predictor function and probability vector q
      #############################################################################
      Kbeta1 <- K%*%beta1 
      eta <- beta0 + Kbeta1 
      q <- 1/(1+exp(eta)) 
      q[q==0] <- 1e-5
      
      #############################################################################
      ## Algorithm 2nd step: Metropolis Hastings-based r update
      #############################################################################
      rnew <- msm::rtnorm(1, r, sqrt(s), lower=0) # Proposal
      ratio <- sum(dnbinom(y, rnew ,q, log=T), na.rm = T) - 
        sum(dnbinom(y, r, q, log=T), na.rm = T) +      # Likelihood difference, diffuse prior for r
        msm::dtnorm(rnew, r, sqrt(s), 0, log=T) - msm::dtnorm(r, rnew, sqrt(s), 0, log=T)  # Proposal not symmetric
      if (is.na(ratio) == F & log(runif(1))<ratio) {
        r <- rnew                                        # Update r
        Acc <- Acc+1
      }
      
      #############################################################################
      ## Algorithm 3rd step: Polya-gamma weights and latent response update
      #############################################################################
      w<-BayesLogit::rpg(N,y+r,eta)                   
      z<-(y-r)/(2*w)                                  
      
      
      #########################################################
      ## Step 4: ICAR precisions for beta0 and beta1  
      #########################################################
      
      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + w +  nug_sig1
      
      # K^T W^(1/2), needed to construct K^TWK for B1_mat                     
      sqrt_w_K <-spam::diag(sqrt(w)) %*% K
      
      B1_mat <- tau1 * Qicar 
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2
      
      B1_star_mat <- spam::bdiag.spam(spam::spam(T0, p, p, sparse = T, doDiag=FALSE),  B1_mat)
      B1_star_mat <- (B1_star_mat) + t(sqrt_w_K) %*% sqrt_w_K
      
      #############################################################################
      ## Step 5: sample beta0 and beta1 
      #############################################################################
      beta0 <- as.vector(spam::rmvnorm.canonical(1, b = w * (z - as.vector(K %*% beta1)), Q = B0_mat))
      b1 <- as.vector(t(sqrt_w_K) %*% (sqrt(w) * (z - beta0)))
      beta1 <- as.vector(spam::rmvnorm.canonical(1, b = b1, Q = B1_star_mat))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta  <- c(beta0, beta1)
      
      #############################################################################
      ## Step 6: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   # f' Q f
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      beta1_spatial <- beta1[-c(1:p)]                   # extracting spatial slopes: beta1(s)'s and not the covariates
      quad1 <- as.numeric(t(beta1_spatial) %*% (Qicar %*% beta1_spatial))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        R[j] <- r
        Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }
  
  if(verbose == "TRUE"){print(paste0("SpaceBF: NB model fitting completed!"))}
  
  return(NB = list(Beta = Beta, R = R, Tau = Tau, Deviance = Deviance)) 
}
