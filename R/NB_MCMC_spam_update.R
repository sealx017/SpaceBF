#' @title Function for the NB fused model
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


NB_model <-function(y1, y2, X = NULL, G, nIter = 5000, beta_thres = 10, 
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
  G[lower.tri(G)] = 0
  sorted_indices <- which(G == 1, arr.ind = T)          # indices of non-zero elements of the G matrix
  sorted_indices <- sorted_indices[order(sorted_indices[,1], sorted_indices[,2]), ]  
  
  # MCMC features
  thin <- 1				            # thinning interval
  burn <- nIter/2		          # burn-in
  lastit <- (nIter-burn)/thin	# last stored value
  
  # Storage objects
  Beta <- matrix(0, lastit, (2*N))                   # stores betas
  if(p > 0){Beta <- matrix(0, lastit, (2*N + p))}          
  Tau <- matrix(0, lastit, 2)                        # stores global shrinkage parameters
  Lambda <- matrix(0, lastit, (2*half_p_mst))        # stores local shrinkage parameters, we are not worried about their convergence
  Deviance <- matrix(0, lastit, 1)                   # stores deviance
  R <- matrix(0,lastit,1)                            # stores over-dispersion parameter

  # Initialize
  beta0 <- beta1 <- rep(1, N)      # beta0's and beta1's
  if(p != 0){beta_cov = rep(0, p)  # slopes for extra covariates if any
  beta1 <- c(beta_cov, beta1)      # putting the slopes with beta1's
  beta <- c(beta0, beta1)
  }else{beta<-c(beta0, beta1)}
  
  r <- 1                           # NB over-dispersion parameter
  
  # Initialize precision matrices using spam
  B0_mat <- B1_mat <- spam::spam(0, nrow = N, ncol = N)  
  
  if(p == 1){
  B1_star_mat = Matrix::bdiag(1,  B1_mat)
  }else if(p > 1){B1_star_mat = spam::bdiag(spam::spam(0, p, p, sparse = T, doDiag=FALSE),  B1_mat)}
  
  zeta40 = zeta41 = rep(0.1, half_p_mst)             # local shrinkage parameters
  gamma40 = gamma41 = rep(0.1, half_p_mst)           # local shrinkage latent parameters
  tau240 = tau241 = 0.5                              # global shrinkage parameters    
  epsilon20 = epsilon21 =  1                         # global shrinkage latent parameters  
  T0 = 0.1                                           # prior precision of the additional covariates
 
 
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
        
        
        #############################################################################
        ## Algorithm 4th step: update the precision matrices of the beta0 and beta1's
        #############################################################################
        off_diag_val0 <- 1/zeta40/tau240
        off_diag_val1 <- 1/zeta41/tau241
        
        B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0   # adjust only non-zero 
        B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   # elements of precision matrices
        spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0

        spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1             # nug_sig1 and nug_sig2 are
        spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2             # additional precision terms, i.e., 

                                                                                # extra normal priors on the 
                                                                                # beta's, helping to stabilize 
                                                                                # the estimates
        # Add diagonal weight matrix to B0_mat using spam
        B0_mat <- B0_mat + spam::diag.spam(w)
        
        # K^T W^(1/2), needed to construct K^TWK for B1)mat                     
        sqrt_w_K <- K
        sqrt_w_K@entries <- sqrt(w) * sqrt_w_K@entries
        B1_mat <- B1_mat + t(sqrt_w_K) %*% sqrt_w_K
        
        ##################################################
        ## Algorithm 5th step: sample beta's
        ##################################################
        
        beta0 <- as.vector(spam::rmvnorm.canonical(1, 
                 w*(z - Kbeta1), (B0_mat)))

        beta1 <- as.vector(spam::rmvnorm.canonical(1, 
                  t(sqrt_w_K)%*%(sqrt(w)*(z - beta0)), (B1_mat)))
        
        beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
        beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
        
        beta <- c(beta0, beta1)

        ##################################################
        ## Algorithm 6th step: update shrinkage parameters
        ##################################################
        beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
        beta_diff1 <- (beta1[sorted_indices[,1]] - beta1[sorted_indices[,2]])^2 + nugget
        
        # Vectorized calculation for zeta40, gamma40, zeta41, gamma41
        zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
        gamma40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
        
        # Vectorized calculation for zeta41, gamma41, others
        zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
        gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
        
        tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
        epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
        
        tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
        epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
        
        
        if (i> burn & i%%thin==0) {
          j <- (i-burn)/thin
          Beta[j, ] <- beta
          R[j] <- r
          Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
          Tau[j, ] <- c(tau240, tau241)
          Lambda[j, ] <- c(zeta40, zeta41)
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
      
      
      #############################################################################
      ## Algorithm 4th step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      off_diag_val0 <- 1/zeta40/tau240
      off_diag_val1 <- 1/zeta41/tau241
      
      B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0  
      B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   
      spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0
      
      spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1             # nug_sig1 and nug_sig2 are
      spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2             # additional precision terms, i.e., 
      

      # Add diagonal weight matrix to B0_mat using spam
      B0_mat <- B0_mat + spam::diag.spam(w)
      
      # K^T W^(1/2), needed to construct K^TWK for B1)mat                     
      sqrt_w_K <- K
      sqrt_w_K@entries <- sqrt(w) * sqrt_w_K@entries

      B1_star_mat <- spam::bdiag(spam::spam(T0, p, p, sparse = T, doDiag=FALSE),  B1_mat)
      B1_star_mat <- (B1_star_mat) + t(sqrt_w_K) %*% sqrt_w_K
      
      ##################################################
      ## Algorithm 5th step: sample beta's
      ##################################################
      beta0 = as.vector(spam::rmvnorm.canonical(1, 
              w*(z - Kbeta1), B0_mat))

      beta1 = as.vector(spam::rmvnorm.canonical(1, 
              t(sqrt_w_K)%*%(sqrt(w)*(z - beta0)), B1_star_mat))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta  <- c(beta0, beta1)
      
      ##################################################
      ## Algorithm 6th step: update shrinkage parameters
      ##################################################
      beta1_spatial <- beta1[-c(1:p)]                   # extracting spatial slopes: beta1(s)'s and not the covariates
      beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
      beta_diff1 <- (beta1_spatial[sorted_indices[,1]] - beta1_spatial[sorted_indices[,2]])^2 + nugget
      
      zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
      gamma40 <-  MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
      
      zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
      gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
      
      tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
      epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
      
      tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
      epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        R[j] <- r
        Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
        Tau[j, ] <- c(tau240, tau241)
        Lambda[j, ] <- c(zeta40, zeta41)
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
      
      
      #############################################################################
      ## Algorithm 4th step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      off_diag_val0 <- 1/zeta40/tau240
      off_diag_val1 <- 1/zeta41/tau241
      
      B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0   # adjust only non-zero 
      B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   # elements of precision matrices
      spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0
      
      spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1             # nug_sig1 and nug_sig2 are
      spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2             # additional precision terms, i.e., 
      

      # Add diagonal weight matrix to B0_mat using spam
      B0_mat <- B0_mat + spam::diag.spam(w)
      
      # K^T W^(1/2), needed to construct K^TWK for B1)mat                     
      sqrt_w_K <- K
      sqrt_w_K@entries <- sqrt(w) * sqrt_w_K@entries
      B1_mat <- B1_mat + t(sqrt_w_K) %*% sqrt_w_K
      
      ##################################################
      ## Algorithm 5th step: sample beta's
      ##################################################
      
      beta0 <- as.vector(spam::rmvnorm.canonical(1, 
                w*(z - Kbeta1), B0_mat))
      
      beta1 <- as.vector(spam::rmvnorm.canonical(1, 
               t(sqrt_w_K)%*%(sqrt(w)*(z - beta0)), B1_mat))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      
      beta = c(beta0, beta1)
      
      ##################################################
      ## Algorithm 6th step: update shrinkage parameters
      ##################################################
      beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
      beta_diff1 <- (beta1[sorted_indices[,1]] - beta1[sorted_indices[,2]])^2 + nugget
      
      # Vectorized calculation for zeta40, gamma40, zeta41, gamma41
      zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
      gamma40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
      
      # Vectorized calculation for zeta41, gamma41, others
      zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
      gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
      
      tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
      epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
      
      tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
      epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        R[j] <- r
        Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
        Tau[j, ] <- c(tau240, tau241)
        Lambda[j, ] <- c(zeta40, zeta41)
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
      
      
      #############################################################################
      ## Algorithm 4th step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      off_diag_val0 = 1/zeta40/tau240
      off_diag_val1 = 1/zeta41/tau241
      
      B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0  
      B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   
      spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0
      
      spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1             # nug_sig1 and nug_sig2 are
      spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2             # additional precision terms, i.e., 
      
      
      # Add diagonal weight matrix to B0_mat using spam
      B0_mat <- B0_mat + spam::diag.spam(w)
      
      # K^T W^(1/2), needed to construct K^TWK for B1)mat                     
      sqrt_w_K <- K
      sqrt_w_K@entries <- sqrt(w) * sqrt_w_K@entries
      
      B1_star_mat <- spam::bdiag(spam::spam(T0, p, p, sparse = T, doDiag=FALSE),  B1_mat)
      B1_star_mat <- (B1_star_mat) + t(sqrt_w_K) %*% sqrt_w_K
      
      
      ##################################################
      ## Algorithm 5th step: sample beta's
      ##################################################
      beta0 <- as.vector(spam::rmvnorm.canonical(1, 
                w*(z - Kbeta1), B0_mat))
      
      beta1 <- as.vector(spam::rmvnorm.canonical(1, 
                t(sqrt_w_K)%*%(sqrt(w)*(z - beta0)), B1_star_mat))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta  <- c(beta0, beta1)
      
      ##################################################
      ## Algorithm 6th step: update shrinkage parameters
      ##################################################
      beta1_spatial <- beta1[-c(1:p)]                   # extracting spatial slopes: beta1(s)'s and not the covariates
      beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
      beta_diff1 <- (beta1_spatial[sorted_indices[,1]] - beta1_spatial[sorted_indices[,2]])^2 + nugget
      
      zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
      gamma40 <-  MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
      
      zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
      gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
      
      tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
      epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
      
      tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
      epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        R[j] <- r
        Deviance[j]<- -2*sum(dnbinom(y,r,q,log=T))
        Tau[j, ] <- c(tau240, tau241)
        Lambda[j, ] <- c(zeta40, zeta41)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }
  
  if(verbose == "TRUE"){print(paste0("SpaceBF: NB model fitting completed!"))}
  
  return(NB = list(Beta = Beta, R = R, Lambda = Lambda, Tau = Tau, Deviance = Deviance)) 
}
