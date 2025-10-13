#' @import Matrix
#' @title Function for the Gaussian fused model
#'
#' @param y1 is the expression vector of gene 1, gene 1 will be regressed on gene 2
#' @param y2 is the expression vector of gene 2
#' @param X is the matrix of additional cell-level covariates
#' @param nIter is the number of MCMC iterations
#' @param beta_thres is a hard threshold on the absolute values of beta's
#' @param nug_sig1 is the additional precision on beta0's
#' @param nug_sig2 is the additional precision on beta1's
#' @param scale_by_sigma is TRUE or FALSE, whether to scale the beta's by variance of y or not 
#' @param verbose is TRUE or FALSE, to show or suppress progression updates
#' @return It returns the estimated beta's, sigma, shrinkage parameters, and deviance
#'
#' @export



Gauss_ICAR_model<-function(y1, y2, X = NULL, G, nIter = 5000, beta_thres = 10, 
                      nug_sig1 = 0.05, nug_sig2 = 0.1, scale_by_sigma = FALSE, 
                      verbose = "TRUE"){
  
  # variable definitions
  y <- y1
  K <- spam::as.spam(diag(scale(log1p(y2))[, 1]))       # design matrix with scaled log1p(y2)'s on the diagonal
  p <- ifelse(!is.null(X), ncol(X), 0)                  # number of covariates
  N <- length(y)                                        # number of samples
  p_mst <- length(which(G == 1))                        # twice the number of connections
  half_p_mst <- p_mst/2
  
  # non-zero indices from the adjacency matrix
  W_sp <- spam::as.spam(G)                              # ICAR adjacency is just G
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
  Sigma2 <- matrix(0,lastit,1)                       # stores variance of y
  
  # Initialize
  beta0 <- beta1 <- rep(1, N)      # beta0's and beta1's
  if(p != 0){beta_cov = rep(0, p)  # slopes for extra covariates if any
  beta1 <- c(beta_cov, beta1)      # putting the slopes with beta1's
  beta <- c(beta0, beta1)
  }else{beta<-c(beta0, beta1)}
  
  sigma2 <- 1                      # variance of y
  
  # Building ICAR precision structure
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
  
  if(p == 1){
    B1_star_mat = Matrix::bdiag(1,  B1_mat)
  }else if(p > 1){B1_star_mat = spam::bdiag(spam::spam(0, p, p, sparse = T, doDiag=FALSE),  B1_mat)}
  

  T0 <- 0.1                                           # prior precision of the additional covariates
  l <- rep(0, N)                                     # latent vector for CRT-based update of r
  a <- b <- 1                                        # gamma prior parameters on r for CRT                      
  Acc <- 0; s <- 0.01                                # counter for MH-based update of r, and proposal SD 
  del1 <- 1; del2 <- 1                               # extra prior parameters on shrinkage parameters, unused
  nugget <- 1e-8
  
  
  if(verbose == "TRUE"){print(paste0("SpaceBF: Gaussian model fitting has started."))}
  
  if(scale_by_sigma == FALSE & p == 0){
    # algorithm without scaling the coefficients by the variance of y
    KTy <- K %*% y                                                # computing and keeping some quantities fixed
    diag_K <- spam::diag(K)                                  # vectorized K
    KTK<-K^2                                                      # K is diagonal, K^2 = K^TK
    ID_mat <- spam::diag(1, N)                               # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the ICAR precision matrices
      #############################################################################

      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + nug_sig1*sigma2 

      B1_mat <- tau1 * Qicar 
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2*sigma2 
      
      B0_mat <- (B0_mat) + ID_mat 
      B1_mat <- (B1_mat) + KTK 
      
      ##################################################
      ## Algorithm 2nd step: sample beta's
      ##################################################
      beta0 <- as.vector(spam::rmvnorm.canonical(1, (y - diag_K*beta1)/sigma2, 
               B0_mat/sigma2))
      
      beta1 <- as.vector(spam::rmvnorm.canonical(1, (KTy - diag_K*beta0)/sigma2, 
               B1_mat/sigma2))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta <- c(beta0, beta1)
      
      ##################################################
      ## Algorithm 3rd step: sample sigma2
      ##################################################
      RSS <- sum((y - beta0 - diag_K*beta1)^2)
      sigma2 <- MCMCpack::rinvgamma(1, shape = (N/2 + 1), scale = (0.5*RSS)) 
      
      #############################################################################
      ## Algorithm 4th step: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   # f' Q f
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      quad1 <- as.numeric(t(beta1) %*% (Qicar %*% beta1))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(scale_by_sigma == FALSE & p > 0){
    # algorithm without scaling the coefficients by the variance of y
    K <- spam::cbind(X, K)                                   # redefine K to include the covariates matrix X
    KTy <- K %*% y                                                # computing and keeping some quantities fixed
    diag_K <- spam::diag(K)                                       # vectorized K
    KTK<- t(K) %*% K                                                    
    ID_mat <- spam::diag(1, N)                                    # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the ICAR precision matrices
      #############################################################################
      
      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + nug_sig1*sigma2 
      
      B1_mat <- tau1 * Qicar 
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2*sigma2 
      
      B0_mat <- (B0_mat) + ID_mat 
      B1_star_mat <- spam::bdiag(spam::spam(T0, p, p),  B1_mat)
      B1_star_mat <- (B1_star_mat) + KTK 
      
      ##################################################
      ## Algorithm 2nd step: sample beta's
      ##################################################
      beta1_spatial <- beta1[-c(1:p)]                   # extracting spatial slopes: beta1(s)'s and not the covariates
      Kbeta1 <- K %*% beta1
      KTbeta0 <- c(crossprod(X, beta0), diag_K*beta0)   # efficient crossprod(K, beta0)
      beta0 <- as.vector(spam::rmvnorm.canonical(1, (y - Kbeta1)/sigma2, 
                                                 B0_mat/sigma2))
      
      beta1 <- as.vector(spam::rmvnorm.canonical(1, (KTy - KTbeta0)/sigma2, 
                                                 B1_star_mat/sigma2))
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta <- c(beta0, beta1)
      
      ##################################################
      ## Algorithm 3rd step: sample sigma2
      ##################################################
      RSS <- sum((y - beta0 - Kbeta1)^2)
      sigma2 <- MCMCpack::rinvgamma(1, shape = (N/2 + 1), scale = (0.5*RSS)) 
      
      #############################################################################
      ## Algorithm 4th step: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   # f' Q f
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      quad1 <- as.numeric(t(beta1_spatial) %*% (Qicar %*% beta1_spatial))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(scale_by_sigma == TRUE & p == 0){
    # algorithm without scaling the coefficients by the variance of y
    KTy <- K %*% y                                                # computing and keeping some quantities fixed
    diag_K <- spam::diag(K)                                  # vectorized K
    KTK<-K^2                                                      # K is diagonal, K^2 = K^TK
    ID_mat <- spam::diag(1, N)                               # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the ICAR precision matrices
      #############################################################################
      
      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + nug_sig1
      
      B1_mat <- tau1 * Qicar 
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2
      
      B0_mat_noK <- B0_mat
      B1_mat_noK <- B1_mat
      
      B0_mat <- (B0_mat) + ID_mat 
      B1_mat <- (B1_mat) + KTK 

      ##################################################
      ## Algorithm 2nd step: sample beta's
      ##################################################
      beta0 <- as.vector(spam::rmvnorm.canonical(1, (y - diag_K*beta1)/sigma2, 
                                                 B0_mat/sigma2))
      
      beta1 <- as.vector(spam::rmvnorm.canonical(1, (KTy - diag_K*beta0)/sigma2, 
                                                 B1_mat/sigma2))
      
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta <- c(beta0, beta1)
      
      ##################################################
      ## Algorithm 3rd step: sample sigma2
      ##################################################
      RSS <- sum((y - beta0 - diag_K*beta1)^2)
      betaBinvbeta <- sum(beta0 * crossprod(B0_mat_noK, beta0)) +
        sum(beta1 * crossprod(B1_mat_noK, beta1)) 
      
      sigma2 <- (MCMCpack::rinvgamma(1, shape = (3*N/2 + 1),
                 scale = (0.5*RSS + 0.5*betaBinvbeta))) 
      
      #############################################################################
      ## Algorithm 4th step: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   # f' Q f
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      quad1 <- as.numeric(t(beta1) %*% (Qicar %*% beta1))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(scale_by_sigma == TRUE & p > 0){
    # same algorithm as above but with additional covariate adjustment
    
    # algorithm without scaling the coefficients by the variance of y
    K <- spam::cbind(X, K)                                   # redefine K to include the covariates matrix X
    KTy <- t(K) %*% y                                  # computing and keeping some quantities fixed
    diag_K <- spam::diag(K)                       # vectorized K
    KTK<- K %*% K                                      # K is diagonal, K^2 = K^TK
    ID_mat <- spam::diag(1, N)                    # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      #############################################################################
      ## Algorithm 1st step: update the ICAR precision matrices
      #############################################################################
      
      ## beta0 precision: tau0*Qicar + diag(w) + ridge
      B0_mat <- tau0 * Qicar
      spam::diag(B0_mat) <- spam::diag(B0_mat) + nug_sig1
      
      B1_mat <- tau1 * Qicar 
      spam::diag(B1_mat) <- spam::diag(B1_mat) + nug_sig2
      
      B0_mat_noK <- B0_mat
      B1_mat_noK <- B1_mat
      
      B0_mat <- (B0_mat) + ID_mat 
      B1_star_mat <- spam::bdiag(spam::spam(T0, p, p),  B1_mat)
      B1_star_mat <- (B1_star_mat) + KTK 
      
      
      ##################################################
      ## Algorithm 2nd step: sample beta's
      ##################################################
      Kbeta1 <- K%*%beta1
      KTbeta0 <- c(crossprod(X, beta0), diag_K*beta0)   # efficient crossprod(K, beta0)
      beta0 <- as.vector(spam::rmvnorm.canonical(1, (y - Kbeta1)/sigma2, 
                                                 B0_mat/sigma2))
      
      beta1 <- as.vector(spam::rmvnorm.canonical(1, (KTy - crossprod(K, beta0))/sigma2, 
                                                 B1_star_mat/sigma2))
      beta0 <- pmax(pmin(beta0, beta_thres), -beta_thres)
      beta1 <- pmax(pmin(beta1, beta_thres), -beta_thres)
      beta <- c(beta0, beta1)
      beta1_spatial <- beta1[-c(1:p)]                   # extracting spatial slopes: beta1(s)'s and not the covariates
      
      ##################################################
      ## Algorithm 3rd step: sample sigma2
      ##################################################
      RSS <- sum((y - beta0 - Kbeta1)^2)
      betaBinvbeta <- sum(beta0 * crossprod(B0_mat_noK, beta0)) +
        sum(beta1_spatial*crossprod(B1_mat_noK,  beta1_spatial)) 
      
      sigma2 <- (MCMCpack::rinvgamma(1, shape = (3*N/2 + 1),
                  scale = (0.5*RSS + 0.5*betaBinvbeta))) 
      
      #############################################################################
      ## Algorithm 4th step: update ICAR scales tau0 and tau1 (Gamma)
      #############################################################################
      quad0 <- as.numeric(t(beta0) %*% (Qicar %*% beta0))   # f' Q f
      tau0  <- rgamma(1, shape = a_tau0 + 0.5*(N - n_comp), rate = b_tau0 + 0.5*quad0)
      
      quad1 <- as.numeric(t(beta1_spatial) %*% (Qicar %*% beta1_spatial))
      tau1  <- rgamma(1, shape = a_tau1 + 0.5*(N - n_comp), rate = b_tau1 + 0.5*quad1)
      
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau0, tau1)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }
  if(verbose == "TRUE"){print(paste0("SpaceBF: Gaussian model fitting completed!"))}
  return(Gauss = list(Beta = Beta, Sigma2 = Sigma2, Lambda = Lambda, Tau = Tau,
                      Deviance = Deviance))  
}

