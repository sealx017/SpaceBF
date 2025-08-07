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



Gauss_model<-function(y1, y2, X = NULL, G, nIter = 5000, beta_thres = 10, 
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
  

  if(verbose == "TRUE"){print(paste0("SpaceBF: Gaussian model fitting has started."))}

  if(scale_by_sigma == FALSE & p == 0){
    # algorithm without scaling the coefficients by the variance of y
    KTy <- K %*% y                                                # computing and keeping some quantities fixed
    diag_K <- spam::diag_spam(K)                                  # vectorized K
    KTK<-K^2                                                      # K is diagonal, K^2 = K^TK
    ID_mat <- spam::diag_spam(1, N)                               # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      off_diag_val0 <- 1/zeta40/tau240
      off_diag_val1 <- 1/zeta41/tau241
      
      B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0   # adjust only non-zero 
      B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   # elements of precision matrices
      spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0
      
      spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1*sigma2             # nug_sig1 and nug_sig2 are
      spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2*sigma2             # additional precision terms, i.e., 
      
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

      ##################################################
      ## Algorithm 4th step: update shrinkage parameters
      ##################################################
      beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
      beta_diff1 <- (beta1[sorted_indices[,1]] - beta1[sorted_indices[,2]])^2 + nugget
      
      # Vectorized calculation for zeta40, gamma40, zeta41, gamma41
      zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
      gamma40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
      
      zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
      gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
      
      # Single step calculations for tau240, epsilon20, tau241, epsilon21
      tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
      epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
      
      tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
      epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau240, tau241)
        Lambda[j, ] <- c(zeta40, zeta41)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(scale_by_sigma == FALSE & p > 0){
    # algorithm without scaling the coefficients by the variance of y
    KTy <- K %*% y                                                # computing and keeping some quantities fixed
    diag_K <- spam::diag_spam(K)                                  # vectorized K
    KTK<-K^2                                                      # K is diagonal, K^2 = K^TK
    ID_mat <- spam::diag_spam(1, N)                               # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      off_diag_val0 <- 1/zeta40/tau240
      off_diag_val1 <- 1/zeta41/tau241
      
      B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0   # adjust only non-zero 
      B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   # elements of precision matrices
      spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0
      
      spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1*sigma2             # nug_sig1 and nug_sig2 are
      spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2*sigma2             # additional precision terms, i.e., 
      
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
      
      ##################################################
      ## Algorithm 4th step: update shrinkage parameters
      ##################################################
      beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
      beta_diff1 <- (beta1_spatial[sorted_indices[,1]] - beta1_spatial[sorted_indices[,2]])^2 + nugget
      
      zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
      gamma40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
      
      zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
      gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
      
      tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
      epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
      
      tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
      epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau240, tau241)
        Lambda[j, ] <- c(zeta40, zeta41)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(scale_by_sigma == TRUE & p == 0){
    # algorithm without scaling the coefficients by the variance of y
    KTy <- K %*% y                                                # computing and keeping some quantities fixed
    diag_K <- spam::diag_spam(K)                                  # vectorized K
    KTK<-K^2                                                      # K is diagonal, K^2 = K^TK
    ID_mat <- spam::diag_spam(1, N)                               # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      off_diag_val0 <- 1/zeta40/tau240
      off_diag_val1 <- 1/zeta41/tau241
      
      B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0   # adjust only non-zero 
      B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   # elements of precision matrices
      spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0
      
      spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1            # nug_sig1 and nug_sig2 are
      spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2            # additional precision terms, i.e., 
      
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
      
      ##################################################
      ## Algorithm 4th step: update shrinkage parameters
      ##################################################
      beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
      beta_diff1 <- (beta1[sorted_indices[,1]] - beta1[sorted_indices[,2]])^2 + nugget
      
      zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
      gamma40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
      
      zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
      gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
      
      tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
      epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
      
      tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
      epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau240, tau241)
        Lambda[j, ] <- c(zeta40, zeta41)
      }
      if(verbose == "TRUE"){
        svMisc::progress(i, nIter, progress.bar = FALSE)
      }else{if(i%%500 == 0){print(paste0(i, " / ", nIter))}
      }
    }
  }else if(scale_by_sigma == TRUE & p > 0){
    # same algorithm as above but with additional covariate adjustment

    # algorithm without scaling the coefficients by the variance of y
    K <- cbind(X, K)                                   # redefine K to include the covariates matrix X
    KTy <- t(K) %*% y                                  # computing and keeping some quantities fixed
    diag_K <- spam::diag_spam(K)                       # vectorized K
    KTK<- K %*% K                                      # K is diagonal, K^2 = K^TK
    ID_mat <- spam::diag_spam(1, N)                    # identity matrix
    
    for (i in 1:nIter){
      #############################################################################
      ## Algorithm 1st step: update the precision matrices of the beta0 and beta1's
      #############################################################################
      off_diag_val0 <- 1/zeta40/tau240
      off_diag_val1 <- 1/zeta41/tau241
      
      B0_mat[sorted_indices]  <- B0_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val0   # adjust only non-zero 
      B1_mat[sorted_indices]  <- B1_mat[sorted_indices[, c(2, 1)]] <- -off_diag_val1   # elements of precision matrices
      spam::diag(B0_mat)  <-  spam::diag(B1_mat)  <- 0
      
      spam::diag(B0_mat) <- abs(spam::rowSums(B0_mat)) + nug_sig1            # nug_sig1 and nug_sig2 are
      spam::diag(B1_mat) <- abs(spam::rowSums(B1_mat)) + nug_sig2             # additional precision terms, i.e., 
      
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
      
      ##################################################
      ## Algorithm 4th step: update shrinkage parameters
      ##################################################
      beta_diff0 <- (beta0[sorted_indices[,1]] - beta0[sorted_indices[,2]])^2 + nugget
      beta_diff1 <- (beta1_spatial[sorted_indices[,1]] - beta1_spatial[sorted_indices[,2]])^2 + nugget
      
      zeta40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma40 + beta_diff0/(2 * tau240))
      gamma40 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta40)
      
      zeta41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1/gamma41 + beta_diff1/(2 * tau241))
      gamma41 <- MCMCpack::rinvgamma(half_p_mst, 1, scale = 1 + 1/zeta41)
      
      tau240 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon20 + sum(beta_diff0 / zeta40) / 2)
      epsilon20 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau240)
      
      tau241 <- MCMCpack::rinvgamma(1, (half_p_mst + 1) / 2, scale = 1/epsilon21 + sum(beta_diff1 / zeta41) / 2)
      epsilon21 <- MCMCpack::rinvgamma(1, 1, scale = 1 + 1/tau241)
      
      if (i> burn & i%%thin==0) {
        j <- (i-burn)/thin
        Beta[j, ] <- beta
        Sigma2[j] <- sigma2
        Deviance[j,]<- -2*(-RSS/2/sigma2 - N/2*log(sigma2))
        Tau[j, ] <- c(tau240, tau241)
        Lambda[j, ] <- c(zeta40, zeta41)
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

