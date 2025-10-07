#' @title Main parent function for fitting different models
#'
#' @param y1 is the expression vector of gene 1, gene 1 will be regressed on gene 2
#' @param y2 is the expression vector of gene 2
#' @param X is the matrix of additional cell-level covariates
#' @param which.model denotes the model to be used, "NB" or "Gaussian"
#' @param nIter is the number of MCMC iterations
#' @param beta_thres is a hard threshold on the absolute values of beta's
#' @param nug_sig is the additional precision parameters on beta's, either NULL or a vector of size two
#' @param which.r.sampler is the choice of sampling r in the NB model, either CRT or MH
#' @param scale_by_sigma is TRUE or FALSE for the Gaussian model, whether to scale the beta's by variance of y or not 
#' @param verbose is TRUE or FALSE, to show or suppress progression updates
#' @return It returns the estimated beta's, r or sigma, shrinkage parameters, and deviance
#'
#' @export


SpaceBF <-function(y1, y2, X = NULL, G, which.model = "NB", nIter = 5000, 
          beta_thres = 10, nug_sig = NULL, which.r.sampler = "CRT", 
          scale_by_sigma = FALSE, verbose = "TRUE"){
  
  # dimension checks
  if(!is.null(dim(y1))){
    y1 = c(y1[, 1])
    print("y1 is a matrix, first column is considered. ")
  }
  
  if(!is.null(dim(y2))){
    y2 = c(y2[, 1])
    print("y2 is a matrix, first column is considered. ")
  }
  
  if(!is.null(X)){
    if(is.null(dim(X))){
      X = as.matrix(X)
      print("X is a vector, making it a matrix. ")
    }
  }
  
  # feasibility checks
  if(length(y1) != length(y2)){print("Stop! Length of y1 and y2 does not match."); break;}
  if(!is.null(X)){if(length(y1) != nrow(X)){print("Stop! Length of y1 and number of rows of X do not match."); break;}}
  if(Matrix::isSymmetric(G) == "FALSE" | max(G) > 1){print("Stop! G needs to be symmetric and binary."); break;}
  if(length(y1) != nrow(G)){print("Stop! Length of y1 and number of rows of G do not match."); break;}
  if(verbose == "TRUE" & length(y1) > 2000){print(paste0("Number of cells is more than 2000, you might want to reduce the # MCMC iterations from ", nIter, " to a smaller value."))}
  if(verbose == "TRUE" & length(table(y1)) > length(y1)*(3/4) & which.model == "NB"){print(paste0("The data is possibly continuous and the Gaussian model might be ideal."))}
  
  if(is.null(nug_sig)){
  if(mean(y1 == 0) > 0.5 | mean(y2 == 0) > 0.75){   # if there are a lot of 0 values either in y1 or y2,
                                                    # amplifying the non-identifiability issue,                                              
                                                    # put extra regularization on the beta's
    nug_sig1 <- 0.2
    nug_sig2 <- 0.2   
  }else{
    nug_sig1 <- 0.1
    nug_sig2 <- 0.1
  }}else if(is.vector(nug_sig)){
    if(length(nug_sig) > 1){
    nug_sig1 <- nug_sig[1]
    nug_sig2 <- nug_sig[2]
    }else{nug_sig1 <- nug_sig2 <- nug_sig}
  }else{print("Stop! Provide precision parameters as a vector of two or keep it NULL."); break;}
  
  if(which.model == "NB"){
    
    fitted_model <- NB_model(y1, y2, X, G, nIter, beta_thres, nug_sig1, 
                    nug_sig2, which.r.sampler, verbose)
  
  }else{
  
    fitted_model <- Gauss_model(y1, y2, X, G, nIter, beta_thres, nug_sig1, 
                    nug_sig2, scale_by_sigma, verbose)
  }
  
  return(fitted_model)
}
  
  