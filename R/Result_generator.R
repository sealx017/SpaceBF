#' @title Local-level or cell-specific credible interval
#'
#' @param Beta is a dataframe of dimension (nIter x N) with MCMC beta samples, for N cells/locations
#' @param p.ths is the size of the credible interval
#' @return a heatmap object from the pheatmap package.
#'
#' @export
#'

local_credible <-function(Beta, p.ths = 0.95){
  CI_of_betas = t(apply(Beta, 2, 
  function(x) unlist(bayestestR::ci(x, ci = p.ths, method = "ETI"))[-1]))
  zero_finder = apply(CI_of_betas, 1, 
  function(x){ifelse(x[1] <= 0 & x[2] >= 0, 0, 1)})
  return(zero_finder)
}

#' @title Perform global significance test
#'
#' @param Beta is a dataframe of dimension (nIter x N) with MCMC beta samples, for N cells/locations
#' @param summary is the type of posterior summary measure to use, either "median" or "mean"
#' @param local.p.ths is the size of the credible interval
#' @return a heatmap object from the pheatmap package.
#'
#' @export
#'

global_res <- function(Beta, summary = "median", local.p.ths = 0.95){
  
  # finding cells with non-zero local estimate based on credible interval
  CI <- local_credible(Beta, p.ths = local.p.ths)
  
  # ROPE p-value based on the value of the posterior mean/median of beta's
  if(summary == "median"){
  global_slope_dir <- apply(Beta, 1, median)         # median of the local-level beta's
  global_pd <- min(1, 2*(1 - bayestestR::p_direction(global_slope_dir,
               method = "KernSmooth")$pd))           # global p-value    
  global_slope_dir_summ <- median(global_slope_dir)  # overall summary or global association index
  }else{
  global_slope_dir <- apply(Beta, 1, mean)           # mean of the local-level beta's
  global_pd <- min(1, 2*(1 - bayestestR::p_direction(global_slope_dir,
                 method = "KernSmooth")$pd))         # global p-value 
  global_slope_dir_summ <- mean(global_slope_dir)    # overall summary or global association index
  }
  
  # create result output
  global_res = data.frame(global.beta = global_slope_dir_summ, 
  global.pval = global_pd, num.nonzero.local = sum(CI), 
  Geweke.global.beta = abs(coda::geweke.diag(global_slope_dir, 0.1, 0.5)$z))

  return(global_res)
}


#' @title Perform local significance test
#'
#' @param Beta is a dataframe of dimension (nIter x N) with MCMC beta samples, for N cells/locations
#' @param summary is the type of posterior summary measure to use, either "median" or "mean"
#' @param local.p.ths is the size of the credible interval
#' @return a heatmap object from the pheatmap package.
#'
#' @export
#'

local_res <- function(Beta, summary = "median", local.p.ths = 0.95){
  
  # compute posterior mean/median of beta's
  if(summary == "median"){
    local_slope_dir <- apply(Beta, 2, median)     # median of the local-level beta's
  }else{
    local_slope_dir <- apply(Beta, 1, mean)       # mean of the local-level beta's
  }
  # compute local p-values
  local_pd <- apply(Beta, 2, function(x){min(1, 2*(1 - bayestestR::p_direction(x,
                method = "KernSmooth")$pd))})     
  # check Geweke convergence
  Geweke_local =  apply(Beta, 2, function(x) {abs(coda::geweke.diag(x, 0.1, 0.5)$z)})
  
  # create result output
  local_res = data.frame(local.beta = local_slope_dir, 
                          local.pval = local_pd,  
                          Geweke.local.beta = Geweke_local)
  
  return(local_res)
}

#' @import RColorBrewer
#' @import ggplot2
#' @import circlize
#' @import latex2exp
#' @import patchwork
#' @import sf
#' @title Generates plots of the estimated slope surface along with the expression profiles
#'
#' @param Beta is a dataframe of dimension (nIter x N) with MCMC beta samples, for N cells/locations
#' @param y1 is the expression vector of gene 1
#' @param y2 is the expression vector of gene 2
#' @param coords is the 2D matrix of xy coordinates
#' @param which.model denotes the model to be used, "NB" or "Gaussian"
#' @param beta_ths is the threshold to hide extreme values of betas
#' @param exp_ths_min is the minimum value of the expression profiles
#' @param exp_ths_max is the maximum value of the expression profiles
#' @param local.p.ths is the size of the credible interval
#' @return a ggplot object with 4 sub-figures, including expression of y1 and y2, 
#' and estimated "z-score"-type beta surface and its truncated version based on local significant
#' @export
#'



plot_estimates <- function(Beta, y1, y2, coords, which.model = "NB", 
                  beta_ths = 1.9, exp_ths_min = 0, exp_ths_max = 2,
                  local.p.ths = 0.9)
{
  
  if(which.model == "NB"){
  location_data_with_beta <- data.frame(coords, log(y1+1), log(y2+1))
  }else{
  location_data_with_beta <- data.frame(coords, y1, y2)
  }
  colnames(location_data_with_beta) <- c("x", "y", "y1", "y2")
  
  location_data_with_beta$beta <- colMeans(Beta)
  sigfinder <- local_credible(Beta, p.ths = local.p.ths)
  
  p.sf<- sf::st_as_sf(location_data_with_beta, coords = c("x", "y"))
  
  mean_beta <- round(mean(p.sf$beta), 3)                    # mean of the beta estimates
  SD_beta <- round(sd(p.sf$beta), 3)                        # SD of the beta estimates
   
  p.sf$beta <- scale(p.sf$beta, scale = T)[,1]              # plot the z-scores for interpretability
  p.sf$beta <- pmax(pmin(p.sf$beta, beta_ths), -beta_ths)   # trimming outliers
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  sc_LH <- scale_colour_gradientn(colours = myPalette(10), 
  limits = c(-beta_ths, beta_ths), breaks = seq(-beta_ths, beta_ths, length.out = 5))
  sc_LH$breaks <- round(sc_LH$breaks, digits = 2)
  
  p.sf$y1 <- pmax(pmin(p.sf$y1, exp_ths_max), -exp_ths_min)
  p.sf$y2 <- pmax(pmin(p.sf$y2, exp_ths_max), -exp_ths_min)
  
  betas <- ggplot(data =  p.sf) +
    geom_sf(aes(color = beta), size = 0.5) +  
    labs(color = TeX(sprintf("$\\beta_1(s) - \\bar{\\beta}_1$"))) + 
    ggtitle(TeX(sprintf("$\\bar{\\beta}_1 = %.2f$, $\\sigma_{\\beta} = %.2f$", mean_beta, SD_beta))) +            #TeX(sprintf("$\\sigma_B = %f $", SD_beta))) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="none", legend.text=element_text(size=9), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
    ) + sc_LH
  
  
  p.sf$beta[sigfinder == 0] <- -50

  betas2 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = beta), size = 0.5) +  
    labs(color = TeX(sprintf("$\\frac{\\beta_1(s) - \\bar{\\beta}_1}{\\sigma_{\\beta}}$"))) + 
    ggtitle(TeX(sprintf("$\\bar{\\beta}_1 = %f $", mean_beta))) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="right", legend.text=element_text(size=9), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12)) + sc_LH
  
  
  p1 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = y1), size = 0.5) +   
    scale_color_gradientn(colours = (viridis::viridis(10)), 
    limits = c(exp_ths_min, exp_ths_max)) +
    labs(color = ifelse(which.model == "NB", "log1p(expr)", "expr")) + 
    ggtitle(paste0("Gene: 1")) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), 
      legend.position="left", legend.text=element_text(size=9), 
      legend.title = element_text(size=10),
      legend.key.size = unit(0.5, 'cm'),
      plot.title = element_text(hjust = 0.5, size = 12))
  
  p2 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = y2), size = 0.5) +   
    scale_color_gradientn(colours = (viridis::viridis(10)), 
    limits = c(exp_ths_min, exp_ths_max)) +
    labs(color = ifelse(which.model == "NB", "log1p(expr)", "expr")) + 
    ggtitle(paste0("Gene: 2")) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
      legend.position="none", legend.text=element_text(size=8), 
      legend.title = element_text(size=10),
      legend.key.size = unit(0.5, 'cm'),
      plot.title = element_text(hjust = 0.5, size = 12))
  
  all_plots = wrap_plots(list(p1, p2, betas, betas2))
  return(all_plots)
}


