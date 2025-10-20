#' @title MCMC convergence plot
#'
#' @param x is a vector of MCMC samples, of dimension nIter x 1
#' @return a plot of MCMC samples with posterior mean in color green
#'
#' @export
#'

MCMC_plot <-function(x){
  plot(x, type = "l", ylim = c(-1*abs(max(x)), 1.2*abs(max(x))))
  abline(h = mean(x), col = "green")
}

#' @title Local-level or cell-specific credible interval (CI)
#'
#' @param Beta is a dataframe of dimension (nIter x N) with MCMC beta samples, for N cells/locations
#' @param local.p.ths is the size of the CI for each location
#' @return a N x 1 vector of 0/1, depending on the CI for a location contains 0 or not.
#'
#' @export
#'

local_credible <-function(Beta, local.p.ths = 0.9){
  CI_of_betas = t(apply(Beta, 2, 
  function(x) unlist(bayestestR::ci(x, ci = local.p.ths, method = "ETI"))[-1]))
  zero_finder = apply(CI_of_betas, 1, 
  function(x){ifelse(x[1] <= 0 & x[2] >= 0, 0, 1)})
  return(zero_finder)
}

#' @title Perform global significance test
#'
#' @param Beta is a dataframe of dimension (nIter x N) with MCMC beta samples, for N cells/locations
#' @param summary is the type of posterior summary measure to use, either "median" or "mean"
#' @param local.p.ths is the size of the credible interval for location-specific tests
#' @return a list of global beta estimate (mean or median across locations),
#' global p-value, number of locations where local CI does not contain 0 (i.e., significant), 
#' and Geweke convergence-diagnostic statistic for the global beta
#'
#' @export
#'

global_res <- function(Beta, summary = "median", local.p.ths = 0.95){
  
  # finding cells with non-zero local estimate based on credible interval
  CI <- local_credible(Beta, local.p.ths = local.p.ths)
  
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
#' @return a data.frame of local beta estimates (N rows for N locations), local p-values, and 
#' Geweke convergence-diagnostic statistic for each local beta.
#'
#' @export
#'

local_res <- function(Beta, summary = "median"){
  
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
  sigfinder <- local_credible(Beta, local.p.ths)
  
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
    geom_sf(aes(color = beta), size = 1) +  
    labs(color = TeX(sprintf("$\\beta_1(s) - \\bar{\\beta}_1$"))) + 
    ggtitle(TeX(sprintf("$\\bar{\\beta}_1 = %.2f$, $\\sigma_{\\beta} = %.2f$", mean_beta, SD_beta))) +            #TeX(sprintf("$\\sigma_B = %f $", SD_beta))) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="none", legend.text=element_text(size=9), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
    ) + sc_LH
  
  
  p.sf$beta[sigfinder == 0] <- -50 # arbitrary value to grey out insignificant results

  betas2 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = beta), size = 1) +  
    labs(color = TeX(sprintf("$\\frac{\\beta_1(s) - \\bar{\\beta}_1}{\\sigma_{\\beta}}$"))) + 
    ggtitle(TeX(sprintf("$\\bar{\\beta}_1 = %f $", mean_beta))) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="right", legend.text=element_text(size=9), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12)) + sc_LH
  
  
  p1 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = y1), size = 1) +   
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
    geom_sf(aes(color = y2), size = 1) +   
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


#' @title Function for comparing HS vs, ICAR local coefficients and p-values
#'
#' @param L_results_first is the output of local_res() function on the slope coefficients from the first model
#' @param L_results_second is the output of local_res() function on the slope coefficients from the second model
#' @return The kernel covariance matrix of dimension N x N, between N locations
#'
#' @export

plot_local_comparisons <- function(L_results_first = L_results, L_results_second = L_results_ICAR, point_alpha = 0.6, point_size = 0.5) {
  
  df <- data.frame(
    beta      = as.numeric(L_results_first$local.beta),
    beta_icar = as.numeric(L_results_second$local.beta),
    mlp       = -log10(pmax(as.numeric(L_results_first$local.pval), 1e-10)),
    mlp_icar  = -log10(pmax(as.numeric(L_results_second$local.pval), 1e-10))
  )
  
  r1 <- cor(df$beta, df$beta_icar, use = "complete.obs")
  r2 <- cor(df$mlp,  df$mlp_icar,  use = "complete.obs")
  
  library(patchwork)
  library(ggplot2)
  
  # common theme to normalize spacing
  thm <- theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title.position = "plot",
      plot.title = element_text(margin = margin(b = 4)),
      plot.subtitle = element_text(margin = margin(b = 6)),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  p1 <- ggplot(df, aes(x = beta, y = beta_icar)) +
    geom_point(alpha = point_alpha, size = point_size, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal() +
    labs(
      title = "Local β: First vs Second Method",
      subtitle = sprintf("r = %.2f", r1),
      x = "Local β (first method)",
      y = "Local β (second method)"
    ) +
    thm
  
  p2 <- ggplot(df, aes(x = mlp, y = mlp_icar)) +
    geom_point(alpha = point_alpha, size = point_size, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal() +
    labs(
      title = expression(paste(-log[10], " p: First vs Second Method")),
      subtitle = sprintf("r = %.2f", r2),
      x = expression(-log[10](p)~"(First)"),
      y = expression(-log[10](p)~"(Second)")
    ) +
    thm
  
  # arrange with equal widths
  wrapped <- p1 + p2 + plot_layout(nrow = 1, widths = c(1, 1))
  return(wrapped)
}


