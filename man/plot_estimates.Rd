% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Result_generator.R
\name{plot_estimates}
\alias{plot_estimates}
\title{Generates plots of the estimated slope surface along with the expression profiles}
\usage{
plot_estimates(
  Beta,
  y1,
  y2,
  coords,
  which.model = "NB",
  beta_ths = 1.9,
  exp_ths_min = 0,
  exp_ths_max = 2,
  local.p.ths = 0.9
)
}
\arguments{
\item{Beta}{is a dataframe of dimension (nIter x N) with MCMC beta samples, for N cells/locations}

\item{y1}{is the expression vector of gene 1}

\item{y2}{is the expression vector of gene 2}

\item{coords}{is the 2D matrix of xy coordinates}

\item{which.model}{denotes the model to be used, "NB" or "Gaussian"}

\item{beta_ths}{is the threshold to hide extreme values of betas}

\item{exp_ths_min}{is the minimum value of the expression profiles}

\item{exp_ths_max}{is the maximum value of the expression profiles}

\item{local.p.ths}{is the size of the credible interval}
}
\value{
a ggplot object with 4 sub-figures, including expression of y1 and y2, 
and estimated "z-score"-type beta surface and its truncated version based on local significant
}
\description{
Generates plots of the estimated slope surface along with the expression profiles
}
