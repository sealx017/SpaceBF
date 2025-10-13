#' @title Function for MST construction
#'
#' @param x is a dist object or an Euclidean distance matrix
#' @param type of the MST, random or direct
#' @return It returns the MST inside a list
#'
#' @export

MST_construction <- function(x, type = "random"){
  if(type == "random"){
  # construct random MST
  x_d <- as.matrix(x)
  x_d[col(x_d) < row(x_d)] <- 0
  distance_values = x_d[col(x_d) > row(x_d)]
  
  set.seed(2024)
  x_d[col(x_d) > row(x_d)] <- x_d[col(x_d) > row(x_d)] + 
    runif(length(x_d[col(x_d) > row(x_d)]), 0, 0.5)  # adding arbitrary, small random number to avoid ties (although fossil package already does something similar)
  x_d <- x_d + t(x_d)
  
  mstobj <- fossil::dino.mst(x_d)
  }else{
  # construct MST directly
  mstobj <- fossil::dino.mst(x)
  }
  return(list(mstobj, type))
}

#' @import igraph
#' @import ggraph
#' @import ggplot2
#' @title Function for plotting the MST 
#'
#' @param mst_out is the list output from function MST_construction
#' @param coords is the 2D matrix of coordinates
#' @param title.hjust denotes the position of the plot title, default is centered
#' @param title.size  denotes the size of the plot title
#' @param node.size denotes the size of each node/location of the MST
#' @return The MST plot
#'
#' @export

print_mst <- function(mst_out, coords, title.hjust = 0.5, title.size = 10, node.size = 1){
  adj.mat <- mst_out[[1]]
  type <- mst_out[[2]]
  
  graph = graph_from_adjacency_matrix(adj.mat, mode = "undirected")
  vertex_attr(graph, name = 'x') <- coords[,1]
  vertex_attr(graph, name = 'y') <- coords[,2]
  lay = create_layout(graph, layout = coords)
  
  mstgraph = ggraph(lay) + 
    geom_edge_link(edge_colour = "grey40") + 
    geom_node_point(size = node.size) +
    ggtitle(ifelse(type == "random", "Random MST", "MST")) +
    theme(plot.title = element_text(hjust = title.hjust, 
          size = title.size))
  
  print(mstgraph)
}

#' @title Function for kernel matrix construction based on the distance matrix 
#'
#' @param x is a dist object or an N x N Euclidean distance matrix between N locations
#' @param l is the lengthscale, controlling the exponential decay
#' @param type denotes the type of kernel function, exponential or radial basis function (RBF)
#' @return The kernel covariance matrix of dimension N x N, between N locations
#'
#' @export

kernel_mat <- function(x, l, type = "exponential"){
  if(type == "exponential"){
    # exponential
    kernel_matrix = exp(-as.matrix(x)/l) 
  }else{
    # RBF
    kernel_matrix = exp(-as.matrix(x)^2/2/l^2)
  }
  return(kernel_matrix)
}

#' @title Function for comparing HS vs, ICAR local coefficients and p-values
#'
#' @param L_results is the output of local_res() function on the slope coefficients from the HS model
#' @param L_results_ICAR is the output of local_res() function on the slope coefficients from the ICAR model
#' @return The kernel covariance matrix of dimension N x N, between N locations
#'
#' @export

plot_local_comparisons <- function(L_results, L_results_ICAR, point_alpha = 0.6, point_size = 1.6) {

    df <- data.frame(
      beta      = as.numeric(L_results$local.beta),
      beta_icar = as.numeric(L_results_ICAR$local.beta),
      mlp       = -log10(pmax(as.numeric(L_results$local.pval), 1e-300)),
      mlp_icar  = -log10(pmax(as.numeric(L_results_ICAR$local.pval), 1e-300))
    )
  
  r1 <- cor(df$beta, df$beta_icar, use = "complete.obs")
  r2 <- cor(df$mlp,  df$mlp_icar,  use = "complete.obs")
  
  library(ggplot2)
  
  p1 <- ggplot(df, aes(x = beta, y = beta_icar)) +
    geom_point(alpha = point_alpha, size = point_size, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal() +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("r = %.2f", r1), hjust = -0.1, vjust = 1.1) +
    labs(title = "Local β: HS vs ICAR",
         x = "Local β (HS)",
         y = "Local β (ICAR)") +
    theme_minimal(base_size = 12)
  
  p2 <- ggplot(df, aes(x = mlp, y = mlp_icar)) +
    geom_point(alpha = point_alpha, size = point_size, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    coord_equal() +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("r = %.2f", r2), hjust = -0.1, vjust = 1.1) +
    labs(title = expression(paste(-log[10], " p: HS vs ICAR")),
         x = expression(-log[10](p)~"(HS)"),
         y = expression(-log[10](p)~"(ICAR)")) +
    theme_minimal(base_size = 12)
  
  
    return(patchwork::wrap_plots(p1, p2, ncol = 2))
}


