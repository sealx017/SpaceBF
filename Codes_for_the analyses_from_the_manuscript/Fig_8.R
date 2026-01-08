# =================================
# SpaceBF installation and loading
# =================================
devtools::install_github('sealx017/SpaceBF', force = T)
library(SpaceBF)

# ==================================================
# loading the melanoma coordinates data frame
# ==================================================

data(mela_coords) # three columns corresponding to xy coordinates and celltypes
coords <- mela_coords

# compute Euclidean distance matrix and construct the MST
dist_mat <- dist(coords[, 1:2]) # consider xy coordinates only
N <- nrow(coords)               # number of locations
mst.out <- MST_construction(dist_mat, type = "random") 
mst.adjmat <- mst.out[[1]]      # MST adjacency matrix, to be used as G in the main function later
print_mst(mst.out, coords[, 1:2])

# construct a mutual 3-NN network
kNN_obj <- SpaceBF::kNN_construction(coords[, 1:2], k = 3)
SpaceBF::print_mst(kNN_obj, coords[, 1:2])
w_knn <- as.matrix(kNN_obj[[1]]) # kNN adjacency matrix

# construct Delaunay network 
w <- MERINGUE::getSpatialNeighbors(coords[, 1:2], filterDist = 5) 

# ==============================================
# exponential kernel covariance matrix function
# ==============================================

kernel_mat <- function(x, r){         # x is the distance matrix or dist object, r is the lengthscale
  return(exp(-as.matrix(x)/r))
}

# ==============================================
# specify and create storage directory
# ==============================================

path <- "/Users/sealso/Research/Development/Results/RDS/Sim_3/"   # specify the folder address on your system
dir.create(file.path(path), showWarnings = FALSE)

# ==============================================
# main simulation function
# ==============================================

sim_runner<-function(iter){
  k <- 1
  for(rho in c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)){
    fin_res <- NULL
    r_count <- 1
    for(r in c(1, 2)){
      num <- r*10
      set.seed(2025*num*iter)
      
      # construct the kernel covariance matrix for y2
      K <- kernel_mat(dist_mat, 7.2)  
      
      NB_q2 <- 0.5; r1 <- r2 <- 1;
      U2 <- pnorm(mvnfast::rmvn(1, rep(0, N), K)) 
      y2 <- c(qnbinom(U2, r2, NB_q2))
      # Loop until the variance of y2 is non-zero
      while (var(y2) == 0 | mean(y2 == 0) > 0.8) {
        U2 <-pnorm(mvnfast::rmvn(1, rep(0, N), K)) 
        y2 <- c(qnbinom(U2, r2, NB_q2))  # Generate new values for y2
      }
      
      mu <- c(mvnfast::rmvn(1, rep(0, N), 0.5*K)) # spatially correlated beta0
      rho_vec <- rep(0, N) # spatially varying beta1
      
      # circular partitioning
      rx <- ry <- r
      sd_in <- 0.3
      sd_out <- 0.6
      
      cx <- median(coords[, 1]);  cy <- median(coords[, 2])
      inside <- ((coords[,1] - cx)^2 / rx^2 + (coords[,2] - cy)^2 / ry^2) <= 1
      # piecewise-constant coefficient beta(s) with a jump across boundary
      beta_piece <- ifelse(inside, rho, -rho)
      # heteroskedasticity by region 
      beta_noise <- rnorm(N, 0, ifelse(inside, sd_in, sd_out))
      
      # final spatially varying coefficient beta1(s)
      rho_vec <- beta_piece + beta_noise
      
      eta1 <- rho_vec*scale(log(y2+1))[,1] + mu 
      
      NB_q1 <- 1/(1+exp(eta1)) # dnbinom fn uses q=1-psi
      y1 <- rnbinom(N, prob = NB_q1, size = r1)
      y1[y1 > 1500] <- 1500 # truncate extreme values
      
      require(sdmTMB)
      snow <- data.frame(y1 = y1, y2 = scale(log1p(y2)), coords[, 1:2])
      mesh_snow1 <-sdmTMB::make_mesh(snow, xy_cols = c("x", "y"), cutoff = 0.5)
      #plot(mesh_snow1$mesh)                      # triangles + nodes
      #points(snow$x, snow$y, pch = 16, cex = .5, col = "red")  # data locations
      
      fit_owl1 <- sdmTMB(
        y1 ~ 0,
        spatial_varying = ~ 0 + y2,
        data = snow,
        mesh = mesh_snow1,
        family = nbinom2(link = "log"),
        spatial = "on",
        silent = FALSE
      )
      zeta_s1 <- predict(fit_owl1, newdata = snow)
      snow$estimate1 <-  zeta_s1$zeta_s_y2
      
      mesh_snow2 <-sdmTMB::make_mesh(snow, xy_cols = c("x", "y"), cutoff = 2)
      #plot(mesh_snow2$mesh)                      # triangles + nodes
      #points(snow$x, snow$y, pch = 16, cex = .5, col = "red")  # data locations
      
      fit_owl2 <- sdmTMB(
        y1 ~ 0, #+ y2,
        spatial_varying = ~ 0 + y2,
        data = snow,
        mesh = mesh_snow2,
        family = nbinom2(link = "log"),
        spatial = "on",
        silent = FALSE
      )
      zeta_s2 <- predict(fit_owl2, newdata = snow)
      snow$estimate2 <- zeta_s2$zeta_s_y2
      
      # SpaceBF with different priors
      nIter <- 5000
      
      # MST horseshoe
      MM <- SpaceBF::SpaceBF(y1, y2, X = NULL, nug_sig = c(0.01, 0.01),  beta_thres = 1e+2,
                             G = mst.adjmat, which.model = "NB", nIter = nIter, verbose = F) 
      # Delaunay ICAR
      MM_ICAR2 <- SpaceBF::SpaceBF(y1, y2, X = NULL, nug_sig = c(0.01, 0.01), beta_thres = 1e+2,
                                   G = w, which.model = "NB", which.prior = "ICAR", nIter = nIter, verbose = F) 
      # kNN ICAR
      MM_ICAR3 <- SpaceBF::SpaceBF(y1, y2, X = NULL, nug_sig = c(0.01, 0.01), beta_thres = 1e+2,
                                   G = w_knn, which.model = "NB", which.prior = "ICAR", nIter = nIter, verbose = F)
      # Delaunay horseshoe
      MM_ICAR4 <- SpaceBF::SpaceBF(y1, y2, X = NULL, nug_sig = c(0.01, 0.01), beta_thres = 1e+2,
                                   G = w, which.model = "NB", which.prior = "HS", nIter = nIter, verbose = F)
      # kNN horseshoe
      MM_ICAR5 <- SpaceBF::SpaceBF(y1, y2, X = NULL, nug_sig = c(0.01, 0.01), beta_thres = 1e+2,
                                   G = w_knn, which.model = "NB", which.prior = "HS",
                                   nIter = nIter, verbose = F)
      
      fin_res[[r_count]]  <- list(SpaceBF = MM, SpaceBF_ICAR1 =  MM_ICAR2, 
                                  SpaceBF_ICAR2 =  MM_ICAR3, SpaceBF_HS1 =  MM_ICAR4, 
                                  SpaceBF_HS2 =  MM_ICAR5, sdmTMB = snow[, c("estimate1", "estimate2")], 
                                  y1, y2)
      r_count <- r_count + 1
    }
    saveRDS(fin_res, file = paste0(path,"/Prior_comparison_MSE_Simulation_design_3_circular_iter_", iter, "_rho_", rho, ".rds")) 
  }
}


numCores <- parallel::detectCores()
numCores

cl <- parallel::makeCluster(numCores-1)
parallel::clusterEvalQ(cl, 2 + 2)
parallel::clusterExport(cl, varlist = c("sim_runner", "kernel_mat", "coords", "N", "mst.adjmat",
                                        "w", "w_knn", "dist_mat", "path"))

library(pbapply)

# (optional) nicer progress bar in console
old_pb <- pboptions(type = "timer")  # or "txt", "win", etc.

# parallelize across lengthscale values with progress
system.time({
  res <- pbapply::pblapply(
    X  = as.list(1:100), 
    FUN = sim_runner,
    cl  = cl                      # use your existing PSOCK/cluster
  )
})

parallel::stopCluster(cl)




