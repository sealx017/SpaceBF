# =================================
# SpaceBF installation and loading
# =================================

devtools::install_github('sealx017/SpaceBF', quiet = TRUE)
suppressWarnings(library(SpaceBF))

# ==================================================
# loading the melanoma coordinates data frame
# ==================================================

data("mela_coords") # three columns corresponding to xy coordinates and celltypes
coords <- mela_coords

# compute Euclidean distance matrix and construct the MST
dist_mat <- dist(coords[, 1:2]) # consider xy coordinates only
N <- nrow(coords)               # number of locations
mst.out <- MST_construction(dist_mat, type = "random") 
mst.adjmat <- mst.out[[1]]      # MST adjacency matrix, to be used as G in the main function later
print_mst(mst.out, coords[, 1:2])

# Delauney network for running MERINGUE
w <- MERINGUE::getSpatialNeighbors(coords[, 1:2], filterDist = 5) # constructing Delauney network using  MERINGUE
MERINGUE::plotNetwork(coords[, 1:2], w)

# k-NN network for Lee's L
k1 <- spdep::knn2nb(spdep::knearneigh(coords[, 1:2], k = 1))      # can vary the # neighbors, k
all.linked <- max(unlist(spdep::nbdists(k1, coords[, 1:2])))      # can be relaxed
col.nb.0.all <- spdep::dnearneigh(coords[, 1:2], 0, all.linked)
spdep::plot.nb(col.nb.0.all, coords[, 1:2])
nb.w <- spdep::nb2listw(col.nb.0.all, style="W", zero.policy=T)

# ==============================================
# exponential kernel covariance matrix function
# ==============================================

kernel_mat <- function(x, r){         # x is the distance matrix or dist object, r is the lengthscale
  return(exp(-as.matrix(x)/r))
}

# ==============================================
# specify and create storage directory
# ==============================================

path <- "/Users/sealso/Research/Development/Results/RDS/Sim_1/"   # specify the folder address on your system
dir.create(file.path(path), showWarnings = FALSE)

# ==============================================
# main simulation function
# ==============================================

sim_runner<-function(r){
  k <- 1
  num <- r*10
  for(rho in c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)){
    fin_res <- NULL 
    for(iter in 1:100){
      set.seed(2024*num*iter*k)
      
      # construct the kernel covariance matrix
      K <- kernel_mat(dist_mat, r)  # r denotes the lengthscale
      
      # copula based NB-distributed predictor variable simulation
      r1 <- 1; r2 <- 1              # NB over-dispersion parameters of the two variables
      NB_q2 <- 0.5                
      U2 <- pnorm(mvnfast::rmvn(1, rep(0, N), K)) 
      y2 <- c(qnbinom(U2, r2, NB_q2))
      
      while (var(y2) == 0) {                        # Loop until the variance of y2 is non-zero
        U2 <- pnorm(mvnfast::rmvn(1, rep(0, N), K)) 
        y2 <- c(qnbinom(U2, r2, NB_q2))             # Generate new values for y2
      }
      
      # direct NB-distributed outcome variable simulation
      mu <- c(mvnfast::rmvn(1, rep(0, N), 0.5*K))   # spatially correlated intercept
      eta1 <- rho*scale(log(y2+1))[,1] + mu         # fixed slope rho
      NB_q1 <- 1/(1+exp(eta1)) 
      y1 <- rnbinom(N, prob = NB_q1, size = r1)     
      
      # SpaceBF
      MM <- SpaceBF::SpaceBF(y1, y2, X = NULL, G = mst.adjmat, which.model = "NB", nIter = 5000) 
      
      # MERINGUE
      names(y1) <- names(y2) <- rownames(w)
      MER_stat <- MERINGUE::spatialCrossCor(log(y1+1), log(y2+1), w)
      MER_pval <- MERINGUE::spatialCrossCorTest(log(y1+1), log(y2+1), w, n = 2000, plot=F)
      
      # Lee's L
      res_xy <- spdep::lee.mc(log(y1 + 1), log(y2 + 1), nb.w, nsim = 2000, alternative = "two.sided")
      
      # SpaGene
      st_meta <- cbind.data.frame(as.character(1:N), coords[, 1:2])
      colnames(st_meta) <- c("cell", "x", "y")
      st_data <- rbind.data.frame(log(y1+1), log(y2+1))
      rownames(st_data) <- c("y1", "y2")
      colnames(st_data) <- as.character(1:N)
      st_data <- as.matrix(st_data)
      
      mob_lr <- SpaGene::SpaGene_LR(st_data, st_meta[,-1],  
              LRpair = matrix(c("y1","y2"), 1, 2), knn = 4, normalize = F) 

      fin_res[[iter]]  = list(SpaceBF = MM, Lees = c(res_xy$statistic, res_xy$p.value),
                              MER = c(MER_stat, MER_pval),
                              SpaGene = mob_lr[1,3:4], 
                              PearsonCorr = cor.test(y1, y2)$p.value,
                              y1, y2)
    }
    
    saveRDS(fin_res, file = paste0(path,"Simulation_design_1_with_num_", 
                                   num, "_eff_size_rho_", rho, ".rds")) # num = r*10, r is the lengthscale
    
  }
  k = k+1  
}

# ==============================================
# run the simulation across lengthscales
# ==============================================

numCores <- parallel::detectCores()
numCores

cl <- parallel::makeCluster(numCores-1)
parallel::clusterEvalQ(cl, 2 + 2)
parallel::clusterExport(cl, varlist = c("sim_runner", "kernel_mat", "coords", "N", 
                                        "w", "nb.w", "dist_mat", "mst.adjmat", "path"))

# parallelize across lengthscale values
system.time(MERINGUE_res <- parallel::parLapply(cl,  as.list(c(3.6, 7.2, 18)), sim_runner))
parallel::stopCluster(cl)

# ==============================================
# extract results
# ==============================================

r1 <- r2 <- 1
NB_q1 <- NB_q2 <- 0.5
r <- c(3.6, 7.2, 18)
all_power <- NULL
num_all <- r*10

for(num in num_all){
  for(rho in c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)){
    power_beta1 <- gew_all <- NULL
    Beta1_res <- Beta1_meds <- NULL

    fin_res <- readRDS(file = paste0(path,"Simulation_design_1_with_num_", 
                      num, "_eff_size_rho_", rho, ".rds"))

    for(iter in 1:100){
        MM <- fin_res[[iter]][[1]]
        Beta1_meds <- rbind.data.frame(Beta1_meds, c(global_res(MM$Beta[, -c(1:N)], summary = "median", 
                      local.p.ths = 0.9)))
    }
    Beta1_res <- rbind.data.frame(Beta1_res, c(mean(Beta1_meds[, 1]), 
                mean(ifelse(Beta1_meds[, 2] < 0.05, 1, 0))))
    Beta1_res <- as.data.frame(Beta1_res)

    BV <- MER <- Spagene <- cortest <- NULL
    for(iter in 1:100){
      BV <- c(BV, ifelse(fin_res[[iter]][[2]][2]<0.05, 1, 0))
      MER <- c(MER, ifelse(fin_res[[iter]][[3]][2]<0.05, 1, 0))
      Spagene <- c(Spagene, ifelse(fin_res[[iter]][[4]][2]<0.05, 1, 0))
      cortest <- c(cortest, ifelse(fin_res[[iter]][[5]]<0.05, 1, 0)) 
    }
    
    power_beta1 <- rbind.data.frame(power_beta1, c(Beta1_res[, 2],
                   mean(BV), mean(MER), mean(Spagene), mean(cortest)))
    power_beta1 <- data.frame(power_beta1)
    colnames(power_beta1) <- c("SpaceBF","Lee's_L", "MERINGUE", "SpaGene", "PearsonCorr")
    power_beta1$rho <- rho 
    power_beta1$lengthscale <- num/10 
    
    all_power <- rbind.data.frame(all_power, power_beta1)
  }
}

print(all_power) # SpatialDM, not shown here, was run on Python using the original package


