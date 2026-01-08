# =================================
# SpaceBF installation and loading
# =================================

devtools::install_github('sealx017/SpaceBF', quiet = TRUE)
suppressWarnings(library(SpaceBF))

# ==================================================
# loading the three data frames with gene expression, 
# xy-coordinates, and LR pairs, respectively
# ==================================================

data("SpaceBF_melanoma")
coords_data_mela <- SpaceBF_melanoma$coords_data_mela
gene_data_mela <- SpaceBF_melanoma$gene_data_mela
LR_pairs <- SpaceBF_melanoma$LR_pairs
# gene_data_mela is p x N matrix, where p denotes the number of genes and N denotes the number of locations
# coords_data_mela is N x 2 matrix, 2 columns denoting xy-coordinates
# LR_pairs is 161 x 2 matrix, with 161 LR pairs as described in the manuscript

# =========================
# construct the random MST 
# =========================

N <- nrow(coords_data_mela)           # number of locations
dist_mat <- dist(coords_data_mela)  

mst.out <- MST_construction(dist_mat, type = "random") 
mst.adjmat <- mst.out[[1]]            # MST adjacency matrix, to be used as G in the main function later
print_mst(mst.out, coords_data_mela[, 1:2])

# ========================================
# function to run SpaceBF across LR pairs
# ========================================

Each_LR <- function(x){
  y1 <- gene_data_mela[LR_pairs[x, 2], ]
  y2 <- gene_data_mela[LR_pairs[x, 1], ]
  beta0s_beta1s <- SpaceBF::SpaceBF(y1, y2, X = NULL, G = mst.adjmat, 
                                    which.model = "NB", nIter = 5000)  # decrease nIter for quicker results
  return(beta0s_beta1s)
}

# ========================================
# parallelized run across cores
# ========================================

set.seed(2025*2025)
numCores <- parallel::detectCores()
numCores

cl <- parallel::makeCluster(numCores - 1)
parallel::clusterEvalQ(cl, 2 + 2)
parallel::clusterExport(cl, varlist = c("gene_data_mela", "mst.adjmat",
                                        "LR_pairs", "Each_LR"))
system.time(temp2 <- parallel::parLapply(cl, (1:nrow(LR_pairs)), Each_LR))
parallel::stopCluster(cl)

# ========================================================================
# MCMC convergence for two LR pairs, four randomly selected coefficients
# ========================================================================

LR_pairs_name <- paste0(LR_pairs[, 1], "_", LR_pairs[, 2])
MM <- temp2[[which(LR_pairs_name == "IGF2_IGF1R")]]      
ran_sam <- sample(1:N, 4)
par(mfrow = c(2, 2))
MCMC_plot(MM$Beta[, (N + ran_sam[1])])
MCMC_plot(MM$Beta[, (N + ran_sam[2])])
MCMC_plot(MM$Beta[, (N + ran_sam[3])])
MCMC_plot(MM$Beta[, (N + ran_sam[4])])

MM <- temp2[[which(LR_pairs_name == "SPP1_CD44")]]      
ran_sam <- sample(1:N, 4)
par(mfrow = c(2, 2))
MCMC_plot(MM$Beta[, (N + ran_sam[1])])
MCMC_plot(MM$Beta[, (N + ran_sam[2])])
MCMC_plot(MM$Beta[, (N + ran_sam[3])])
MCMC_plot(MM$Beta[, (N + ran_sam[4])])

