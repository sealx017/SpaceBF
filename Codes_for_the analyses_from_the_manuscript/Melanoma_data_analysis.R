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

ck = colMeans(beta0s_beta1s$Beta[, -c(1:N)])

# ========================================
# parallelized run across cores
# ========================================
numCores <- parallel::detectCores()
numCores

cl <- parallel::makeCluster(numCores - 1)
parallel::clusterEvalQ(cl, 2 + 2)
parallel::clusterExport(cl, varlist = c("gene_data_mela", "mst.adjmat",
                                        "LR_pairs", "Each_LR"))
system.time(temp2 <- parallel::parLapply(cl, (1:nrow(LR_pairs)), Each_LR))
parallel::stopCluster(cl)

# unparallelized run (not recommended)
# sapply(1:nrow(LR_pairs), Each_LR)

# ===================================
# MCMC convergence for a random run
# ===================================
set.seed(2025*2025)

MM <- temp2[[sample((1:nrow(LR_pairs)), 1)]]       # choosing one LR pair randomly
ran_sam <- sample(1:N, 4)
par(mfrow = c(2, 2))
MCMC_plot(MM$Beta[, (N + ran_sam[1])])
MCMC_plot(MM$Beta[, (N + ran_sam[2])])
MCMC_plot(MM$Beta[, (N + ran_sam[3])])
MCMC_plot(MM$Beta[, (N + ran_sam[4])])

# ===================================
# global and local results
# ===================================
all_global_res <- NULL
for(i in (1:nrow(LR_pairs))){
  MM <- temp2[[i]]
  all_global_res <- rbind(all_global_res,
                    global_res(MM$Beta[, -c(1:N)], summary = "median", 
                    local.p.ths = 0.9))
}
all_global_res <- cbind.data.frame(all_global_res, LR_pairs)
all_global_res <- all_global_res[order(all_global_res$global.pval), ]

print(all_global_res)  

#L_results <- local_res(MM$Beta[, -c(1:N)],  local.p.ths = 0.9)
#print(L_results[1:5, ]) 

# ===================================
# plot local estimates
# ===================================
lowest_pval_LR_pair <- all_global_res[1, c(5, 6)]              # plotting the local estimates for the LR pair with lowest p-value
lowest_pval_index <- which(LR_pairs[, 1] == as.character(lowest_pval_LR_pair[1]) & 
                           LR_pairs[, 2] == as.character(lowest_pval_LR_pair[2]))

MM <- temp2[[lowest_pval_index]]
plot_estimates(MM$Beta[, -c(1:N)], y1 = gene_data_mela[LR_pairs[lowest_pval_index, 2], ], 
               y2 =  gene_data_mela[LR_pairs[lowest_pval_index, 1], ], coords_data_mela)



