# =================================
# SpaceBF installation and loading
# =================================
devtools::install_github('sealx017/SpaceBF', quiet = TRUE)
suppressWarnings(library(SpaceBF))

# =========================================================================
# download the proteomics dataset from: https://zenodo.org/records/15866928
# and store inside a "path" directory
# =========================================================================
path <- "/Users/sealso/Research/Collaboration/Peggi/Data/"
peptide_exp <- readRDS(paste0(path, "Peptide_expression_in_spots.rds"))
coords_data <- readRDS(paste0(path, "Spot_coordinates.rds"))
peptide_info <- readRDS(paste0(path,"Peptide_description.rds"))
peptides <- peptide_info[, 1]

# =========================
# construct the random MST 
# =========================
N <- nrow(coords_data)           # number of locations
dist_mat <- dist(coords_data)  

mst.out <- MST_construction(dist_mat, type = "random") # takes a lot of time when N is large, igraph package has one implementaion of MST as well
mst.adjmat <- mst.out[[1]]                             # MST adjacency matrix, to be used as G in the main function later
print_mst(mst.out, coords_data)

# ============================================
# function to run SpaceBF across peptide pairs
# ============================================
peptide_pair <- matrix(0, (length(peptides)*(length(peptides) - 1)/2), 2)
k <- 1
for(i in 1:length(peptides)){
  if(i < length(peptides)){
    for(j in (i+1):length(peptides)){
      peptide_pair[k, ] <- c(i, j)
      k <- k+1
    }}
}    


Each_pair <- function(x){
  y1 <- peptide_exp[, peptide_pair[x, 1]]
  y2 <- peptide_exp[, peptide_pair[x, 2]]
  beta0s_beta1s <- SpaceBF::SpaceBF(y1, y2, X = NULL, G = mst.adjmat, 
                                    which.model = "Gauss", nIter = 2)  # decrease nIter for quicker results
  return(beta0s_beta1s)
}


# ========================================
# parallelized run across cores
# ========================================
numCores <- parallel::detectCores()
numCores

cl <- parallel::makeCluster(numCores - 1)
parallel::clusterEvalQ(cl, 2 + 2)
parallel::clusterExport(cl, varlist = c("peptide_exp", "mst.adjmat",
                                        "peptide_pair", "Each_pair"))
system.time(temp2 <- parallel::parLapply(cl, (1:nrow(peptide_pair)), Each_pair))
parallel::stopCluster(cl)

# unparallelized run (not recommended)
# sapply(1:nrow(peptide_pair), Each_pair)

# ===================================
# MCMC convergence for a random run
# ===================================
set.seed(2025*2025)

MM <- temp2[[sample((1:nrow(peptide_pair)), 1)]]       # choosing one LR pair randomly
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
for(i in (1:nrow(peptide_pair))){
  MM <- temp2[[i]]
  all_global_res <- rbind(all_global_res,
                          global_res(MM$Beta[, -c(1:N)], summary = "median", 
                          local.p.ths = 0.9))
}
all_global_res <- cbind.data.frame(all_global_res, peptide_pair)
all_global_res <- all_global_res[order(all_global_res$global.pval), ]

print(all_global_res)  

#L_results <- local_res(MM$Beta[, -c(1:N)],  local.p.ths = 0.9)
#print(L_results[1:5, ]) 

# ===================================
# plot local estimates
# ===================================
lowest_pval_peptide_pair <- all_global_res[1, c(5, 6)]              # plotting the local estimates for the LR pair with lowest p-value
lowest_pval_index <- which(peptide_pairs[, 1] == as.character(lowest_pval_peptide_pair[1]) & 
                     peptide_pairs[, 2] == as.character(lowest_pval_peptide_pair[2]))

MM <- temp2[[lowest_pval_index]]
plot_estimates(MM$Beta[, -c(1:N)], y1 = peptide_exp[LR_pairs[lowest_pval_index, 2], ], 
               y2 =  gene_data_mela[LR_pairs[lowest_pval_index, 1], ], coords_data)



