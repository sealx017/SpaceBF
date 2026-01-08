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
# construct the kNN graph
# =========================

N <- nrow(coords_data_mela)           # number of locations
kNN_obj <- SpaceBF::kNN_construction(coords_data_mela, k = 4)
w_knn <- as.matrix(kNN_obj[[1]]) 


# ========================================
# function to run SpaceBF across LR pairs
# ========================================

Each_LR <- function(x){
  y1 <- gene_data_mela[LR_pairs[x, 2], ]
  y2 <- gene_data_mela[LR_pairs[x, 1], ]
  beta0s_beta1s <- SpaceBF::SpaceBF(y1, y2, X = NULL, G = w_knn, 
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
parallel::clusterExport(cl, varlist = c("gene_data_mela", "w_knn",
                                        "LR_pairs", "Each_LR"))
system.time(temp2 <- parallel::parLapply(cl, (1:nrow(LR_pairs)), Each_LR))
parallel::stopCluster(cl)

# ===================================
# MCMC convergence for a random run
# ===================================

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

# ===================================
# plot local estimates
# ===================================

require(RColorBrewer)
require(latex2exp)

expression_min <- 0
expression_max <- 3
LH_max <- 1.6
LH_min <- -1.6 
LH_plots_n <- list()

for(l in 1:nrow(LR_pairs)){
  
  location_data_with_beta <- as.data.frame(cbind(coords_data_mela, 
                             log(gene_data_mela[LR_pairs[l, 2], ] + 1)))
  location_data_with_beta$y2 <- log(gene_data_mela[LR_pairs[l, 1], ] + 1)
  colnames(location_data_with_beta)[1:3] <- c("x", "y", "y1")
  
  ck <- temp2[[l]]$Beta[, -c(1:(N))]
  sigfinder <- SpaceBF::local_credible(ck)
  
  location_data_with_beta$LH <- robustbase::colMedians(ck) #ck
  p.sf <- sf::st_as_sf(location_data_with_beta, coords = c("x", "y"))
  mean_beta <- round(mean(p.sf$LH), 3)
  SD_beta <- round(sd(p.sf$LH), 3)
  
  p.sf$LH <- scale(p.sf$LH, scale = T)[,1]
  p.sf$LH <- pmin(pmax(p.sf$LH, LH_min), LH_max)
  
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  sc_LH <- scale_colour_gradientn(colours = myPalette(10), 
           limits=c(LH_min, LH_max), breaks = seq(LH_min, LH_max, length.out = 5))
  sc_LH$breaks <- round(sc_LH$breaks, digits = 2)
  
  p.sf$y1 <- pmin(pmax(p.sf$y1, expression_min), expression_max)
  p.sf$y2 <- pmin(pmax(p.sf$y2, expression_min), expression_max)

  
  betas <- ggplot(data =  p.sf) +
    geom_sf(aes(color = LH), size = 1.5) +  
    labs(color = TeX(sprintf("$\\beta_1(s) - \\bar{\\beta}_1$"))) + 
    ggtitle(TeX(sprintf("$\\bar{\\beta}_1 = %.2f$, $\\sigma_{\\beta} = %.2f$", mean_beta, SD_beta))) +            #TeX(sprintf("$\\sigma_B = %f $", SD_beta))) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="none", legend.text=element_text(size=9), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
    ) + sc_LH
  
  
  p.sf$LH[sigfinder == 0] = -5
  
  betas2 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = LH), size = 1.5) +  
    labs(color = TeX(sprintf("$\\frac{\\beta_1(s) - \\bar{\\beta}_1}{\\sigma_{\\beta}}$"))) + 
    ggtitle(TeX(sprintf("$\\bar{\\beta}_1 = %f $", mean_beta))) + 
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="right", legend.text=element_text(size=9), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12)) + sc_LH
  
  
  p1 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = y1), size = 1.5) +   
    scale_color_gradientn(colours = (viridis::viridis(10, direction = -1)), 
    limits=c(expression_min, expression_max)) +
    labs(color = "log1p(expr)") + 
    ggtitle(paste0("Receptor: ", LR_pairs[l, 2]))+
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="left", legend.text=element_text(size=9), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12))
  
  p2 <- ggplot(data =  p.sf) +
    geom_sf(aes(color = y2), size = 1.5) +   
    scale_color_gradientn(colours = (viridis::viridis(100, direction = -1)), 
    limits=c(expression_min, expression_max)) +
    labs(color = "Expression") + ggtitle(paste0("Ligand: ", LR_pairs[l, 1]))+
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position="none", legend.text=element_text(size=8), 
          legend.title = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, size = 12))
  
  LH_plots_n[[l]] =  patchwork::wrap_plots(list(p1, p2, betas, betas2), nrow = 1)
  print(l)
}  

# =====================================================================
# LR pair names and extract the slope plot of three demonstrated pairs
# =====================================================================

LR_pairs_name <- paste0(LR_pairs[, 1], "_", LR_pairs[, 2])
LH_plots_n[[which(LR_pairs_name == "IGF2_IGF1R")]]
LH_plots_n[[which(LR_pairs_name == "PTPRC_CD22")]]
LH_plots_n[[which(LR_pairs_name == "SPP1_CD44")]]


# ====================================
# Finding clusters of scaled slopes
# ====================================

ths <- 0.05
k <- 1
local_beta_estimates_top <- NULL
for(l in which((all_global_res$global.pval < ths))){
  local_beta_estimates_top[[k]] <- robustbase::colMedians((temp2[[l]]$Beta)[, -c(1:N)])
  local_beta_estimates_top[[k]] <- scale(local_beta_estimates_top[[k]], center = T)[,1]  # centering to capture only the relative structure between spots
  k <- k + 1
}

local_beta_estimates_top <- do.call(rbind, local_beta_estimates_top)
spot_clustering <- hclust(dist(t(local_beta_estimates_top)), "ward.D2")
location_data_with_cluster <- as.data.frame(coords_data_mela)
location_data_with_cluster$cluster <- as.factor(cutree(spot_clustering, k = 3))    # 3 clusters
p.sf <- sf::st_as_sf(location_data_with_cluster,  coords = c("x", "y"))

ggplot(data =  p.sf) + geom_sf(aes(color = cluster), size = 1) +
  labs(color = "Cluster") + ggtitle("Clusters based on slope surfaces") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.position="bottom", legend.text=element_text(size=5), 
        legend.title = element_text(size=5),
        plot.title = element_text(hjust = 0.5, size = 5))  +
  guides(colour = guide_legend(override.aes = list(size=2)))


