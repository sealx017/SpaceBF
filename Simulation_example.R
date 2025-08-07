library(anndata)
library(ggplot2)
library(MERINGUE)

mela <- read_h5ad("/Users/sealso/Research/Development/Data/SpaDM_melanoma.h5ad")
mela_ligand<-read.table("/Users/sealso/Research/Development/Data/SpaDM_melanoma_ligand.csv", header = T, sep = ",")
mela_receptor<-read.table("/Users/sealso/Research/Development/Data/SpaDM_melanoma_receptor.csv",  header = T, sep = ",")
LR_pairs <- paste0(mela_ligand$Ligand0, "_", mela_receptor$Receptor0)
LR_Overlap<-intersect(mela_receptor$X, LR_pairs)
LR_Overlap<-LR_pairs

LR_names <- reshape2::colsplit(LR_Overlap, "_", c("ligand", "receptor"))
raw_data_mela = t(as.matrix(as.matrix(mela$raw)))
coords = mela$obsm$spatial
colnames(coords) = c("x", "y")
rownames(coords) = 1:293

# RCTD cell/spot types
celltypes = as.matrix(colnames(mela$obs)[apply(mela$obs, 1, which.max)])
rownames(celltypes) = 1:nrow(coords)


saveRDS(cbind.data.frame(coords, celltypes), 
        "/Users/sealso/Documents/GitHub/SpaceBF/Data/Melanoma_coords.rds")
mela_coords <- cbind.data.frame(coords, celltypes)
save(mela_coords, file = "/Users/sealso/Documents/GitHub/SpaceBF/Data/mela_coords.rda")


# raw count expression data
raw_data_mela = t(as.matrix(mela$raw))

# coordinates of the cells/spots (293)
coords = mela$obsm$spatial
colnames(coords) = c("x", "y")
rownames(coords) = 1:nrow(coords)



x <- dist(coords)
x_d <- as.matrix(x)
x_d[col(x_d) < row(x_d)] <- 0
x_d[col(x_d) > row(x_d)] <- x_d[col(x_d) > row(x_d)] + 
  runif(length(x_d[col(x_d) > row(x_d)]), 0, 0.1)
x_d <- x_d + t(x_d)

mstobj <- fossil::dino.mst(x_d)

w <- MERINGUE::getSpatialNeighbors(coords, filterDist = 3)
plotNetwork(coords, w)


k1 <- spdep::knn2nb(spdep::knearneigh(coords))
#> Warning: neighbour object has 13 sub-graphs
all.linked <- max(unlist(spdep::nbdists(k1, coords)))
col.nb.0.all <- spdep::dnearneigh(coords, 0, all.linked)
spdep::plot.nb(col.nb.0.all, coords, add=TRUE)
nb.w <- spdep::nb2listw(col.nb.0.all, style="W", zero.policy=T)

#nb.w <- spdep::nb2listw(col.tri.nb, style="B", zero.policy=F)

### find a minimum spanning tree using dino
mstobj <- fossil::dino.mst(x)
G = mstobj
N = 293
#-----------------------------------------

set.seed(2025*2)
K <- kernel_mat(x, l = 7.2, type = "exponential") 

# generate bivariate Gaussian data across space: (y1, y2)^T, y1, y2 are N x 1
nu = 0.5 # correlation term 
y_latent_both <- c(MASS::mvrnorm(1, rep(0, 2*N), 
                                 kronecker(matrix(c(1, nu, nu, 1), nrow = 2), K))) # kronecker product sturcture
y1_latent <- y_latent_both[1:N]; y2_latent <- y_latent_both[-c(1:N)]

U1 = pnorm(y1_latent); U2 = pnorm(y2_latent) # transformation using the CDF of standard normal
NB_q1 <- NB_q2 <- 0.2                        # NB probabilities 
r1 <- r2 <- 1                                # NB dispersion parameters

y1 <- qnbinom(U1, r1, NB_q1)                 # inverse of NB CDF for both variables
y2 <- qnbinom(U2, r2, NB_q2)

#NBM <- NB_model(y1, y2, X = NULL, G = mst.adjmat, nIter = 5000, verbose = "TRUE")
#GM <- Gauss_model(y1, y2, X = NULL, G = mst.adjmat, nIter = 5000, verbose = "TRUE")

MM <- SpaceBF(y1, y2, X = NULL, G = mst.adjmat, which.model = "NB")
G_results <- global_res(MM$Beta[, -c(1:N)], local.p.ths = 0.9)
L_results <- local_res(MM$Beta[, -c(1:N)], local.p.ths = 0.9)

########################

library(anndata)
library(ggplot2)
library(MERINGUE)




# Gaussian --------------------------------
N = 5000
sigma2 = 10
ck = ck2 = NULL
X_1 = scale(runif(N, 0, 1000))
X_2 = scale(runif(N, 0, 1000))

x <- dist(cbind(X_1, X_2))
### find a minimum spanning tree using dino
mstobj <- fossil::dino.mst(x)
kernel_mat <- function(x, l){
  return(exp(-as.matrix(x)/l))
}

#-----------------------------------------


set.seed(2025*2)
K <- kernel_mat(x, l = 7.2) 

# generate bivariate Gaussian data across space: (y1, y2)^T, y1, y2 are N x 1
nu = 0.5 # correlation term 
y_latent_both <- c(MASS::mvrnorm(1, rep(0, 2*N), 
                   kronecker(matrix(c(1, nu, nu, 1), nrow = 2), K))) # kronecker product sturcture
y1_latent <- y_latent_both[1:N]; y2_latent <- y_latent_both[-c(1:N)]

U1 = pnorm(y1_latent); U2 = pnorm(y2_latent) # transformation using the CDF of standard normal
NB_q1 <- NB_q2 <- 0.2                        # NB probabilities 
r1 <- r2 <- 1                                # NB dispersion parameters

y1 <- qnbinom(U1, r1, NB_q1)                 # inverse of NB CDF for both variables
y2 <- qnbinom(U2, r2, NB_q2)

#NBM <- NB_model(y1, y2, X = NULL, G = mst.adjmat, nIter = 5000, verbose = "TRUE")
#GM <- Gauss_model(y1, y2, X = NULL, G = mst.adjmat, nIter = 5000, verbose = "TRUE")
start_time <- proc.time()
MM <- SpaceBF(y1, y2, X = NULL, G = mstobj, which.model = "NB", nIter = 1000)
print(proc.time() - start_time)

G_results <- global_res(MM$Beta[, -c(1:N)], local.p.ths = 0.9)
L_results <- local_res(MM$Beta[, -c(1:N)], local.p.ths = 0.9)

start_time <- proc.time()
ck <- NB_model(y1, y2, X = NULL, G = mstobj, nIter = 1000, nug_sig1 = 0.2,nug_sig2 = 0.2)
print(proc.time() - start_time)
