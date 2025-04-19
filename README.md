
<table style="width:100%">
<tr>
<td style="text-align:left; vertical-align:middle;">
<h1>
SpaceBF
</h1>
<p>
<strong>Authors:</strong> Souvik Seal and Brian Neelon
</p>
</td>
<td style="text-align:right; padding-left:400px;">
<img src="SpaceBF_logo.png" width="140" height="150"/>
</td>
</tr>
</table>

This *R* package implements the models proposed in the paper “SpaceBF:
Spatial Coexpression Analysis Using Bayesian Fused Approaches in Spatial
Omics Datasets.” It enables the analysis of spatially varying
co-expression between pairs of molecules, such as ligand-receptor
interactions, in spatial transcriptomics (ST) and mass spectrometry
imaging (MSI) datasets. It is applicable to broader spatial datasets as
well, to study spatially varying interaction between any two variables.

## Loading the package

We install and load the developmental version of SpaceBF from GitHub.

``` r

#suppressWarnings(devtools::install_github('sealx017/SpaceBF', quiet = TRUE))
#library(SpaceBF)

# packages providing functional help
library(Matrix)
library(igraph)
library(MERINGUE)

# packages providing plotting help
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(ggraph)
library(latex2exp)
library(patchwork)
source("R/Additional_functions.R")
source("R/Gauss_MCMC.R")
source("R/NB_MCMC.R")
source("R/SpaceBF_main.R")
source("R/Result_generator.R")
```

## Load the coordinates of the Melanoma dataset and compute the MST

``` r
coords <- readRDS("/Users/sealso/Documents/GitHub/SpaceBF/Data/Melanoma_coords.rds")
data("mela_coords") # three columns corresponding to xy coordinates and celltypes
# Warning in data("mela_coords"): data set 'mela_coords' not found

# compute Euclidean distance matrix and construct MST
x <- dist(coords[, 1:2]) # consider xy coordinates only
N <- nrow(coords)        # number of locations
mst.out <- MST_construction(x, type = "random")
mst.adjmat <- mst.out[[1]]
print_mst(mst.out, coords[, 1:2])
```

<img src="README_files/figure-gfm/loading the coords-1.png" width="50%" />

``` r

# optional network to compute MERINGUE cross-correlation
#w <- getSpatialNeighbors(coords[, 1:2], filterDist = 5) # constructing Delauney netwrok 
                                                        # using  MERINGUE
#plotNetwork(coords[, 1:2], w)
```

## Construct the kernel covariance matrix

This function can be used to perform the entire analysis at one go,
starting from summary function estimation to the univariate and
multivariate FANOVA-based association analyses.

``` r

# spatial kernel covariance matrix with exponential kernel function
K <- kernel_mat(x, l = 7.2, type = "exponential") 

# visualizing the kernel function to get an idea about spatial structure
palettes <-  sort(hcl.pals("divergingx"))
col_fun <- colorRamp2(c(0, 0.5, 1), hcl_palette = palettes[15], reverse = T) # spectral palette
ComplexHeatmap::Heatmap(K, cluster_rows = F, cluster_columns = F, col = 
  col_fun, show_row_names = FALSE,  show_column_names = FALSE, 
  heatmap_legend_param = list(title = "value"))
```

<img src="README_files/figure-gfm/construct kernel matrix-1.png" width="50%" />

## Simulate count expression data using bivariate Gaussian copula

``` r
set.seed(2025*2)

# generate bivariate Gaussian data across space: (y1, y2)^T, y1, y2 are N x 1
nu = 0.5                                          # correlation term 
y_latent_both <- c(MASS::mvrnorm(1, rep(0, 2*N), 
kronecker(matrix(c(1, nu, nu, 1), nrow = 2), K))) #kronecker product sturcture

y1_latent <- y_latent_both[1:N]; y2_latent <- y_latent_both[-c(1:N)]
U1 = pnorm(y1_latent); U2 = pnorm(y2_latent) # transformation using the CDF of standard normal

NB_q1 <- NB_q2 <- 0.2              # NB probabilities 
r1 <- r2 <- 1                      # NB dispersion parameters

y1 <- qnbinom(U1, r1, NB_q1)       # inverse of NB CDF for both variables
y2 <- qnbinom(U2, r2, NB_q2)
```

The function takes several input parameters out of which the key ones
are displayed below. “fixed_r” corresponds to the grid values of radius
r. “Summary_function” corresponds to which summary function to be used:
g, K, or L. “Hard_ths” corresponds to a QC threshold in terms of cell
counts. If in an image, a particular cell type ‘A’ has less than
“Hard_ths” many cells, the image will be dropped from any pairwise
comparisons involving cell type ‘A’. “homogeneous” specifies whether
homogeneous poisson point process (PPP) (and corresponding g function)
or inhomogeneous PPP (and corresponding g function) is to be used.
“interaction_adjustment” specifies whether to consider the interaction
adjusted version of SpaceANOVA Mult. “perm” is used to employ the
permutation envelope-based adjustment to the summary functions in
presence of holes in images and “nPerm” denotes the number of
permutations to be used. “cores” is the number of cores to be used, if
more than 1, mclapply (parallel package) is called to construct the
permutation envelope.

``` r

MM <- SpaceBF(y1, y2, X = NULL, G = mst.adjmat, which.model = "NB", nIter = 50)
# [1] "SpaceBF: NB model fitting has started."
# Progress:  1 on 50  Progress:  2 on 50  Progress:  3 on 50  Progress:  4 on 50  Progress:  5 on 50  Progress:  6 on 50  Progress:  7 on 50  Progress:  8 on 50  Progress:  9 on 50  Progress: 10 on 50  Progress: 11 on 50  Progress: 12 on 50  Progress: 13 on 50  Progress: 14 on 50  Progress: 15 on 50  Progress: 16 on 50  Progress: 17 on 50  Progress: 18 on 50  Progress: 19 on 50  Progress: 20 on 50  Progress: 21 on 50  Progress: 22 on 50  Progress: 23 on 50  Progress: 24 on 50  Progress: 25 on 50  Progress: 26 on 50  Progress: 27 on 50  Progress: 28 on 50  Progress: 29 on 50  Progress: 30 on 50  Progress: 31 on 50  Progress: 32 on 50  Progress: 33 on 50  Progress: 34 on 50  Progress: 35 on 50  Progress: 36 on 50  Progress: 37 on 50  Progress: 38 on 50  Progress: 39 on 50  Progress: 40 on 50  Progress: 41 on 50  Progress: 42 on 50  Progress: 43 on 50  Progress: 44 on 50  Progress: 45 on 50  Progress: 46 on 50  Progress: 47 on 50  Progress: 48 on 50  Progress: 49 on 50  Progress: 50 on 50  [1] "SpaceBF: NB model fitting completed!"
G_results <- global_res(MM$Beta[, -c(1:N)], local.p.ths = 0.9)
L_results <- local_res(MM$Beta[, -c(1:N)], local.p.ths = 0.9)
```

## Extracting p-values

Extract the p-values of the univariate and multivariate tests:
SpaceANOVA Univ. and SpaceANOVA Mult. and print.

``` r
# p_res = p_extract(Final_result)
# Univ_p = p_res[[1]]
# Mult_p = p_res[[2]]
# print(Univ_p)
```

## Heatmap of -log10 of p-values

Display the -log10 of the p-values using heatmap.

``` r
# Plot.heatmap(Univ_p, main = "SpaceANOVA Univ.")
# Plot.heatmap(Mult_p, main = "SpaceANOVA Mult.")
```

## Visualizing summary functions

Next, we can check the summary functions corresponding to those pairs of
cell types whose p-values turned out to be significant. Here, we
visualize the g-functions of pairs: (alpha, delta) and (beta, Th). For
every pair, the first panel shows subject-level mean functions (averaged
across the images of a subject), while the second panel shows
image-level summary functions. We also display the point-wise F-values
which might help to understand which particular values of the radius r
are most influential, i.e., where the maximum difference between the
group-level summary functions are observed.

### Pair: (alpha, delta)

We notice that the spatial co-occurrence of alpha and delta cells is
positive in both the groups but higher in the “Onset” group.

``` r
# Pointwise_Fvals = Final_result[[1]][[2]]
# Functional_results = Final_result[[2]]
# 
# # Pair: (alpha, delta)
# Plot.functions(Functional_results, pair = c("alpha", "delta"), Fvals = Pointwise_Fvals)
```

### Pair: (beta, Th)

We notice that the spatial co-occurrence of (beta, Th) is negative in
both the groups meaning that the cell types avoid each other. But the
degree of avoidance decreases in the “Onset” group.

``` r
# # Pair: (beta, Tc)
# Plot.functions(Functional_results, pair = c("beta", "Tc"), Fvals = Pointwise_Fvals)
```

## Visualizing cellular organization

We can visualize cellular organization in different images of every
subject. Here, 4 images each are shown for two subjects, “6126” from
group “Non-diabetic” and “6414” from group “Onset”.

``` r
# table(IMC_T1DM$cellType)
# palette = c("darkorchid1","red", "cyan", "grey", "blue", "green") #assign colors to cell types 
# 
# Plot.cellTypes(data = IMC_T1DM, ID = "6126", palette = palette)
# Plot.cellTypes(data = IMC_T1DM, ID = "6414", palette = palette)
```

## References

1\) Seal, S., & Neelon, B. (2025). SpaceBF: Spatial coexpression
analysis using Bayesian Fused approaches in spatial omics datasets.
bioRxiv, 2025-03.

2\) K. Thrane, H. Eriksson, J. Maaskola, J. Hansson, and J. Lundeberg
(2018). Spatially resolved transcriptomics enables dissection of genetic
heterogeneity in stage III cutaneous malignant melanoma. Cancer
research, 78(20):5970–5979.
