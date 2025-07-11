#' mela_coords dataset
#'
#' The location coordinates of the melanoma dataset
#'
#' @format A data frame with N = 293 rows (spots) and 3 columns (xy coordinates and celltypes):
#' \describe{
#'   \item{x}{x co-ordinate of the cell}
#'   \item{y}{y co-ordinate of the cell}
#'   \item{celltypes}{Celltypes detected by the RCTD software based on the gene expression profile}
#' }
#' @source \url{https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables}
#' @usage data(mela_coords)
"mela_coords"

#' SpaceBF_melanoma dataset
#'
#' The full melanoma dataset, consisting gene expression, xy-coordinates, and LR pairs
#'
#' @format A data frame with N = 292 rows (spots, 1 removed based on QC) and 3 columns (xy coordinates and celltypes):
#' @source \url{https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables}
#' @usage data(mela_coords)
"SpaceBF_melanoma"
