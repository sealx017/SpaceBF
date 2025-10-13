#' mela_coords dataset
#'
#' The location coordinates of the melanoma dataset
#'
#' @format A data frame with N = 293 rows (spots) and 3 columns (xy coordinates and celltypes):
#' \describe{
#'   \item{x}{x co-ordinate of the cell}
#'   \item{y}{y co-ordinate of the cell}
#'   \item{celltypes}{Celltypes detected by RCTD}
#' }
#' @source \url{https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables}
#' @usage data(mela_coords)
#' @keywords datasets
"mela_coords"

#' SpaceBF_melanoma dataset
#'
#' The full melanoma dataset: gene expression, xy-coordinates, and LR pairs
#'
#' @format A list with components:
#' \describe{
#'   \item{gene_data_mela}{gene expression data}
#'   \item{coords_data_mela}{xy co-ordinates}
#'   \item{LR_pairs}{161 ligand–receptor pairs}
#' }
#' @source \url{https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables}
#' @usage data(SpaceBF_melanoma)
#' @keywords datasets
"SpaceBF_melanoma"

#' SpaceBF_cSCC dataset
#'
#' The full cSCC dataset: gene expression, xy-coordinates, and keratin pairs
#'
#' @format A list of three data frames
#' @source \url{https://www.cell.com/cell/fulltext/S0092-8674(20)30672-3}
#' @usage data(SpaceBF_cSCC)
#' @keywords datasets
"SpaceBF_cSCC"
