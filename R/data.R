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
#' @export
"mela_coords"

#' SpaceBF_melanoma dataset
#'
#' The full melanoma dataset, consisting gene expression, xy-coordinates, and LR pairs
#'
#' @format A list of three data frames, with gene expression, co-ordinates, and LR pairs
#' \describe{
#'   \item{gene_data_mela}{gene expression data}
#'   \item{coords_data_mela}{xy co-ordinates}
#'   \item{LR_pairs}{161 LR pairs}
#' }
#' @source \url{https://aacrjournals.org/cancerres/article/78/20/5970/631815/Spatially-Resolved-Transcriptomics-Enables}
#' @usage data(SpaceBF_melanoma)
#' @export
"SpaceBF_melanoma"

#' SpaceBF_cSCC dataset
#'
#' The full cSCC dataset, consisting gene expression, xy-coordinates, and keratin pairs
#'
#' @format A list of three data frames, with gene expression, co-ordinates, and keratin pairs
#' @source \url{https://www.cell.com/cell/fulltext/S0092-8674(20)30672-3}
#' @usage data(SpaceBF_cSCC)
#' @export
"SpaceBF_cSCC"