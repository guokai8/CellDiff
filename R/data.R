#' Example CellChat objects dataset
#'
#' A list containing two example CellChat objects representing different conditions.
#' These objects have been preprocessed and contain the necessary slots for cell-cell
#' communication analysis using the ContriDiff function and other CellDiff package functions.
#'
#' @format A named list with 2 CellChat objects:
#' \describe{
#'   \item{NL}{CellChat object representing the Normal Lung condition}
#'   \item{LS}{CellChat object representing the Lung Scleroderma condition}
#' }
#'
#' @details
#' Each CellChat object contains various slots including:
#' \itemize{
#'   \item{idents}{Cell type identities}
#'   \item{net}{Communication network information}
#'   \item{LR}{Ligand-receptor pair information}
#'   \item{DB}{Database of interactions}
#' }
#'
#' These CellChat objects have been processed with the standard CellChat workflow:
#' 1. Created from single-cell expression data
#' 2. Identification of overexpressed ligands/receptors
#' 3. Inference of cell-cell communication network
#' 4. Network analysis
#'
#' @source Lung single-cell RNA sequencing data from normal and scleroderma patients.
#'
#' @examples
#' # Load the example data
#' data(celllist)
#'
#' # View structure of the list
#' str(celllist, max.level = 1)
#'
#' # Access the individual CellChat objects
#' nl_cellchat <- celllist$NL
#' ls_cellchat <- celllist$LS
#'
#' # Compare ligand-receptor contributions
#' ContriDiff(
#'   object.list = celllist,
#'   signaling = "TNF"
#' )
"celllist"
