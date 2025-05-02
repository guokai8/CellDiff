#' Example CellChat objects dataset with multiple conditions
#'
#' A list containing three example CellChat objects representing different genetic conditions.
#' These objects have been preprocessed and contain the necessary slots for cell-cell
#' communication analysis using CellDiff package functions, particularly for multi-condition comparison.
#'
#' @format A named list with 3 CellChat objects:
#' \describe{
#'   \item{WT}{CellChat object representing the Wild Type condition (reference)}
#'   \item{KO}{CellChat object representing the Knockout condition}
#'   \item{DKO}{CellChat object representing the Double Knockout condition}
#' }
#'
#' @details
#' Each CellChat object contains various slots including:
#' \itemize{
#'   \item{idents}{Cell type identities}
#'   \item{net}{Communication network information}
#'   \item{netP}{Pathway-level communication network information}
#'   \item{LR}{Ligand-receptor pair information}
#'   \item{DB}{Database of interactions}
#' }
#'
#' These CellChat objects have been processed with the standard CellChat workflow:
#' 1. Created from single-cell expression data
#' 2. Identification of overexpressed ligands/receptors
#' 3. Inference of cell-cell communication network
#' 4. Computation of network centrality scores
#' 5. Aggregation of network at the pathway level
#'
#' The WT condition serves as the reference for multi-condition comparisons.
#'
#' @note The dataset contains synthetic data generated for demonstration purposes.
#'
#' @source Single-cell RNA sequencing data from a mouse model comparing wild-type,
#'   knockout, and double knockout conditions.
#'
#' @examples
#' # Load the example data
#' data(cellchatlist)
#'
#' # View structure of the list
#' str(cellchatlist, max.level = 1)
#'
#' # Access the individual CellChat objects
#' wt_cellchat <- cellchatlist$WT
#' ko_cellchat <- cellchatlist$KO
#' dko_cellchat <- cellchatlist$DKO
#'
#' # Multi-condition comparisons with WT as reference
#' # Compare pathway differences
#' pathway_results <- rankDiffM(
#'   object.list = cellchatlist,
#'   comparison_method = "all_vs_ref",
#'   reference = "WT"
#' )
#'
#' # Create multi-condition heatmap
#' heatDiffM(
#'   object.list = cellchatlist,
#'   comparison = c("WT", "KO", "DKO"),
#'   reference = "WT",
#'   measure = "sender"
#' )
#'
#' # Compare sender-receiver roles
#' scatterDiff2DM(
#'   object.list = cellchatlist,
#'   comparison_method = "all_vs_ref",
#'   reference = "WT"
#' )
"cellchatlist"

#' Original example CellChat objects dataset with two conditions
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
#'   \item{netP}{Pathway-level communication network information}
#'   \item{LR}{Ligand-receptor pair information}
#'   \item{DB}{Database of interactions}
#' }
#'
#' These CellChat objects have been processed with the standard CellChat workflow:
#' 1. Created from single-cell expression data
#' 2. Identification of overexpressed ligands/receptors
#' 3. Inference of cell-cell communication network
#' 4. Computation of network centrality scores
#' 5. Aggregation of network at the pathway level
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
