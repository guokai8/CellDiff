% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatDiffM.R
\name{heatDiffM}
\alias{heatDiffM}
\title{Compare signaling patterns across multiple conditions with a heatmap visualization}
\usage{
heatDiffM(
  object.list,
  comparison = NULL,
  reference = NULL,
  measure = c("both", "sender", "receiver"),
  signaling = NULL,
  color.heatmap = "custom",
  color.use = NULL,
  title = NULL,
  width = 10,
  height = 8,
  font.size = 8,
  font.size.title = 10,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  show_values = FALSE,
  value_digits = 2,
  signif_mark = FALSE,
  border_color = "white",
  show_heatmap_border = TRUE,
  border_width = 1,
  use_log2fc = FALSE,
  transpose = TRUE,
  return_data = FALSE,
  save_plot = FALSE,
  save_name = "heatmap_diff.pdf",
  filter_zeros = FALSE,
  filter_thresh = 0,
  slot.name = "netP",
  cluster_conditions = FALSE,
  split_heatmap = FALSE,
  highlight_reference = TRUE,
  show_condition_names = TRUE,
  show_annotation = c("both", "row", "column", "none"),
  comparison_method = c("all_vs_ref", "all_vs_all"),
  big_heatmap = FALSE
)
}
\arguments{
\item{object.list}{A list of CellChat objects to compare.}

\item{comparison}{A numeric vector specifying the indices of the objects to compare.
If NULL (default), all objects will be compared.}

\item{reference}{Integer or character, index or name of the reference object to compare against.
If NULL (default), the first object in comparison will be used as reference.}

\item{measure}{Character, specifies the signaling role to analyze: "sender" (outgoing signaling),
"receiver" (incoming signaling), or "both" (overall signaling).}

\item{signaling}{Character vector of pathway names to include in the analysis. If NULL (default),
all pathways present in any condition will be used.}

\item{color.heatmap}{Either a string specifying a predefined color scheme ("custom", "RdBu", "BuOr")
or a vector of colors to create a custom gradient. Default is "custom" (blue-white-red).}

\item{color.use}{Vector of colors for cell type annotations. If NULL, colors will be generated automatically.}

\item{title}{Character, custom title for the heatmap. If NULL, a default title will be generated.}

\item{width}{Numeric, width of the heatmap in centimeters.}

\item{height}{Numeric, height of the heatmap in centimeters.}

\item{font.size}{Numeric, font size for row and column labels.}

\item{font.size.title}{Numeric, font size for the heatmap title.}

\item{cluster_rows}{Logical, whether to cluster rows (cell types in default orientation).}

\item{cluster_cols}{Logical, whether to cluster columns (pathways in default orientation).}

\item{show_row_dend}{Logical, whether to display the row dendrogram when clustering.}

\item{show_column_dend}{Logical, whether to display the column dendrogram when clustering.}

\item{show_rownames}{Logical, whether to display row names.}

\item{show_colnames}{Logical, whether to display column names.}

\item{show_values}{Logical, whether to display the difference values in cells.}

\item{value_digits}{Integer, number of decimal places to display for cell values.}

\item{signif_mark}{Logical, whether to mark significant differences.}

\item{border_color}{Character, color for cell borders.}

\item{show_heatmap_border}{Logical, whether to display cell borders.}

\item{border_width}{Numeric, width of cell borders.}

\item{use_log2fc}{Logical, whether to use log2 fold change instead of simple difference.}

\item{transpose}{Logical, whether to transpose the heatmap (pathways as rows, cell types as columns).}

\item{return_data}{Logical, whether to return a list with the heatmap and data matrices (TRUE)
or just the heatmap object (FALSE).}

\item{save_plot}{Logical, whether to save the heatmap to a file.}

\item{save_name}{Character, filename for saving the plot (if save_plot=TRUE).}

\item{filter_zeros}{Logical, whether to filter out pathways with no signal.}

\item{filter_thresh}{Numeric, threshold for filtering pathways based on signal strength.}

\item{slot.name}{Character, name of the slot to extract data from (default: "netP").}

\item{cluster_conditions}{Logical, whether to cluster conditions (default: FALSE).}

\item{split_heatmap}{Logical, whether to split the heatmap by measure (useful when comparing multiple measures).}

\item{highlight_reference}{Logical, whether to highlight the reference condition (default: TRUE).}

\item{show_condition_names}{Logical, whether to display condition names (default: TRUE).}

\item{show_annotation}{Character, type of annotation to show: "both", "row", "column", or "none".}

\item{comparison_method}{Character, method for comparing conditions: "all_vs_ref" (all vs. reference) or
"all_vs_all" (pairwise comparisons between all conditions).}

\item{big_heatmap}{Logical, whether to create a single big heatmap with side-by-side comparisons (default: FALSE).}
}
\value{
If return_data=TRUE, returns a list containing:
  \itemize{
    \item heatmap: The ComplexHeatmap object
    \item data_list: List of matrices for each condition
    \item normalized_list: List of normalized matrices for each condition
    \item diff_matrix: Difference matrix used for visualization (normalized)
    \item comparison: The comparison indices used
    \item reference: The reference index used
    \item condition_names: Names of the conditions compared
    \item measure: The measure type used
    \item pathways: Names of pathways included in the analysis
    \item cell_types: Names of cell types included in the analysis
  }
  If return_data=FALSE, returns only the ComplexHeatmap object.
}
\description{
This function creates a differential heatmap visualization comparing signaling patterns
across multiple conditions. It calculates normalized signaling strength for each cell type
and pathway, then visualizes the differences between conditions.
}
\examples{
# Basic usage with default parameters
heatDiffM(cellchat.list, measure = "sender")

# Advanced usage with customization
heatDiffM(cellchat.list,
         measure = "both",
         reference = "Normal",
         color.heatmap = c("purple", "white", "green"),
         use_log2fc = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_values = TRUE,
         filter_thresh = 0.1)

# Big side-by-side heatmap
heatDiffM(cellchat.list,
         measure = "sender",
         big_heatmap = TRUE)

}
