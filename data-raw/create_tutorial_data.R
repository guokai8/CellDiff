# Create Tutorial Seurat Object for CellDiff Testing
# This script generates a synthetic Seurat object with multiple conditions

library(Seurat)
library(Matrix)

set.seed(123)

# Parameters
n_genes <- 2000
n_cells_per_condition <- 500
conditions <- c("WT", "KO", "DKO")
cell_types <- c("T_cells", "B_cells", "Monocytes", "NK_cells", "Dendritic")

# Create synthetic expression data
create_synthetic_data <- function(n_genes, n_cells, cell_type, condition) {
  # Base expression
  base_expr <- matrix(
    rnbinom(n_genes * n_cells, mu = 0.5, size = 0.1),
    nrow = n_genes,
    ncol = n_cells
  )

  # Add cell type-specific markers
  if (cell_type == "T_cells") {
    marker_genes <- 1:50
    base_expr[marker_genes, ] <- base_expr[marker_genes, ] +
      rnbinom(50 * n_cells, mu = 5, size = 0.5)
  } else if (cell_type == "B_cells") {
    marker_genes <- 51:100
    base_expr[marker_genes, ] <- base_expr[marker_genes, ] +
      rnbinom(50 * n_cells, mu = 5, size = 0.5)
  } else if (cell_type == "Monocytes") {
    marker_genes <- 101:150
    base_expr[marker_genes, ] <- base_expr[marker_genes, ] +
      rnbinom(50 * n_cells, mu = 5, size = 0.5)
  } else if (cell_type == "NK_cells") {
    marker_genes <- 151:200
    base_expr[marker_genes, ] <- base_expr[marker_genes, ] +
      rnbinom(50 * n_cells, mu = 5, size = 0.5)
  } else if (cell_type == "Dendritic") {
    marker_genes <- 201:250
    base_expr[marker_genes, ] <- base_expr[marker_genes, ] +
      rnbinom(50 * n_cells, mu = 5, size = 0.5)
  }

  # Add condition-specific effects (simulate differential communication)
  if (condition == "KO") {
    # Reduce some ligand expression in KO
    ligand_genes <- 300:320
    base_expr[ligand_genes, ] <- base_expr[ligand_genes, ] * 0.3
  } else if (condition == "DKO") {
    # Further reduce ligand expression in DKO
    ligand_genes <- 300:320
    base_expr[ligand_genes, ] <- base_expr[ligand_genes, ] * 0.1
    # Also reduce receptor expression
    receptor_genes <- 350:370
    base_expr[receptor_genes, ] <- base_expr[receptor_genes, ] * 0.4
  }

  return(base_expr)
}

# Generate data for each condition and cell type
all_data <- list()
all_metadata <- list()
counter <- 1

for (cond in conditions) {
  for (ct in cell_types) {
    n_cells_this_type <- round(n_cells_per_condition / length(cell_types))

    # Create expression data
    expr_data <- create_synthetic_data(n_genes, n_cells_this_type, ct, cond)

    # Create cell names
    cell_names <- paste0(cond, "_", ct, "_", 1:n_cells_this_type)
    colnames(expr_data) <- cell_names

    # Create metadata
    metadata <- data.frame(
      cell_id = cell_names,
      condition = cond,
      cell_type = ct,
      nCount_RNA = colSums(expr_data),
      nFeature_RNA = colSums(expr_data > 0),
      row.names = cell_names,
      stringsAsFactors = FALSE
    )

    all_data[[counter]] <- expr_data
    all_metadata[[counter]] <- metadata
    counter <- counter + 1
  }
}

# Combine all data
combined_expr <- do.call(cbind, all_data)
combined_metadata <- do.call(rbind, all_metadata)

# Create gene names (include some ligand-receptor gene names)
gene_names <- c(
  # Ligands (300-349)
  paste0("LIGAND_", 1:50),
  # Receptors (350-399)
  paste0("RECEPTOR_", 1:50),
  # Cell type markers
  paste0("CD", 1:200),
  # Other genes
  paste0("GENE_", 1:(n_genes - 300))
)
rownames(combined_expr) <- gene_names

# Create Seurat object
pbmc_tutorial <- CreateSeuratObject(
  counts = combined_expr,
  meta.data = combined_metadata,
  project = "CellDiff_Tutorial"
)

# Add some common genes that CellChat might recognize
# This is simplified - real data would have actual gene names
pbmc_tutorial <- NormalizeData(pbmc_tutorial)
pbmc_tutorial <- FindVariableFeatures(pbmc_tutorial, nfeatures = 500)
pbmc_tutorial <- ScaleData(pbmc_tutorial)
pbmc_tutorial <- RunPCA(pbmc_tutorial, npcs = 20, verbose = FALSE)
pbmc_tutorial <- RunUMAP(pbmc_tutorial, dims = 1:20, verbose = FALSE)

# Print summary
cat("\n=======================================================\n")
cat("Tutorial Seurat Object Created\n")
cat("=======================================================\n\n")
cat("Conditions:", paste(unique(pbmc_tutorial$condition), collapse = ", "), "\n")
cat("Cell types:", paste(unique(pbmc_tutorial$cell_type), collapse = ", "), "\n")
cat("Total cells:", ncol(pbmc_tutorial), "\n")
cat("Total genes:", nrow(pbmc_tutorial), "\n\n")

cat("Cells per condition:\n")
print(table(pbmc_tutorial$condition))
cat("\nCells per cell type:\n")
print(table(pbmc_tutorial$cell_type))
cat("\nCells per condition and cell type:\n")
print(table(pbmc_tutorial$condition, pbmc_tutorial$cell_type))

# Save the object
save(pbmc_tutorial, file = "data/pbmc_tutorial.rda", compress = "xz")

cat("\n\nSaved to: data/pbmc_tutorial.rda\n")
cat("\nUsage:\n")
cat("  data(pbmc_tutorial)\n")
cat("  cellchat_list <- runCellChat(pbmc_tutorial, group.by = 'condition', species = 'human')\n")
cat("  results <- compareCell(cellchat_list, reference = 'WT')\n\n")
