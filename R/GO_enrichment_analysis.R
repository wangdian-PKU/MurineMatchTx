.get_symbol_map <- function(species) {
  if (species == "hsa") {
    return(toTable(get("org.Hs.egSYMBOL", envir = asNamespace("org.Hs.eg.db"))))
  }
  if (species == "mmu") {
    return(toTable(get("org.Mm.egSYMBOL", envir = asNamespace("org.Mm.eg.db"))))
  }
  toTable(get("org.Rn.egSYMBOL", envir = asNamespace("org.Rn.eg.db")))
}

#' Gene Ontology (GO) Enrichment Analysis
#'
#' This function performs GO enrichment analysis (BP, MF, CC)
#' for normalized differential expression results.
#'
#' @description
#' - Extracts significant DEGs with the shared inclusive log2 fold-change and
#'   exclusive adjusted-p thresholds.
#' - Reads input data automatically from .tsv or .csv files.
#' - Converts gene IDs based on the selected species (hsa, mmu, rno).
#' - Supports user-defined ontology (BP, MF, CC).
#' - Performs GO enrichment analysis and saves the plot to a PDF file.
#'
#' @details
#' - Input Data: The function accepts explicit edgeR, DESeq2, limma, or custom
#'   DEG formats.
#' - Supported File Formats: .tsv (tab-separated) and .csv (comma-separated).
#' - Supported Species: hsa (human), mmu (mouse), rno (rat).
#' - Ontology Options (ont):
#'   - "BP": Biological Process (default)
#'   - "MF": Molecular Function
#'   - "CC": Cellular Component
#' - Visualization: GO term similarity is visualized using simplifyGO().
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom simplifyEnrichment simplifyGO GO_similarity
#' @importFrom dplyr filter mutate
#' @importFrom AnnotationDbi toTable
#' @importFrom tools file_ext
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom org.Rn.eg.db org.Rn.eg.db
#'
#' @param deg_file Character. Path to the edgeR DEG results file (.tsv or .csv).
#' @param species Character. One of hsa, mmu, or rno.
#' @param org_db Character. The corresponding `OrgDb` annotation database (e.g., "org.Hs.eg.db").
#' @param ont Character. Ontology type for GO enrichment ("BP", "MF", "CC"). Default is "BP".
#' @param column_title Character. The title used for GO similarity visualization.
#' @param width Integer. The width of PDF (default: 10.36).
#' @param height Integer. The height of PDF (default: 6.35).
#' @param output_path Character. Directory to save the PDF file (default: "./enrichment").
#' @param deg_method One of `"edgeR"`, `"DESeq2"`, `"limma"`, or `"custom"`.
#' @param lfc_col,padj_col,gene_col Optional explicit column names. Custom
#'   inputs require `lfc_col` and `padj_col`.
#' @param deg_thresholds A list with `log2fc` and `fdr` shared with the other
#'   DEG modules.
#'
#' @return Saves the GO enrichment analysis plot as a .pdf file.
#'
#' @examples
#' \dontrun{
#' # Example usage for human (TCGA data)
#' GO_enrichment_analysis(
#'   deg_file = "./data/TCGA_DEG.tsv",
#'   species = "hsa",
#'   org_db = "org.Hs.eg.db",
#'   ont = "BP",
#'   column_title = "TCGA GO BP terms"
#' )
#' }
#'
#' @export
GO_enrichment_analysis <- function(deg_file,
                                   species,
                                   org_db,
                                   ont = "BP",
                                   column_title,
                                   width = 10.36,
                                   height = 6.35,
                                   output_path = "./enrichment",
                                   deg_method = "edgeR",
                                   lfc_col = NULL,
                                   padj_col = NULL,
                                   gene_col = NULL,
                                   deg_thresholds = NULL) {
  # Check if the "species" parameter is valid
  species_list <- c("hsa", "mmu", "rno")
  if (!species %in% species_list) {
    stop("Error: `species` must be one of 'hsa', 'mmu', or 'rno'.")
  }

  # Check if the "ont" parameter is valid
  if (!ont %in% c("BP", "MF", "CC")) {
    stop("Error: `ont` must be one of 'BP', 'MF', or 'CC'.")
  }

  thresholds <- .validate_deg_thresholds(deg_thresholds)
  deg_data <- .normalize_deg_file(
    deg_file,
    "DEG",
    deg_method,
    lfc_col,
    padj_col,
    gene_col
  )

  # Load species-specific gene ID mapping
  EG2Symbol <- .get_symbol_map(species)

  # Filter significant DEGs using the shared rule.
  geneLists <- deg_data[.classify_deg(deg_data, thresholds) != "Stable", , drop = FALSE]
  geneLists$symbol <- geneLists$gene

  # Convert gene IDs
  results <- merge(geneLists, EG2Symbol, by = "symbol", all.x = TRUE)
  results <- na.omit(results)
  id <- results$gene_id # Extract gene ID column

  # Perform GO enrichment analysis
  GO_results <- enrichGO(
    OrgDb = get(org_db), gene = id,
    ont = ont, pvalueCutoff = 0.05, readable = TRUE
  )

  GO_result_df <- as.data.frame(GO_results)

  go_terms <- na.omit(GO_result_df$ID)
  go_terms <- go_terms[grepl("^GO:\\d{7}$", go_terms)]

  if (length(go_terms) < 2) {
    stop("Error: The number of effective GO terms is too small (< 2) to conduct go similarity analysis.
         Please check whether the input differential expression matrix is too small or the annotation fails.")
  }

  mat <- GO_similarity(go_terms, ont = ont)


  # Create output directory
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  # Define output file name
  pdf_file <- file.path(output_path, paste0(column_title, ".pdf"))

  # Save GO enrichment analysis result to PDF
  grDevices::pdf(pdf_file, width = width, height = height)
  simplifyGO(mat, column_title = column_title)
  grDevices::dev.off()

  message("GO enrichment analysis saved to: ", pdf_file)
}
