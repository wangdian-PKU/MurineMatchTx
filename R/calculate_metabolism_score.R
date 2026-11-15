#' Calculate Metabolism Signature Score
#'
#' This function computes scores for the complete IOBR metabolism signature
#' collection. The recommended `"relative_zscore"` method standardizes each
#' gene across samples within one dataset and then averages the standardized
#' expression of the genes in each pathway. Run it separately for TCGA and for
#' each mouse dataset before combining the returned score tables.
#'
#' This function and `plot_TME_boxplot()` independently consume the TPM matrix
#' returned by `plot_TME_barplot()`. The metabolism scores, rather than the TME
#' boxplots, are the input to `plot_metabolism_heatmap()`.
#'
#' When current human-reference metabolism signatures are applied to mouse samples, the
#' resulting scores are exploratory, method-specific descriptions. They do not enter cross-species similarity
#' scoring or model-selection decisions and do not establish cross-species validation without
#' mouse-specific signatures or orthogonal metabolic evidence.
#'
#' @importFrom IOBR calculate_sig_score
#'
#' @param eset_tpm Data frame or matrix. TPM-transformed RNA-seq expression matrix,
#'   output from `plot_TME_barplot()`.
#' @param input_source Character. Recorded TPM source: `"provided_tpm"` or
#'   `"count2tpm"`. Required for external TPM data without the trusted source
#'   attached by `plot_TME_barplot()`. ComBat-corrected input is rejected. This
#'   0.99.0 pre-release contract avoids treating unlabelled external
#'   transformations as TPM data.
#' @param method Character. Method for score calculation. The recommended
#'   default is `"relative_zscore"`; the original IOBR methods `"pca"`,
#'   `"ssgsea"`, `"zscore"`, and `"integration"` remain available.
#' @param mini_gene_count Integer. Minimum number of available genes required
#'   per pathway. The default of 1 retains single-gene pathways. Pathways below
#'   the threshold are kept as `NA` rather than silently removed.
#'
#' @return A data frame containing metabolism scores per sample.
#'
#' @examples
#' \dontrun{
#' sig_meta <- calculate_metabolism_score(eset_tpm, method = "relative_zscore",
#'                                         mini_gene_count = 1,
#'                                         input_source = "provided_tpm")
#' }
#'
#' @seealso [plot_TME_barplot()], [plot_TME_boxplot()],
#'   [plot_metabolism_heatmap()]
#'
#' @export
calculate_metabolism_score <- function(eset_tpm,
                                       method = "relative_zscore",
                                       mini_gene_count = 1,
                                       input_source = NULL) {
  # Input validation
  if (!is.data.frame(eset_tpm) && !is.matrix(eset_tpm)) {
    stop("Error: 'eset_tpm' must be a data frame or matrix.")
  }
  input_info <- .validate_expression_input(
    data_matrix = eset_tpm,
    input_type = "TPM",
    input_source = input_source
  )
  valid_methods <- c("relative_zscore", "pca", "ssgsea", "zscore", "integration")
  if (!is.character(method) || length(method) != 1L || !method %in% valid_methods) {
    stop(
      "Error: 'method' must be one of 'relative_zscore', 'pca', ",
      "'ssgsea', 'zscore', or 'integration'."
    )
  }
  if (!is.numeric(mini_gene_count) || length(mini_gene_count) != 1L ||
      is.na(mini_gene_count) || mini_gene_count < 1 || mini_gene_count %% 1 != 0) {
    stop("Error: 'mini_gene_count' must be a positive integer.")
  }

  expression_matrix <- as.matrix(eset_tpm)
  if (is.null(rownames(expression_matrix)) || anyDuplicated(rownames(expression_matrix))) {
    stop("Error: 'eset_tpm' must have unique gene symbols as row names.")
  }
  if (is.null(colnames(expression_matrix)) || anyDuplicated(colnames(expression_matrix))) {
    stop("Error: 'eset_tpm' must have unique sample IDs as column names.")
  }

  signatures <- .get_metabolism_signatures()
  pathway_gene_counts <- vapply(
    signatures,
    function(genes) length(intersect(genes, rownames(expression_matrix))),
    integer(1)
  )

  if (method == "relative_zscore") {
    if (ncol(expression_matrix) < 2L) {
      stop("Error: 'relative_zscore' requires at least two samples.")
    }

    log_expression <- log2(expression_matrix + 1)
    gene_means <- rowMeans(log_expression)
    gene_sds <- apply(log_expression, 1, stats::sd)
    standardized <- sweep(log_expression, 1, gene_means, FUN = "-")
    standardized <- sweep(standardized, 1, gene_sds, FUN = "/")
    standardized[!is.finite(standardized)] <- 0

    score_matrix <- matrix(
      NA_real_,
      nrow = ncol(standardized),
      ncol = length(signatures),
      dimnames = list(colnames(standardized), names(signatures))
    )

    eligible <- names(signatures)[pathway_gene_counts >= mini_gene_count]
    for (pathway in eligible) {
      genes <- intersect(signatures[[pathway]], rownames(standardized))
      score_matrix[, pathway] <- colMeans(standardized[genes, , drop = FALSE])
    }

    sig_meta <- data.frame(
      ID = rownames(score_matrix),
      score_matrix,
      check.names = FALSE,
      row.names = NULL
    )
  } else {
    sig_meta <- calculate_sig_score(
      pdata = NULL,
      eset = as.data.frame(expression_matrix),
      signature = signatures,
      method = method,
      mini_gene_count = mini_gene_count
    )

    missing_pathways <- setdiff(names(signatures), colnames(sig_meta))
    sig_meta[missing_pathways] <- NA_real_
    sig_meta <- as.data.frame(sig_meta[, c("ID", names(signatures)), drop = FALSE],
                              check.names = FALSE)
  }

  attr(sig_meta, "MurineMatchTx_tpm_source") <- input_info$input_source
  attr(sig_meta, "MurineMatchTx_metabolism_method") <- method
  attr(sig_meta, "MurineMatchTx_pathway_gene_counts") <- pathway_gene_counts

  return(sig_meta)
}
