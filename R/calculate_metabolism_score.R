#' Calculate Metabolism Signature Score
#'
#' This function computes the metabolism-related gene signature score using different methods.
#'
#' @importFrom IOBR calculate_sig_score
#'
#' @param eset_tpm Data frame or matrix. TPM-transformed RNA-seq expression matrix,
#'   output from `plot_TME_barplot()`.
#' @param input_source Character. Recorded TPM source: `"provided_tpm"` or
#'   `"count2tpm"`. Defaults to the source attached by `plot_TME_barplot()` or
#'   `"provided_tpm"` for directly supplied TPM data. ComBat-corrected input is
#'   rejected. An external transformation without trusted provenance cannot be
#'   identified from its values alone; declare `"combat_corrected"` when
#'   applicable.
#' @param method Character. Method for score calculation ("pca", "ssgsea", "zscore", "integration"). Default: "pca".
#' @param mini_gene_count Integer. Minimum gene count required per signature. Default: 2.
#'
#' @return A data frame containing metabolism scores per sample.
#'
#' @examples
#' \dontrun{
#' sig_meta <- calculate_metabolism_score(eset_tpm, method = "pca", mini_gene_count = 2)
#' }
#'
#' @export
calculate_metabolism_score <- function(eset_tpm,
                                       method = "pca",
                                       mini_gene_count = 2,
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
  if (!method %in% c("pca", "ssgsea", "zscore", "integration")) {
    stop("Error: 'method' must be one of 'pca', 'ssgsea', 'zscore', or 'integration'.")
  }

  # Calculate signature score using IOBR::calculate_sig_score

  sig_meta <- calculate_sig_score(
    pdata = NULL,
    eset = as.data.frame(eset_tpm),
    signature = signature_metabolism,
    method = method, # It can also be ssgsea or Z-score method
    mini_gene_count = mini_gene_count
  )

  attr(sig_meta, "MurineMatchTx_tpm_source") <- input_info$input_source

  return(sig_meta)
}
