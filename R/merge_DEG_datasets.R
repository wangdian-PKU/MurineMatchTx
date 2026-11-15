#' Merge Differential Expression Gene (DEG) Datasets
#'
#' Read human and mouse DEG result files through one explicit input contract,
#' classify genes with shared thresholds, and combine the normalized results.
#'
#' @param human_file_path Path to the human DEG result file.
#' @param mouse_files Named list of paths to mouse DEG result files.
#' @param logFC_value Deprecated-compatible fold-change threshold used only when
#'   `deg_thresholds` is `NULL`.
#' @param FDR_value Deprecated-compatible adjusted-p threshold used only when
#'   `deg_thresholds` is `NULL`.
#' @param deg_method One of `"edgeR"`, `"DESeq2"`, `"limma"`, or `"custom"`.
#'   A named character vector may specify a method for `TCGA` and each mouse
#'   dataset separately.
#' @param lfc_col,padj_col,gene_col Optional explicit column names. Each may be
#'   scalar or named by dataset. Custom inputs require `lfc_col` and `padj_col`.
#' @param deg_thresholds A list with `log2fc` and `fdr`. The classification rule
#'   is `logFC >= log2fc` or `logFC <= -log2fc`, together with `FDR < fdr`.
#'
#' @return A data frame with `gene`, `logFC`, `FDR`, `contrast`, and `change`.
#'   The validated thresholds are stored in the `deg_thresholds` attribute.
#'
#' @importFrom dplyr bind_rows
#' @export
merge_DEG_datasets <- function(human_file_path,
                               mouse_files,
                               logFC_value = 1,
                               FDR_value = 0.05,
                               deg_method = "edgeR",
                               lfc_col = NULL,
                               padj_col = NULL,
                               gene_col = NULL,
                               deg_thresholds = NULL) {
  if (!is.list(mouse_files) || is.null(names(mouse_files)) || any(!nzchar(names(mouse_files)))) {
    stop("'mouse_files' must be a named list.", call. = FALSE)
  }
  thresholds <- .validate_deg_thresholds(deg_thresholds, logFC_value, FDR_value)

  human <- .normalize_deg_file(
    human_file_path, "TCGA", deg_method, lfc_col, padj_col, gene_col
  )
  human$contrast <- "TCGA"

  mouse <- lapply(names(mouse_files), function(dataset) {
    normalized <- .normalize_deg_file(
      mouse_files[[dataset]], dataset, deg_method, lfc_col, padj_col, gene_col
    )
    normalized$contrast <- dataset
    normalized
  })

  merged <- dplyr::bind_rows(c(list(human), mouse))
  merged$change <- .classify_deg(merged, thresholds)
  merged$contrast <- factor(merged$contrast, levels = c("TCGA", names(mouse_files)))
  attr(merged, "deg_thresholds") <- thresholds
  merged
}
