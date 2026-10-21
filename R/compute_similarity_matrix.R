.legacy_trom_metadata <- function(expression_matrix,
                                  group_all,
                                  tcga_sample_size,
                                  tcga_normal_count) {
  sample_ids <- colnames(expression_matrix)
  if (is.null(group_all) || is.null(tcga_sample_size) || is.null(tcga_normal_count)) {
    stop(
      "Supply 'sample_metadata'; the legacy path requires group_all, tcga_sample_size, and tcga_normal_count.",
      call. = FALSE
    )
  }
  if (length(group_all) != length(sample_ids) || tcga_sample_size < 1 ||
      tcga_sample_size > length(sample_ids) || tcga_normal_count < 0 ||
      tcga_normal_count > tcga_sample_size) {
    stop("Legacy sample counts do not match the expression matrix.", call. = FALSE)
  }
  tissue <- ifelse(grepl("normal", group_all, ignore.case = TRUE), "normal", "tumor")
  dataset <- sub("_(normal|tumor).*$", "", group_all, ignore.case = TRUE)
  data.frame(
    sample_id = sample_ids,
    species = c(rep("human", tcga_sample_size), rep("mouse", length(sample_ids) - tcga_sample_size)),
    dataset = dataset,
    model = dataset,
    tissue = tissue,
    stringsAsFactors = FALSE
  )
}

#' Compute Similarity Matrix Between Human and Mouse Samples
#'
#' Split an expression matrix by explicit species metadata, compute TROM scores,
#' and retain only tissue-matched normal-normal and tumor-tumor comparisons.
#'
#' @param expression_matrix Numeric gene-by-sample expression matrix.
#' @param group_all,tcga_sample_size,tcga_normal_count Deprecated-compatible
#'   positional inputs used only when `sample_metadata` is `NULL`.
#' @param merge_clinic_CNV Optional human clinical/CNV table keyed by
#'   `PATIENT_ID`. Missing patients are retained with missing annotations.
#' @param cnv_genes,clinical_vars Requested annotation columns.
#' @param sample_metadata Data frame containing at least `sample_id`, `species`,
#'   `dataset`, `model`, and `tissue`.
#'
#' @return A list with raw and tissue-filtered similarity matrices, the
#'   annotated human-sample table, and metadata in expression-column order.
#'
#' @importFrom TROM ws.trom
#' @export
compute_similarity_matrix <- function(expression_matrix,
                                      group_all = NULL,
                                      tcga_sample_size = NULL,
                                      tcga_normal_count = NULL,
                                      merge_clinic_CNV = NULL,
                                      cnv_genes = character(),
                                      clinical_vars = character(),
                                      sample_metadata = NULL) {
  expression_matrix <- as.matrix(expression_matrix)
  if (!is.numeric(expression_matrix) || is.null(colnames(expression_matrix)) || is.null(rownames(expression_matrix))) {
    stop("'expression_matrix' must be a numeric matrix with gene and sample names.", call. = FALSE)
  }
  if (anyNA(expression_matrix) || any(!is.finite(expression_matrix))) {
    stop("'expression_matrix' must contain only finite, non-missing values.", call. = FALSE)
  }
  if (anyDuplicated(rownames(expression_matrix))) {
    stop("'expression_matrix' has duplicate gene identifiers.", call. = FALSE)
  }

  if (is.null(sample_metadata)) {
    warning(
      "Count- and position-based TROM inputs are deprecated; supply explicit 'sample_metadata'.",
      call. = FALSE
    )
    sample_metadata <- .legacy_trom_metadata(
      expression_matrix,
      group_all,
      tcga_sample_size,
      tcga_normal_count
    )
  }
  sample_metadata <- .validate_sample_metadata(sample_metadata, colnames(expression_matrix))

  human_metadata <- sample_metadata[sample_metadata$species == "human", , drop = FALSE]
  mouse_metadata <- sample_metadata[sample_metadata$species == "mouse", , drop = FALSE]
  if (nrow(human_metadata) == 0 || nrow(mouse_metadata) == 0) {
    stop("TROM requires at least one human and one mouse sample.", call. = FALSE)
  }

  make_trom_input <- function(sample_ids) {
    data.frame(
      gene_id = rownames(expression_matrix),
      expression_matrix[, sample_ids, drop = FALSE],
      check.names = FALSE
    )
  }
  human_data <- make_trom_input(human_metadata$sample_id)
  mouse_data <- make_trom_input(mouse_metadata$sample_id)

  similarity_matrix <- local({
    old_directory <- setwd(tempdir())
    on.exit(setwd(old_directory), add = TRUE)
    withCallingHandlers(
      TROM::ws.trom(
        sp_gene_expr = human_data,
        single = FALSE,
        sp_gene_expr2 = mouse_data,
        provide = FALSE,
        save_overlap_genes = FALSE
      ),
      warning = function(warning_condition) {
        if (grepl("Please use 'rowNames' instead of 'row.names'", conditionMessage(warning_condition), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  })
  similarity_matrix <- as.matrix(similarity_matrix)
  similarity_matrix <- similarity_matrix[
    human_metadata$sample_id,
    mouse_metadata$sample_id,
    drop = FALSE
  ]

  same_tissue <- outer(human_metadata$tissue, mouse_metadata$tissue, FUN = "==")
  filtered_similarity_matrix <- similarity_matrix
  filtered_similarity_matrix[!same_tissue] <- NA_real_

  requested_annotations <- c(clinical_vars, cnv_genes)
  if (is.null(merge_clinic_CNV)) {
    merge_clinic_CNV <- data.frame(PATIENT_ID = character(), stringsAsFactors = FALSE)
  }
  if (!is.data.frame(merge_clinic_CNV) || !"PATIENT_ID" %in% colnames(merge_clinic_CNV)) {
    stop("'merge_clinic_CNV' must be a data frame containing 'PATIENT_ID'.", call. = FALSE)
  }
  duplicate_patients <- unique(merge_clinic_CNV$PATIENT_ID[duplicated(merge_clinic_CNV$PATIENT_ID)])
  if (length(duplicate_patients) > 0) {
    stop(
      "'merge_clinic_CNV' has duplicate PATIENT_ID values: ",
      paste(utils::head(duplicate_patients, 5), collapse = ", "),
      call. = FALSE
    )
  }
  missing_annotations <- setdiff(requested_annotations, colnames(merge_clinic_CNV))
  if (length(missing_annotations) > 0) {
    stop("Clinical/CNV table is missing requested columns: ", paste(missing_annotations, collapse = ", "), call. = FALSE)
  }

  patient_ids <- gsub("-", ".", substr(human_metadata$sample_id, 1, 12), fixed = TRUE)
  clinical_index <- match(patient_ids, as.character(merge_clinic_CNV$PATIENT_ID))
  clinical_data <- merge_clinic_CNV[clinical_index, requested_annotations, drop = FALSE]
  rownames(clinical_data) <- NULL

  annotated <- data.frame(
    PATIENT_ID = patient_ids,
    sample_id = human_metadata$sample_id,
    filtered_similarity_matrix,
    tissue = human_metadata$tissue,
    dataset = human_metadata$dataset,
    model = human_metadata$model,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (length(requested_annotations) > 0) {
    annotated <- cbind(annotated, clinical_data)
  }

  list(
    similarity_matrix = similarity_matrix,
    filtered_similarity_matrix = filtered_similarity_matrix,
    dm_trom_del_backgroud_clinic = annotated,
    sample_metadata = sample_metadata
  )
}
