.validate_deg_thresholds <- function(deg_thresholds = NULL, log2fc = 1, fdr = 0.05) {
  if (is.null(deg_thresholds)) {
    deg_thresholds <- list(log2fc = log2fc, fdr = fdr)
  }
  if (!is.list(deg_thresholds) || !all(c("log2fc", "fdr") %in% names(deg_thresholds))) {
    stop("'deg_thresholds' must be a list with 'log2fc' and 'fdr'.", call. = FALSE)
  }

  log2fc <- deg_thresholds$log2fc
  fdr <- deg_thresholds$fdr
  if (!is.numeric(log2fc) || length(log2fc) != 1 || is.na(log2fc) || !is.finite(log2fc) || log2fc < 0) {
    stop("'deg_thresholds$log2fc' must be one non-negative finite number.", call. = FALSE)
  }
  if (!is.numeric(fdr) || length(fdr) != 1 || is.na(fdr) || !is.finite(fdr) || fdr < 0 || fdr > 1) {
    stop("'deg_thresholds$fdr' must be one finite number between 0 and 1.", call. = FALSE)
  }
  list(log2fc = as.numeric(log2fc), fdr = as.numeric(fdr))
}

.classify_deg <- function(deg_data, deg_thresholds) {
  thresholds <- .validate_deg_thresholds(deg_thresholds)
  if (!all(c("logFC", "FDR") %in% colnames(deg_data))) {
    stop("Normalized DEG data must contain 'logFC' and 'FDR'.", call. = FALSE)
  }
  ifelse(
    deg_data$FDR < thresholds$fdr & deg_data$logFC >= thresholds$log2fc,
    "Up",
    ifelse(
      deg_data$FDR < thresholds$fdr & deg_data$logFC <= -thresholds$log2fc,
      "Down",
      "Stable"
    )
  )
}

.normalize_deg_input <- function(deg_data,
                                 deg_method = c("edgeR", "DESeq2", "limma", "custom"),
                                 lfc_col = NULL,
                                 padj_col = NULL,
                                 gene_col = NULL) {
  deg_method <- match.arg(deg_method)
  if (!is.data.frame(deg_data)) {
    stop("DEG input must be a data frame.", call. = FALSE)
  }

  defaults <- switch(
    deg_method,
    edgeR = c(lfc = "logFC", padj = "FDR"),
    DESeq2 = c(lfc = "log2FoldChange", padj = "padj"),
    limma = c(lfc = "logFC", padj = "adj.P.Val"),
    custom = c(lfc = NA_character_, padj = NA_character_)
  )
  if (deg_method == "custom" && (is.null(lfc_col) || is.null(padj_col))) {
    stop("'custom' DEG input requires explicit 'lfc_col' and 'padj_col'.", call. = FALSE)
  }
  lfc_col <- if (is.null(lfc_col)) unname(defaults["lfc"]) else lfc_col
  padj_col <- if (is.null(padj_col)) unname(defaults["padj"]) else padj_col

  requested_columns <- c(lfc_col, padj_col, gene_col)
  requested_columns <- requested_columns[!is.null(requested_columns) & !is.na(requested_columns)]
  missing_columns <- setdiff(requested_columns, colnames(deg_data))
  if (length(missing_columns) > 0) {
    stop("DEG input is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  if (is.null(gene_col)) {
    genes <- rownames(deg_data)
    automatic_rows <- identical(genes, as.character(seq_len(nrow(deg_data))))
    if (nrow(deg_data) > 0 && (is.null(genes) || automatic_rows)) {
      stop("Gene identifiers must be row names or supplied with 'gene_col'.", call. = FALSE)
    }
  } else {
    genes <- deg_data[[gene_col]]
  }

  log_fc <- deg_data[[lfc_col]]
  adjusted_p <- deg_data[[padj_col]]
  if (!is.numeric(log_fc) || !is.numeric(adjusted_p)) {
    stop("DEG fold-change and adjusted-p columns must be numeric.", call. = FALSE)
  }

  genes <- as.character(genes)
  if (anyNA(genes) || any(!nzchar(genes)) || anyNA(log_fc) || anyNA(adjusted_p)) {
    stop("DEG input contains missing gene, fold-change, or adjusted-p values.", call. = FALSE)
  }
  if (any(!is.finite(log_fc)) || any(!is.finite(adjusted_p))) {
    stop("DEG fold-change and adjusted-p values must be finite.", call. = FALSE)
  }
  if (any(adjusted_p < 0 | adjusted_p > 1)) {
    stop("Adjusted-p values must be between 0 and 1.", call. = FALSE)
  }

  duplicate_genes <- unique(genes[duplicated(genes)])
  if (length(duplicate_genes) > 0) {
    stop("DEG input has duplicate genes: ", paste(utils::head(duplicate_genes, 5), collapse = ", "), call. = FALSE)
  }

  data.frame(
    gene = genes,
    logFC = as.numeric(log_fc),
    FDR = as.numeric(adjusted_p),
    stringsAsFactors = FALSE
  )
}

.read_deg_file <- function(file_path, gene_col = NULL) {
  extension <- tolower(tools::file_ext(file_path))
  if (!extension %in% c("tsv", "csv")) {
    stop("Unsupported DEG file format: only .tsv and .csv are allowed.", call. = FALSE)
  }

  reader <- if (extension == "tsv") {
    function(...) utils::read.table(..., sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    function(...) utils::read.csv(..., header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (is.null(gene_col)) {
    reader(file_path, row.names = 1)
  } else {
    reader(file_path)
  }
}

.resolve_dataset_value <- function(value, dataset, argument, allow_null = TRUE) {
  if (is.null(value)) {
    if (allow_null) return(NULL)
    stop("'", argument, "' must be supplied.", call. = FALSE)
  }
  if (is.list(value)) {
    if (!is.null(names(value)) && dataset %in% names(value)) return(value[[dataset]])
    if (length(value) == 1) return(value[[1]])
  } else {
    if (!is.null(names(value)) && any(nzchar(names(value)))) {
      if (dataset %in% names(value)) return(unname(value[[dataset]]))
    } else if (length(value) == 1) {
      return(value[[1]])
    }
  }
  stop("'", argument, "' must be scalar or named for dataset '", dataset, "'.", call. = FALSE)
}

.normalize_deg_file <- function(file_path, dataset, deg_method, lfc_col, padj_col, gene_col) {
  dataset_gene_col <- .resolve_dataset_value(gene_col, dataset, "gene_col")
  data <- .read_deg_file(file_path, dataset_gene_col)
  .normalize_deg_input(
    data,
    deg_method = .resolve_dataset_value(deg_method, dataset, "deg_method", allow_null = FALSE),
    lfc_col = .resolve_dataset_value(lfc_col, dataset, "lfc_col"),
    padj_col = .resolve_dataset_value(padj_col, dataset, "padj_col"),
    gene_col = dataset_gene_col
  )
}

.strip_ensembl_version <- function(gene_ids) {
  sub("^((?:ENS[A-Z]*G)[0-9]+)\\.[0-9]+$", "\\1", as.character(gene_ids), perl = TRUE)
}

.mapping_confidence_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x == 1)
  tolower(as.character(x)) %in% c("1", "true", "yes")
}

.audit_ortholog_mapping <- function(ortholog_map) {
  required <- c(
    "human_gene_id", "mouse_gene_id",
    "human_gene_id_version", "mouse_gene_id_version",
    "homology_type", "is_high_confidence", "release"
  )
  if (!is.data.frame(ortholog_map)) {
    stop("'ortholog_map' must be a data frame.", call. = FALSE)
  }
  missing_columns <- setdiff(required, colnames(ortholog_map))
  if (length(missing_columns) > 0) {
    stop("'ortholog_map' is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  mapping <- as.data.frame(ortholog_map, stringsAsFactors = FALSE)
  mapping$human_gene_id <- as.character(mapping$human_gene_id)
  mapping$mouse_gene_id <- as.character(mapping$mouse_gene_id)
  mapping$human_gene_id_version <- as.character(mapping$human_gene_id_version)
  mapping$mouse_gene_id_version <- as.character(mapping$mouse_gene_id_version)

  complete <- stats::complete.cases(mapping[, required]) &
    apply(mapping[, required, drop = FALSE], 1, function(row) all(nzchar(as.character(row))))
  complete_mapping <- mapping[complete, , drop = FALSE]
  inconsistent <-
    .strip_ensembl_version(complete_mapping$human_gene_id_version) != complete_mapping$human_gene_id |
    .strip_ensembl_version(complete_mapping$mouse_gene_id_version) != complete_mapping$mouse_gene_id
  if (any(inconsistent)) {
    stop("versioned ortholog IDs do not match their stable IDs.", call. = FALSE)
  }

  audit <- data.frame(
    step = "input",
    input_rows = nrow(mapping),
    removed_rows = 0L,
    retained_rows = nrow(mapping),
    flagged_rows = 0L,
    stringsAsFactors = FALSE
  )
  add_step <- function(step, before, after) {
    data.frame(
      step = step,
      input_rows = nrow(before),
      removed_rows = nrow(before) - nrow(after),
      retained_rows = nrow(after),
      flagged_rows = 0L,
      stringsAsFactors = FALSE
    )
  }
  add_symbol_audit <- function(step, data, column) {
    flagged <- 0L
    if (column %in% colnames(data)) {
      symbols <- as.character(data[[column]])
      valid <- !is.na(symbols) & nzchar(symbols)
      flagged <- sum(valid & (duplicated(symbols) | duplicated(symbols, fromLast = TRUE)))
    }
    data.frame(
      step = step,
      input_rows = nrow(data),
      removed_rows = 0L,
      retained_rows = nrow(data),
      flagged_rows = flagged,
      stringsAsFactors = FALSE
    )
  }

  audit <- rbind(audit, add_step("complete_identifiers", mapping, complete_mapping))
  audit <- rbind(
    audit,
    add_symbol_audit("duplicate_human_symbol_labels", complete_mapping, "human_gene_symbol"),
    add_symbol_audit("duplicate_mouse_symbol_labels", complete_mapping, "mouse_gene_symbol")
  )
  without_one_to_many <- complete_mapping[
    !complete_mapping$homology_type %in% c("ortholog_one2many", "ortholog_many2one"),
    ,
    drop = FALSE
  ]
  audit <- rbind(audit, add_step("exclude_one_to_many", complete_mapping, without_one_to_many))
  without_many_to_many <- without_one_to_many[
    without_one_to_many$homology_type != "ortholog_many2many",
    ,
    drop = FALSE
  ]
  audit <- rbind(audit, add_step("exclude_many_to_many", without_one_to_many, without_many_to_many))
  one_to_one <- without_many_to_many[
    without_many_to_many$homology_type == "ortholog_one2one",
    ,
    drop = FALSE
  ]
  audit <- rbind(audit, add_step("exclude_other_homology", without_many_to_many, one_to_one))
  high_confidence <- one_to_one[.mapping_confidence_flag(one_to_one$is_high_confidence), , drop = FALSE]
  audit <- rbind(audit, add_step("exclude_low_confidence", one_to_one, high_confidence))

  ambiguous <- duplicated(high_confidence$human_gene_id) |
    duplicated(high_confidence$human_gene_id, fromLast = TRUE) |
    duplicated(high_confidence$mouse_gene_id) |
    duplicated(high_confidence$mouse_gene_id, fromLast = TRUE)
  unique_mapping <- high_confidence[!ambiguous, , drop = FALSE]
  unique_mapping <- unique_mapping[order(unique_mapping$human_gene_id, unique_mapping$mouse_gene_id), , drop = FALSE]
  rownames(unique_mapping) <- NULL
  audit <- rbind(audit, add_step("exclude_duplicate_ids", high_confidence, unique_mapping))
  rownames(audit) <- NULL

  list(mapping = unique_mapping, audit = audit)
}
