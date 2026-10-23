#' Prepare DEG Data for Heatmap Visualization
#'
#' Normalize DEG inputs, apply one shared threshold object, and align human and
#' mouse genes through an audited ortholog mapping.
#'
#' @param tcga_file Path to the human DEG result file.
#' @param mouse_files Named list of mouse DEG result files.
#' @param logFC,FDR Deprecated-compatible thresholds used only when
#'   `deg_thresholds` is `NULL`.
#' @param inTax,outTax Taxonomy IDs used only for the explicit HomoloGene
#'   compatibility path.
#' @param deg_method One of `"edgeR"`, `"DESeq2"`, `"limma"`, or `"custom"`;
#'   may be a named vector by dataset.
#' @param lfc_col,padj_col,gene_col Optional explicit column names, scalar or
#'   named by dataset.
#' @param deg_thresholds A list with `log2fc` and `fdr`.
#' @param ortholog_map Standardized Compara mapping with stable and versioned
#'   human/mouse IDs, `homology_type`, `is_high_confidence`, and `release`.
#' @param mapping_source `"compara"` for the primary audited mapping or
#'   `"homologene"` for an explicit compatibility/sensitivity analysis.
#'
#' @return A list containing normalized human data, aligned mouse data, the
#'   human gene list, mapping audit, per-dataset gene-loss audit, retained
#'   mapping, and thresholds.
#'
#' @export
prepare_DEG_heatmap_data <- function(tcga_file,
                                     mouse_files,
                                     logFC = 1,
                                     FDR = 0.05,
                                     inTax = 9606,
                                     outTax = 10090,
                                     deg_method = "edgeR",
                                     lfc_col = NULL,
                                     padj_col = NULL,
                                     gene_col = NULL,
                                     deg_thresholds = NULL,
                                     ortholog_map = NULL,
                                     mapping_source = c("compara", "homologene")) {
  if (!is.list(mouse_files) || is.null(names(mouse_files)) || any(!nzchar(names(mouse_files)))) {
    stop("'mouse_files' must be a named list.", call. = FALSE)
  }
  mapping_source <- match.arg(mapping_source)
  thresholds <- .validate_deg_thresholds(deg_thresholds, logFC, FDR)

  human_all <- .normalize_deg_file(
    tcga_file, "TCGA", deg_method, lfc_col, padj_col, gene_col
  )
  mouse <- lapply(names(mouse_files), function(dataset) {
    .normalize_deg_file(
      mouse_files[[dataset]], dataset, deg_method, lfc_col, padj_col, gene_col
    )
  })
  names(mouse) <- names(mouse_files)

  human_all$change <- .classify_deg(human_all, thresholds)
  human <- human_all[human_all$change != "Stable", , drop = FALSE]
  human_significant_count <- nrow(human)
  human$stable_gene_id <- .strip_ensembl_version(human$gene)
  duplicate_human <- unique(human$stable_gene_id[duplicated(human$stable_gene_id)])
  if (length(duplicate_human) > 0) {
    stop(
      "Human DEG IDs are duplicate after removing Ensembl versions: ",
      paste(utils::head(duplicate_human, 5), collapse = ", "),
      call. = FALSE
    )
  }

  if (mapping_source == "compara") {
    if (is.null(ortholog_map)) {
      stop("Primary analysis requires an explicit 'ortholog_map'.", call. = FALSE)
    }
  } else {
    warning(
      "HomoloGene is an explicit sensitivity/compatibility mapping, not the primary high-confidence mapping.",
      call. = FALSE
    )
    homologs <- homologene::homologene(human$stable_gene_id, inTax = inTax, outTax = outTax)
    if (nrow(homologs) == 0) {
      ortholog_map <- data.frame(
        human_gene_id = character(), mouse_gene_id = character(),
        human_gene_id_version = character(), mouse_gene_id_version = character(),
        homology_type = character(), is_high_confidence = logical(), release = character()
      )
    } else {
      human_column <- as.character(inTax)
      mouse_column <- as.character(outTax)
      ortholog_map <- data.frame(
        human_gene_id = as.character(homologs[[human_column]]),
        mouse_gene_id = as.character(homologs[[mouse_column]]),
        human_gene_id_version = as.character(homologs[[human_column]]),
        mouse_gene_id_version = as.character(homologs[[mouse_column]]),
        homology_type = "ortholog_one2one",
        is_high_confidence = TRUE,
        release = "HomoloGene build 68",
        stringsAsFactors = FALSE
      )
    }
  }

  mapping_result <- .audit_ortholog_mapping(ortholog_map)
  retained_mapping <- mapping_result$mapping
  human <- human[order(-human$logFC, human$stable_gene_id), , drop = FALSE]
  mapping_index <- match(human$stable_gene_id, retained_mapping$human_gene_id)
  keep <- !is.na(mapping_index)
  human <- human[keep, , drop = FALSE]
  mapping_index <- mapping_index[keep]
  human$mouse_gene_id <- retained_mapping$mouse_gene_id[mapping_index]

  processed_tcga <- human[, c("logFC", "FDR"), drop = FALSE]
  rownames(processed_tcga) <- human$stable_gene_id

  processed_mouse_result <- lapply(names(mouse), function(dataset) {
    mouse_deg <- mouse[[dataset]]
    mouse_deg$stable_gene_id <- .strip_ensembl_version(mouse_deg$gene)
    duplicate_mouse <- unique(mouse_deg$stable_gene_id[duplicated(mouse_deg$stable_gene_id)])
    if (length(duplicate_mouse) > 0) {
      stop(
        "Mouse DEG IDs are duplicate after removing Ensembl versions: ",
        paste(utils::head(duplicate_mouse, 5), collapse = ", "),
        call. = FALSE
      )
    }
    mouse_index <- match(human$mouse_gene_id, mouse_deg$stable_gene_id)
    mouse_change <- rep(NA_character_, nrow(human))
    present <- !is.na(mouse_index)
    mouse_change[present] <- .classify_deg(mouse_deg[mouse_index[present], , drop = FALSE], thresholds)
    concordant <- ifelse(
      human$change == "Up" & mouse_change == "Up",
      1,
      ifelse(human$change == "Down" & mouse_change == "Down", -1, NA_real_)
    )
    adjusted_p <- rep(NA_real_, nrow(human))
    adjusted_p[present] <- mouse_deg$FDR[mouse_index[present]]
    output <- data.frame(change = concordant, FDR = adjusted_p)
    rownames(output) <- human$stable_gene_id
    mapped_input <- mouse_deg$stable_gene_id %in% retained_mapping$mouse_gene_id
    list(
      data = output,
      audit = data.frame(
        dataset = dataset,
        step = c("normalized_input", "primary_mapping", "selected_human_orthologs_present"),
        input_genes = c(nrow(mouse_deg), nrow(mouse_deg), nrow(human)),
        removed_genes = c(0L, sum(!mapped_input), sum(!present)),
        retained_genes = c(nrow(mouse_deg), sum(mapped_input), sum(present)),
        stringsAsFactors = FALSE
      )
    )
  })
  names(processed_mouse_result) <- names(mouse)
  processed_mouse <- lapply(processed_mouse_result, `[[`, "data")

  gene_loss_audit <- rbind(
    data.frame(
      dataset = "TCGA",
      step = c("normalized_input", "significant_genes", "primary_mapping"),
      input_genes = c(nrow(human_all), nrow(human_all), human_significant_count),
      removed_genes = c(
        0L,
        nrow(human_all) - human_significant_count,
        human_significant_count - nrow(human)
      ),
      retained_genes = c(nrow(human_all), human_significant_count, nrow(human)),
      stringsAsFactors = FALSE
    ),
    do.call(rbind, lapply(processed_mouse_result, `[[`, "audit"))
  )
  rownames(gene_loss_audit) <- NULL

  mapping_audit <- rbind(
    mapping_result$audit,
    data.frame(
      step = "human_significant_genes_mapped",
      input_rows = sum(human_all$change != "Stable"),
      removed_rows = sum(human_all$change != "Stable") - nrow(human),
      retained_rows = nrow(human),
      flagged_rows = 0L,
      stringsAsFactors = FALSE
    )
  )

  list(
    processed_tcga = processed_tcga,
    processed_mouse = processed_mouse,
    tcga_gene_list = human$stable_gene_id,
    mapping_audit = mapping_audit,
    gene_loss_audit = gene_loss_audit,
    ortholog_map = retained_mapping,
    deg_thresholds = thresholds
  )
}
