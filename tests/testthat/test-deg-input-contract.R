make_method_deg <- function(method) {
  values <- list(
    gene = c("g_up", "g_down", "g_fdr_boundary", "g_stable"),
    lfc = c(1, -1, 2, 0.5),
    padj = c(0.01, 0.01, 0.05, 0.01)
  )
  if (method == "edgeR") {
    return(data.frame(logFC = values$lfc, FDR = values$padj, row.names = values$gene))
  }
  if (method == "DESeq2") {
    return(data.frame(log2FoldChange = values$lfc, padj = values$padj, row.names = values$gene))
  }
  data.frame(logFC = values$lfc, adj.P.Val = values$padj, row.names = values$gene)
}

test_that("edgeR, DESeq2, and limma fixtures normalize identically", {
  normalized <- lapply(c("edgeR", "DESeq2", "limma"), function(method) {
    human <- write_deg_fixture(make_method_deg(method))
    mouse <- write_deg_fixture(make_method_deg(method))
    result <- merge_DEG_datasets(
      human,
      list(ModelA = mouse),
      deg_method = method,
      deg_thresholds = list(log2fc = 1, fdr = 0.05)
    )
    result$contrast <- as.character(result$contrast)
    result[, c("gene", "logFC", "FDR", "contrast", "change")]
  })

  expect_equal(normalized[[1]], normalized[[2]])
  expect_equal(normalized[[1]], normalized[[3]])
  expect_identical(
    normalized[[1]]$change[normalized[[1]]$contrast == "TCGA"],
    c("Up", "Down", "Stable", "Stable")
  )
})

test_that("DESeq2 normalization ignores column order and unrelated columns", {
  input <- data.frame(
    padj = c(0.01, 0.02),
    unrelated = c("a", "b"),
    log2FoldChange = c(1.5, -1.5),
    row.names = c("g1", "g2")
  )

  expect_equal(
    MurineMatchTx:::.normalize_deg_input(input, deg_method = "DESeq2"),
    data.frame(gene = c("g1", "g2"), logFC = c(1.5, -1.5), FDR = c(0.01, 0.02))
  )
})

test_that("DEG heatmap selects the named change column", {
  thresholds <- list(log2fc = 1, fdr = 0.05)
  deg_data <- list(
    processed_tcga = data.frame(logFC = c(1, -1), FDR = c(0.01, 0.02)),
    processed_mouse = list(
      ModelA = data.frame(FDR = c(0.01, 0.02), change = c(1, -1))
    ),
    tcga_gene_list = c("g1", "g2"),
    deg_thresholds = thresholds
  )
  output_path <- paste0(tempfile("deg-heatmap-"), "/")
  dir.create(output_path)

  result <- plot_DEG_heatmap(
    deg_data,
    mouse_files = list(ModelA = "unused"),
    color_palette = c("1" = "#ff7676", "-1" = "#66d4ff", "0.01" = "grey", "0.02" = "grey", "NA" = "white"),
    output_path = output_path
  )

  expect_identical(as.numeric(result[, "ModelA"]), c(1, -1))
})

test_that("custom DEG columns require explicit mappings", {
  custom <- data.frame(
    symbol = c("g1", "g2"),
    effect = c(1.5, -1.5),
    q_value = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  expect_error(
    MurineMatchTx:::.normalize_deg_input(custom, deg_method = "custom"),
    "lfc_col.*padj_col"
  )
  expect_equal(
    MurineMatchTx:::.normalize_deg_input(
      custom,
      deg_method = "custom",
      lfc_col = "effect",
      padj_col = "q_value",
      gene_col = "symbol"
    ),
    data.frame(gene = c("g1", "g2"), logFC = c(1.5, -1.5), FDR = c(0.01, 0.02))
  )
})

test_that("DEG normalization rejects duplicates, missing values, and non-numeric columns", {
  expect_error(
    MurineMatchTx:::.normalize_deg_input(
      data.frame(gene = c("g1", "g1"), logFC = c(1, 2), FDR = c(0.01, 0.02)),
      "edgeR",
      gene_col = "gene"
    ),
    "duplicate.*g1"
  )
  expect_error(
    MurineMatchTx:::.normalize_deg_input(
      data.frame(gene = "g1", logFC = NA_real_, FDR = 0.01),
      "edgeR",
      gene_col = "gene"
    ),
    "missing"
  )
  expect_error(
    MurineMatchTx:::.normalize_deg_input(
      data.frame(gene = "g1", logFC = "high", FDR = 0.01),
      "edgeR",
      gene_col = "gene"
    ),
    "numeric"
  )

  empty <- MurineMatchTx:::.normalize_deg_input(
    data.frame(gene = character(), logFC = numeric(), FDR = numeric()),
    "edgeR",
    gene_col = "gene"
  )
  expect_equal(nrow(empty), 0)
  expect_identical(colnames(empty), c("gene", "logFC", "FDR"))
})

test_that("DEG threshold object applies inclusive fold-change and exclusive FDR boundaries", {
  deg <- data.frame(
    gene = c("up", "down", "fdr_boundary", "custom_only"),
    logFC = c(1, -1, 2, 0.75),
    FDR = c(0.049, 0.049, 0.05, 0.049)
  )
  default_change <- MurineMatchTx:::.classify_deg(
    deg,
    list(log2fc = 1, fdr = 0.05)
  )
  custom_change <- MurineMatchTx:::.classify_deg(
    deg,
    list(log2fc = 0.5, fdr = 0.1)
  )

  expect_identical(default_change, c("Up", "Down", "Stable", "Stable"))
  expect_identical(custom_change, c("Up", "Down", "Up", "Up"))
  expect_true(all(c("deg_method", "deg_thresholds") %in% names(formals(GO_enrichment_analysis))))
  expect_true("deg_thresholds" %in% names(formals(plot_DEG_heatmap)))
  expect_true("deg_thresholds" %in% names(formals(plot_DEG_barplot)))
  expect_true("deg_thresholds" %in% names(formals(plot_DEG_volcano)))
})
