make_mapping_fixture <- function() {
  data.frame(
    human_gene_id = c("ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004", "ENSG000004", "ENSG000006"),
    mouse_gene_id = c("ENSMUSG000001", "ENSMUSG000002", "ENSMUSG000003", "ENSMUSG000004", "ENSMUSG000005", "ENSMUSG000006"),
    human_gene_id_version = c("ENSG000001.1", "ENSG000002.1", "ENSG000003.1", "ENSG000004.1", "ENSG000004.1", "ENSG000006.2"),
    mouse_gene_id_version = c("ENSMUSG000001.1", "ENSMUSG000002.1", "ENSMUSG000003.1", "ENSMUSG000004.1", "ENSMUSG000005.1", "ENSMUSG000006.2"),
    homology_type = c("ortholog_one2one", "ortholog_one2many", "ortholog_one2one", "ortholog_one2one", "ortholog_one2one", "ortholog_one2one"),
    is_high_confidence = c(1, 1, 0, 1, 1, 1),
    release = rep(116L, 6),
    human_gene_symbol = c("shared", "B", "C", "D", "D", "shared"),
    mouse_gene_symbol = c("a", "b", "c", "d", "e", "f"),
    stringsAsFactors = FALSE
  )
}

test_that("ortholog audit conservatively retains high-confidence unique one-to-one pairs", {
  mapping <- make_mapping_fixture()
  audit_a <- MurineMatchTx:::.audit_ortholog_mapping(mapping)
  audit_b <- MurineMatchTx:::.audit_ortholog_mapping(mapping[c(6, 4, 2, 5, 1, 3), ])

  expect_identical(
    audit_a$mapping[, c("human_gene_id", "mouse_gene_id")],
    audit_b$mapping[, c("human_gene_id", "mouse_gene_id")]
  )
  expect_identical(audit_a$mapping$human_gene_id, c("ENSG000001", "ENSG000006"))
  expect_true(all(c("step", "input_rows", "removed_rows", "retained_rows") %in% colnames(audit_a$audit)))
  expect_equal(audit_a$audit$removed_rows[audit_a$audit$step == "exclude_one_to_many"], 1)
  expect_equal(audit_a$audit$removed_rows[audit_a$audit$step == "exclude_low_confidence"], 1)
  expect_equal(audit_a$audit$removed_rows[audit_a$audit$step == "exclude_duplicate_ids"], 2)
  expect_equal(audit_a$audit$flagged_rows[audit_a$audit$step == "duplicate_human_symbol_labels"], 4)
  expect_equal(tail(audit_a$audit$retained_rows, 1), 2)
})

test_that("ortholog audit rejects inconsistent stable and versioned IDs", {
  mapping <- make_mapping_fixture()[1, ]
  mapping$human_gene_id_version <- "ENSG999999.1"

  expect_error(
    MurineMatchTx:::.audit_ortholog_mapping(mapping),
    "versioned.*stable"
  )
})

test_that("DEG heatmap preparation uses mapped mouse IDs and shared thresholds", {
  human <- data.frame(
    logFC = c(1, -1, 2),
    FDR = c(0.01, 0.01, 0.05),
    row.names = c("ENSG000001.1", "ENSG000006.2", "ENSG000007.1")
  )
  mouse <- data.frame(
    logFC = c(1, -2, 5),
    FDR = c(0.01, 0.01, 0.01),
    row.names = c("ENSMUSG000001.1", "ENSMUSG000006.2", "ENSMUSG000007.1")
  )
  human_file <- write_deg_fixture(human)
  mouse_file <- write_deg_fixture(mouse)
  mapping <- make_mapping_fixture()[c(1, 6), ]

  result <- prepare_DEG_heatmap_data(
    human_file,
    list(ModelA = mouse_file),
    deg_method = "edgeR",
    deg_thresholds = list(log2fc = 1, fdr = 0.05),
    ortholog_map = mapping
  )

  expect_identical(result$tcga_gene_list, c("ENSG000001", "ENSG000006"))
  expect_identical(result$processed_mouse$ModelA[[1]], c(1, -1))
  expect_equal(tail(result$mapping_audit$retained_rows, 1), 2)
  expect_identical(result$deg_thresholds, list(log2fc = 1, fdr = 0.05))
})

test_that("Compara mapping is required by default and HomoloGene is explicit compatibility input", {
  human <- write_deg_fixture(data.frame(logFC = 1.5, FDR = 0.01, row.names = "TP53"))
  mouse <- write_deg_fixture(data.frame(logFC = 1.5, FDR = 0.01, row.names = "Trp53"))

  expect_error(
    prepare_DEG_heatmap_data(human, list(ModelA = mouse), deg_method = "edgeR"),
    "ortholog_map"
  )

  expect_warning(
    compatibility <- prepare_DEG_heatmap_data(
      human,
      list(ModelA = mouse),
      deg_method = "edgeR",
      mapping_source = "homologene"
    ),
    "sensitivity|compatibility"
  )
  expect_identical(compatibility$tcga_gene_list, "TP53")
})
