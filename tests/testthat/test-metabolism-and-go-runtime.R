test_that("calculate_metabolism_score uses bundled metabolism signatures", {
  eset_tpm <- data.frame(
    s1 = c(1, 2, 3, 4, 5, 6, 7, 8),
    s2 = c(2, 3, 4, 5, 6, 7, 8, 9),
    row.names = c("HK1", "LDHA", "PFKM", "ENO1", "NDUFA1", "COX4I1", "ATP5F1A", "ATP5F1B")
  )

  sig_meta <- calculate_metabolism_score(eset_tpm, method = "zscore", mini_gene_count = 1)

  expect_s3_class(sig_meta, "data.frame")
  expect_true("ID" %in% colnames(sig_meta))
  expect_true(any(colnames(sig_meta) != "ID"))
})

test_that("plot_metabolism_heatmap requires explicit grouping and writes output", {
  sig_meta <- data.frame(
    ID = c("s1", "s2", "s3", "s4"),
    Glycolysis = c(1, 2, 3, 4),
    OXPHOS = c(4, 3, 2, 1)
  )

  expect_error(
    plot_metabolism_heatmap(sig_meta, group_as = c("A", "A", "B"), output_path = tempdir()),
    "same length"
  )

  output_dir <- tempfile(pattern = "metabolism-plot-")
  heatmap_obj <- plot_metabolism_heatmap(
    sig_meta,
    group_as = c("A", "A", "B", "B"),
    output_path = output_dir,
    height = 6
  )

  expect_true(file.exists(file.path(output_dir, "heatmap.pdf")))
  expect_s4_class(heatmap_obj, "Heatmap")
})

test_that("GO symbol map lookup works for supported species", {
  human_map <- MurineMatchTx:::.get_symbol_map("hsa")
  mouse_map <- MurineMatchTx:::.get_symbol_map("mmu")

  expect_s3_class(human_map, "data.frame")
  expect_s3_class(mouse_map, "data.frame")
  expect_true(all(c("gene_id", "symbol") %in% colnames(human_map)))
  expect_true(all(c("gene_id", "symbol") %in% colnames(mouse_map)))
})
