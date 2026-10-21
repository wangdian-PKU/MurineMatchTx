test_that("plot_TME_barplot rejects invalid count inputs before conversion", {
  continuous_matrix <- matrix(
    c(10.25, 20.5, 30.75, 40.125),
    nrow = 2,
    dimnames = list(c("GENE1", "GENE2"), c("sample1", "sample2"))
  )
  negative_matrix <- matrix(
    c(10, -1, 30, 40),
    nrow = 2,
    dimnames = list(c("GENE1", "GENE2"), c("sample1", "sample2"))
  )

  expect_error(
    plot_TME_barplot(
      continuous_matrix,
      input_type = "count",
      input_source = "raw_count",
      method = "cibersort",
      group_as = c("A", "B")
    ),
    "whole numbers"
  )
  expect_error(
    plot_TME_barplot(
      negative_matrix,
      input_type = "count",
      input_source = "raw_count",
      method = "cibersort",
      group_as = c("A", "B")
    ),
    "non-negative"
  )
})

test_that("plot_TME_barplot rejects explicitly labelled ComBat input", {
  count_matrix <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("GENE1", "GENE2"), c("sample1", "sample2"))
  )

  expect_error(
    plot_TME_barplot(
      count_matrix,
      input_type = "count",
      input_source = "combat_corrected",
      method = "cibersort",
      group_as = c("A", "B")
    ),
    "ComBat"
  )
})

test_that("expression input validation requires source for external matrices", {
  count_matrix <- matrix(c(10, 20, 30, 40), nrow = 2)
  tpm_matrix <- matrix(c(0, 12.5, 30, 45.5), nrow = 2)

  expect_error(
    MurineMatchTx:::.validate_expression_input(count_matrix, "count"),
    "explicit 'input_source'"
  )
  expect_error(
    MurineMatchTx:::.validate_expression_input(tpm_matrix, "TPM"),
    "explicit 'input_source'"
  )
  expect_identical(
    MurineMatchTx:::.validate_expression_input(
      count_matrix,
      "count",
      input_source = "raw_count"
    )$input_source,
    "raw_count"
  )
  expect_identical(
    MurineMatchTx:::.validate_expression_input(
      tpm_matrix,
      "TPM",
      input_source = "provided_tpm"
    )$input_source,
    "provided_tpm"
  )
  expect_identical(
    MurineMatchTx:::.validate_expression_input(
      tpm_matrix,
      "TPM",
      input_source = "count2tpm"
    )$input_source,
    "count2tpm"
  )
  expect_error(
    MurineMatchTx:::.validate_expression_input(
      count_matrix,
      "count",
      input_source = "combat_corrected"
    ),
    "ComBat"
  )
  expect_error(
    MurineMatchTx:::.validate_expression_input(
      -tpm_matrix,
      "TPM",
      input_source = "provided_tpm"
    ),
    "non-negative"
  )
})

test_that("calculate_metabolism_score rejects negative TPM before scoring", {
  eset_tpm <- data.frame(
    s1 = c(1, -2, 3),
    s2 = c(4, 5, 6),
    row.names = c("HK1", "LDHA", "PFKM")
  )

  expect_error(
    calculate_metabolism_score(
      eset_tpm,
      method = "zscore",
      mini_gene_count = 1,
      input_source = "provided_tpm"
    ),
    "non-negative"
  )
})

test_that("calculate_metabolism_score rejects explicitly labelled ComBat TPM", {
  eset_tpm <- data.frame(
    s1 = c(1, 2, 3),
    s2 = c(4, 5, 6),
    row.names = c("HK1", "LDHA", "PFKM")
  )

  expect_error(
    calculate_metabolism_score(
      eset_tpm,
      method = "zscore",
      mini_gene_count = 1,
      input_source = "combat_corrected"
    ),
    "ComBat"
  )
})

test_that("calculate_metabolism_score accepts TPM matrices from the TME workflow", {
  eset_tpm <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8,
      2, 3, 4, 5, 6, 7, 8, 9),
    ncol = 2,
    dimnames = list(
      c("HK1", "LDHA", "PFKM", "ENO1", "NDUFA1", "COX4I1", "ATP5F1A", "ATP5F1B"),
      c("s1", "s2")
    )
  )
  attr(eset_tpm, "MurineMatchTx_tpm_source") <- "count2tpm"

  sig_meta <- calculate_metabolism_score(eset_tpm, method = "zscore", mini_gene_count = 1)

  expect_s3_class(sig_meta, "data.frame")
  expect_identical(attr(sig_meta, "MurineMatchTx_tpm_source"), "count2tpm")
})

test_that("input source extension preserves existing positional arguments", {
  expect_identical(
    names(formals(plot_TME_barplot))[1:4],
    c("data_matrix", "input_type", "method", "group_as")
  )
  expect_identical(
    names(formals(calculate_metabolism_score))[1:3],
    c("eset_tpm", "method", "mini_gene_count")
  )
})
