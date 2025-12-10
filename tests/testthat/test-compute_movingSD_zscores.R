# Key changes vs your original test:
#   We now build a full DAList-like object with slots: data, annotation, metadata, design, eBayes_fit, results, tags, and set class(raw) <- c("DAList", "list").
#   We include tags$contrast_vector, which should satisfy the “Cannot determine contrast names…” logic if it’s still used internally.
#   We keep your original assertions about movingSDs and logFC_z_scores.
#
# No extra arguments to compute_movingSD_zscores() other than DAList and binsize.
# We set raw$filtered_proteins_per_contrast so the function can grab contrast_names from there
# devtools::test(filter = "compute_movingSD_zscores")

test_that("compute_moving_sd_and_zscores works with valid input", {
  set.seed(123)
  n_proteins <- 200
  n_samples  <- 5
  
  # Simulate expression data
  expr <- matrix(
    rnorm(n_proteins * n_samples),
    nrow = n_proteins,
    ncol = n_samples,
    dimnames = list(
      paste0("P", seq_len(n_proteins)),
      paste0("S", seq_len(n_samples))
    )
  )
  
  # Annotation: uniprot_id must be unique
  annotation <- data.frame(
    uniprot_id = rownames(expr),
    stringsAsFactors = FALSE
  )
  
  # Metadata with a simple two-group design
  metadata <- data.frame(
    sample = colnames(expr),
    group  = rep(c("control", "treat"), length.out = n_samples),
    stringsAsFactors = TRUE
  )
  rownames(metadata) <- metadata$sample
  
  # Design: simple ~ group model
  design_matrix <- model.matrix(~ group, data = metadata)
  design <- list(
    design_matrix  = design_matrix,
    design_formula = ~ group,
    random_factor  = NULL
  )
  
  # Simulate per-contrast logFC results
  contrast_names <- c("treat_vs_control", "drug_vs_control")
  results <- list(
    treat_vs_control = data.frame(
      logFC = rnorm(n_proteins),
      row.names = rownames(expr)
    ),
    drug_vs_control = data.frame(
      logFC = rnorm(n_proteins),
      row.names = rownames(expr)
    )
  )
  
  # Minimal eBayes_fit list with correct class
  eBayes_fit <- list(
    treat_vs_control = structure(list(), class = "MArrayLM"),
    drug_vs_control  = structure(list(), class = "MArrayLM")
  )
  
  # filtered_proteins_per_contrast: just give all proteins for each contrast
  filtered_proteins_per_contrast <- setNames(
    replicate(length(contrast_names),
              rownames(expr),
              simplify = FALSE),
    contrast_names
  )
  
  # Assemble minimal valid DAList
  raw <- list(
    data                        = expr,
    annotation                  = annotation,
    metadata                    = metadata,
    design                      = design,
    eBayes_fit                  = eBayes_fit,
    results                     = results,
    filtered_proteins_per_contrast = filtered_proteins_per_contrast,
    tags                        = list()
  )
  class(raw) <- c("DAList", "list")
  
  binsize <- 50
  
  raw <- compute_movingSD_zscores(
    DAList  = raw,
    binsize = binsize
  )
  
  # Check output structure
  expect_true("movingSDs" %in% names(raw$tags))
  expect_true("logFC_z_scores" %in% names(raw$tags))
  
  # Check each comparison has expected output
  for (comp in contrast_names) {
    expect_equal(length(raw$tags$movingSDs[[comp]]), n_proteins)
    expect_equal(length(raw$tags$logFC_z_scores[[comp]]), n_proteins)
    expect_false(any(is.na(raw$tags$movingSDs[[comp]])))
    expect_false(any(is.na(raw$tags$logFC_z_scores[[comp]])))
  }
})
