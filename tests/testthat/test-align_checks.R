# tests/testthat/test-align_check.R

test_that("align_data_and_metadata groups by condition and natural-sorts within group", {
  # data columns intentionally out of natural order
  dat <- data.frame(
    A1  = c(1, 2, 3),
    A10 = c(4, 5, 6),
    A2  = c(7, 8, 9),
    B1  = c(10, 11, 12)
  )
  rownames(dat) <- paste0("Prot", 1:3)
  
  meta <- data.frame(
    sample = c("A1", "A2", "A10", "B1"),
    group  = c("A",  "A",  "A",   "B"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample
  
  out <- align_data_and_metadata(
    data = dat,
    metadata = meta,
    sample_col = "sample",
    group_col = "group",
    prefer_group_blocks = TRUE,
    strict = FALSE
  )
  
  expect_true(is.list(out))
  expect_true(all(c("data", "metadata", "changes") %in% names(out)))
  
  # Expect groups ordered A then B; within A: A1, A2, A10 (natural order), then B1
  expect_identical(colnames(out$data), c("A1", "A2", "A10", "B1"))
  expect_identical(rownames(out$metadata), c("A1", "A2", "A10", "B1"))
  
  # changes log shape
  expect_s3_class(out$changes, "data.frame")
  expect_true(all(c("position", "sample", "old_pos_in_metadata", "old_pos_in_data") %in% names(out$changes)))
})

test_that("align_data_and_metadata preserves data order when prefer_group_blocks=FALSE (metadata reorders only)", {
  dat <- data.frame(
    S1  = 1:3,
    S10 = 4:6,
    S2  = 7:9
  )
  rownames(dat) <- paste0("Prot", 1:3)
  
  meta <- data.frame(
    sample = c("S1", "S2", "S10"),
    group  = c("A",  "A",  "A"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample
  
  out <- align_data_and_metadata(
    data = dat,
    metadata = meta,
    sample_col = "sample",
    group_col = "group",
    prefer_group_blocks = FALSE,   # do not regroup; just match data's order
    strict = FALSE
  )
  
  # Data columns unchanged
  expect_identical(colnames(out$data), colnames(dat))
  # Metadata reordered to match data columns exactly
  expect_identical(rownames(out$metadata), colnames(dat))
})

test_that("align_data_and_metadata errors when strict=TRUE and order differs", {
  dat <- data.frame(
    X1 = 1:3,
    X2 = 4:6
  )
  rownames(dat) <- paste0("P", 1:3)
  
  meta <- data.frame(
    sample = c("X2", "X1"),
    group  = c("G",  "G"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample
  
  expect_error(
    align_data_and_metadata(
      data = dat,
      metadata = meta,
      sample_col = "sample",
      group_col = "group",
      prefer_group_blocks = FALSE,
      strict = TRUE         # should error due to order mismatch
    ),
    "Order mismatch detected"
  )
})

test_that("align_data_and_metadata errors on mismatched sample sets", {
  dat <- data.frame(
    S1 = 1:3,
    S2 = 4:6
  )
  rownames(dat) <- paste0("P", 1:3)
  
  meta <- data.frame(
    sample = c("S1", "S3"),  # S3 not in data; S2 missing
    group  = c("A", "A"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample
  
  expect_error(
    align_data_and_metadata(
      data = dat,
      metadata = meta,
      sample_col = "sample",
      group_col = "group",
      prefer_group_blocks = TRUE,
      strict = FALSE
    ),
    "Sample ID mismatch"
  )
})

test_that("align_data_and_metadata errors on duplicate sample IDs (data or metadata)", {
  # Duplicate in data
  dat_dup <- data.frame(
    D1 = 1:3,
    D1 = 4:6,   # duplicate name in columns
    check.names = FALSE
  )
  rownames(dat_dup) <- paste0("P", 1:3)
  
  meta <- data.frame(
    sample = c("D1", "D2"),
    group  = c("G",  "G"),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample
  
  expect_error(
    align_data_and_metadata(
      data = dat_dup,
      metadata = meta,
      sample_col = "sample"
    ),
    "Duplicate sample IDs in `data`"
  )
  
  # Duplicate in metadata
  dat <- data.frame(
    D1 = 1:3,
    D2 = 4:6
  )
  rownames(dat) <- paste0("P", 1:3)
  
  meta_dup <- data.frame(
    sample = c("D1", "D1"),   # duplicate in metadata
    group  = c("G",  "G"),
    stringsAsFactors = FALSE
  )
  rownames(meta_dup) <- meta_dup$sample
  
  expect_error(
    align_data_and_metadata(
      data = dat,
      metadata = meta_dup,
      sample_col = "sample"
    ),
    "Duplicate sample IDs in `metadata`"
  )
})

test_that("sample_col is respected; rownames not required when sample_col provided", {
  dat <- data.frame(
    A1 = 1:3,
    A2 = 4:6
  )
  rownames(dat) <- paste0("Prot", 1:3)
  
  meta <- data.frame(
    sample = c("A2", "A1"),
    group  = c("X",  "X"),
    stringsAsFactors = FALSE
  )
  # intentionally no rownames(meta)
  
  out <- align_data_and_metadata(
    data = dat,
    metadata = meta,
    sample_col = "sample",
    group_col = "group",
    prefer_group_blocks = FALSE
  )
  
  expect_identical(colnames(out$data), c("A1", "A2"))
  expect_identical(rownames(out$metadata), c("A1", "A2"))
})

test_that("check_data_metadata_match passes when aligned and fails when not", {
  dat <- data.frame(S1 = 1:3, S2 = 4:6)
  rownames(dat) <- paste0("P", 1:3)
  
  meta_ok <- data.frame(group = c("A", "B"), row.names = c("S1", "S2"))
  expect_invisible(check_data_metadata_match(dat, meta_ok))
  
  meta_bad <- data.frame(group = c("A", "B"), row.names = c("S1", "S3"))
  expect_error(check_data_metadata_match(dat, meta_bad), "Sample ID mismatch")
})
