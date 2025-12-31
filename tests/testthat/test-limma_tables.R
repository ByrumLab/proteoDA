# A lot of these functions are tricky to test,
# the base code could probably use a re-factor.
# But, tried to test as much as possible.

test_that("make_excel_hyperlinks adds hyperlinks to col", {
  df <- data.frame(x = c("A",
                         "B",
                         "C",
                         NA))
  expect <- data.frame(x = c("HYPERLINK(\"test.com/A\", \"A\")",
                             "HYPERLINK(\"test.com/B\", \"B\")",
                             "HYPERLINK(\"test.com/C\", \"C\")",
                             NA))
  class(expect[,"x"])  <- "formula"

  expect_equal(make_excel_hyperlinks(df, url.col = 1, url = "test.com/"), expect)

})


test_that("write_limma_excel returns consistent output", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink("temp.xlsx"), add = T)

  suppressMessages(
      output <- write_limma_excel(filename = "temp.xlsx",
                        statlist = input$results,
                        annotation = input$annotation,
                        data = input$data,
                        norm.method = input$tags$norm_method,
                        pval_thresh = input$tags$DA_criteria$pval_thresh,
                        lfc_thresh = input$tags$DA_criteria$lfc_thresh,
                        add_filter = T)
  )

  # Track the the worksheet output, should be consistent
  expect_snapshot(output$wb$worksheets)

})

test_that("create_combined_results returns expected output", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  a <- input$results$control
  colnames(a) <- paste(colnames(a), "control", sep = "_")
  b <- input$results$treatment
  colnames(b) <- paste(colnames(b), "treatment", sep = "_")

  expected_output <- cbind(input$annotation,
                           input$data,
                           a,
                           b)

  expect_equal(create_combined_results(annotation = input$annotation,
                                       data = input$data,
                                       statlist = input$results),
              expected_output)
})


test_that("summarize_contrast_DA returns expected output", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))
  # Add in a mix of all significance types in the fake data
  input$results$control$sig.PVal <- c(rep(-1, 4),
                                      rep(0, 3),
                                      rep(1, 3))
  input$results$control$sig.FDR <- c(rep(-1, 3),
                                      rep(0, 4),
                                      rep(1, 3))
  input$results$treatment$sig.PVal <- c(rep(-1, 4),
                                        rep(0, 3),
                                        rep(1, 3))
  input$results$treatment$sig.FDR <- c(rep(-1, 3),
                                       rep(0, 4),
                                       rep(1, 2),
                                       NA)

  control_expect <- data.frame(contrast = "control",
                               type = c("down", "nonsig", "up"),
                               sig.PVal = as.character(c(4,3,3)),
                               sig.FDR = as.character(c(3,4,3)))
  treatment_expect <- data.frame(contrast = "treatment",
                               type = c("down", "nonsig", "up"),
                               sig.PVal = as.character(c(4,3,3)),
                               sig.FDR = as.character(c(3,4,2)))

  expect_equal(summarize_contrast_DA("control", input$results), control_expect)
  expect_equal(summarize_contrast_DA("treatment", input$results), treatment_expect)

})


test_that("write_limma_tables gives error when there are no results", {

  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))


  expect_error(write_limma_tables(input, overwrite = T), "does not have results")
})



test_that("write_limma_tables informs user of output dir and filenames when not supplied", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("per_contrast_results/control.csv",
                   "per_contrast_results/treatment.csv",
                   "per_contrast_results",
                   "DA_summary.csv",
                   "combined_results.csv",
                   "results.xlsx"), recursive = T), add = T)
  # A message for the directory and one for each file
  suppressMessages(expect_message(write_limma_tables(input, overwrite = T), "current working directory:"))
  suppressMessages(expect_message(write_limma_tables(input, overwrite = T), "DA_summary.csv"))
  suppressMessages(expect_message(write_limma_tables(input, overwrite = T), "per_contrast_results"))
  suppressMessages(expect_message(write_limma_tables(input, overwrite = T), "combined_results.csv"))
  suppressMessages(expect_message(write_limma_tables(input, overwrite = T), "results.xlsx"))

})



test_that("write_limma_tables checks validity of user-suppplied filenames", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("per_contrast_results/control.csv",
                   "per_contrast_results/treatment.csv",
                   "DA_summary.csv",
                   "combined_results.csv",
                   "results.xlsx"), recursive = T), add = T)

  suppressMessages(
    expect_error(
      write_limma_tables(input, overwrite = T, summary_csv = "xxx.xxx"),
      "Invalid file name"
    )
  )
  suppressMessages(
    expect_error(
      write_limma_tables(input, overwrite = T, combined_file_csv = "xxx.xxx"),
      "Invalid file name"
    )
  )
  suppressMessages(
    expect_error(
      write_limma_tables(input, overwrite = T, spreadsheet_xlsx = "xxx.xxx"),
    "Invalid file name"
    )
  )
})





test_that("write_limma_tables won't overwrite when false", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("per_contrast_results/control.csv",
                   "per_contrast_results/treatment.csv",
                   "DA_summary.csv",
                   "combined_results.csv",
                   "results.xlsx"), recursive = T), add = T)

  suppressMessages(write_limma_tables(input))
  suppressMessages(
    expect_error(
      write_limma_tables(input,
                         overwrite = F),
      "Change \\`output_dir\\`"
    )
  )
})


test_that("write_limma_tables messages user when overwriting", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("per_contrast_results/control.csv",
                   "per_contrast_results/treatment.csv",
                   "DA_summary.csv",
                   "combined_results.csv",
                   "results.xlsx",
                   "per_contrast_results"),
                 recursive = T), add = T)

  suppressMessages(write_limma_tables(input))
  suppressMessages(
    expect_message(
      write_limma_tables(input,
                         overwrite = T),
      "Overwriting"
    )
  )
})

test_that("write_limma_tables creates output files", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  expected_files_default <-  c("per_contrast_results/control.csv",
                               "per_contrast_results/treatment.csv",
                               "DA_summary.csv",
                               "combined_results.csv",
                               "results.xlsx")

  on.exit(unlink(c(expected_files_default,
                   "per_contrast_results"),
                 recursive = T), add = T)

  suppressMessages(write_limma_tables(input, overwrite = F))
  expect_true(all(file.exists(expected_files_default)))

  # Non-standard directory and files
  expected_files_custom <-   c("tempfortest/per_contrast_results_custom/control.csv",
                               "tempfortest/per_contrast_results_custom/treatment.csv",
                               "tempfortest/sum.csv",
                               "tempfortest/combo.csv",
                               "tempfortest/res.xlsx")
  on.exit(unlink(c(expected_files_custom,
                   "tempfortest"),
                 recursive = T), add = T)
  suppressMessages(
    write_limma_tables(input,
                       output_dir = "tempfortest",
                       overwrite = F,
                       contrasts_subdir = "per_contrast_results_custom",
                       summary_csv = "sum.csv",
                       combined_file_csv = "combo.csv",
                       spreadsheet_xlsx = "res.xlsx")
  )
  expect_true(all(file.exists(expected_files_custom)))
})


test_that("write_limma_tables creates output files when more than 12 stat results", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))


  input$results <- rep(input$results, 7)
  expected_files_default <-  c("per_contrast_results/control.csv",
                               "per_contrast_results/treatment.csv",
                               "DA_summary.csv",
                               "combined_results.csv",
                               "results.xlsx")

  on.exit(unlink(c(expected_files_default,
                   "per_contrast_results"),
                 recursive = T), add = T)

  suppressMessages(write_limma_tables(input, overwrite = F))
  expect_true(all(file.exists(expected_files_default)))

})


test_that("write_limma_tables overwrite cleans per-contrast dir and uses data_per_contrast intensities", {
  
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))
  
  # --- Create a synthetic 'normalized' data_per_contrast that differs from input$data ---
  # Make it obvious (add a constant), and also make row sets differ per contrast.
  stopifnot(all(c("control", "treatment") %in% names(input$results)))
  
  d_raw <- input$data
  
  # pick two overlapping subsets to mimic your real scenario
  n <- nrow(d_raw)
  idx1 <- seq_len(max(1, floor(n * 0.6)))
  idx2 <- seq.int(max(1, floor(n * 0.4)), n)
  
  d_control   <- d_raw[idx1, , drop = FALSE] + 1000
  d_treatment <- d_raw[idx2, , drop = FALSE] + 2000
  
  # --- NEW: make the overlap identical so get_export_data does not warn ---
  overlap <- intersect(rownames(d_control), rownames(d_treatment))
  if (length(overlap) > 0) {
    d_treatment[overlap, ] <- d_control[overlap, , drop = FALSE]
  }
  
  input$data_per_contrast <- list(
    control   = d_control,
    treatment = d_treatment
  )
  
  # Ensure norm.method is scalar to avoid grep warning
  if (length(input$tags$norm_method) > 1) input$tags$norm_method <- input$tags$norm_method[1]
  
  # --- Run in a temp output directory ---
  outdir <- file.path(tempdir(), "limma_tables_overwrite_regression")
  on.exit(unlink(outdir, recursive = TRUE, force = TRUE), add = TRUE)
  
  # 1) First run creates outputs
  suppressMessages(
    write_limma_tables(input, output_dir = outdir, overwrite = TRUE)
  )
  
  # Create a stale file that should NOT survive overwrite, and should NOT become an Excel sheet
  stale_csv <- file.path(outdir, "per_contrast_results", "stale.csv")
  dir.create(dirname(stale_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data.frame(a = 1), stale_csv, row.names = FALSE)
  expect_true(file.exists(stale_csv))
  
  # 2) Second run with overwrite should wipe stale artifacts
  suppressMessages(
    write_limma_tables(input, output_dir = outdir, overwrite = TRUE)
  )
  
  # stale CSV must be gone (this is the core regression fix)
  expect_false(file.exists(stale_csv))
  
  # Excel should not have a "stale" worksheet
  xlsx_file <- file.path(outdir, "results.xlsx")
  sheets <- openxlsx::getSheetNames(xlsx_file)
  expect_false("stale" %in% sheets)
  
  # Basic expected tabs should exist
  expect_true(all(c("Overview", "Protein Results") %in% sheets))
  
  # --- Check exported intensities come from data_per_contrast (normalized) not raw ---
  # Pick a protein that exists ONLY in control subset to avoid "first non-NA wins" ambiguity.
  prot_only_control <- setdiff(rownames(d_control), rownames(d_treatment))
  if (length(prot_only_control) == 0) {
    # fallback: choose any protein in control; still usually ok if it filled from control first
    prot_only_control <- rownames(d_control)[1]
  } else {
    prot_only_control <- prot_only_control[1]
  }
  
  combined_csv <- file.path(outdir, "combined_results.csv")
  combo <- utils::read.csv(combined_csv, check.names = FALSE)
  
  # Find the sample columns (they should match colnames(input$data))
  sample_cols <- intersect(colnames(d_raw), colnames(combo))
  expect_true(length(sample_cols) > 0)
  
  # Locate row for the chosen protein.
  # create_combined_results uses rownames(data) as ordering, and the first column in combo
  # isn't guaranteed to be the rowname. So we match by annotation key if possible.
  # Prefer uniprot_id if present, else use row position match via annotation rownames.
  if ("uniprot_id" %in% colnames(input$annotation) && "uniprot_id" %in% colnames(combo)) {
    id_val <- input$annotation[prot_only_control, "uniprot_id"]
    row_i <- which(combo$uniprot_id == id_val)[1]
  } else {
    # fallback: use rownames alignment assumption (works if annotation rownames were preserved)
    row_i <- match(prot_only_control, rownames(input$annotation))
  }
  
  expect_true(!is.na(row_i))
  
  # Compare one value (first sample) to ensure it's the normalized (control) value.
  s1 <- sample_cols[1]
  exported_val  <- as.numeric(combo[row_i, s1])
  expected_norm <- as.numeric(d_control[prot_only_control, s1])
  
  expect_equal(exported_val, expected_norm)
})




