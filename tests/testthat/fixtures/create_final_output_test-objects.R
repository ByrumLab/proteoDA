

# An object with the completed pipeline -----------------------------------
test_metadata <- data.frame(sample_ID = paste0("sample", 1:6),
                            group = rep(c("A", "B"), 3),
                            sex = rep(c("M", "F"), 3),
                            treatment = c(rep("control", 3),
                                          rep("treatment", 3)))
rownames(test_metadata) <- test_metadata$sample_ID
test_data <- as.data.frame(matrix(data = 1:60, nrow = 10))
colnames(test_data) <- test_metadata$sample_ID
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10))
input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata) |>
  normalize_data("log2") |>
  add_design(~ 0 + treatment) |>
  fit_limma_model() |>
  extract_DA_results()

saveRDS(object = input, file = "tests/testthat/fixtures/final_output_input.rds")
