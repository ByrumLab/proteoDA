# Input object ------------------------------------------------------------

# An input DAList object for testing adding contrasts
# Data and annotations portions are simple,
# metadata is what matters here
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
  add_design(~ 0 + treatment)
saveRDS(object = input, file = "tests/testthat/fixtures/add_contrasts_input.rds")




# Create some good and bad contrast CSV files -----------------------------
# Good contrast
data.frame(
  x = c("Treatment_vs_Control= treatment - control")
) |>
  write.table(file = "tests/testthat/fixtures/good_contrast.csv", col.names = FALSE, row.names = F, sep = ",")

# Multiple column
data.frame(
  x = c("Treatment_vs_Control= treatment - control"),
  y = c("second column")
) |>
  write.table(file = "tests/testthat/fixtures/bad_contrast_mult_col.csv", col.names = F, row.names = F, sep = ",")

# Group not present column
data.frame(
  x = c("Treatment_vs_Control= treatment - XXX")
) |>
  write.table(file = "tests/testthat/fixtures/bad_contrast_wrong_group.csv", col.names = F, row.names = F, sep = ",")


