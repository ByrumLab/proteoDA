

# Input object ------------------------------------------------------------

# An input DAList object for testing adding design
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
test_annotation <- data.frame(protein_ID = paste0("protein", 1:10))
input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)
saveRDS(object = input, file = "tests/testthat/fixtures/add_design_input.rds")


# Output objects ----------------------------------------------------------

# One term intercept
out_treat_intercept <- input
design_matrix <- as.matrix(
  data.frame(
    `Intercept` = c(1,1,1,1,1,1),
    treatment = c(0,0,0,1,1,1),
    row.names = paste0("sample", 1:6)
    )
)
attr(design_matrix, "assign") <- c(0,1)
attr(design_matrix, "contrasts") <- list("treatment" = "contr.treatment")

out_treat_intercept$design <- list(design_formula = "~treatment",
                                         design_matrix = design_matrix)

saveRDS(object = out_treat_intercept, file = "tests/testthat/fixtures/add_design_output_treat_intercept.rds")


# One term no intercept
out_treat_nointercept <- input
design_matrix <- as.matrix(
  data.frame(
    control = c(1,1,1,0,0,0),
    treatment = c(0,0,0,1,1,1),
    row.names = paste0("sample", 1:6)
  )
)
attr(design_matrix, "assign") <- c(1,1)
attr(design_matrix, "contrasts") <- list("treatment" = "contr.treatment")

out_treat_intercept$design <- list(design_formula = "~0 + treatment",
                                   design_matrix = design_matrix)

saveRDS(object = out_treat_intercept, file = "tests/testthat/fixtures/add_design_output_treat_nointercept.rds")

# interaction intercept
out_interaction_intercept <- input
design_matrix <- as.matrix(
  data.frame(
    `Intercept` = c(1,1,1,1,1,1),
    treatment = c(0,0,0,1,1,1),
    M = c(1,0,1,0,1,0),
    treatment.sexM = c(0,0,0,0,1,0),
    row.names = paste0("sample", 1:6)
  )
)
attr(design_matrix, "assign") <- c(0:3)
attr(design_matrix, "contrasts") <- list("treatment" = "contr.treatment",
                                         "sex" = "contr.treatment")

out_interaction_intercept$design <- list(design_formula = "~treatment * sex",
                                   design_matrix = design_matrix)

saveRDS(object = out_interaction_intercept, file = "tests/testthat/fixtures/add_design_output_interact_intercept.rds")


# interaction no intercept
out_interaction_nointercept <- input
design_matrix <- as.matrix(
  data.frame(
    control = c(1,1,1,0,0,0),
    treatment = c(0,0,0,1,1,1),
    M = c(1,0,1,0,1,0),
    treatment.sexM = c(0,0,0,0,1,0),
    row.names = paste0("sample", 1:6)
  )
)
attr(design_matrix, "assign") <- c(1, 1:3)
attr(design_matrix, "contrasts") <- list("treatment" = "contr.treatment",
                                         "sex" = "contr.treatment")

out_interaction_nointercept$design <- list(design_formula = "~0 + treatment * sex",
                                         design_matrix = design_matrix)

saveRDS(object = out_interaction_nointercept, file = "tests/testthat/fixtures/add_design_output_interact_nointercept.rds")

# Random effect
out_treat_mixed <- input
design_matrix <- as.matrix(
  data.frame(
    `Intercept` = c(1,1,1,1,1,1),
    treatment = c(0,0,0,1,1,1),
    row.names = paste0("sample", 1:6)
  )
)
attr(design_matrix, "assign") <- c(0,1)
attr(design_matrix, "contrasts") <- list("treatment" = "contr.treatment")

out_treat_mixed$design <- list(design_formula = "~treatment + (1 | group)",
                                   design_matrix = design_matrix,
                                   random_factor = "group")

saveRDS(object = out_treat_mixed, file = "tests/testthat/fixtures/add_design_output_treat_mixed.rds")


