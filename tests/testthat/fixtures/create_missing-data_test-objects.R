
# Missing to zero objects -------------------------------------------------
test_metadata <- data.frame(sample_ID = paste0("sample", 1:10))
rownames(test_metadata) <- test_metadata$sample_ID
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10))

test_data <- matrix(data = c(1:100), nrow = 10)
diag(test_data) <- c(rep(NA, 5), rep(-9, 5))
test_data <- as.data.frame(test_data)
colnames(test_data) <- test_metadata$sample_ID

input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)

saveRDS(object = input, file = "tests/testthat/fixtures/missing_to_zero_input.rds")

# Set up output data
out_data_NA <- test_data
out_data_NA[is.na(out_data_NA)] <- 0

out_data_neg9 <- test_data
out_data_neg9[out_data_neg9 == -9] <- 0

out_data_both <- test_data
out_data_both[is.na(out_data_both)] <- 0
out_data_both[out_data_both == -9] <- 0

output_NA <- DAList(data = out_data_NA,
                    annotation = test_annotation,
                    metadata = test_metadata)

output_neg9 <- DAList(data = out_data_neg9,
                    annotation = test_annotation,
                    metadata = test_metadata)

output_both <- DAList(data = out_data_both,
                    annotation = test_annotation,
                    metadata = test_metadata)

saveRDS(object = output_NA, file = "tests/testthat/fixtures/missing_to_zero_output_NA.rds")
saveRDS(object = output_neg9, file = "tests/testthat/fixtures/missing_to_zero_output_neg9.rds")
saveRDS(object = output_both, file = "tests/testthat/fixtures/missing_to_zero_output_both.rds")


# Zero to missing objects -------------------------------------------------
test_metadata <- data.frame(sample_ID = paste0("sample", 1:10))
rownames(test_metadata) <- test_metadata$sample_ID
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10))

test_data <- matrix(data = c(1:100), nrow = 10)
diag(test_data) <- 0
test_data <- as.data.frame(test_data)
colnames(test_data) <- test_metadata$sample_ID

input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)

saveRDS(object = input, file = "tests/testthat/fixtures/zero_to_missing_input.rds")

out_data <- test_data
out_data[out_data == 0] <- NA

output <- DAList(data = out_data,
                 annotation = test_annotation,
                 metadata = test_metadata)

saveRDS(object = output, file = "tests/testthat/fixtures/zero_to_missing_output.rds")
