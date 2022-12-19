

# Sample filtering --------------------------------------------------------

# An input DAList object for testing sample-level filtering
# Data and annotations portions are simple,
# metadata is what matters here
# Set up a few different columns that all
# filter out the last five samples
test_metadata <- data.frame(sample_ID = paste0("sample", 1:10),
                            keeper = c(T,T,T,T,T,F,F,F,F,F),
                            numeric = 1:10)
rownames(test_metadata) <- test_metadata$sample_ID
test_data <- as.data.frame(matrix(data = 1:100, nrow = 10))
colnames(test_data) <- test_metadata$sample_ID
test_annotation <- data.frame(protein_ID = paste0("protein", 1:10))
input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)
saveRDS(object = input, file = "tests/testthat/fixtures/filter_samples_input.rds")

# And a dataframe of how it should be after filtering
filtered <- input
filtered$data <- filtered$data[,1:5]
filtered$metadata <- filtered$metadata[1:5,]
saveRDS(object = filtered, file = "tests/testthat/fixtures/filter_samples_output.rds")


# Protein filter, contaminants --------------------------------------------






# Protein filter,  annotation ---------------------------------------------






