

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

rm(list = ls())

# Protein filter, contaminants --------------------------------------------
test_metadata <- data.frame(sample_ID = paste0("sample", 1:10))
rownames(test_metadata) <- test_metadata$sample_ID
test_data <- as.data.frame(matrix(data = 1:100, nrow = 10))
colnames(test_data) <- test_metadata$sample_ID
test_annotation <- data.frame(protein_ID = paste0("protein", 1:10),
                              Protein.Name = c(rep("ok", 5),
                                               "DECOY",
                                               "DECOY",
                                               "Group of Something",
                                               "Group of Something",
                                               "Group of Something"))
input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)
saveRDS(object = input, file = "tests/testthat/fixtures/filter_protein_contam_input.rds")

# And a dataframe of how it should be after filtering
filtered <- input
filtered$data <- filtered$data[1:5,]
filtered$annotation <- filtered$annotation[1:5,]
# Set tags manually
filtered$tags$filter_proteins_by_annotation <- list(condition = substitute(!(stringr::str_detect(Protein.Name, "DECOY"))),
                                                    condition = substitute(!(stringr::str_detect(Protein.Name, "Group of"))))

saveRDS(object = filtered, file = "tests/testthat/fixtures/filter_protein_contam_output.rds")

rm(list = ls())


# Protein filter,  annotation ---------------------------------------------
test_metadata <- data.frame(sample_ID = paste0("sample", 1:10))
rownames(test_metadata) <- test_metadata$sample_ID
test_data <- as.data.frame(matrix(data = 1:100, nrow = 10))
colnames(test_data) <- test_metadata$sample_ID
test_annotation <- data.frame(protein_ID = paste0("protein", 1:10),
                              Protein.Name = c(rep("ok", 9),
                                               "keratin"),
                              molecular_weight = c(1:10),
                              extra = c(1:9, NA))
input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)
saveRDS(object = input, file = "tests/testthat/fixtures/filter_protein_annotation_input.rds")

# Some dataframes of how it should look after filtering
filtered <- input
filtered$data <- filtered$data[1:9,]
filtered$annotation <- filtered$annotation[1:9,]

# set tags manually on separate objects
filtered_ID <- filtered
filtered_ID$tags$filter_proteins_by_annotation <- list(condition = substitute(protein_ID != "protein10"))

filtered_name <- filtered
filtered_name$tags$filter_proteins_by_annotation <- list(condition = substitute(!stringr::str_detect(Protein.Name, "keratin")))

filtered_MW <- filtered
filtered_MW$tags$filter_proteins_by_annotation <- list(condition = substitute(molecular_weight < 10))


saveRDS(object = filtered_ID, file = "tests/testthat/fixtures/filter_protein_annotation_output_ID.rds")
saveRDS(object = filtered_name, file = "tests/testthat/fixtures/filter_protein_annotation_output_name.rds")
saveRDS(object = filtered_MW, file = "tests/testthat/fixtures/filter_protein_annotation_output_MW.rds")

rm(list = ls())




