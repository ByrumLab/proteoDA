

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
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10))
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
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10),
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
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10),
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
filtered_ID$tags$filter_proteins_by_annotation <- list(condition = substitute(uniprot_id != "protein10"))

filtered_name <- filtered
filtered_name$tags$filter_proteins_by_annotation <- list(condition = substitute(!stringr::str_detect(Protein.Name, "keratin")))

filtered_MW <- filtered
filtered_MW$tags$filter_proteins_by_annotation <- list(condition = substitute(molecular_weight < 10))


saveRDS(object = filtered_ID, file = "tests/testthat/fixtures/filter_protein_annotation_output_ID.rds")
saveRDS(object = filtered_name, file = "tests/testthat/fixtures/filter_protein_annotation_output_name.rds")
saveRDS(object = filtered_MW, file = "tests/testthat/fixtures/filter_protein_annotation_output_MW.rds")

rm(list = ls())


# Protein filter, by group number non-missing ----------------------------------
test_metadata <- data.frame(sample_ID = paste0("sample", 1:12),
                            group = c(rep("A", 4),
                                      rep("B", 4),
                                      rep("C", 4)))
rownames(test_metadata) <- test_metadata$sample_ID
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10))

test_data <- as.data.frame(matrix(data = 1:120, nrow = 10))
colnames(test_data) <- test_metadata$sample_ID

# Set proteins to different levels of missingness
# Protein 1 not present in first group at all
# protein 2 only present in 1 sample in first two groups
# protein 3 only present in 1 sample in all groups
# Make last sample all 0s. This shouldn't affect things
test_data[1,1:4] <- NA
test_data[2,c(1,5)] <- NA
test_data[3, c(1:3, 5:7, 9:11)] <- NA
test_data[,12] <- 0

input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)
saveRDS(object = input, file = "tests/testthat/fixtures/filter_proteins_by_group_input.rds")


# min_reps = 2, min_groups = 1
# should filter out protein 3
output_21 <- DAList(data = test_data[-3,],
                    annotation = test_annotation[-3,,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_group = list(min_reps = 2,
                                                                min_groups = 1,
                                                                grouping_column = "group")))
saveRDS(object = output_21, file = "tests/testthat/fixtures/filter_proteins_by_group_output21.rds")

# min_reps = 1, min_groups = 3
# should filter out protein 1
output_13 <- DAList(data = test_data[-1,],
                    annotation = test_annotation[-1,,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_group = list(min_reps = 1,
                                                                min_groups = 3,
                                                                grouping_column = "group")))
saveRDS(object = output_13, file = "tests/testthat/fixtures/filter_proteins_by_group_output13.rds")


# min_reps = 2, min_groups = 3
# should filter out protein 1 and 3
output_23 <- DAList(data = test_data[c(-1,-3),],
                    annotation = test_annotation[c(-1,-3),,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_group = list(min_reps = 2,
                                                                min_groups = 3,
                                                                grouping_column = "group")))
saveRDS(object = output_23, file = "tests/testthat/fixtures/filter_proteins_by_group_output23.rds")



# min_reps = 3, min_groups = 3
# should filter out proteins 1, and 3
output_33 <- DAList(data = test_data[c(-1, -3),],
                    annotation = test_annotation[c(-1, -3),,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_group = list(min_reps = 3,
                                                                min_groups = 3,
                                                                grouping_column = "group")))
saveRDS(object = output_33, file = "tests/testthat/fixtures/filter_proteins_by_group_output33.rds")

# min_reps = 4, min_groups = 3
# should filter out proteins 1, 2 and 3
output_43 <- DAList(data = test_data[c(-1, -2, -3),],
                    annotation = test_annotation[c(-1, -2, -3),,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_group = list(min_reps = 4,
                                                                min_groups = 3,
                                                                grouping_column = "group")))
saveRDS(object = output_43, file = "tests/testthat/fixtures/filter_proteins_by_group_output43.rds")

# Protein filter, by group proportion non-missing ----------------------------------
test_metadata <- data.frame(sample_ID = paste0("sample", 1:12),
                            group = c(rep("A", 4),
                                      rep("B", 4),
                                      rep("C", 4)))
rownames(test_metadata) <- test_metadata$sample_ID
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:10))

test_data <- as.data.frame(matrix(data = 1:120, nrow = 10))
colnames(test_data) <- test_metadata$sample_ID

# Set proteins to different levels of missingness
# Thresholds: 0, 25, 50, 75, 100
# Protein 1 not present in first group at all
# protein 2
# protein 3 only present in 1 sample in all groups
# Make last sample all 0s. This shouldn't affect things

test_data[1,1:4] <- NA
test_data[2,2:4] <- NA
test_data[3,3:4] <- NA
test_data[4,4] <- NA
test_data[,12] <- 0

input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)
saveRDS(object = input, file = "tests/testthat/fixtures/filter_proteins_by_proportion_input.rds")

# min_prop = 0
# should remove nothing
output_0 <- DAList(data = test_data[,],
                    annotation = test_annotation[,,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_proportion = list(min_prop = 0,
                                                                     grouping_column = "group")))
saveRDS(object = output_0, file = "tests/testthat/fixtures/filter_proteins_by_proportion_output0.rds")

# min_prop = 25
# should remove protein 1
output_25 <- DAList(data = test_data[c(-1),],
                   annotation = test_annotation[c(-1),,drop = F],
                   metadata = test_metadata,
                   tags = list(filter_proteins_by_proportion = list(min_prop = 0.25,
                                                                    grouping_column = "group")))
saveRDS(object = output_25, file = "tests/testthat/fixtures/filter_proteins_by_proportion_output25.rds")

# min_prop = 50
# should remove protein 1 and 2
output_50 <- DAList(data = test_data[c(-1,-2),],
                    annotation = test_annotation[c(-1,-2),,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_proportion = list(min_prop = 0.50,
                                                                     grouping_column = "group")))
saveRDS(object = output_50, file = "tests/testthat/fixtures/filter_proteins_by_proportion_output50.rds")

# min_prop = 75
# should remove protein 1 and 2 and 3
output_75 <- DAList(data = test_data[c(-1,-2, -3),],
                    annotation = test_annotation[c(-1,-2, -3),,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_proportion = list(min_prop = 0.75,
                                                                     grouping_column = "group")))
saveRDS(object = output_75, file = "tests/testthat/fixtures/filter_proteins_by_proportion_output75.rds")

# min_prop = 1
# should remove protein 1 and 2 and 3 and 4
output_1 <- DAList(data = test_data[c(-1,-2, -3, -4),],
                    annotation = test_annotation[c(-1,-2, -3, -4),,drop = F],
                    metadata = test_metadata,
                    tags = list(filter_proteins_by_proportion = list(min_prop = 1,
                                                                     grouping_column = "group")))
saveRDS(object = output_1, file = "tests/testthat/fixtures/filter_proteins_by_proportion_output1.rds")
