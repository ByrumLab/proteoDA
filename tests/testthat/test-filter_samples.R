


test_that("filter_samples removes sample data as expected", {

  # An initial, ugly test to start getting things going.
  # But, should probably move this to a separate script
  # within the test quite where I make a few different objects for testing
  # Maybe do this after I create an user function for making an initial DAList?


  # Assemble a data frame to test filtering on
  test_metadata <- data.frame(sample_ID = paste0("sample", 1:10),
                              keeper = c(T,T,T,T,T,F,F,F,F,F))
  rownames(test_metadata) <- test_metadata$sample_ID
  test_data <- as.data.frame(matrix(data = 1:100, nrow = 10))
  colnames(test_data) <- test_metadata$sample_ID
  test_annotation <- data.frame(protein_ID = paste0("protein", 1:10))
  to_filter <- new_DAList(x = list(data = test_data,
                                    annotation = test_annotation,
                                    metadata = test_metadata

  ))

  rm(test_data, test_annotation, test_metadata)

  # And a dataframe of how it should be after filtering
  filtered_metadata <- data.frame(sample_ID = paste0("sample", 1:5),
                                  keeper = c(T,T,T,T,T))
  rownames(filtered_metadata) <- filtered_metadata$sample_ID
  filtered_data <- as.data.frame(matrix(data = 1:100, nrow = 10))[,1:5]
  colnames(filtered_data) <- filtered_metadata$sample_ID
  filtered_annotation <- data.frame(protein_ID = paste0("protein", 1:10))
  filtered <- new_DAList(x = list(data = filtered_data,
                                   annotation = filtered_annotation,
                                   metadata = filtered_metadata

  ))

  expect_equal(filter_samples(to_filter, keeper == T), filtered)
  expect_equal(filter_samples(to_filter, sample_ID %notin% paste0("sample", 6:10)), filtered)
})
