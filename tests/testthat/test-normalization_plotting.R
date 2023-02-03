test_that("eval_pn_metric_for_plot output has proper structure", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))


  normList <- apply_all_normalizations(input$data)

  groupings <- alist(sex = input$metadata$sex,
                     group =  input$metadata$group,
                     treatment = input$metadata$treatment)

  # Expected structure for variability stats
  expected_rows <- c(sex = length(normList)*length(unique(input$metadata$sex)),
                     group = length(normList)*length(unique(input$metadata$group)),
                     treatment = length(normList)*length(unique(input$metadata$treatment)))

  for (grouping in names(groupings)) {
    for (metric in c("PCV", "PMAD", "PEV")) {
      output <- eval_pn_metric_for_plot(normList, groupings[[grouping]], metric)
      # proper colnames
      expect_equal(colnames(output), c("method", "group", "value"))
      # proper groups
      expect_equal(unique(output$group), unique(eval(groupings[[grouping]])))
      # Proper normalizations
      expect_equal(unique(as.character(output$method)), names(normList))
      # proper row numbers
      expect_equal(nrow(output), expected_rows[[grouping]])
    }
  }

  # Expected structure for COR
  expected_rows <- c(sex = length(normList)*(6 +1),
                     group = length(normList)*3,
                     treatment = length(normList)*(3 + 3))

  # Expected structure for COR
  for (grouping in names(groupings)) {
    for (metric in c("COR")) {
      output <- eval_pn_metric_for_plot(normList, groupings[[grouping]], metric)
      # proper colnames
      expect_equal(colnames(output), c("method", "group", "value"))
      # Proper normalizations
      expect_equal(unique(as.character(output$method)), names(normList))
      # proper row numbers
      expect_equal(nrow(output), expected_rows[[grouping]])
    }
  }

  # Expected structure for log2ratio
  expected_rows <- c(sex = length(normList)*nrow(input$data)*1,
                     group = length(normList)*nrow(input$data)*3,
                     treatment = length(normList)*nrow(input$data)*1)

  # Expected structure for log2ratio
  for (grouping in names(groupings)) {
    for (metric in c("log2ratio")) {
      output <- eval_pn_metric_for_plot(normList, groupings[[grouping]], metric)
      # proper colnames
      expect_equal(colnames(output), c("method", "value"))
      # Proper normalizations
      expect_equal(unique(as.character(output$method)), names(normList))
      # proper row numbers
      expect_equal(nrow(output), expected_rows[[grouping]])
    }
  }
})


# At the moment, not going to do testing of the
# various plotting functions, which should all output ggplot objects
# and testing all of them for correctness would be tricky.
