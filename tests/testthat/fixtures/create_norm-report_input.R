# Input object ------------------------------------------------------------

# An input DAList object for testing adding design
# Data and annotations portions are simple,
# metadata is what matters here
test_metadata <- data.frame(sample_ID = paste0("sample", 1:6),
                            group = rep(c("A", "B", "C"), 2),
                            sex = c(rep(c("M", "F"), 2), "M", "M"),
                            treatment = c(rep("control", 3),
                                          rep("treatment", 3)))
rownames(test_metadata) <- test_metadata$sample_ID
set.seed(42)
test_data <- data.frame(a = as.integer(runif(n = 100, min = 10000, max = 20000)),
                        b = as.integer(runif(n = 100, min = 10000, max = 20000)),
                        c = as.integer(runif(n = 100, min = 10000, max = 20000)),
                        d = as.integer(runif(n = 100, min = 10000, max = 20000)),
                        e = as.integer(runif(n = 100, min = 10000, max = 20000)),
                        f = as.integer(runif(n = 100, min = 10000, max = 20000)))
colnames(test_data) <- test_metadata$sample_ID
test_annotation <- data.frame(uniprot_id = paste0("protein", 1:100))
input <- DAList(data = test_data,
                annotation = test_annotation,
                metadata = test_metadata)
saveRDS(object = input, file = "tests/testthat/fixtures/norm_report_input.rds")

