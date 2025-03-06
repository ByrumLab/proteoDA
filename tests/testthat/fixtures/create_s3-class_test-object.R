

# Data frames of input ----------------------------------------------------
raw_data <- data.frame(sampleA = c(10000, 15000, 12500),
                       sampleB = c(20000, 30000, 25000),
                       sampleC = c(18000, 23000, 20500),
                       sampleD = c(36000, 46000, 41000))
protein_info <- data.frame(uniprot_id = c("A0A023ZSD2", "M9NH73", "H9AZR4"),
                           gene = c("Mcf", "shh", "WntA"))
sample_info <- data.frame(sample_id = c("sampleA", "sampleB",
                                        "sampleC", "sampleD"),
                          group = c("control",   "control",
                                    "treatment", "treatment"))
rownames(sample_info) <- colnames(raw_data)

input <- list(data = raw_data,
     annotation = protein_info,
     metadata = sample_info)
saveRDS(object = input, file = "tests/testthat/fixtures/s3-class_input.rds")

