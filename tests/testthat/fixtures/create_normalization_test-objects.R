set.seed(42)
input <- as.matrix(data.frame(V1 = runif(500, 0, 100),
                            V2 =  runif(500, 0, 100),
                            V3 =  runif(500, 0, 100),
                            V4 =  runif(500, 0, 100)))

# Add in some 0s and NAs
input[1:5,2:4] <- 0
input[6:10,1:2] <- NA
saveRDS(object = input, file = "tests/testthat/fixtures/normalization_external_fxn_input.rds")
