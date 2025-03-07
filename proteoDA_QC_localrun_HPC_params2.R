#library(devtools) # nocov start
library("proteoDA")
library(tidyverse)
library(yaml)

#devtools::load_all() # run local functions without installing

## proteoDA_params.R should be in same directory as this file
## test this 
source("proteoDA_params.R")

# DiaNN quan output
setwd(working_dir)

input_data <- read.csv(input_quan)
sample_metadata <- read.csv(metadata)

head(input_data)
head(sample_metadata)

intensity_data <- input_data[,sample_start:ncol(input_data)] # select columns 5 to 21
annotation_data <- input_data[,anno_start:anno_end] # select columns 1 to 4

rownames(sample_metadata)
rownames(sample_metadata) <- sample_metadata$sample
#rownames of metadata must match column names of intensity data
rownames(sample_metadata)
colnames(intensity_data)

##########
### Create DAList object
##########

raw <- DAList(data = intensity_data,
              annotation = annotation_data,
              metadata = sample_metadata,
              design = NULL,
              eBayes_fit = NULL,
              results = NULL,
              tags = NULL)

######

# Convert missing values
filtered_samples <- zero_to_missing(raw)

### in general require 2/3 of biological replicates to have a value in 1 sample group
filtered_proteins <- filter_proteins_by_group(filtered_samples,
                                              min_reps = filt_min_reps,
                                              min_groups = filt_min_groups,
                                              grouping_column = group)
###### OUTPUT
#Keeping only protein entries with non-missing intensity in at least 2 samples in at least 1 group
#Filtered 35 entries from the dataset leaving 8358 entries for analysis

# if using the Rshiny raw norm values from DiaNN
#normlog2 <- normalize_data(filtered_proteins,
#                           norm_method = "log2")

# using the log2 norm output from diann_quan

#norm <- filtered_proteins
norm <- filtered_samples
norm$tags$norm_method <- "diann_quan"


###############
## ADD QC PLOTS HERE
##################
# either copy script "proteoDA_QC_plots.R" or use Rscript to run it after this one
##
# WRITE QC REPORT
# Violin plot, PCA, dendrograms, missing value heatmap
###

# write_qc_report(norm,
#                 color_column = group,
#                 label_column = NULL,
#                 output_dir = QC_dir,
#                 filename = QC_report,
#                 overwrite = T,
#                 top_proteins = pca_top_prot,        # number of most variable proteins
#                 standardize = TRUE,
#                 pca_axes = c(1,2),          # first 2 PCs
#                 dist_metric = "euclidean",  # stats::dist for options
#                 clust_method = "complete",  # stats::hclust for options
#                 show_all_proteins = F)      # only proteins with missing data)


# limma model setup
no_intercept <- add_design(norm,
                            design_formula = design)

# add sample group comparisons
# defining in code
#no_intercept <- add_contrasts(no_intercept,
 #                             contrasts_vector = "CHP212_Trt_vs_CHP212_ctrl = CHP212_Trt - CHP212_Ctrl")

no_intercept <- add_contrasts(no_intercept,
                              contrasts_file = contrasts)

############
## RUN LIMMA
########

# First, fit the model
fit <- fit_limma_model(no_intercept)

# Extract results ---------------------------------------------------------

# Define the significance thresholds for Volcano and MD plots
# pvalue threshold for both unadjusted and the adjusted p-value - default 0.05
# log2 fold change cutoff - default is set to 2-fold change
results <- extract_DA_results(fit,
                              pval_thresh = p.val,
                              lfc_thresh = logFC,
                              adj_method = "BH",
                              extract_intercept = F) #For models with intercept terms

###########
## Standard deviation for biological replicates

##############
# WRITE DATA TABLES TO FILE
# includes an overall Excel file with all comparisons added
####

write_limma_tables(results,
                   output_dir = DA_dir,
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = T)


##############
# WRITE PLOTS
# includes an html for each comparison that is interactive and stand alone
# open html in web browser to explore the data
#########

## table_columns must match the protein annotation information. Displayed in the html Table
## title_column is displayed on the protein intensity plot (can be changed to Genes)

results <- missing_to_zero(results)

write_limma_plots(results,
                  grouping_column = group,
                  table_columns = DA_table_cols,
                  title_column = DA_title_col,
                  output_dir = DA_dir,
                  tmp_subdir = tmp_subdir,
                  overwrite = T,
                  height = 1000,
                  width = 1000)


# ##### FINAL DAList()
# save for QC plots
saveRDS(results, "results.rds")

###########
## convert 02_DA/static_plots/.pdfs to .pngs

####
#dir.create(QC_dir)
if (!dir.exists(QC_dir)) {dir.create(QC_dir)} 

library(pdftools)

# Define the output directory
output_dir <- file.path(DA_dir, "static_plots")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Convert Adjusted p-value Volcano plots from PDF to PNG
for (i in 1:length(results$results)) {

  # Correctly construct the input PDF file path
  pdf_file <- file.path(DA_dir, "static_plots", paste0(names(results$results)[i], "-volcano-adjusted-pval.pdf"))

  # Ensure the PDF file exists before converting
  if (file.exists(pdf_file)) {

    # Convert to PNG with proper filename pattern
    pdf_convert(
      pdf_file,
      format = "png",
      pages = NULL,
      filenames = file.path(output_dir, paste0(names(results$results)[i], "-volcano-adjusted-pval_%d.png")),
      dpi = 600
    )

  } else {
    message("File not found: ", pdf_file)
  }
}


for (i in 1:length(results$results)) {

  # Correctly construct the input PDF file path
  pdf_file <- file.path(DA_dir, "static_plots", paste0(names(results$results)[i], "-volcano-raw-pval.pdf"))

  # Ensure the PDF file exists before converting
  if (file.exists(pdf_file)) {

    # Convert to PNG with proper filename pattern
    pdf_convert(
      pdf_file,
      format = "png",
      pages = NULL,
      filenames = file.path(output_dir, paste0(names(results$results)[i], "-volcano-raw-pval_%d.png")),
      dpi = 600
    )

  } else {
    message("File not found: ", pdf_file)
  }
}

# Define the output directory
output_dir <- file.path(DA_dir, "static_plots")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Convert Adjusted p-value Volcano plots from PDF to PNG
for (i in 1:length(results$results)) {

  # Correctly construct the input PDF file path
  pdf_file <- file.path(DA_dir, "static_plots", paste0(names(results$results)[i], "-MD-adjusted-pval.pdf"))

  # Ensure the PDF file exists before converting
  if (file.exists(pdf_file)) {

    # Convert to PNG with proper filename pattern
    pdf_convert(
      pdf_file,
      format = "png",
      pages = NULL,
      filenames = file.path(output_dir, paste0(names(results$results)[i], "-MD-adjusted-pval_%d.png")),
      dpi = 600
    )

  } else {
    message("File not found: ", pdf_file)
  }
}


for (i in 1:length(results$results)) {

  # Correctly construct the input PDF file path
  pdf_file <- file.path(DA_dir, "static_plots", paste0(names(results$results)[i], "-MD-raw-pval.pdf"))

  # Ensure the PDF file exists before converting
  if (file.exists(pdf_file)) {

    # Convert to PNG with proper filename pattern
    pdf_convert(
      pdf_file,
      format = "png",
      pages = NULL,
      filenames = file.path(output_dir, paste0(names(results$results)[i], "-MD-raw-pval_%d.png")),
      dpi = 600
    )

  } else {
    message("File not found: ", pdf_file)
  }
}


### UPDATE _variables.yml with params for Project Name and authors.

variables <- yaml.load_file("_variables.yml")
variables$title <- project_name
variables$author <- author
variables$author2 <- author2

# You can also add new variables
# variables$new_var <- 42

write_yaml(variables, "_variables.yml")
