#library(devtools) # nocov start
library("proteoDAstjude")
#library(tidyverse)
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

########
## USE THIS TO FILTER PROTEINS GLOBALLY ACROSS ALL GROUPS
### in general require 2/3 of biological replicates to have a value in 1 sample group
filtered_proteins <- filter_proteins_by_group(filtered_samples,
                                              min_reps = filt_min_reps,
                                              min_groups = filt_min_groups,
                                              grouping_column = group)

########
## USE THIS TO FILTER PROTEINS PER CONTRAST SEPARATELY
# creates multiple DAList objects for analysis.

filtered_DALists <- filter_proteins_per_contrast(
  DAList = filtered_samples,
  contrasts_file = "data/contrasts.csv",
  min_reps = 3,
  require_both_groups = FALSE,
  grouping_column = "group"
)
summary_df <- attr(filtered_DALists, "retention_summary")
write.csv(summary_df, "filtered_protein_counts.csv", row.names = FALSE)


# if using the raw norm values from DiaNN and proteoDAs normalization function
#normlog2 <- normalize_data(filtered_proteins,
#                           norm_method = "log2")

# using the log2 norm output from diann_quan

### set the DAList to the filter proteins per contrast
norm <- filtered_DALists

### set the DAList to the filter proteins globally
#norm <- filtered_proteins
## set the DAList to remove samples 
# norm <- filtered_samples 

## set the tag so it doesn't complain about not running normalization function.
## not needed if using proteoDA normalization
norm$tags$norm_method <- "diann_quan"



############
## RUN LIMMA
########

# Define the significance thresholds for Volcano and MD plots
# pvalue threshold for both unadjusted and the adjusted p-value - default 0.05
# log2 fold change cutoff - default is set to 2-fold change


###### Calculate Rolling Standard Deviation
# if passes G-test, then fill in missing values if significant
## plot Volcano significance using moving_sd parameter 
# performed within the new run_filtered_limma_analysis function 
#results <- compute_movingSD_zscores(results, binsize = 1000)

results <- run_filtered_limma_analysis(
  DAList = norm,
  design_formula = design,
  pval_thresh = 0.05,
  lfc_thresh = 0,
  adj_method = "BH",
  contrasts_file  = contrasts 
)

##############
# WRITE DATA TABLES TO FILE
# includes an overall Excel file with all comparisons added
####

# Build statlist with full protein list across all contrasts
statlist <- build_statlist(
  DAList = results,
  stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSDs", "logFC_z_scores", "sig.PVal", "sig.FDR"))
# annotation_cols = c("uniprot_id", "Genes", "Accession.Number", "Identified.Peptides", "Protein.Description")
#)

write_limma_tables(results,
                   output_dir = DA_dir,
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = T,
                   color_palette = NULL,
                   add_contrast_sheets = TRUE,
                   statlist = statlist)


##############
# WRITE PLOTS
# includes an html for each comparison that is interactive and stand alone
# open html in web browser to explore the data
#########

## table_columns must match the protein annotation information. Displayed in the html Table
## title_column is displayed on the protein intensity plot (can be changed to Genes)

saveRDS(results, "results_NA.rds") # for QC reports 

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
variables$project <- project_name

# You can also add new variables
# variables$new_var <- 42

write_yaml(variables, "_variables.yml")

sink("sessionInfo.txt")
sessionInfo()
sink()
