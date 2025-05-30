# updated proteoDAstjude package May 27, 2025 - fixed logFC_zscore
library(proteoDAstjude)
library(yaml)

# Load required packages to convert DIANN quan to proteoDA format
library(readr)
library(stringr)

#######
# load project parameters 
source("proteoDA_params.R")

############
## File format conversion - DIANN quan to proteoDA
###########

# Read the CSV files for testing
# input_data <- read.csv("data/uni_prot_quan_rmNA_norm.csv")

# argument passed from bash 
args <- commandArgs(trailingOnly = TRUE)
diann_quan <- args[1]
input_data <- read.csv(diann_quan)

sample_metadata <- read.csv(metadata)


# Step 1: Create 'uniprot_id' by extracting the first ID from 'Protein.Group'
# Create 'uniprot_id' as the first column
uniprot_split <- strsplit(input_data$Protein.Group, ";")
uniprot_id <- sapply(uniprot_split, `[`, 1)
input_data <- cbind(uniprot_id = uniprot_id, input_data)

# Step 2: Replace NA and blank entries safely, preserving numeric types
for (col in names(input_data)) {
  if (is.numeric(input_data[[col]])) {
    input_data[[col]][is.na(input_data[[col]])] <- 0
  } else {
    input_data[[col]][is.na(input_data[[col]]) | input_data[[col]] == ""] <- "0"
  }
}

# # === Step 3: Rename sample columns based on sample_metadata ===
# # Create a named vector: names = file column, values = sample column
# col_rename_map <- setNames(sample_metadata$sample, sample_metadata$file)
# names(input_data)
# # Find which of those file names are actual column names in input_data
# matching_cols <- intersect(names(col_rename_map), names(input_data))
# 
# # Rename the matching columns
# names(input_data)[names(input_data) %in% matching_cols] <- col_rename_map[matching_cols]

# === Step 4: Select relevant columns ===
#avg_intensity_index <- which(names(input_data) == "Average.Intensity")
#prefix_cols <- names(input_data)[1:avg_intensity_index]
#sample_cols <- sample_metadata$sample
#final_cols <- c(prefix_cols, intersect(sample_cols, names(input_data)))

# Subset the data frame
#reformatted_df <- input_data[, final_cols]

# Save the reformatted data
#write_csv(reformatted_df, "uni_prot_quan_rmNA_norm_proteoDA.csv")
#########################

#input_data <- reformatted_df

###########################
# Step 0: create anno before removing for the intensity data. 
annotation_data <- input_data[,anno_start:anno_end] # select protein annotation 

# Step 1: Create mapping from file name to sample name
col_rename_map <- setNames(sample_metadata$sample, sample_metadata$file)

# Step 2: Subset the map to only those files present in input_data
matching_files <- intersect(names(col_rename_map), colnames(input_data))

# Step 3: Rename the columns in input_data
colnames(input_data)[match(matching_files, colnames(input_data))] <- col_rename_map[matching_files]

# Step 4: Reorder columns to match the sample_metadata$sample order
ordered_samples <- sample_metadata$sample
# Keep only those sample names that are now in input_data
ordered_samples <- ordered_samples[ordered_samples %in% colnames(input_data)]
intensity_data <- input_data[, ordered_samples]

head(intensity_data)
head(sample_metadata)

#intensity_data <- input_data[,sample_start:ncol(input_data)] # select sample columns 


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
table(filtered_samples$metadata$group)

###############
##  Filtering proteins with missing values
## Option 1: filter_proteins_by_group: filters globally using all samples in the data
## Option 2: filter_proteins_per_contrast: filters proteins for each sample group comparison separately
################

# ### OPTION 1 
# in general require 2/3 of biological replicates to have a value in 1 sample group
# filtered_proteins <- filter_proteins_by_group(filtered_samples,
#                                               min_reps = filt_min_reps,
#                                               min_groups = filt_min_groups,
#                                               grouping_column = group)

# ### OPTION 2
# Filter proteins separately for each contrast (both groups must meet threshold if require_both_groups = TRUE)
filtered_DALists <- filter_proteins_per_contrast(
  DAList = filtered_samples,
  contrasts_file = contrasts,
  min_reps = filt_min_reps,
  require_both_groups = require_both_groups,
  grouping_column = group
)
# Create summary_df from filtered_proteins_per_contrast
summary_df <- data.frame(
  contrast = names(filtered_DALists$filtered_proteins_per_contrast),
  filtered_proteins = sapply(filtered_DALists$filtered_proteins_per_contrast, length)
)

filtered_DALists2 <- filter_proteins_per_contrast(
  DAList = filtered_samples,
  contrasts_file = contrasts,
  min_reps = 2,
  require_both_groups = FALSE,
  grouping_column = group
)
# Create summary_df from filtered_proteins_per_contrast
summary_df <- data.frame(
  contrast = names(filtered_DALists$filtered_proteins_per_contrast),
  filtered_proteins = sapply(filtered_DALists$filtered_proteins_per_contrast, length)
)
write.csv(summary_df, "filtered_protein_counts.csv", row.names = FALSE)

# GTest Imputation using min_val
source("C:/Users/sbyrum/OneDrive - St. Jude Children's Research Hospital/Documents/Development/proteoDAstjude/HPC_diann_proteoda2/R/Gtest_impute_v3.R")
filtered_DAList_Gtest2 <- impute_missing_by_gtest(filtered_DALists, 
                                                 contrast = NULL, 
                                                 grouping_column = "group",
                                                 p_threshold = p.val)

# test imput function by filtering so 2 reps have a value in 1 group (per contrast)
filtered_DAList_Gtest3 <- impute_missing_by_gtest(filtered_DALists2, 
                                                  contrast = NULL, 
                                                  grouping_column = "group",
                                                  p_threshold = p.val)


filtered_DAList_Gtest4 <- impute_missing_by_gtest(filtered_samples, 
                                                  contrast = contrasts, 
                                                  grouping_column = "group",
                                                  p_threshold = p.val)

stopifnot(all(sapply(names(filtered_DAList_Gtest2$data_per_contrast), function(ct) {
  identical(rownames(filtered_DAList_Gtest2$data_per_contrast[[ct]]), filtered_DAList_Gtest2$filtered_proteins_per_contrast[[ct]])
})))

# if using the Rshiny raw norm values from DiaNN
#normlog2 <- normalize_data(filtered_proteins,
#                           norm_method = "log2")



# using the log2 norm output from diann_quan

norm <- filtered_DAList_Gtest
norm <- filtered_DAList_Gtest2 # fixed data in imputation
#norm <- filtered_proteins
#norm <- filtered_samples
norm$tags$norm_method <- "diann_quan"

############
## RUN LIMMA
########

# Run limma model for each contrast, computes SDs and z-scores
results <- run_filtered_limma_analysis(
  DAList = norm,
  design_formula = design,
  pval_thresh = p.val,
  lfc_thresh = logFC,
  adj_method = "BH",
  binsize = bin_size,
  plot_movingSD = TRUE
  # binsize = "auto",
  # binsize_range = c(100, 250, 500, 1000, 1500)
)


###########
## Standard deviation for biological replicates
###### Calculate Rolling Standard Deviation
# if passes G-test, then fill in missing values if significant


##############
# WRITE DATA TABLES TO FILE
# includes an overall Excel file with all comparisons added
####

# Build statlist with full protein list across all contrasts
# Fills in missing stats with blanks or NAs in order to build a full data table
statlist <- build_statlist(
  DAList = results,
  stat_cols = stat_cols)
 

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

saveRDS(results, "results_NA.rds")

results <- missing_to_zero(results)

write_limma_plots(results,
                  grouping_column = group,
                  table_columns = DA_table_cols,
                  title_column = DA_title_col,
                  output_dir = DA_dir,
                  tmp_subdir = tmp_subdir,
                  overwrite = T,
                  height = 1000,     # sets the container size in the html file
                  width = 1000)      # sets the container size in the html file


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
