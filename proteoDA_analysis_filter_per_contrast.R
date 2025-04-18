library(devtools)
devtools::install_github("ByrumLab/proteoDAstjude",
                         auth_token="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", force=TRUE)


#library(devtools) # nocov start
library(proteoDAstjude)
#library(tidyverse)
library(yaml)
#library(tidyverse)

source("proteoDA_params.R")
source("R/compute_movingSD_zscores.R")
source("R/filter_proteins.R")
source("R/limma_analysis.R")
source("R/limma_tables.R")
source("R/s3_class_v2.R")
source("R/utils.R")
source("R/SD_within_group.R")
source("R/missing_data.R")


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
library(ggplot2)

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

## subset into different comparisons in order to filter proteins by comparisons

#filtered_samples <- filter_samples(filtered_samples, group != "K48ac_6mon") # removes K48ac_24mon
#filtered_samples <- filter_samples(filtered_samples, group != "K48ac_24mon") # 
#filtered_samples <- filter_samples(filtered_samples, group != "K48ac_30mon") # 
#filtered_samples <- filter_samples(filtered_samples, group != "Ub_6mon") # 
#filtered_samples <- filter_samples(filtered_samples, group != "Ub_24mon") # 
#filtered_samples <- filter_samples(filtered_samples, group != "Ub_30mon") # 
#filtered_samples <- filter_samples(filtered_samples, group != "IgG_24mon") # 

table(filtered_samples$metadata$group)

# ### in general require 2/3 of biological replicates to have a value in 1 sample group
# filtered_proteins <- filter_proteins_by_group(filtered_samples,
#                                               min_reps = filt_min_reps,
#                                               min_groups = filt_min_groups,
#                                               grouping_column = group)


# Filter proteins separately for each contrast (both groups must meet threshold)
# creates multiple DAList objects for analysis. 

filtered_DALists <- filter_proteins_per_contrast(
  DAList = filtered_samples,
  contrasts_file = contrasts,
  min_reps = 3,
  require_both_groups = TRUE,
  grouping_column = "group"
)
summary_df <- attr(filtered_DALists, "retention_summary")
write.csv(summary_df, "filtered_protein_counts.csv", row.names = FALSE)

# You now have:
# filtered_DALists[["GroupA_vs_GroupB"]]  # A filtered DAList for this contrast
# filtered_DALists[["GroupC_vs_GroupD"]]  # A separate filtered DAList for another contrast


# if using the Rshiny raw norm values from DiaNN
#normlog2 <- normalize_data(filtered_proteins,
#                           norm_method = "log2")

# using the log2 norm output from diann_quan

norm <- filtered_DALists
#norm <- filtered_proteins
#norm <- filtered_samples
norm$tags$norm_method <- "diann_quan"

############
## RUN LIMMA
########
# Load or prepare your filtered DAList with filtered_proteins_per_contrast already included
# Typically from something like:
# filtered_DAList <- filter_proteins_per_contrast(DAList, contrasts_file, ...)

# Run filtered limma model for each contrast, computes SDs and z-scores
source("R/add_design.R")
source("R/add_contrasts.R")
results <- run_filtered_limma_analysis(
  DAList = norm,
  design_formula = design,
  pval_thresh = 0.05,
  lfc_thresh = 1,
  adj_method = "BH"
)


###########
## Standard deviation for biological replicates
###### Calculate Rolling Standard Deviation
# if passes G-test, then fill in missing values if significant

## plot Volcano significance using moving_sd parameter 
#source("R/compute_movingSD_zscores.R")


# source("R/compute_movingSD_zscores_v3.R")
# results_SDs <- compute_movingSD_zscores(results, 
#                                     binsize = "auto",
#                                     binsize_range = c(50, 100, 200, 400, 500), 
#                                     plot = TRUE)

### Calculate SD_within_group
# Example inputs
#count_data <- your_result$data      # From DAList
#metadata <- your_result$metadata    # Also from DAList

# Run function
#sd_report <- analyze_within_group_SDs(count_data, metadata)

#####in development
source("R/SD_within_group.R")
sd_report <- analyze_within_group_SDs(results$data, 
                                      results$metadata,
                                      group_column = "group",
                                      output_dir = QC_dir,
                                      outlier_method = c("IQR", "z-score"),
                                      z_thresh = 2.5)

# needs work
#source("R/ggplot_binsize_cv.R")
#cv_df <- plot_binsize_cv_curve(raw, binsize_range = c(50, 100, 200, 400, 500, 800))
#cv_df

##############
# WRITE DATA TABLES TO FILE
# includes an overall Excel file with all comparisons added
####
#source("R/s3_class.R")
#source("R/utils.R")

source("R/limma_tables_v5.R", echo = TRUE, max.deparse.length = Inf)
#exists("summarize_contrast_DA")

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

saveRDS(results, "results_NA.rds")
source("R/s3_class_v2.R")

results <- missing_to_zero(results)

### for testing
### limma report is in relation to proteoDA and not proteoDAstjude due to extra flag: line 255 in limma_report.R
source("R/limma_report.R")
#library(ggplot2)
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
