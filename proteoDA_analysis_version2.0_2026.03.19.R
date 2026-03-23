## proteoDA 2.0.0 Recommended Analysis Workflow
# Stephanie D. Byrum
# March 19, 2026

#This workflow
# 1. inputs raw MS2 intensities from DIANN search results
# 2. Adds the required "uniprot_id" column to the data
# 3. finds any blank, NA, NaN sample intentisties and converts to zero 
#      so we know specifically what values we are using for missing data. 
#      Then we can convert to NA using zero_to_missing() function
# 3. splits the input file into annotation and sample data 
# 4. renames the mass spec file names from the DIANN search result sample column names to
#    the sample names in the Sample_metadata.csv so they are less than 30 characters and are R friendly
# 5. check alignment of data sample columns and the Sample_metadata.csv using the new align_data_metadata() function
# 6. Create the raw DAList object that will be used for the rest of the pipeline
# 7. Filter samples and/or proteins
# 8. Evaluate Normalization methods and Normalize
# 9. Optional run missing value imputation 
# 10. Run limma analysis and extract results
# 11. save rds objects for the QC_report.md and PowerPoint reports

###==========
## INSTALLATION
## install version 2.0 in a separate R library
##============
dir.create("~/R/proteoDA_test_lib", showWarnings = FALSE, recursive = TRUE)
.libPaths(c("~/R/proteoDA_test_lib", .libPaths()))
.libPaths()
remotes::install_local(".", force = TRUE, build_vignettes = TRUE)
.rs.restartR()
##==============

##==============
## Load libraries 
##=============
library(readr)
library(stringr)

.libPaths(c("~/R/proteoDA_test_lib", .libPaths()))
library(proteoDA)
# check package version 
packageVersion("proteoDA")
# where is it installed
find.package("proteoDA")
# package citation 
citation("proteoDA")
# run vignettes
browseVignettes("proteoDA")

##===============
# load project parameters 
##===============
source("proteoDA_params.R")

####===========================
## INPUT DATA from bash or parameter file
## fix uniprot_id required column
## replace blanks or Na or NaN with 0
###===========================

# load from parameter file
input_data <- read.csv(input_quan)

# argument passed from bash 
 # args <- commandArgs(trailingOnly = TRUE)
 # diann_quan <- args[1]
 # input_data <- read.csv(diann_quan)

##========================
# ----Step 1: Create 'uniprot_id' by extracting the first ID from 'Protein.Group'------
# Create 'uniprot_id' as the first column
##======================
uniprot_split <- strsplit(input_data$Protein.Group, ";")
uniprot_id <- sapply(uniprot_split, `[`, 1)
input_data <- cbind(uniprot_id = uniprot_id, input_data)

##========================
# ----Step 2: Replace NA and blank entries safely, preserving numeric types-------
##========================
for (col in names(input_data)) {
  if (is.numeric(input_data[[col]])) {
    input_data[[col]][is.na(input_data[[col]])] <- 0
  } else {
    input_data[[col]][is.na(input_data[[col]]) | input_data[[col]] == ""] <- "0"
  }
}

##========================
## Load Sample metadata
## Subset annotation from input file
##========================

# ----Step 3: Load sample metadata sheet  
sample_metadata <- read.csv(metadata)

# Step 4: create anno from input data. 
annotation_data <- input_data[,anno_start:anno_end] # select protein annotation 

#########
## --------check the sample metadata and data have the same samples in the same order
##
## If your columns are currently misassigned because of lexicographic ordering (1, 10, 11, 2…), 
##       the function uses name matching, not position, to realign—so group labels will no longer be misapplied.
## If any sample IDs are missing/duplicated in either table, 
##           it aborts with a clean diagnostic list of the offenders.
## 
#### WATCH ----------
# If some files in metadata$file don’t appear as columns in input_data,
##        the renamer silently ignores them (we subset with intersect). You’ll still get a hard error later if the sets don’t match (from the helper).

## If you prefer metadata order over grouping, set prefer_group_blocks = FALSE 
##            (the helper will then only reorder metadata to match data columns).

## You can optionally call check_data_metadata_match() inside validate_DAList() for light-weight runtime checks without reordering.
############

##========================
# --- Rename data columns: file -> sample (keep this) -------------------------
# Requires metadata columns: 'file' and 'sample'
##========================
rename_map <- setNames(sample_metadata$sample, sample_metadata$file)
# check samples in the input_data
present_files <- intersect(names(rename_map), colnames(input_data))
colnames(input_data)[match(present_files, colnames(input_data))] <- rename_map[present_files]

# Build intensity matrix after renaming; keep your existing sample subset if needed
# (Optional) If you only want samples listed in metadata:
# after renaming file -> sample
dat_ids <- colnames(input_data)
keep_samples <- dat_ids[dat_ids %in% sample_metadata$sample]   # preserves current order
intensity_data <- input_data[, keep_samples, drop = FALSE]
head(intensity_data)
# trim metadata to those samples but do NOT reorder yet
sample_metadata <- subset(sample_metadata, sample %in% keep_samples)

# Row names for metadata should be sample IDs
rownames(sample_metadata) <- sample_metadata$sample

head(intensity_data)
sample_metadata

##========================
# --- Auto-align data & metadata (replaces Step 4) ----------------------------
##========================
aligned <- align_data_and_metadata(
  data                = intensity_data,
  metadata            = sample_metadata,
  sample_col          = "sample",   # your metadata column with sample IDs
  group_col           = "group",    # optional but recommended
  prefer_group_blocks = FALSE,       # keep replicates together + natural sort
  strict              = FALSE
)

# Verify what changed
aligned$changes

# save the aligned data and metadata for loading into the DAList
intensity_data  <- aligned$data
sample_metadata <- aligned$metadata
# verify the data columns and sample metadata names match
head(intensity_data)
head(sample_metadata)

##========================
### Create DAList object
##========================

raw <- DAList(data = intensity_data,
              annotation = annotation_data,
              metadata = sample_metadata,
              design = NULL,
              eBayes_fit = NULL,
              results = NULL,
              tags = NULL)

##========================
# Convert missing values
##========================
filtered_samples <- zero_to_missing(raw)
# verify the number of replicates per group
table(filtered_samples$metadata$group)

##========================
##  Filtering proteins with missing values
## Option 1: filter_proteins_by_group: filters globally using all samples in the data
## Option 2: filter_proteins_per_contrast: filters proteins for each sample group comparison separately
##========================

# ### OPTION 1 
# in general require 2/3 of biological replicates to have a value in 1 sample group
# filtered_proteins <- filter_proteins_by_group(filtered_samples,
#                                               min_reps = filt_min_reps,
#                                               min_groups = filt_min_groups,
#                                               grouping_column = group)

# ### OPTION 2 - NEW
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

write.csv(summary_df, "filtered_protein_counts.csv", row.names = FALSE)

##========================
# GTest Imputation using min_val
## NEW option for imputation of missing values
##========================
# filtered_DAList_Gtest <- impute_missing_by_gtest(filtered_DALists, 
#                                                  contrast = NULL, 
#                                                  grouping_column = "group",
#                                                  p_threshold = p.val)
# #######
## if skipping filter proteins per contrasts then use a contrasts.csv file
# filtered_DAList_Gtest4 <- impute_missing_by_gtest(filtered_samples, 
#                                                   contrast = contrasts, 
#                                                   grouping_column = "group",
#                                                   p_threshold = p.val)

# stopifnot(all(sapply(names(filtered_DAList_Gtest2$data_per_contrast), function(ct) {
#   identical(rownames(filtered_DAList_Gtest2$data_per_contrast[[ct]]), filtered_DAList_Gtest2$filtered_proteins_per_contrast[[ct]])
# })))

##========================
## NORMALIZATION
## methods = "log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi"
##========================

# write normalization report to choose the best method
write_norm_report(
  filtered_DALists,
  grouping_column   = group,
  output_dir        = file.path(QC_dir, "norm_evaluation"),
  filename          = "normalization_MD.pdf",
  overwrite         = TRUE,
  suppress_zoom_legend = FALSE,
  use_ggrastr       = FALSE,
  input_is_log2     = FALSE,    # NEW
  per_contrast      = TRUE,     # NEW
  contrasts         = NULL,     #c("M1Y1_vs_Ref", "M2Y1_vs_Ref"), NEW optional- restrict to selected contrasts
  sample_id_col     = "sample",  # NEW <— set this to the metadata column that equals colnames()
  metrics_csv     = file.path(QC_dir, "norm_evaluation/metrics_PRE.csv"),
  include_MD_plots  = TRUE       # NEW <- turn off MA/MD page
  )

# # legacy function version 1.0.0
# write_norm_report(
#     filtered_DALists,
#     grouping_column = group,
#     output_dir = "QC_report",
#     filename = "normalization.pdf",
#     overwrite = TRUE
# )


# Cyclic loess (limma::normalizeCyclicLoess); 
# if groups supplied, normalize within each group independently; 
# otherwise normalize globally
norm <- normalize_data(
  filtered_DALists,
  norm_method   = "cycloess",
  input_is_log2 = FALSE,           # NEW set TRUE if your per-contrast data are already log2
  groups        = filtered_DALists$metadata$group # NEW for cyclic loess only will perform within groups
)

## evaluate the normalized data - 
# useful for cyclic loess within group evaluation. 
# The write_norm_report evaluates global normalization methods
write_norm_eval_report(
                      norm,            # already normalized per contrast
  norm_label       = "cycloess",
  grouping_column  = group,
  sample_id_col    = "sample",       # this matched earlier
  output_dir       = file.path(QC_dir, "norm_evaluation"),
  filename         = "cyclic_loess_eval.pdf",
  overwrite        = TRUE,
  metrics_csv     = file.path("QC_report", "metrics_POST_cycloess.csv")
)

# Interpretation:
#   Cyclic loess effectively stabilized within-group intensities (lower PMAD/PEV) without flattening the global intensity shift (the bimodality remains).
# COR remains high, so normalization didn’t distort relative rankings.
# This confirms cyclic loess normalization worked as intended for your mixed-species ratio design.

# proteoDA version 1.0 legacy code. 
# norm <- normalize_data(filtered_proteins,
#                            norm_method = "cycloess")
# 

##========================
# ----- If the data is already normalized, just add a tag for the method used
##========================
#norm <- filtered_DALists
#norm$tags$norm_method <- "diann_quan"

##========================
# Perseus IMPUTATION normal distribution
# 1) Start with log2 data (with NAs)
##========================

imputed <- perseus_impute(norm,
                          shift = 1.8, width = 0.3,
                          robust = TRUE,
                          save_before_after = TRUE,
                          seed = 1)

# Generate Plot
## plot currently only works with filtered_proteins_per_contrast data slots
plots <- write_perseus_imputation_plots(
  DAList    = imputed,
  out_dir   = file.path(QC_dir,"PerseusPlots"),  # set to NULL to not save files
  bins      = 30,
  facet_ncol   = 4,
  overlay      = TRUE,
  width        = 7,
  height       = 7,
  dpi          = 300,
  device       = "png"
)

##========================
## RUN LIMMA
# Run limma model for each contrast, computes SDs and z-scores
# runs add_design, add_contrast, lmFit, eBayes
##========================

results <- run_filtered_limma_analysis(
  DAList = norm,
  design_formula = design,
  pval_thresh = p.val,
  lfc_thresh = logFC,
  adj_method = "BH",
  binsize = bin_size,  # new auto detect size based on N
  plot_movingSD = FALSE,
  contrasts_file = contrasts
)

# analysis imputed data
results2 <- run_filtered_limma_analysis(
  DAList = imputed,
  design_formula = design,
  pval_thresh = p.val,
  lfc_thresh = logFC,
  adj_method = "BH",
  binsize = bin_size,
  plot_movingSD = TRUE,
  contrasts_file = contrasts
)

##========================
#---- write movingSD bin size report per contrasts
##========================
results <- write_movingSD_report(
  DAList         = results,
  out_dir        = file.path(DA_dir,"movingSD"),
  binsize        = "auto",    # can use what was selected above or let it choose it again
  contrasts_file = contrasts,
  device         = "png"  # pdf, png, or both
)

# imputed data
results2 <- write_movingSD_report(
  DAList         = results2,
  out_dir        = file.path(DA_dir,"imputed/movingSD"),
  binsize        = "auto",    # can use what was selected above or let it choose it again
  contrasts_file = contrasts,
  device         = "png"  # pdf, png, or both
)

##========================
# WRITE DATA TABLES TO FILE
# includes an overall Excel file with all comparisons added
# Build statlist with full protein list across all contrasts
# Fills in missing stats with blanks or NAs in order to build a full data table
# statistical columns are defined in proteoDA_params.R stat_cols parameter
##========================

statlist <- build_statlist(
  DAList = results,
  stat_cols = stat_cols)

statlist2 <- build_statlist(
  DAList = results2,
  stat_cols = stat_cols) 

##========================
## Generate Excel and CSV tables
##========================
write_limma_tables(results,
                   output_dir = DA_dir,
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,           # summary of significance, used in PowerPoint
                   combined_file_csv = NULL,   # useful for downstream analysis
                   spreadsheet_xlsx = NULL,    # Generate Excel file
                   add_filter = T,             # Excel filter on columns
                   color_palette = NULL,       # change color scheme
                   add_contrast_sheets = TRUE, # for large number of contrasts, turn off
                   statlist = statlist,
                   annot_cols = DA_table_cols   # from params.R or c("uniprot_id", "Genes")
                   )

# imputed data 
write_limma_tables(results2,
                   output_dir = file.path(DA_dir,"imputed"),
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = T,
                   color_palette = NULL,
                   add_contrast_sheets = TRUE,
                   statlist = statlist2,
                   annot_cols = DA_table_cols)

##========================
# WRITE PLOTS
# includes an html for each comparison that is interactive and stand alone
# open html in web browser to explore the data
## table_columns must match the protein annotation information. Displayed in the html Table
## title_column is displayed on the protein intensity plot (can be changed to Genes)
##========================

# ##### FINAL DAList()
# save rds object for QC_report.Rmd
# need data with NA and one with 0 for missing value heatmap and PCA
saveRDS(results, "results_NA.rds")
# convert NAs to 0
results0 <- missing_to_zero(results)
# save for QC plots
saveRDS(results0, "results.rds")

write_limma_plots(
    DAList = results0,
    grouping_column = group,    
    table_columns = DA_table_cols,
    title_column = DA_title_col,
    height = 1000,
    width = 1000,
    output_dir = DA_dir,
    overwrite = T,
    control_proteins = ctrl_proteins, #NEW - highlight volcano static plots - ex: c("P58004", "Q12766", "MDM2")
    highlight_by = "uniprot_id",
    image_formats = c("pdf","png"))

# imputed data 
write_limma_plots(
  DAList = results2,
  grouping_column = group,
  table_columns = DA_table_cols,
  title_column = DA_title_col,
  height = 1000,
  width = 1000,
  output_dir = file.path(DA_dir,"imputed"),
  overwrite = T,
  control_proteins = ctrl_proteins, #c("P58004", "Q12766", "MDM2")
  highlight_by = "uniprot_id",
  image_formats = c("pdf","png"))


##========================
### ---- UPDATE _variables.yml with params for Project Name and authors.
###=======================
library(yaml)
variables <- yaml.load_file("_variables.yml")
variables$title <- project_name
variables$author <- author
variables$author2 <- author2
variables$project <- project_name

# You can also add new variables
# variables$new_var <- 42

write_yaml(variables, "_variables.yml")

# save R session information
sink("sessionInfo.txt")
sessionInfo()
sink()
