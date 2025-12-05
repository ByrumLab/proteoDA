# updated proteoDAstjude package May 27, 2025 - fixed logFC_zscore
library(proteoDAstjude)
library(yaml)

# Load required packages to convert DIANN quan to proteoDA format
library(readr)
library(stringr)

# for testing locally 
devtools::load_all()
devtools::build_vignettes()
devtools::check(document = FALSE) # checks vignettes for R CMD check
#######
# load project parameters 
source("proteoDA_params.R")

############
## INPUT DATA from bash or parameter file
## fix uniprot_id required column
## replace blanks or Na or NaN with 0
###########

# Read the CSV files for testing
#input_data <- read.csv("data/uni_prot_quan_rmNA_norm.csv")
input_data <- read.csv(input_quan)

# argument passed from bash 
 args <- commandArgs(trailingOnly = TRUE)
 diann_quan <- args[1]
 input_data <- read.csv(diann_quan)

# ----Step 1: Create 'uniprot_id' by extracting the first ID from 'Protein.Group'------
# Create 'uniprot_id' as the first column
uniprot_split <- strsplit(input_data$Protein.Group, ";")
uniprot_id <- sapply(uniprot_split, `[`, 1)
input_data <- cbind(uniprot_id = uniprot_id, input_data)

# ----Step 2: Replace NA and blank entries safely, preserving numeric types-------
for (col in names(input_data)) {
  if (is.numeric(input_data[[col]])) {
    input_data[[col]][is.na(input_data[[col]])] <- 0
  } else {
    input_data[[col]][is.na(input_data[[col]]) | input_data[[col]] == ""] <- "0"
  }
}

###########################
## Load Sample metadata
## Subset annotation from input file
###########################

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

#source("R/align_checks.R")
# --- Rename data columns: file -> sample (keep this) -------------------------
# Requires metadata columns: 'file' and 'sample'
rename_map <- setNames(sample_metadata$sample, sample_metadata$file)

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
head(sample_metadata)

# --- Auto-align data & metadata (replaces Step 4) ----------------------------
aligned <- align_data_and_metadata(
  data                = intensity_data,
  metadata            = sample_metadata,
  sample_col          = "sample",   # your metadata column with sample IDs
  group_col           = "group",    # optional but recommended
  prefer_group_blocks = FALSE,       # keep replicates together + natural sort
  strict              = FALSE
)

# (Optional) peek at what changed
#  head(aligned$changes)
aligned$changes

intensity_data  <- aligned$data
sample_metadata <- aligned$metadata
# aligned$changes  # optional: inspect what moved

head(intensity_data)
head(sample_metadata)


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

write.csv(summary_df, "filtered_protein_counts.csv", row.names = FALSE)

# GTest Imputation using min_val
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

# if using the Rshiny raw norm values from DiaNN

##############
## NORMALIZATION
## methods = "log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi"
###############

# write normalization report to choose the best method
write_norm_report(
  filtered_DALists,
  grouping_column   = group,
  output_dir        = "QC_report",
  filename          = "normalization_new.pdf",
  overwrite         = TRUE,
  suppress_zoom_legend = FALSE,
  use_ggrastr       = FALSE,
  input_is_log2     = FALSE,    # NEW
  contrasts         = NULL,     # NEW — restrict to selected contrasts
  sample_id_col     = "sample",  # NEW <— set this to the metadata column that equals colnames()
  metrics_csv     = file.path("QC_report", "metrics_PRE.csv"))

######
## Evaluation ----
# That double-hump pattern fits exactly what you’d expect from a complex mixture 
# where two proteomes (human + yeast) contribute in different proportions across your samples:
#   One peak corresponds mostly to yeast-dominant proteins,
#   the other to human-dominant proteins,
#   and their relative amplitudes shift according to your mixture ratios.
#   
# So here the bimodality isn’t a normalization artifact — it’s a true compositional signature of your biological design.
# You want to preserve this shape through normalization, not flatten it out.
# 
# What to do:
# 1. Prefer within-group cyclic loess (groups=) or VSN, which stabilize variance but keep global intensity differences.
# 2. Avoid quantile or global median/mean normalizations 
#   if your goal is to analyze the relative abundance of human vs yeast proteins,
#   because those would force identical sample distributions and erase the biological mixture pattern.
# 3. After normalization, check per-group density plots: each group should remain unimodal with shifted centers,
#   and the pooled bimodal curve should persist.

# # legacy function
# write_norm_report(
#     filtered_DALists,
#     grouping_column = group,
#     output_dir = "QC_report",
#     filename = "normalization.pdf",
#     overwrite = TRUE
# )

# Normalize each contrast independently; leave DAList$data untouched
# Anywhere your report previously pulled from DAList$data, switch to:
#   If per-contrast path is used: use DAList$data_per_contrast[[contrast]] for the main normalized matrix, 
#   and pull before/after summaries from DAList$normalization_per_contrast[[contrast]]$diagnostics.
#   Else: keep using DAList$data and DAList$tags$normalized.

# Double-hump = distinct underlying protein populations (here, human vs yeast mixtures)
# → keep it! It represents your intended biology, not a normalization issue.


# Cyclic loess (limma::normalizeCyclicLoess); if groups supplied, normalize within each group independently; otherwise normalize globally
norm <- normalize_data(
  filtered_DALists,
  norm_method   = "cycloess",
  input_is_log2 = FALSE,                    # set TRUE if your per-contrast data are already log2
  groups        = filtered_DALists$metadata$group # for cyclic loess only will perform within groups
)


## evaluate the normalized data - not really necessary since already shows results in first evaluation.
## maybe change it to show by pair-wise samples instead of overlay? 
write_norm_eval_report(
  norm,            # already normalized per contrast
  norm_label       = "cycloess",
  grouping_column  = group,
  sample_id_col    = "sample",       # this matched earlier
  output_dir       = "QC_report",
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
#                            norm_method = "log2")
# 

##################
# ----- If the data is already normalized, just add a tag for the method used
###########
#norm <- filtered_DALists
#norm$tags$norm_method <- "diann_quan"

###################
# Perseus IMPUTATION normal distribution
#####
#### Version normalization_v3 ----------------------
# 1) Start with log2 data (with NAs)

#source("R/normalization_v2.R")
#source("s3_class.R")

# use DAList
imputed <- perseus_impute(norm,
                          shift = 1.8, width = 0.3,
                          robust = TRUE,
                          save_before_after = TRUE,
                          seed = 1)


# after: imputed <- perseus_impute(norm, save_before_after = TRUE, seed = 1)
plots <- write_perseus_imputation_plots(
  DAList    = imputed,
  out_dir   = "PerseusPlots",  # set to NULL to not save files
  bins      = 30,
  facet_ncol   = 4,
  overlay      = TRUE,
  width        = 7,
  height       = 7,
  dpi          = 300,
  device       = "png"
)

# To preview one in RStudio:
# print(plots[["ZR_fus_vs_RELA"]])

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
 # binsize = bin_size,
  plot_movingSD = FALSE,
  contrasts_file = contrasts,
  binsize = "auto"
 # binsize_range = c(100, 250, 500, 1000, 1500, 5000)  # removed this from parameters since now chooses based on N
)

#########
#---- write movingSD bin size report per contrasts
#########
results <- write_movingSD_report(
  DAList         = results,
  out_dir        = "vignettes/images/movingSD", # "QC_report/movingSD",
  binsize        = "auto",    # can use what was selected abover or let it choose it again
  contrasts_file = contrasts,
  device         = "png"  # pdf, png, or both
)
###########

results2 <- run_filtered_limma_analysis(
  DAList = imputed,
  design_formula = design,
  pval_thresh = p.val,
  lfc_thresh = logFC,
  adj_method = "BH",
  binsize = bin_size,
  plot_movingSD = TRUE,
  contrasts_file = contrasts
  # binsize = "auto",
  # binsize_range = c(100, 250, 500, 1000, 1500)
)

##############
# WRITE DATA TABLES TO FILE
# includes an overall Excel file with all comparisons added
####

# Build statlist with full protein list across all contrasts
# Fills in missing stats with blanks or NAs in order to build a full data table
statlist <- build_statlist(
  DAList = results,
  stat_cols = stat_cols)

statlist2 <- build_statlist(
  DAList = results2,
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

write_limma_tables(results2,
                   output_dir = "impute_out",
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = T,
                   color_palette = NULL,
                   add_contrast_sheets = TRUE,
                   statlist = statlist2)

##############
# WRITE PLOTS
# includes an html for each comparison that is interactive and stand alone
# open html in web browser to explore the data
#########

## table_columns must match the protein annotation information. Displayed in the html Table
## title_column is displayed on the protein intensity plot (can be changed to Genes)

saveRDS(results, "results_NA.rds")

results0 <- missing_to_zero(results)

write_limma_plots(
    DAList = results0,
    grouping_column = group,
    table_columns = DA_table_cols,
    title_column = DA_title_col,
    height = 1000,
    width = 1000,
    output_dir = DA_dir,
    overwrite = T,
    control_proteins = ctrl_proteins, #c("P58004", "Q12766", "MDM2")
    highlight_by = "uniprot_id",
    image_formats = c("pdf","png"))

write_limma_plots(
  DAList = results2,
  grouping_column = group,
  table_columns = DA_table_cols,
  title_column = DA_title_col,
  height = 1000,
  width = 1000,
  output_dir = "impute_out",
  overwrite = T,
  control_proteins = ctrl_proteins, #c("P58004", "Q12766", "MDM2")
  highlight_by = "uniprot_id",
  image_formats = c("pdf","png"))

# ##### FINAL DAList()
# save for QC plots
saveRDS(results0, "results.rds")

##########
### ---- UPDATE _variables.yml with params for Project Name and authors.
###########

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
