# updated proteoDAstjude package May 27, 2025 - fixed logFC_zscore
library(proteoDAstjude)
library(yaml)

# Load required packages to convert DIANN quan to proteoDA format
library(readr)
library(stringr)

# for testing locally 
devtools::load_all()
#######
# load project parameters 
source("proteoDA_params.R")

############
## File format conversion - DIANN quan to proteoDA
###########

# Read the CSV files for testing
#input_data <- read.csv("data/uni_prot_quan_rmNA_norm.csv") # DIANN
input_data <- read.csv("data/20250910_Durbigrp_071825_Protein_proteoDA_input.csv") # spectronaut with all samples

# argument passed from bash 
# args <- commandArgs(trailingOnly = TRUE)
# diann_quan <- args[1]
# input_data <- read.csv(diann_quan)

sample_metadata <- read.csv(metadata)

#################
# Step 1: Create 'uniprot_id' by extracting the first ID from 'Protein.Group'
# Create 'uniprot_id' as the first column
##################
#uniprot_split <- strsplit(input_data$Protein.Group, ";")
#uniprot_id <- sapply(uniprot_split, `[`, 1)
#input_data <- cbind(uniprot_id = uniprot_id, input_data)

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

# === Step 4: Select relevant columns === WHEN the AVerage.Intensity column is present
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
#### NO Avg Int column to find the start column of the Samples in the data
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
filtered_proteins <- filter_proteins_by_group(filtered_samples,
                                              min_reps = filt_min_reps,
                                              min_groups = filt_min_groups,
                                              grouping_column = group)

# Keeping only protein entries with non-missing intensity in at least 2 samples in at least 1 group
# Filtered 174 entries from the dataset leaving 1874 entries for analysis

# ### OPTION 2
# Filter proteins separately for each contrast (both groups must meet threshold if require_both_groups = TRUE)
# filtered_DALists <- filter_proteins_per_contrast(
#   DAList = filtered_samples,
#   contrasts_file = contrasts,
#   min_reps = filt_min_reps,
#   require_both_groups = require_both_groups,
#   grouping_column = group
# )
# # Create summary_df from filtered_proteins_per_contrast
# summary_df <- data.frame(
#   contrast = names(filtered_DALists$filtered_proteins_per_contrast),
#   filtered_proteins = sapply(filtered_DALists$filtered_proteins_per_contrast, length)
# )
# 
# write.csv(summary_df, "filtered_protein_counts.csv", row.names = FALSE)

# GTest Imputation using min_val
# filtered_DAList_Gtest <- impute_missing_by_gtest(filtered_DALists, 
#                                                  contrast = NULL, 
#                                                  grouping_column = "group",
#                                                  p_threshold = p.val)
# #######
# ## if skipping filter proteins per contrasts then use a contrasts.csv file
# filtered_DAList_Gtest4 <- impute_missing_by_gtest(filtered_samples, 
#                                                   contrast = contrasts, 
#                                                   grouping_column = "group",
#                                                   p_threshold = p.val)

# stopifnot(all(sapply(names(filtered_DAList_Gtest2$data_per_contrast), function(ct) {
#   identical(rownames(filtered_DAList_Gtest2$data_per_contrast[[ct]]), filtered_DAList_Gtest2$filtered_proteins_per_contrast[[ct]])
# })))

# if using the Rshiny raw norm values from DiaNN
#normlog2 <- normalize_data(filtered_proteins,
#                           norm_method = "log2")



#### using the log2 norm output from diann_quan
#norm <- filtered_DAList_Gtest
#norm <- filtered_proteins
#norm <- filtered_samples
#norm$tags$norm_method <- "diann_quan"

############
## NORMALIZATION
############
write_norm_report(filtered_proteins,
                  grouping_column = "group",
                  output_dir = "01_QC_report",
                  filename = "normalization.pdf",
                  overwrite = T,
                  suppress_zoom_legend = FALSE,
                  use_ggrastr = FALSE)

# methods = log2, median, mean, vsn, quantile, cycloess, rlr, gi
norm <- normalize_data(filtered_proteins,
                             norm_method = "cycloess")

norm_med <- normalize_data(filtered_proteins, norm_method = "median")
norm_vsn <- normalize_data(filtered_proteins, norm_method = "vsn")

################
## Perseus normal distribution imputation Method
## Perseus requires log2 transformed data
## added function to cycloess to run as independent groups and not globally
## The min_obs_per_sample guard makes the function stable for sparse samples by borrowing pooled stats.
## This matches the standard Perseus logic without any condition-dependent behavior.
####
#source("R/normalization_v3.R")

# Example
impute_data <- filtered_proteins

# norm_ind$data <- log2Norm(norm_ind$data)                 # your existing helper
# norm_ind$data <- perseus_impute(norm_ind$data, 
#                               robust = TRUE, # If you prefer classic mean/SD instead of robust median/MAD, set robust = FALSE.
#                               seed=1)   # impute into low tail
# 
# norm_ind_cyc <- normalize_data(norm_ind,
#                          norm_method = "cycloess",
#                          groups = norm_ind$metadata$group)

#### Version normalization_v3 ----------------------
# 1) Start with log2 data (with NAs)
X_log2 <- log2Norm(impute_data$data)

# 2) Keep a copy for plotting comparison
X_before <- X_log2

# 3) Impute
X_after <- perseus_impute(X_log2, shift = 1.8, width = 0.3, robust = TRUE, seed = 1)

write.csv(X_after,"Log2_Perseus_Imputed_Values.csv")

# 4) Visualize (all samples, or pick a subset by name or index)
p <- plot_perseus_imputation(X_before, X_after, bins = 60, facet_ncol = 3)
print(p)

ggsave(file = file.path(QC_dir,"Perseus_Histogram_Imputation.png"), # use .extension such as .png, .pdf, .jpeg
       plot = p, 
       path = NULL, 
       width = 10, 
       height = 10, 
       units = c("in"), 
       dpi = 600)

# 5) Proceed to within-group cyclic loess, etc.
impute_data$data <- X_after

#source("R/s3_class.R")

# norm_ind <- normalize_data(impute_data, norm_method = "cycloess",
#                            input_is_log2 = TRUE,                     # <— prevent double log
#                            groups = impute_data$metadata$group)

# norm_cyc <- filtered_proteins
# norm_cyc$data <- X_after

norm_cyc <- normalize_data(impute_data, norm_method = "cycloess",
                           input_is_log2 = TRUE,                     # <— prevent double log
                           groups = NULL)

# norm_med <- filtered_proteins
# norm_med$data <- X_after
# 
# norm_med <- normalize_data(norm_med, norm_method = "median",
#                            input_is_log2 = TRUE)                      # <— prevent double log
#                            
# norm_vsn <- filtered_proteins
# norm_vsn$data <- X_after   # perseus imputed
# norm_vsn <- normalize_data(norm_vsn, norm_method = "vsn", 
#                            input_is_log2 = TRUE)


### OPTION B -- use cycloessNorm directly
# norm_ind$data <- cycloessNorm(
#   logDat = as.matrix(X_after),              # already log2
#   groups = norm_ind$metadata$group
# )
# norm_ind$tags$normalized <- TRUE
# norm_ind$tags$norm_method <- "cycloess"


####----------------------------------
# write_qc_report(norm_ind,
#                 color_column = "group",
#                 label_column = NULL,
#                 output_dir = "01_QC_report",
#                 filename = "QC_report_Perseus_cycloess.pdf",
#                 overwrite = T,
#                 top_proteins = 500,        # number of most variable proteins
#                 standardize = TRUE,
#                 pca_axes = c(1,2),          # first 2 PCs
#                 dist_metric = "euclidean",  # stats::dist for options
#                 clust_method = "complete",  # stats::hclust for options
#                 show_all_proteins = F)      # only proteins with missing data)

write_qc_report(norm_cyc,
                color_column = "group",
                label_column = NULL,
                output_dir = "01_QC_report",
                filename = "QC_report_Perseus_cycloess.pdf",
                overwrite = T,
                top_proteins = 500,        # number of most variable proteins
                standardize = TRUE,
                pca_axes = c(1,2),          # first 2 PCs
                dist_metric = "euclidean",  # stats::dist for options
                clust_method = "complete",  # stats::hclust for options
                show_all_proteins = F)      # only proteins with missing data)

# write_qc_report(norm_med,
#                 color_column = "group",
#                 label_column = NULL,
#                 output_dir = "01_QC_report",
#                 filename = "QC_report_Perseus_median.pdf",
#                 overwrite = T,
#                 top_proteins = 500,        # number of most variable proteins
#                 standardize = TRUE,
#                 pca_axes = c(1,2),          # first 2 PCs
#                 dist_metric = "euclidean",  # stats::dist for options
#                 clust_method = "complete",  # stats::hclust for options
#                 show_all_proteins = F)      # only proteins with missing data)
# 
# write_qc_report(norm_vsn,
#                 color_column = "group",
#                 label_column = NULL,
#                 output_dir = "01_QC_report",
#                 filename = "QC_report_Perseus_vsn.pdf",
#                 overwrite = T,
#                 top_proteins = 500,        # number of most variable proteins
#                 standardize = TRUE,
#                 pca_axes = c(1,2),          # first 2 PCs
#                 dist_metric = "euclidean",  # stats::dist for options
#                 clust_method = "complete",  # stats::hclust for options
#                 show_all_proteins = F)      # only proteins with missing data)
# 
############
## RUN LIMMA
########
# Run limma model for each contrast, computes SDs and z-scores
# results_ind <- run_filtered_limma_analysis(
#   DAList = norm_ind,
#   design_formula = design,
#   contrasts_file = contrasts,
#   pval_thresh = p.val,
#   lfc_thresh = logFC,
#   adj_method = "BH",
#   binsize = bin_size,
#   plot_movingSD = TRUE
#   # binsize = "auto",
#   # binsize_range = c(100, 250, 500, 1000, 1500)
# )

results_cyc <- run_filtered_limma_analysis(
  DAList = norm_cyc,
  design_formula = design,
  contrasts_file = contrasts,
  pval_thresh = p.val,
  lfc_thresh = logFC,
  adj_method = "BH",
  binsize = bin_size,
  plot_movingSD = FALSE
  # binsize = "auto",
  # binsize_range = c(100, 250, 500, 1000, 1500)
)


#source("R/factorial_design.R")
# --- C) Rescale the unscaled average to a true average & rename ---
# (keeps t/P/q/B unchanged; updates contrast matrix & tags if you use the upgraded helper)

lab <- "Avg_Treat_Bio_vs_DMSO_unscaled"

# results_group_ind <- rescale_contrast_logFC(
#   x = results_ind,
#   label = lab,
#   factor = 0.5,
#   new_label = "Avg_Treat_Bio_vs_DMSO"
# )

results_cyc <- rescale_contrast_logFC(
  x = results_cyc,
  label = lab,
  factor = 0.5,
  new_label = "Avg_Treat_Bio_vs_DMSO"
)

# --- D) Build statlist from the SAME object you will write ---
# statlist <- build_statlist(
#   DAList = results_group_ind,
#   stat_cols = stat_cols)

statlist <- build_statlist(
  DAList = results_cyc,
  stat_cols = stat_cols)

# --- E) Write tables and plots (note grouping_column is a STRING) ---
# write_limma_tables(results_group_ind,
#                    output_dir = "Results_Perseus_GroupMeans",
#                    overwrite = T,
#                    contrasts_subdir = NULL,
#                    summary_csv=NULL,
#                    combined_file_csv = NULL,
#                    spreadsheet_xlsx = NULL,
#                    add_filter = T,
#                    color_palette = NULL,
#                    add_contrast_sheets = TRUE,
#                    statlist = statlist)
# 
# write_limma_plots(
#   DAList = results_group_ind,
#   grouping_column = group,
#   table_columns = DA_table_cols,
#   title_column = DA_title_col,
#   height = 1000,
#   width = 1000,
#   output_dir = "Results_Perseus_GroupMeans",
#   overwrite = T,
#   control_proteins = ctrl_proteins, #c("P58004", "Q12766", "MDM2")
#   highlight_by = "uniprot_id")

write_limma_tables(results_cyc,
                   output_dir = "Results_Perseus_CyclicLoess_GroupMeans",
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = T,
                   color_palette = NULL,
                   add_contrast_sheets = TRUE,
                   statlist = statlist)

write_limma_plots(
  DAList = results_cyc,
  grouping_column = group,
  table_columns = DA_table_cols,
  title_column = DA_title_col,
  height = 1000,
  width = 1000,
  output_dir = "Results_Perseus_CyclicLoess_GroupMeans",
  overwrite = T,
  control_proteins = ctrl_proteins, #c("P58004", "Q12766", "MDM2")
  highlight_by = "uniprot_id")

# --- F) Interpretation + classification (uses group-means labels) ---
contrast_map_gm <- c(
  Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",      # scaled
  Interaction_CHLA_vs_SKNF1  = "Interaction_CHLA_vs_SKNF1",
  TreatEffect_CHLA90         = "TreatEffect_CHLA90",
  TreatEffect_SKNF1          = "TreatEffect_SKNF1",
  CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
  CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO"
)

ifi <- interpret_protein_factorial(
  DA_results  = results_cyc,
  protein     = "Q16666",
  protein_col = NULL,
  contrast_map= contrast_map_gm
)

ifi$table
ifi$summary

factab_groupMeans <- classify_factorial_batch(
  DA_results   = results_cyc,
  contrast_map = contrast_map_gm,
  p_thresh     = 0.05,
  lfc_thresh   = 0.5
)
#write.csv(factab_groupMeans, "Results_GroupMeans/classify_Perseus_groupMeans_results_lfcThresh1.csv", row.names = FALSE)
write.csv(factab_groupMeans, "Results_Perseus_CyclicLoess_GroupMeans/classify_Perseus_cyclic_groupMeans_results_lfcThresh0.5.csv", row.names = FALSE)

saveRDS(results_cyc, "Results_Perseus_CyclicLoess_GroupMeans/results_NA_Perseus_Cyc.rds")

results0 <- missing_to_zero(results_group)
saveRDS(results0, "results.rds")

########################
# STOP HERE -------------------------------------------
######
## Extract significant proteins based on avg and interaction models. 

ifi <- interpret_protein_factorial(
  DA_results  = results_group,
  protein     = "Q16666",         # or "IFI16"
  protein_col = NULL,              # we want rowname/annotation resolution
  contrast_map = c(
    Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
    Interaction_CHLA_vs_SKNF1  = "Interaction_CHLA_vs_SKNF1",
    TreatEffect_CHLA90         = "TreatEffect_CHLA90",
    TreatEffect_SKNF1          = "TreatEffect_SKNF1",
    CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
    CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO"))

ifi$summary
ifi$table

###############----------
# assuming your DAList with limma results is `results0`

factab_groupMeans <- classify_factorial_batch(results_group,
                                              contrast_map = c(
                                                Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
                                                Interaction_CHLA_vs_SKNF1  = "Interaction_CHLA_vs_SKNF1",
                                                CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
                                                CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO" ),
                                              p_thresh = 0.05,
                                              lfc_thresh = 1)
subset(factab_groupMeans, id == "Q16666")  # IFI16 should show WTPreferential_Tx

write.csv(factab_groupMeans,"classify_groupMeans_results.csv")


######### ----- Group-means model ----------Testing

design_group = ~0 + group

### Build contrasts with Group-means model.design =~0 + group
writeLines(c(
  
  "TreatEffect_SKNF1= SKNF1_Bio-SKNF1_DMSO",
  "TreatEffect_CHLA90= CHLA90_Bio-CHLA90_DMSO",
  "CHLA90_Bio_vs_SKNF1_Bio= CHLA90_Bio-SKNF1_Bio",
  "CHLA90_DMSO_vs_SKNF1_DMSO= CHLA90_DMSO-SKNF1_DMSO",
  "Interaction_CHLA_vs_SKNF1 = (CHLA90_Bio - CHLA90_DMSO) - (SKNF1_Bio - SKNF1_DMSO)",
  "Avg_Treat_Bio_vs_DMSO_unscaled = (CHLA90_Bio - CHLA90_DMSO) + (SKNF1_Bio - SKNF1_DMSO)"
), "data/contrasts_groupmeans.csv")


# Run limma model for each contrast, computes SDs and z-scores
results_group <- run_filtered_limma_analysis(
  DAList = norm,
  design_formula = design_group,
  contrasts_file = "data/contrasts_groupmeans.csv",
  pval_thresh = p.val,
  lfc_thresh = logFC,
  adj_method = "BH",
  binsize = bin_size,
  plot_movingSD = TRUE
  # binsize = "auto",
  # binsize_range = c(100, 250, 500, 1000, 1500)
)

source("R/factorial_design.R")

lab <- "Avg_Treat_Bio_vs_DMSO_unscaled"
results_group <- rescale_contrast_logFC(
  x = results_group,
  label = lab,
  factor = 0.5,
  new_label = "Avg_Treat_Bio_vs_DMSO"
)

statlist <- build_statlist(
  DAList = results,
  stat_cols = stat_cols)

write_limma_tables(results_group,
                   output_dir = "Results_GroupMeans",
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = T,
                   color_palette = NULL,
                   add_contrast_sheets = TRUE,
                   statlist = statlist)

write_limma_plots(
  DAList = results_group,
  grouping_column = group,
  table_columns = DA_table_cols,
  title_column = DA_title_col,
  height = 1000,
  width = 1000,
  output_dir = "Results_GroupMeans",
  overwrite = T,
  control_proteins = ctrl_proteins, #c("P58004", "Q12766", "MDM2")
  highlight_by = "uniprot_id")

ifi <- interpret_protein_factorial(
  DA_results  = results_group,
  protein     = "Q16666",         # or "IFI16"
  protein_col = NULL,              # we want rowname/annotation resolution
  contrast_map = c(
    Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
    Interaction_CHLA_vs_SKNF1  = "Interaction_CHLA_vs_SKNF1",
    TreatEffect_CHLA90         = "TreatEffect_CHLA90",
    TreatEffect_SKNF1          = "TreatEffect_SKNF1",
    CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
    CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO"))

ifi$summary
ifi$table

###############----------
# assuming your DAList with limma results is `results0`

factab_groupMeans <- classify_factorial_batch(results_group,
                                   contrast_map = c(
                                     Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
                                     Interaction_CHLA_vs_SKNF1  = "Interaction_CHLA_vs_SKNF1",
                                     CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
                                     CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO" ),
                                   p_thresh = 0.05,
                                   lfc_thresh = 1)
subset(factab_groupMeans, id == "Q16666")  # IFI16 should show WTPreferential_Tx

write.csv(factab_groupMeans,"classify_groupMeans_results.csv")

######### - Full FACTORIAL MODEL 
### 2.  Build contrasts with Full Factorial Model. design =~cell * treatment
# these rely on the coefficient parameterization, 
# which can be unintuitive when some groups are all NA for a protein. 
# The group-mean design avoids that: every contrast is just one group mean minus another.

writeLines(c(
  # (CHLA90_Bio) - (SKNF1_Bio)
  # = (Intercept + Bio) - (Intercept + SKNF1 + Bio + SKNF1.treatmentBio)
  "CHLA90_Bio_vs_SKNF1_Bio = (Intercept + Bio) - (Intercept + SKNF1 + Bio + SKNF1.treatmentBio)",
  
  # (CHLA90_DMSO) - (SKNF1_DMSO)
  # = (Intercept) - (Intercept + SKNF1)
  "CHLA90_DMSO_vs_SKNF1_DMSO = Intercept - (Intercept + SKNF1)",
  "TreatEffect_CHLA90 = Bio",
  "TreatEffect_SKNF1  = Bio + SKNF1.treatmentBio",
  
  # (Bio−DMSO in CHLA90) − (Bio−DMSO in SKNF1)
  # = (Bio) − (Bio + SKNF1.treatmentBio)
  "Diff_in_effect_CHLA_SKNF1 = Bio - (Bio + SKNF1.treatmentBio)",
  
  # Average treatment effect across cells: ((CHLA90_Bio - CHLA90_DMSO) + (SKNF1_Bio - SKNF1_DMS0))/2
  # CHLA90: Bio
  # SKNF1:  Bio + SKNF1.treatmentBio
  # mean = Bio + 0.5*SKNF1.treatmentBio
  # Numerator of the average: (CHLA90_Bio + SKNF1_Bio) - (CHLA90_DMSO + SKNF1_DMSO)
  # Written in terms of your design columns: Intercept, SKNF1, Bio, SKNF1.treatmentBio
  "Avg_Treat_Bio_vs_DMSO_unscaled = (Intercept + Bio) + (Intercept + SKNF1 + Bio + SKNF1.treatmentBio) - (Intercept) - (Intercept + SKNF1)"
    # evaluates to 2*Bio + 1*SKNF1.treatmentBio, which is exactly 2 × the desired average (Bio + 0.5*SKNF1.treatmentBio)
), "data/contrasts_factorial.csv")

# sanity check:
readLines("data/contrasts_factorial.csv")
# # should print the four clean lines above, with spaces and backticks
# x <- readLines("data/contrasts_factorial.csv", warn = FALSE)
# # remove BOM and non-breaking spaces / non-ASCII
# x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
# x <- gsub("[\uFEFF\u00A0]", " ", x)   # BOM / NBSP -> space
# x <- trimws(x)
# drop blank lines
# x <- x[nzchar(x)]
# writeLines(x, "data/contrasts_factorial.csv")
# readLines("data/contrasts_factorial.csv")
#1) Build a stable factorial design
#Set baselines explicitly so the columns are predictable. 
#Use DMSO as the treatment baseline so the treatment coefficient is “Bio vs DMSO”.

meta <- norm$metadata
meta$cell      <- relevel(factor(meta$cell), ref = "CHLA90")
meta$treatment <- relevel(factor(meta$treatment), ref = "DMSO")
# Optional: add a group column for plotting
meta$group <- interaction(meta$cell, meta$treatment, sep = "_")
# Replace metadata inside DAList
norm$metadata <- meta

# 2) Full factorial with intercept, no intercept version is not full rank 
design_factorial <- ~ cell * treatment      # equivalent to ~ 1 + cell + treatment + cell:treatment


results_full <- run_filtered_limma_analysis(
  DAList = norm,
  design_formula = design_factorial,
  contrasts_file = "data/contrasts_factorial.csv",
  pval_thresh = p.val,
  lfc_thresh = logFC,
  adj_method = "BH",
  binsize = bin_size,
  plot_movingSD = TRUE
  # binsize = "auto",
  # binsize_range = c(100, 250, 500, 1000, 1500)
)
#sanity check
#mm <- model.matrix(~0 + cell:treatment, data = norm$metadata)
#colnames(mm)
# expect things like: "cellCHLA90:treatmentBio", "cellSKNF1:treatmentBio", ...

### Fix scale of Average expression 
#  Only logFC needs scaling. Limma’s t-stats and p-values are invariant to multiplying a contrast by a constant (SE scales by the same constant), so you don’t need to touch them.
#  If you also stash any per-contrast z-scores derived from logFC (e.g., your “movingSD z-scores”), leave them as-is—they’re already standardized.

source("R/factorial_design.R")

lab <- "Avg_Treat_Bio_vs_DMSO_unscaled"
results_full <- rescale_contrast_logFC(
  x = results_full,
  label = lab,
  factor = 0.5,
  new_label = "Avg_Treat_Bio_vs_DMSO"
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


write_limma_tables(results_full,
                   output_dir = "Results_FullFactorial",
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

results0 <- missing_to_zero(results)

library("ggplot2")
#######

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
    highlight_by = "uniprot_id")

write_limma_plots(
  DAList = results_full,
  grouping_column = group,
  table_columns = DA_table_cols,
  title_column = DA_title_col,
  height = 1000,
  width = 1000,
  output_dir = "Results_FullFactorial",
  overwrite = T,
  control_proteins = ctrl_proteins, #c("P58004", "Q12766", "MDM2")
  highlight_by = "uniprot_id")

## example usage
# write_limma_plots(
#   DAList = dalist,
#   grouping_column = "treatment",
#   table_columns = c("Gene.names", "description"),
#   title_column = "Gene.names",
#   output_dir = "DA_plots",
#   control_proteins = c("P58004", "Q12766"),
#   highlight_by = "uniprot_id"
# )

# ##### FINAL DAList()
# save for QC plots
saveRDS(results, "results.rds")


############
# FULL FACTORIAL DESIGN

source("R/factorial_design.R")
# Example: results object returned by run_filtered_limma_analysis(...)
# out <- interpret_protein_factorial(
#   DA_results = results_full,     # your DAList-like results holder
#   protein    = "Q16666",         # Uniprot or whatever your row IDs are
#   protein_col = NULL,            # set to "Protein" or similar if IDs are in a column
#   alpha = 0.05
# )
# 
# out$table      # tidy per-contrast stats for this protein
# out$summary    # one-line interpretation

ifi <- interpret_protein_factorial(
  DA_results  = results_full,
  protein     = "Q16666",         # or "IFI16"
  protein_col = NULL,              # we want rowname/annotation resolution
  contrast_map = c(
    Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
    Interaction_CHLA_vs_SKNF1  = "Diff_in_effect_CHLA_SKNF1",
    TreatEffect_CHLA90         = "TreatEffect_CHLA90",
    TreatEffect_SKNF1          = "TreatEffect_SKNF1",
    CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
    CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO"))

ifi$summary
ifi$table

###############----------
# assuming your DAList with limma results is `results0`

factab <- classify_factorial_batch(results_full,
                                   contrast_map = c(
                                     Avg_Treat_Bio_vs_DMSO      = "Avg_Treat_Bio_vs_DMSO",
                                     Interaction_CHLA_vs_SKNF1  = "Diff_in_effect_CHLA_SKNF1",
                                     TreatEffect_CHLA90         = "TreatEffect_CHLA90",
                                     TreatEffect_SKNF1          = "TreatEffect_SKNF1",
                                     CHLA90_Bio_vs_SKNF1_Bio    = "CHLA90_Bio_vs_SKNF1_Bio",
                                     CHLA90_DMSO_vs_SKNF1_DMSO  = "CHLA90_DMSO_vs_SKNF1_DMSO" ),
                                   p_thresh = 0.05,
                                   lfc_thresh = 1)
subset(factab, id == "Q16666")  # IFI16 should show WTPreferential_Tx

write.csv(factab,"classify_factorial_results.csv")

################-------------------------
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
