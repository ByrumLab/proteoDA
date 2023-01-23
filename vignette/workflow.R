devtools::load_all()

# Load in data on a bunch of files --------------------------------

?read_DIA_data
#data <- read_DIA_data("for_testing/Samples Report of AvarittNL_081622.csv")

input <- read.csv("vignette/DIA_data.csv")
data <- input[,5:21]
anno <- input[,1:4]
meta <- read.csv("vignette/metafile.csv")
row.names(meta) <- meta$data_column_name

DA <- DAList(data,
               anno,
               meta,
               design = NULL,
               eBayes_fit = NULL,
               results = NULL,
               tags = NULL)



#data <- add_metadata(data, "for_testing/AvarittNL_120522_metafile_DIA.csv")

#data <- add_metadata(data, "vignette/metafile.csv")


# filter samples --------------------------------------------------------------
sub_data <- filter_samples(DA, group != "Pool")

# filter proteins ---------------------------------------------------------

filtered_data <- filter_proteins_by_group(sub_data,
                                          min_reps = 3,
                                          min_groups = 1,
                                          grouping_column = "group")

# Normalization report ----------------------------------------------------

write_norm_report(filtered_data,
                  grouping_column = "group",
                  output_dir = "01_QC_report",
                  filename = "proteiNorm.pdf",
                  overwrite = T,
                  suppress_zoom_legend = FALSE,
                  use_ggrastr = FALSE)

norm <- normalize_data(filtered_data,
                       norm_method = "cycloess")

write_qc_report(norm,
                color_column = "group",
                label_column = NULL,
                output_dir = "01_QC_report",
                filename = "QC_report.pdf",
                overwrite = T,
                top_proteins = 500,        # number of most variable proteins
                standardize = TRUE,
                pca_axes = c(1,2),          # first 2 PCs
                dist_metric = "euclidean",  # stats::dist for options
                clust_method = "complete",  # stats::hclust for options
                show_all_proteins = F)  # only those with missing data

norm <- add_design(norm,
                   design_formula = "~ 0 + group")


norm <- add_contrasts(norm,
                      contrasts_file = "vignette/contrasts.csv")

# Run the analysis --------------------------------------------------------
# Splitting up functionality
# First, fit the model
fit <- fit_limma_model(norm)

# Extract results ---------------------------------------------------------
results <- extract_DA_results(fit,
                              pval_thresh = 0.055,
                              lfc_thresh = 1,
                              adj_method = "BH",
                              extract_intercept = F)

# Write results -----------------------------------------------------------
write_limma_tables(results,
                   output_dir = "02_DA_results",
                   overwrite = T,
                   contrasts_subdir = NULL,
                   summary_csv=NULL,
                   combined_file_csv = NULL,
                   spreadsheet_xlsx = NULL,
                   add_filter = T)

write_limma_plots(results,
                  grouping_column = "group",
                  table_columns = c("uniprot_id","Protein.0me"),
                  title_column = "uniprot_id",
                  output_dir = "02_DA_results",
                  tmp_subdir = "tmp",
                  height = 1000,
                  width = 1000)
