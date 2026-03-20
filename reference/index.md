# Package index

## All functions

- [`DAList()`](https://byrumlab.github.io/proteoDA/reference/DAList.md)
  : Create a DAList
- [`GTest_impute()`](https://byrumlab.github.io/proteoDA/reference/GTest_impute.md)
  : G-test for Detection Imbalance and Conditional Imputation
- [`Lou_HF_DIANN_uni_prot_quan`](https://byrumlab.github.io/proteoDA/reference/Lou_HF_DIANN_uni_prot_quan.md)
  : Example DIANN protein-level quantification for Lou HF dataset
- [`Lou_HF_sample_metadata`](https://byrumlab.github.io/proteoDA/reference/Lou_HF_sample_metadata.md)
  : Sample metadata for Lou HF dataset
- [`Lou_contrasts`](https://byrumlab.github.io/proteoDA/reference/Lou_contrasts.md)
  : Example contrasts table for Lou HF dataset
- [`add_contrasts()`](https://byrumlab.github.io/proteoDA/reference/add_contrasts.md)
  : Prepare limma contrasts matrix
- [`add_design()`](https://byrumlab.github.io/proteoDA/reference/add_design.md)
  : Prepare limma model design matrix
- [`align_data_and_metadata()`](https://byrumlab.github.io/proteoDA/reference/align_data_and_metadata.md)
  : Align and validate sample correspondence between data and metadata
- [`analyze_within_group_SDs()`](https://byrumlab.github.io/proteoDA/reference/analyze_within_group_SDs.md)
  : Analyze Within-Group Log2FC Standard Deviations
- [`build_statlist()`](https://byrumlab.github.io/proteoDA/reference/build_statlist.md)
  : Build a cleaned statlist for export to Excel
- [`check_DA_perc()`](https://byrumlab.github.io/proteoDA/reference/check_DA_perc.md)
  : Check percentage of DA genes
- [`check_data()`](https://byrumlab.github.io/proteoDA/reference/check_data.md)
  : Check if an object is a data.frame or matrix
- [`check_data_metadata_match()`](https://byrumlab.github.io/proteoDA/reference/check_data_metadata_match.md)
  : Check data–metadata alignment consistency
- [`check_dup()`](https://byrumlab.github.io/proteoDA/reference/check_dup.md)
  : Identify Duplicate Values
- [`check_file()`](https://byrumlab.github.io/proteoDA/reference/check_file.md)
  : Check if a file exists
- [`check_int()`](https://byrumlab.github.io/proteoDA/reference/check_int.md)
  : Identify Positive Integer Values
- [`check_logical()`](https://byrumlab.github.io/proteoDA/reference/check_logical.md)
  : Check if a Value is a TRUE or FALSE Logical of Length 1
- [`check_long()`](https://byrumlab.github.io/proteoDA/reference/check_long.md)
  : Identify Values Above a Specified Character Length
- [`check_num()`](https://byrumlab.github.io/proteoDA/reference/check_num.md)
  : Identify Positive Numeric Values
- [`check_string()`](https://byrumlab.github.io/proteoDA/reference/check_string.md)
  : Validate Single Character String
- [`check_syntax()`](https://byrumlab.github.io/proteoDA/reference/check_syntax.md)
  : Check if a vector follows R syntax rules
- [`check_vec()`](https://byrumlab.github.io/proteoDA/reference/check_vec.md)
  : Check if a Vector Includes a Set of Reference Values
- [`classify_factorial_batch()`](https://byrumlab.github.io/proteoDA/reference/classify_factorial_batch.md)
  : Classify proteins under a factorial design
- [`compute_movingSD_zscores()`](https://byrumlab.github.io/proteoDA/reference/compute_movingSD_zscores.md)
  : Compute Rolling Standard Deviations and Z-scores for logFC
  Comparisons
- [`dummy_statmod()`](https://byrumlab.github.io/proteoDA/reference/dummy_statmod.md)
  : Dummy function for statmod import
- [`extract_DA_results()`](https://byrumlab.github.io/proteoDA/reference/extract_DA_results.md)
  : Extract differential abundance results from a model fit
- [`filter_proteins_by_annotation()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_annotation.md)
  : Remove proteins based on annotation data
- [`filter_proteins_by_group()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_group.md)
  : Filter protein data by number of quantified samples in a group
- [`filter_proteins_by_proportion()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_proportion.md)
  : Filter protein data by proportion of quantified samples in group
- [`filter_proteins_per_contrast()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_per_contrast.md)
  : Filter proteins separately for each contrast using sample group info
- [`filter_samples()`](https://byrumlab.github.io/proteoDA/reference/filter_samples.md)
  : Filter samples from a DAList
- [`find_id_column()`](https://byrumlab.github.io/proteoDA/reference/find_id_column.md)
  : Identify ID Column
- [`fit_limma_model()`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md)
  : Fit the limma differential abundance model
- [`get_colors()`](https://byrumlab.github.io/proteoDA/reference/get_colors.md)
  : Generate colors for Sample Groups
- [`` `%notin%` ``](https://byrumlab.github.io/proteoDA/reference/grapes-notin-grapes.md)
  : Not in operator
- [`impute_missing_by_gtest()`](https://byrumlab.github.io/proteoDA/reference/impute_missing_by_gtest.md)
  : Run G-test-based imputation on DAList data or contrast-specific data
- [`interpret_protein_factorial()`](https://byrumlab.github.io/proteoDA/reference/interpret_protein_factorial.md)
  : Interpret a protein across factorial contrasts
- [`missing_to_zero()`](https://byrumlab.github.io/proteoDA/reference/missing_data.md)
  [`zero_to_missing()`](https://byrumlab.github.io/proteoDA/reference/missing_data.md)
  : Conversion of missing values and 0s
- [`new_DAList()`](https://byrumlab.github.io/proteoDA/reference/new_DAList.md)
  : DAList internal constructor
- [`normalize_data()`](https://byrumlab.github.io/proteoDA/reference/normalize_data.md)
  : Normalize data in a DAList (contrast-aware)
- [`perseus_impute()`](https://byrumlab.github.io/proteoDA/reference/perseus_impute.md)
  : Perseus-style MNAR imputation (left-censoring) on log2 data
- [`plot_perseus_imputation()`](https://byrumlab.github.io/proteoDA/reference/plot_perseus_imputation.md)
  : Plot histograms highlighting imputed values (Perseus-style
  imputation)
- [`qc_boxplot()`](https://byrumlab.github.io/proteoDA/reference/qc_boxplot.md)
  : Generate a QC Boxplot
- [`qc_boxplot_beforeNorm()`](https://byrumlab.github.io/proteoDA/reference/qc_boxplot_beforeNorm.md)
  : QC boxplot for data before normalization
- [`qc_dendrogram_subgroups()`](https://byrumlab.github.io/proteoDA/reference/qc_dendrogram_subgroups.md)
  : Dendrogram for QC Report with Subgroups
- [`qc_pca_plot7()`](https://byrumlab.github.io/proteoDA/reference/qc_pca_plot7.md)
  : PCA Plot
- [`qc_pca_plot7_subgroups()`](https://byrumlab.github.io/proteoDA/reference/qc_pca_plot7_subgroups.md)
  : PCA Plot for Subgroups
- [`qc_pca_scree_plot7()`](https://byrumlab.github.io/proteoDA/reference/qc_pca_scree_plot7.md)
  : Scree Plot for PCA
- [`qc_totInt()`](https://byrumlab.github.io/proteoDA/reference/qc_totInt.md)
  : Generate a Barplot of Total Intensities by Sample Group
- [`qc_totInt_by_group()`](https://byrumlab.github.io/proteoDA/reference/qc_totInt_by_group.md)
  : Generate Total Intensity Barplots by Group
- [`qc_violin()`](https://byrumlab.github.io/proteoDA/reference/qc_violin.md)
  : Generate a Violin Plot for Intensity Data
- [`rename_samples()`](https://byrumlab.github.io/proteoDA/reference/rename_samples.md)
  : Rename samples and/or reorganize DGEList by group levels
- [`rescale_contrast_logFC()`](https://byrumlab.github.io/proteoDA/reference/rescale_contrast_logFC.md)
  : Rescale log-fold changes for a single contrast
- [`run_filtered_limma_analysis()`](https://byrumlab.github.io/proteoDA/reference/run_filtered_limma_analysis.md)
  : Run limma analysis per contrast on filtered proteins
- [`split_by_groups()`](https://byrumlab.github.io/proteoDA/reference/split_by_groups.md)
  : Split a data.frame or matrix according to group membership.
- [`validate_contrasts()`](https://byrumlab.github.io/proteoDA/reference/validate_contrasts.md)
  : Validate contrasts format
- [`validate_formula()`](https://byrumlab.github.io/proteoDA/reference/validate_formula.md)
  : Validate a statistical formula
- [`write_limma_plots()`](https://byrumlab.github.io/proteoDA/reference/write_limma_plots.md)
  : Make interactive reports on differential abundance
- [`write_limma_tables()`](https://byrumlab.github.io/proteoDA/reference/write_limma_tables.md)
  : Write tables of limma results
- [`write_movingSD_report()`](https://byrumlab.github.io/proteoDA/reference/write_movingSD_report.md)
  : Write a movingSD QC report (PDF and/or PNGs)
- [`write_norm_eval_report()`](https://byrumlab.github.io/proteoDA/reference/write_norm_eval_report.md)
  : Write a post-normalization evaluation report for a single method
- [`write_norm_report()`](https://byrumlab.github.io/proteoDA/reference/write_norm_report.md)
  : Create a normalization report (optionally contrast-aware)
- [`write_perseus_imputation_plots()`](https://byrumlab.github.io/proteoDA/reference/write_perseus_imputation_plots.md)
  : Write Perseus imputation histograms for each contrast
- [`write_qc_report()`](https://byrumlab.github.io/proteoDA/reference/write_qc_report.md)
  : Create a quality control report
