

library(devtools) # nocov start
library(tidyverse)



# Pipeline ----------------------------------------------------------------

# Load in data on a bunch of files --------------------------------
# Originally these were loaded with our internal
# functions, (read_DIA_data and add_metadata)
# Saved those as RDS files for use in the test pipeline
# upon removal of the internal functions.
for (dataset in c("higgs", "kaul", "lupashin", "ndu", "reb", "zhan")) {
  assign(x = dataset,
         value = readRDS(paste0("for_testing/Example Data/rds_files/", dataset, ".rds")))
}

# filter samples --------------------------------------------------------------
sub_higgs <- filter_samples(higgs, group != "Pool")

sub_ndu <- filter_samples(ndu, group != "Pool")

sub_lupashin <- filter_samples(lupashin, group != "Pool") |>
  filter_samples(!stringr::str_detect(group, "input"))

sub_zhan <- filter_samples(zhan, group != "Pool")

sub_reb <- filter_samples(reb, group != "Pool")

sub_kaul <- filter_samples(kaul, group != "Pool")

# filter proteins ---------------------------------------------------------

filtered_higgs <- filter_proteins_by_group(sub_higgs, min_reps = 5, min_groups = 3)

filtered_ndu <- filter_proteins_by_group(sub_ndu, min_reps = 3, min_groups = 1)

filtered_lupashin <- filter_proteins_by_group(sub_lupashin, min_reps = 1, min_groups = 1)

filtered_zhan <- filter_proteins_by_group(sub_zhan, min_reps = 13, min_groups = 2)

filtered_zhan_2 <- filter_proteins_by_proportion(sub_zhan, min_prop = 0.66)

filtered_reb  <- filter_proteins_by_group(sub_reb, min_reps = 2, min_groups = 1)

filtered_kaul <- filter_proteins_by_group(sub_kaul, min_reps = 4, min_groups = 2)


# Normalization report ----------------------------------------------------
# Higgs
write_norm_report(filtered_higgs,
                  grouping_column = "group",
                  file = "higgs_update_2.pdf", overwrite = T)


# Ndu
write_norm_report(filtered_ndu,
                  output_dir = "temp",
                        grouping_column = "group",
                        file = "ndu_update_with_points_ggrastr.pdf", overwrite = T, use_ggrastr = T)

# Lupashin
write_norm_report(filtered_lupashin,
                        grouping_column = "group",
                        file = "lupashin_update_with_points.pdf", overwrite = T, suppress_zoom_legend = T)
# Zhan
write_norm_report(filtered_zhan_2,
                        grouping_column = "group",
                        file = "zhan_update_2.pdf", overwrite = T)

# Rebello
write_norm_report(filtered_reb,
                        grouping_column = "group",
                        file = "rebello_update_2.pdf",
                        overwrite = T)

# Normalize data ----------------------------------------------------------
norm_lupashin <- filtered_lupashin |>
  normalize_data("rlr")

norm_ndu <- filtered_ndu |>
  normalize_data("median")

norm_reb <- filtered_reb |>
  normalize_data("vsn")

norm_zhan <- filtered_zhan |>
  normalize_data("cycloess")

# Make QC report ----------------------------------------------------------
full_higgs_chain <- higgs |>
  filter_samples(group != "Pool") |>
  filter_proteins_by_group(min_reps = 4, min_groups = 3) |>
  filter_proteins_by_group(min_reps = 5, min_groups = 3) |>
  filter_proteins_by_proportion(min_prop = 1) |>
  normalize_data(norm_method = "cycloess")


write_qc_report(full_higgs_chain,
                output_dir = "update",
                color_column = "group",
                filename = "higgs_qc.pdf",
                overwrite = T)

# To test longer IDs, reverse the
# sample ID strings so they're still unique after truncation
full_higgs_chain$metadata$sampleIDs_rev <- stringi::stri_reverse(full_higgs_chain$metadata$sampleIDs)
write_qc_report(full_higgs_chain,
                output_dir = "update",
                color_column = "group",
                label_column = "sampleIDs_rev",
                filename = "higgs_qc_samplelabs.pdf",
                overwrite = T)

# Some "normal" ones
write_qc_report(norm_ndu,
                output_dir = "update",
                color_column = "group",
                filename = "ndu_qc.pdf",
                overwrite = T)

write_qc_report(norm_ndu,
                output_dir = "update",
                filename = "ndu_qc_batch.pdf",
                overwrite = T)


write_qc_report(norm_lupashin,
                output_dir = "update",
                color_column = "group",
                filename = "lupashin_qc.pdf",
                overwrite = T)

# Use the Zhan data to test different sizes:
# really big (regular data)
# just over 50, and
# just under 50
norm_zhan$metadata$long_ID <- stringr::str_pad(norm_zhan$metadata$sample, width = 25, side = "right", pad = "X")
norm_zhan_50 <- norm_zhan %>%
  filter_samples(as.numeric(number) %% 2 == 0) %>%
  filter_samples(rep(c(T, F), times = 54))

norm_zhan_49 <- norm_zhan %>%
  filter_samples(as.numeric(number) %% 2 == 0) %>%
  filter_samples(c(rep(c(T, F), times = 49), rep(F, 10)))

write_qc_report(norm_zhan,
                output_dir = "update",
                color_column = "group",
                filename = "zhan_qc.pdf",
                overwrite = T)
write_qc_report(norm_zhan,
                output_dir = "update",
                label_column = "long_ID",
                color_column = "group",
                filename = "zhan_qc_ids.pdf",
                overwrite = T)


write_qc_report(norm_zhan_50,
                output_dir = "update",
                color_column = "group",
                filename = "zhan54_qc.pdf",
                overwrite = T)
write_qc_report(norm_zhan_50,
                output_dir = "update",
                label_column = "long_ID",
                color_column = "group",
                filename = "zhan54_qc_ids.pdf",
                overwrite = T)


write_qc_report(norm_zhan_49,
                output_dir = "update",
                color_column = "group",
                filename = "zhan49_qc.pdf",
                overwrite = T)
write_qc_report(norm_zhan_49,
                output_dir = "update",
                label_column = "long_ID",
                color_column = "group",
                filename = "zhan49_qc_ids.pdf",
                overwrite = T)



write_qc_report(norm_reb,
                output_dir = "update",
                color_column = "group",
                filename = "rebello_qc.pdf",
                overwrite = T, standardize = T, top_proteins = nrow(norm_reb$data),
                pca_axes = c(2, 5))



# Make design -------------------------------------------------------------
norm_ndu$metadata <- norm_ndu$metadata |>
  tidyr::separate(group, into = c("treatment", "tissue"), remove = F)

norm_lupashin <- add_design(norm_lupashin,
                            design_formula = "~ 0 + group")

norm_ndu <- add_design(norm_ndu,
                design_formula = "~ 0 + treatment + (1 | tissue)")

norm_ndu <- add_design(norm_ndu,
                       design_formula = "~ 0 + treatment + tissue")


norm_reb <- add_design(norm_reb,
                       "~ group")

norm_zhan <- add_design(norm_zhan,
                        ~group*sex)

full_higgs_chain <- higgs |>
  filter_samples(group != "Pool") |>
  filter_proteins_by_group(min_reps = 4, min_groups = 3) |>
  filter_proteins_by_group(min_reps = 5, min_groups = 3) |>
  filter_proteins_by_proportion(min_prop = 1) |>
  normalize_data(norm_method = "cycloess") |>
  add_design(~0 +group)

# Make contrasts ----------------------------------------------------------
# Higgs
# Will be good test of later code, can let higgs run the model without contrasts

# Ndu
norm_ndu <- norm_ndu |>
  add_design(~ 0 + group) |>
  add_contrasts(contrasts_file = "for_testing/Example Data/NDu_030822_DIA/input_files/kidney_contrasts.txt")

# Lupashin
# norm_lupashin <- norm_lupashin |>
#   add_contrasts(contrasts_file = "for_testing/Example Data/lupashin_030222/contrasts_bad.csv")
norm_lupashin <- norm_lupashin |>
  add_contrasts(contrasts_file = "for_testing/Example Data/lupashin_030222/contrasts.csv")

# Zhan
norm_zhan <- norm_zhan |>
  add_design(~0 + group) |>
  add_contrasts(contrasts_file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/contrasts.txt")

norm_zhan <- norm_zhan |>
  add_design(~0 + group) |>
  add_contrasts(contrasts_vector = c("test=HiR-LoR"))

# Rebello
norm_reb <- norm_reb |>
  add_design(~0 + group) |>
  add_contrasts(contrasts_file = "for_testing/Example Data/rebello/contrasts.csv")


# Run the analysis --------------------------------------------------------
# Splitting up functionality
# First, fit the model
fit_lupashin <- fit_limma_model(norm_lupashin)

fit_ndu<- fit_limma_model(norm_ndu)

fit_zhan <- fit_limma_model(norm_zhan)

fit_reb <- fit_limma_model(norm_reb)

fit_ndu_random <- norm_ndu |>
  add_design(design_formula = "~ 0 + treatment + (1 | tissue)") |>
  add_contrasts(contrasts_vector = "SLCKO_vs_control = SLCKO  - Control") |>
  fit_limma_model()

full_higgs_chain <- higgs |>
  filter_samples(group != "Pool") |>
  filter_proteins_by_group(min_reps = 4, min_groups = 3) |>
  filter_proteins_by_group(min_reps = 5, min_groups = 3) |>
  filter_proteins_by_proportion(min_prop = 1) |>
  normalize_data(norm_method = "cycloess") |>
  add_design(~group)


fit_higgs <- fit_limma_model(full_higgs_chain)


# Extract results ---------------------------------------------------------
results_lupashin <- extract_DA_results(fit_lupashin)
results_ndu <- extract_DA_results(fit_ndu)
results_ndu_random <- extract_DA_results(fit_ndu_random)
results_zhan <- extract_DA_results(fit_zhan)
results_reb <- extract_DA_results(fit_reb)
results_higgs <- extract_DA_results(fit_higgs, extract_intercept = F)


# Write results -----------------------------------------------------------
write_limma_tables(results_lupashin, overwrite = T)

write_limma_tables(results_ndu,
                   overwrite = T)

write_limma_tables(results_ndu_random,
                   output_dir = "Ndu_random_s3obj",
                   overwrite = T)

write_limma_tables(results_reb,
                   output_dir = "reb_s3obj",
                   overwrite = T)

write_limma_tables(results_zhan,
                   output_dir = "zhan_s3obj",
                   overwrite = T)

write_limma_tables(results_higgs,
                   output_dir = "higgs_s3obj",
                   overwrite = T)

# testing report making ---------------------------------------------------


write_limma_plots(results_reb,
                  grouping_column = "group",
                  output_dir = "reb_s3obj/default_cols",
                  overwrite = T)

write_limma_plots(results_reb,
                  grouping_column = "group",
                  output_dir = "reb_s3obj/title_col",
                  title_column = "Accession.Number", overwrite = T)

write_limma_plots(results_reb,
                  output_dir = "reb_s3obj/title_and_table_cols",
                  grouping_column = "group",
                  title_column = "Protein.Name",
                  table_columns = c("Molecular.Weight", "Gene_name"), overwrite = T)


write_limma_plots(results_lupashin,
                  grouping_column = "group",
                  output_dir = "Lupashin_s3obj", overwrite = T)

write_limma_plots(results_ndu,
                  grouping_column = "group",
                  output_dir = "Ndu_s3obj", overwrite = T)

write_limma_plots(results_ndu_random,
                  grouping_column = "group",
                  output_dir = "Ndu_random_s3obj", overwrite = T)


write_limma_plots(results_zhan,
                  grouping_column = "group",
                  output_dir = "zhan_s3obj",
                  overwrite = T)

write_limma_plots(results_zhan,
                  grouping_column = "group",
                  output_dir = "zhan_wide_s3obj",
                   width = 2000, overwrite = T)

write_limma_plots(results_higgs,
                  grouping_column = "group",
                  output_dir = "higgs_s3obj/default_cols", overwrite = T)

# Is the | symbol in protein name the issue???
results_higgs$annotation$long_accession <- stringr::str_pad(results_higgs$annotation$Accession.Number, width = "25", side = "right", pad = "X")
write_limma_plots(results_higgs,
                  title_column = "Protein.Name",
                  grouping_column = "group",
                  output_dir = "higgs_s3obj/title_col/protein_name",
                  overwrite = T)

write_limma_plots(results_higgs,
                  title_column = "long_accession",
                  grouping_column = "group",
                  output_dir = "higgs_s3obj/title_col/long_accession",
                  overwrite = T)

write_limma_plots(results_higgs,
                  title_column = "Accession.Number",
                  grouping_column = "group",
                  output_dir = "higgs_s3obj/title_col/accession_number",
                  overwrite = T)

# nocov start
