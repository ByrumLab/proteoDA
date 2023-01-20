

library(devtools)
library(tidyverse)



# Pipeline ----------------------------------------------------------------


# Load in data on a bunch of files --------------------------------

#ext_bart <- read_DIA_data("for_testing/Example Data/04_Bartholomew_101520_DIA/Samples Report of Bartholomew_101520.CSV") # Missing exclusivity col

higgs <- read_DIA_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv") # worked

ndu <- read_DIA_data("for_testing/Example Data/NDu_030822_DIA/input_files/Samples Report of Du_030822.csv") # Worked

kinter <- read_DIA_data("for_testing/Example Data/19_Kinter_120720_TMT_DIA_AG/Kinter_120720_DIA/Samples Report of Kinter_DIA_022521.csv") # Worked

lupashin <- read_DIA_data("for_testing/Example Data/lupashin_030222/Samples Report of Lupashin_030222.csv") # Worked

zhan <- read_DIA_data("for_testing/Example Data/Zhan_DIA_217_samples/input_files/Samples Report of Zhan_111821_Experiment.csv") # Worked

reb <- read_DIA_data("for_testing/Example Data/rebello/Samples Report of Rebello_040522.csv") # Worked

kaul <- read_DIA_data("for_testing/Example Data/kaul/Samples Report of Kaul_030922.csv") # Worked

ext_porter <- read_DIA_data("for_testing/Example Data/porter/Samples Report of PorterC_100422.csv")

# Add metadata ------------------------------------------------------------

higgs <- add_metadata(higgs, "for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv") # worked

ndu <- add_metadata(ndu, "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv") # worked

ndu_chain <- read_DIA_data("for_testing/Example Data/NDu_030822_DIA/input_files/Samples Report of Du_030822.csv") |>
  add_metadata("for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv")

ndu_chain2 <- read_DIA_data("for_testing/Example Data/NDu_030822_DIA/input_files/Samples Report of Du_030822.csv") |>
  add_metadata("for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv")

lupashin <- add_metadata(lupashin, "for_testing/Example Data/lupashin_030222/Lupashin_030222_metafile_DIA.csv") # Worked

#zhan <- add_metadata(zhan, "for_testing/Example Data/rebello/Rebello_040522_metafile_DIA.csv") # error as expected

zhan <- add_metadata(zhan, "for_testing/Example Data/Zhan_DIA_217_samples/input_files/Zhan_111821_DIA_metadata.csv") # worked

#reb <- add_metadata(reb, "for_testing/Example Data/Zhan_DIA_217_samples/input_files/Zhan_111821_DIA_metadata.csv") # error as expected

reb <- add_metadata(reb, "for_testing/Example Data/rebello/Rebello_040522_metafile_DIA.csv") # worked

kaul <- add_metadata(kaul, "for_testing/Example Data/kaul/Kaul_030922_metafile_DIA.csv")

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

# Kaul
write_norm_report(filtered_kaul,
                        grouping_column = "group",
                        file = "kaul_update_2.pdf",
                        overwrite = T, suppress_zoom_legend = T)

# Normalize data ----------------------------------------------------------
norm_kaul <- filtered_kaul |>
  normalize_data("quantile")

norm_lupashin <- filtered_lupashin |>
  normalize_data("rlr")

norm_ndu <- filtered_ndu |>
  normalize_data("median")

norm_reb <- filtered_reb |>
  normalize_data("vsn")

norm_zhan <- filtered_zhan |>
  normalize_data("cycloess")

# Make QC report ----------------------------------------------------------
full_higgs_chain <- read_DIA_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv") |>
  add_metadata("for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv") |>
  filter_samples(group != "Pool") |>
  filter_proteins_by_group(min_reps = 4, min_groups = 3) |>
  filter_proteins_by_group(min_reps = 5, min_groups = 3) |>
  filter_proteins_by_proportion(min_prop = 1) |>
  normalize_data(norm_method = "cycloess")

write_qc_report(full_higgs_chain,
                color_column = "group",
                filename = "higgs_qc_update.pdf",
                overwrite = T)

write_qc_report(full_higgs_chain,
                color_column = "group",
                label_column = "sampleIDs",
                filename = "higgs_qc_update_samplelabs.pdf",
                overwrite = T)

write_qc_report(norm_ndu,
                color_column = "group",
                filename = "ndu_qc_update.pdf",
                overwrite = T)

write_qc_report(norm_ndu,
                filename = "ndu_qc_update_batch.pdf",
                overwrite = T)


write_qc_report(norm_lupashin,
                color_column = "group",
                filename = "lupashin_qc_update.pdf",
                overwrite = T)

write_qc_report(norm_zhan,
                color_column = "group",
                filename = "zhan_qc_update.pdf",
                overwrite = T)

write_qc_report(norm_reb,
                color_column = "group",
                filename = "rebello_qc_update_changes.pdf",
                overwrite = T, standardize = T, top_proteins = nrow(norm_reb$data),
                pca_axes = c(2, 5))

write_qc_report(norm_kaul,
                color_column = "group",
                filename = "kaul_qc_update.pdf",
                overwrite = T)


# Make design -------------------------------------------------------------
norm_ndu$metadata <- norm_ndu$metadata |>
  separate(group, into = c("treatment", "tissue"), remove = F)

norm_kaul <- add_design(norm_kaul,
                        ~ 0 + group + gender)

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

full_higgs_chain <- read_DIA_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv") |>
  add_metadata("for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv") |>
  filter_samples(group != "Pool") |>
  filter_proteins_contaminants() |>
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

#kaul
norm_kaul$metadata$gender <- factor(norm_kaul$metadata$gender, levels = c("female", "male"))


# Lupashin
norm_lupashin <- norm_lupashin |>
  add_contrasts(contrasts_file = "for_testing/Example Data/lupashin_030222/contrasts_bad.csv")
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

fit_kaul <- fit_limma_model(norm_kaul)


fit_ndu_random <- norm_ndu |>
  add_design(design_formula = "~ 0 + treatment + (1 | tissue)") |>
  fit_limma_model()

full_higgs_chain <- read_DIA_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv") |>
  add_metadata("for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv") |>
  filter_samples(group != "Pool") |>
  filter_proteins_contaminants() |>
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
results_kaul <- extract_DA_results(fit_kaul)
results_higgs <- extract_DA_results(fit_higgs, extract_intercept = F)
names(results_higgs$results)



# Write results -----------------------------------------------------------
write_limma_tables(results_lupashin,
                   output_dir = "Lupashin_s3obj",
                   overwrite = T)

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

write_limma_tables(results_kaul,
                   output_dir = "kaul_s3obj",
                   overwrite = T)

write_limma_tables(results_higgs,
                   output_dir = "higgs_s3obj",
                   overwrite = T)

# testing report making ---------------------------------------------------
write_limma_plots(results_reb,
                  grouping_column = "group",
                  output_dir = "reb_s3obj", key_column = "Protein.Name")

write_limma_plots(results_kaul,
                  grouping_column = "group",
                  output_dir = "kaul_s3obj/")


write_limma_plots(results_reb,
                  grouping_column = "group",
                  output_dir = "reb_s3obj/")

write_limma_plots(results_lupashin,
                  grouping_column = "group",
                  output_dir = "Lupashin_s3obj")

write_limma_plots(results_ndu,
                  grouping_column = "group",
                  output_dir = "Ndu_s3obj")

write_limma_plots(results_ndu_random,
                  grouping_column = "group",
                  output_dir = "Ndu_random_s3obj")


write_limma_plots(results_zhan,
                  grouping_column = "group",
                  output_dir = "zhan_s3obj")

write_limma_plots(results_zhan,
                  grouping_column = "group",
                  output_dir = "zhan_wide_s3obj",
                   width = 2000)

write_limma_plots(results_higgs,
                  grouping_column = "group",
                  output_dir = "higgs_s3obj")
write_limma_plots(results_higgs,
                  #key_column = "Protein.Name",
                  grouping_column = "group",
                  output_dir = "higgs_s3obj")
