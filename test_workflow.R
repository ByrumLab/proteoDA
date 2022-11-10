
# library(pkgnet)
# library(tictoc)
# library(microbenchmark)
library(devtools)
library(tidyverse)
# library(profvis)



# For figuring out argument names
#lsf.str("package:proteomicsDIA")

CreatePackageReport("proteomicsDIA")

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

# Add metadata ------------------------------------------------------------

higgs <- add_metadata(higgs, "for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv") # worked

ndu <- add_metadata(ndu, "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv") # worked

ndu_chain <- read_DIA_data("for_testing/Example Data/NDu_030822_DIA/input_files/Samples Report of Du_030822.csv") %>%
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

sub_lupashin <- filter_samples(lupashin, group != "Pool") %>%
  filter_samples(!stringr::str_detect(group, "input"))

sub_zhan <- filter_samples(zhan, group != "Pool")

sub_reb <- filter_samples(reb, group != "Pool")

sub_kaul <- filter_samples(kaul, group != "Pool")

# filter proteins ---------------------------------------------------------

filtered_higgs <- filter_proteins_contaminants(sub_higgs) %>%
  filter_proteins_by_group(min_reps = 5, min_groups = 3)

filtered_ndu <- filter_proteins_contaminants(sub_ndu) %>%
  filter_proteins_by_group(min_reps = 3, min_groups = 1)

filtered_lupashin <- filter_proteins_contaminants(sub_lupashin) %>%
  filter_proteins_by_group(min_reps = 1, min_groups = 1)

filtered_zhan <- filter_proteins_contaminants(sub_zhan) %>%
  filter_proteins_by_group(min_reps = 13, min_groups = 2)

filtered_zhan_2 <- filter_proteins_contaminants(sub_zhan) %>%
  filter_proteins_by_proportion(min_prop = 0.66)

filtered_reb  <- filter_proteins_contaminants(sub_reb) %>%
  filter_proteins_by_group(min_reps = 2, min_groups = 1)

filtered_kaul <- filter_proteins_contaminants(sub_kaul) %>%
  filter_proteins_by_group(min_reps = 4, min_groups = 2)


# Normalization report ----------------------------------------------------
# Higgs
write_proteinorm_report(full_higgs_chain,
                        grouping_column = "group",
                        file = "higgs_update_2.pdf", overwrite = T)

# Ndu
write_proteinorm_report(filtered_ndu,
                        grouping_column = "group",
                        file = "ndu_update_2.pdf", overwrite = T)

# Lupashin
write_proteinorm_report(filtered_lupashin,
                        grouping_column = "group",
                        file = "lupashin_update_2.pdf", overwrite = T, suppress_zoom_legend = T)
# Zhan
write_proteinorm_report(filtered_zhan_2,
                        grouping_column = "group",
                        file = "zhan_update_2.pdf", overwrite = T)

# Rebello
write_proteinorm_report(filtered_reb,
                        grouping_column = "group",
                        file = "rebello_update_2.pdf",
                        overwrite = T)

# Kaul
write_proteinorm_report(filtered_kaul,
                        grouping_column = "group",
                        file = "kaul_update_2.pdf",
                        overwrite = T, suppress_zoom_legend = T)

# Normalize data ----------------------------------------------------------
norm_kaul <- filtered_kaul %>%
  normalize_data("quantile")

norm_lupashin <- filtered_lupashin %>%
  normalize_data("rlr")

norm_ndu <- filtered_ndu %>%
  normalize_data("median")

norm_reb <- filtered_reb %>%
  normalize_data("vsn")

norm_zhan <- filtered_zhan %>%
  normalize_data("cycloess")

# Make QC report ----------------------------------------------------------
write_qc_report(full_higgs_chain,
                grouping_column = "group",
                file = "higgs_qc_update.pdf",
                overwrite = T)

write_qc_report(full_higgs_chain,
                grouping_column = "group",
                label_column = "sampleIDs",
                file = "higgs_qc_update_samplelabs.pdf",
                overwrite = T)

write_qc_report(norm_ndu,
                grouping_column = "group",
                file = "ndu_qc_update.pdf",
                overwrite = T)

write_qc_report(norm_ndu,
                file = "ndu_qc_update_batch.pdf",
                overwrite = T)


write_qc_report(norm_lupashin,
                grouping_column = "group",
                file = "lupashin_qc_update.pdf",
                overwrite = T)

write_qc_report(norm_zhan,
                grouping_column = "group",
                file = "zhan_qc_update.pdf",
                overwrite = T)

write_qc_report(norm_reb,
                grouping_column = "group",
                file = "rebello_qc_update.pdf",
                overwrite = T)

write_qc_report(norm_kaul,
                grouping_column = "group",
                file = "kaul_qc_update.pdf",
                overwrite = T)

# Make design -------------------------------------------------------------
norm_ndu$metadata <- norm_ndu$metadata %>%
  separate(group, into = c("treatment", "tissue"), remove = F)


norm_kaul <- add_design(norm_kaul,
                        ~ 0 +group + gender)

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

full_higgs_chain <- read_DIA_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv") %>%
  add_metadata("for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv") %>%
  filter_samples(group != "Pool") %>%
  filter_proteins_contaminants() %>%
  filter_proteins_by_group(min_reps = 4, min_groups = 3) %>%
  filter_proteins_by_group(min_reps = 5, min_groups = 3) %>%
  filter_proteins_by_proportion(min_prop = 1) %>%
  normalize_data(method = "cycloess") %>%
  add_design(~0 +group)

# Make contrasts ----------------------------------------------------------
# Higgs
# Will be good test of later code, can let higgs run the model without contrasts

# Ndu
norm_ndu <- norm_ndu %>%
  add_design(~ 0 + group) %>%
  add_contrasts(contrasts_file = "for_testing/Example Data/NDu_030822_DIA/input_files/kidney_contrasts.txt")

#kaul
norm_kaul$metadata$gender <- factor(norm_kaul$metadata$gender, levels = c("male", "female"))
norm_kaul <- norm_kaul %>%
  add_contrasts(contrasts_file = "for_testing/Example Data/kaul/contrasts_designfomula.csv")


# Lupashin
norm_lupashin <- norm_lupashin %>%
  add_contrasts(contrasts_file = "for_testing/Example Data/lupashin_030222/contrasts_bad.csv")
norm_lupashin <- norm_lupashin %>%
  add_contrasts(contrasts_file = "for_testing/Example Data/lupashin_030222/contrasts.csv")

# Zhan
norm_zhan <- norm_zhan %>%
  add_design(~0 + group) %>%
  add_contrasts(contrasts_file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/contrasts.txt")

norm_zhan <- norm_zhan %>%
  add_design(~0 + group) %>%
  add_contrasts(contrasts_vector = c("test=HiR-LoR"))

# Rebello
norm_reb <- norm_reb %>%
  add_design(~0 + group) %>%
  add_contrasts(contrasts_file = "for_testing/Example Data/rebello/contrasts.csv")


# Run the analysis --------------------------------------------------------
# Splitting up functionality
# First, fit the model
fit_lupashin <- fit_limma_model(norm_lupashin)

fit_ndu<- fit_limma_model(norm_ndu)

fit_zhan <- fit_limma_model(norm_zhan)

fit_reb <- fit_limma_model(norm_reb)

fit_kaul <- fit_limma_model(norm_kaul)


fit_ndu_random <- norm_ndu %>%
  add_design(design_formula = "~ 0 + treatment + (1 | tissue)") %>%
  fit_limma_model()

full_higgs_chain <- read_DIA_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv") %>%
  add_metadata("for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv") %>%
  filter_samples(group != "Pool") %>%
  filter_proteins_contaminants() %>%
  filter_proteins_by_group(min_reps = 4, min_groups = 3) %>%
  filter_proteins_by_group(min_reps = 5, min_groups = 3) %>%
  filter_proteins_by_proportion(min_prop = 1) %>%
  normalize_data(method = "cycloess") %>%
  add_design(~0 +group)


fit_higgs <- fit_limma_model(full_higgs_chain)


# Extract results ---------------------------------------------------------
results_lupashin <- extract_limma_DE_results(limma_fit = fit_lupashin)
results_ndu_brain <- extract_limma_DE_results(limma_fit = fit_ndu_brain)
results_ndu_intestine <- extract_limma_DE_results(limma_fit = fit_ndu_intestine)
results_ndu_kidney <- extract_limma_DE_results(limma_fit = fit_ndu_kidney)
results_zhan <- extract_limma_DE_results(limma_fit = fit_zhan)
results_zhan2 <- extract_limma_DE_results(limma_fit = fit_zhan2)
results_reb <- extract_limma_DE_results(limma_fit = fit_reb)
results_kaul <- extract_limma_DE_results(limma_fit = fit_kaul)

waldo::compare(results_zhan, results_zhan2)


# Write results -----------------------------------------------------------
write_limma_tables(model_results = results_lupashin,
                    norm.method = "vsn",
                    annotation = ext_lupashin$annot,
                    ilab = "Lupashin_82928",
                    overwrite = T)

write_limma_tables(model_results = results_ndu_brain,
                    norm.method = "vsn",
                    annotation = ext_ndu$annot,
                    ilab = "ndu_brain_82928")

write_limma_tables(model_results = results_ndu_intestine,
                    norm.method = "vsn",
                    annotation = ext_ndu$annot,
                    ilab = "ndu_intestine_82928")

write_limma_tables(model_results = results_ndu_kidney,
                    norm.method = "vsn",
                    annotation = ext_ndu$annot,
                    ilab = "ndu_kidney_82928")

write_limma_tables(model_results = results_reb,
                    norm.method = "vsn",
                    annotation = ext_reb$annot,
                    ilab = "Rebello_82928")


write_limma_tables(model_results = results_zhan,
                    norm.method = "vsn",
                    annotation = ext_zhan$annot,
                    ilab = "zhan_982974")

write_limma_tables(model_results = results_kaul,
                    norm.method = "cycloess",
                    annotation = ext_kaul$annot,
                    ilab = "kaul_82921", overwrite=T)


# testing report making ---------------------------------------------------
write_limma_plots(model_results = results_reb,
                   annotation = ext_reb$annot,
                   groups = norm_reb$targets$group,
                   output_dir = "output_rebello")

write_limma_plots(model_results = results_kaul,
                   annotation = ext_kaul$annot,
                   groups = norm_kaul$targets$group,
                   output_dir = "output_kaul")

write_limma_plots(model_results = results_reb,
                   annotation = ext_reb$annot,
                   groups = norm_reb$targets$group,
                   output_dir = "output_rebello_wide",
                   width = 1500,
                   height = 1500)


write_limma_plots(model_results = results_lupashin,
                   annotation = ext_lupashin$annot,
                   groups = norm_lupashin$targets$group,
                   output_dir = "output_lupashin")


write_limma_plots(model_results = results_ndu_brain,
                   annotation = ext_ndu$annot,
                   groups = norm_ndu$targets$group,
                   output_dir = "output_ndu_brain")


write_limma_plots(model_results = results_zhan,
                   annotation = ext_zhan$annot,
                   groups = norm_zhan$targets$group,
                   output_dir = "output_zhan")


write_limma_plots(model_results = results_zhan,
                   annotation = ext_zhan$annot,
                   groups = norm_zhan$targets$group,
                   output_dir = "output_zhan_wide",
                   width = 2000)
