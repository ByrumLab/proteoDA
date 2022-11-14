
# library(pkgnet)
# library(tictoc)
# library(microbenchmark)
library(devtools)
# library(profvis)



# For figuring out argument names
#lsf.str("package:proteomicsDIA")

CreatePackageReport("proteomicsDIA")
# source("bin/functions26_CURRENT_031022.r")

# making some temporary changes
# make a change again

# Pipeline ----------------------------------------------------------------


# extract_data on a bunch of files --------------------------------

ext_bart <- read_DIA_data("for_testing/Example Data/04_Bartholomew_101520_DIA/Samples Report of Bartholomew_101520.CSV") # Missing exclusivity col

ext_higgs <- read_DIA_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv") # Worked

ext_ndu <- read_DIA_data("for_testing/Example Data/NDu_030822_DIA/input_files/Samples Report of Du_030822.csv") # Worked

ext_kinter <- read_DIA_data("for_testing/Example Data/19_Kinter_120720_TMT_DIA_AG/Kinter_120720_DIA/Samples Report of Kinter_DIA_022521.csv") # Worked

ext_lupashin <- read_DIA_data("for_testing/Example Data/lupashin_030222/Samples Report of Lupashin_030222.csv") # Worked

ext_zhan <- read_DIA_data("for_testing/Example Data/Zhan_DIA_217_samples/input_files/Samples Report of Zhan_111821_Experiment.csv")

ext_reb <- read_DIA_data("for_testing/Example Data/rebello/Samples Report of Rebello_040522.csv")

ext_kaul <- read_DIA_data("for_testing/Example Data/kaul/Samples Report of Kaul_030922.csv")

ext_porter <- read_DIA_data("for_testing/Example Data/porter/Samples Report of PorterC_100422.csv")

# Make targets ------------------------------------------------------------

target_higgs_metadata <- make_targets(input_file = "for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv",
                             sample_IDs = colnames(ext_higgs$data)) # worked

target_higgs_IDs_only <- make_targets(#file = "for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv",
                             sample_IDs = colnames(ext_higgs$data))

target_ndu <- make_targets(input_file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                           sample_IDs = colnames(ext_ndu$data)) # worked

target_ndu_IDs_only<- make_targets(#file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                           sample_IDs = colnames(ext_ndu$data))

target_ndu_extra_IDs <- make_targets(input_file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                           sample_IDs = c(colnames(ext_ndu$data), "sample_45")) # Gave warning as expected

target_ndu_extra_meta<- make_targets(input_file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                                     sample_IDs = colnames(ext_higgs$data))

target_lupashin <- make_targets(input_file = "for_testing/Example Data/lupashin_030222/Lupashin_030222_metafile_DIA.csv",
                                sample_IDs = colnames(ext_lupashin$data)) # Worked

target_zhan <- make_targets(input_file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/Zhan_111821_DIA_metadata.csv",
                            sample_IDs = colnames(ext_zhan$data))

target_reb <- make_targets(input_file = "for_testing/Example Data/rebello/Rebello_040522_metafile_DIA.csv",
                           sample_IDs = colnames(ext_reb$data))

target_kaul <- make_targets(input_file = "for_testing/Example Data/kaul/Kaul_030922_metafile_DIA.csv",
                            sample_IDs =  colnames(ext_kaul$data))

target_porter <- make_targets("for_testing/Example Data/porter/metafile_DIA_TJT.csv",
                              sample_IDs = colnames(ext_porter$data))

# Subset targets --------------------------------------------------------------
sub_higgs <- subset_targets(targets = target_higgs_metadata,
                            filter_list = list(group = "Pool"))

sub_ndu <- subset_targets(targets = target_ndu,
                          filter_list = list(group = "Pool"))

sub_lupashin <- subset_targets(targets = target_lupashin,
                               filter_list = list(group = c("Pool", "input")))

sub_zhan <- subset_targets(targets = target_zhan,
                           filter_list = list(group = "Pool"))


sub_reb <- subset_targets(targets = target_reb,
                          filter_list = list(group = "Pool"))

sub_kaul <- subset_targets(targets = target_kaul,
                          filter_list = list(group = "Pool"))

sub_porter <- subset_targets(targets = target_porter,
                             filter_list = list(group = "Pool",
                                                sample = "S2117"))


# Process data ------------------------------------------------------------
norm_higgs <- process_data(data = ext_higgs$data,
                           targets = sub_higgs$targets,
                           min.reps = 5,
                           min.grps = 3)

norm_ndu <- process_data(data = ext_ndu$data,
                         targets = sub_ndu$targets,
                         min.reps = 3,
                         min.grps = 1)

norm_lupashin <- process_data(data = ext_lupashin$data,
                              targets = sub_lupashin$targets,
                              min.reps = 1,
                              min.grps = 1)

# tic()
norm_zhan <- process_data(data = ext_zhan$data,
                          targets = sub_zhan$targets,
                          min.reps = 13,
                          min.grps = 2)
# toc()

norm_reb <- process_data(data = ext_reb$data,
                         targets = sub_reb$targets,
                         min.reps = 2,
                         min.grps = 1)

norm_kaul <- process_data(data = ext_kaul$data,
                          targets = sub_kaul$targets,
                          min.reps = 4,
                          min.grps = 2)

norm_porter <- process_data(data = ext_porter$data,
                          targets = sub_porter$targets,
                          min.reps = 4,
                          min.grps = 2)

# Normalization report ----------------------------------------------------
# Higgs
write_proteinorm_report(processed_data = norm_higgs,
                       grouping_column = "group",
                       file = "higgs_update_2.pdf", overwrite = T)

# Ndu
write_proteinorm_report(processed_data = norm_ndu,
                       grouping_column = "group",
                       file = "ndu_update_2.pdf", overwrite = T)

# Lupashin
write_proteinorm_report(processed_data = norm_lupashin,
                       grouping_column = "group",
                       file = "lupashin_update_2.pdf", overwrite = T, suppress_zoom_legend = T)
# Zhan
# tic()
write_proteinorm_report(processed_data = norm_zhan,
                       grouping_column = "group",
                       file = "zhan_update_2.pdf", overwrite = T)
# toc()

# Rebello
write_proteinorm_report(processed_data = norm_reb,
                       grouping_column = "group",
                       file = "rebello_update_2.pdf",
                       overwrite = T)

# Make QC report ----------------------------------------------------------
write_qc_report(processed_data = norm_higgs,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                file = "higgs_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_higgs,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                label_column = "sampleIDs",
                file = "higgs_qc_update_samplelabs.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_ndu,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                file = "ndu_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_ndu,
                chosen_norm_method = "vsn",
                file = "ndu_qc_update_batch.pdf",
                overwrite = T)


write_qc_report(processed_data = norm_lupashin,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                file = "lupashin_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_zhan,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                file = "zhan_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_reb,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                file = "rebello_qc_update.pdf",
                overwrite = T)


write_qc_report(processed_data = norm_porter,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                file = "porter_qc_update.pdf",
                overwrite = T)

norm_porter$targets$group

# Make design -------------------------------------------------------------
des_higgs <- make_design(targets=norm_higgs$targets,
                        group_column = "group",
                        factor_columns = NULL,
                        paired_column = NULL)

des_higgs2 <- make_design(targets=norm_higgs$targets,
                        group_column = "group",
                        factor_columns = NULL,
                        paired_column = NULL,
                        design_formula = "~0 + group")

waldo::compare(des_higgs, des_higgs2)

des_ndu <- make_design(targets = norm_ndu$targets,
                       group_column = "group",
                       factor_columns = NULL,
                       paired_column =  NULL)

des_ndu2 <- make_design(targets = norm_ndu$targets,
                       group_column = "group",
                       factor_columns = "test",
                       paired_column =  NULL,
                       design_formula = "~ 0 + group")

waldo::compare(des_ndu, des_ndu2) # slight difference in order, should be no big deal

des_lupashin <- make_design(targets = norm_lupashin$targets,
                            group_column = "group",
                            factor_columns = NULL,
                            paired_column = NULL)

des_lupashin2 <- make_design(targets = norm_lupashin$targets,
                            group_column = "group",
                            factor_columns = NULL,
                            paired_column = NULL,
                            design_formula = "~ 0 + group")

waldo::compare(des_lupashin, des_lupashin2) # same: different order, but otherwise should be the same.


des_zhan <- make_design(targets = norm_zhan$targets,
                        group_column = "group",
                        factor_columns = NULL,
                        paired_column = NULL)


des_zhan2 <- make_design(targets = norm_zhan$targets,
                        group_column = "group",
                        factor_columns = NULL,
                        paired_column = NULL,
                        design_formula = "~0 + group")
waldo::compare(des_zhan, des_zhan2) # switched order again, but looks the same


des_reb <- make_design(targets = norm_reb$targets,
                           group_column = "group",
                           factor_columns = NULL,
                           paired_column = NULL)

des_kaul <- make_design(targets = norm_kaul$targets,
                           group_column = "group",
                           factor_columns = NULL,
                           paired_column = NULL,
                           design_formula = "~group+gender")

des_kaul2 <- make_design(targets = norm_kaul$targets,
                        group_column = "group",
                        factor_columns = "gender",
                        paired_column = NULL)
waldo::compare(des_kaul, des_kaul2) # order and levels a little different by defauls, but looks OK


# Make contrasts ----------------------------------------------------------
# Higgs
# No higgs contrast file??

# Ndu
contrasts_ndu_kidney <- make_contrasts(file = "for_testing/Example Data/NDu_030822_DIA/input_files/kidney_contrasts.txt",
                                       design = des_ndu$design)
contrasts_ndu_brain <- make_contrasts(file = "for_testing/Example Data/NDu_030822_DIA/input_files/brain_contrasts.txt",
                                      design = des_ndu$design)
contrasts_ndu_intestine <- make_contrasts(file = "for_testing/Example Data/NDu_030822_DIA/input_files/intestine_contrasts.txt",
                                          design = des_ndu$design)

#kaul
contrasts_kaul <- make_contrasts(file = "for_testing/Example Data/kaul/contrasts_designfomula.csv",
                                          design = des_kaul2$design)
# Lupashin
# contrasts_lupashin <- make_contrasts(file = "for_testing/Example Data/lupashin_030222/contrasts_bad.csv",
#                                      design = des_lupashin$design)
contrasts_lupashin <- make_contrasts(file = "for_testing/Example Data/lupashin_030222/contrasts.csv",
                                     design = des_lupashin$design)


# Zhan
contrasts_zhan <- make_contrasts(file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/contrasts.txt",
                                 design = des_zhan$design)

contrasts_zhan2 <- make_contrasts(file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/contrasts.txt",
                                 design = des_zhan2$design)
waldo::compare(contrasts_zhan, contrasts_zhan2)

# Rebello
contrasts_rebello <- make_contrasts(file = "for_testing/Example Data/rebello/contrasts.csv",
                                    design = des_reb$design)


# Run the analysis --------------------------------------------------------
# Splitting up functionality
# First, fit the model

fit_lupashin <- fit_limma_model(data = norm_lupashin$normList[["vsn"]],
                                design_obj = des_lupashin,
                                contrasts_obj = contrasts_lupashin)

fit_ndu_brain <- fit_limma_model(data = norm_ndu$normList[["vsn"]],
                                 design_obj = des_ndu,
                                 contrasts_obj = contrasts_ndu_brain)

fit_ndu_intestine <- fit_limma_model(data = norm_ndu$normList[["vsn"]],
                                     design_obj = des_ndu,
                                     contrasts_obj = contrasts_ndu_intestine)

fit_ndu_kidney <- fit_limma_model(data = norm_ndu$normList[["vsn"]],
                                  design_obj = des_ndu,
                                  contrasts_obj = contrasts_ndu_kidney)

fit_zhan <- fit_limma_model(data = norm_zhan$normList[["vsn"]],
                            design_obj = des_zhan,
                            contrasts_obj = contrasts_zhan)

fit_zhan2 <- fit_limma_model(data = norm_zhan$normList[["vsn"]],
                            design_obj = des_zhan2,
                            contrasts_obj = contrasts_zhan2)

waldo::compare(fit_zhan, fit_zhan2)


fit_reb <- fit_limma_model(data = norm_reb$normList[["vsn"]],
                           design_obj = des_reb,
                           contrasts_obj = contrasts_rebello)

fit_kaul <- fit_limma_model(data = norm_kaul$normList[["cycloess"]],
                           design_obj = des_kaul,
                           contrasts_obj = contrasts_kaul)
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
