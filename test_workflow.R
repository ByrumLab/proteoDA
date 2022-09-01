
library(pkgnet)
library(tictoc)
library(microbenchmark)
library(devtools)
library(profvis)

CreatePackageReport("proteomicsDIA")
# source("bin/functions26_CURRENT_031022.r")

# making some temporary changes
# make a change again

# Pipeline ----------------------------------------------------------------


# extract_data on a bunch of files --------------------------------

ext_bart <- extract_data("for_testing/Example Data/04_Bartholomew_101520_DIA/Samples Report of Bartholomew_101520.csv",
             pipe = "DIA",
             enrich = "protein") # Missing exclusivity col


ext_higgs <- extract_data("for_testing/Example Data/09_Higgs_072721_DIA_AG/Samples Report of Higgs_072721.csv",
             pipe = "DIA",
             enrich = "protein") # Worked


ext_ndu <- extract_data("for_testing/Example Data/NDu_030822_DIA/input_files/Samples Report of Du_030822.csv",
             pipe = "DIA",
             enrich = "protein") # Worked


ext_kinter <- extract_data(file = "for_testing/Example Data/19_Kinter_120720_TMT_DIA_AG/Kinter_120720_DIA/Samples Report of Kinter_DIA_022521.csv",
                            pipe = "DIA",
                            enrich = "protein") # Worked

ext_lupashin <- extract_data(file = "for_testing/Example Data/lupashin_030222/Samples Report of Lupashin_030222.csv",
                             pipe = "DIA",
                             enrich = "protein") # Worked

ext_zhan <- extract_data(file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/Samples Report of Zhan_111821_Experiment.csv",
                         pipe = "DIA",
                         enrich = "protein")

ext_reb <- extract_data(file = "for_testing/Example Data/rebello/Samples Report of Rebello_040522.csv",
                       pipe = "DIA",
                       enrich = "protein")

# Make targets ------------------------------------------------------------

target_higgs <- make_targets(file = "for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv",
                             sampleIDs = colnames(ext_higgs$data),
                             pipe = "DIA",
                             enrich = "protein") # worked

target_higgs2 <- make_targets(#file = "for_testing/Example Data/09_Higgs_072721_DIA_AG/metadata.csv",
                             sampleIDs = colnames(ext_higgs$data),
                             pipe = "DIA",
                             enrich = "protein")

target_ndu <- make_targets(file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                           sampleIDs = colnames(ext_ndu$data),
                           pipe = "DIA",
                           enrich = "protein") # worked

# target_ndu2 <- make_targets(#file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
#                            sampleIDs = colnames(ext_ndu$data),
#                            pipe = "DIA",
#                            enrich = "protein")
#
# target_ndu3 <- make_targets(file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
#                            sampleIDs = c(colnames(ext_ndu$data), "sample_45"),
#                            pipe = "DIA",
#                            enrich = "protein") # Gave warning as expected
#
# target_ndu4 <- make_targets(file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
#                             sampleIDs = colnames(ext_higgs$data),
#                             pipe = "DIA",
#                             enrich = "protein") # Gave error as expected,
#                             # but need to improve the error message and figure this logic out

target_lupashin <- make_targets(file = "for_testing/Example Data/lupashin_030222/Lupashin_030222_metafile_DIA.csv",
                                sampleIDs = colnames(ext_lupashin$data),
                                pipe = "DIA",
                                enrich = "protein") # Worked

target_zhan <- make_targets(file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/Zhan_111821_DIA_metadata.csv",
                            sampleIDs = colnames(ext_zhan$data),
                            pipe = "DIA",
                            enrich = "protein")

target_reb <- make_targets(file = "for_testing/Example Data/rebello/Rebello_040522_metafile_DIA.csv",
                           sampleIDs = colnames(ext_reb$data),
                           pipe = "DIA",
                           enrich = "protein")

# Subset targets --------------------------------------------------------------
sub_higgs <- subset_targets(targets = target_higgs,
                            filter_list = list(group = "Pool"))

sub_ndu <- subset_targets(targets = target_ndu,
                          filter_list = list(group = "Pool"))

sub_lupashin <- subset_targets(targets = target_lupashin,
                               filter_list = list(group = c("Pool", "input")))

sub_zhan <- subset_targets(targets = target_zhan,
                           filter_list = list(group = "Pool"))


sub_reb <- subset_targets(targets = target_reb,
                          filter_list = list(group = "Pool"))


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

tic()
norm_zhan <- process_data(data = ext_zhan$data,
                          targets = sub_zhan$targets,
                          min.reps = 13,
                          min.grps = 2)
toc()

norm_reb <- process_data(data = ext_reb$data,
                         targets = sub_reb$targets,
                         min.reps = 2,
                         min.grps = 1)

# Normalization report ----------------------------------------------------
# Higgs
make_proteinorm_report(normList = norm_higgs$normList,
                       groups = norm_higgs$targets$group,
                       file = "higgs_update.pdf", overwrite = T)

# Ndu
make_proteinorm_report(normList = norm_ndu$normList,
                       groups = norm_ndu$targets$group,
                       file = "ndu_update.pdf", overwrite = T)

# Lupashin
make_proteinorm_report(normList=norm_lupashin$normList,
                       groups = norm_lupashin$targets$group,
                       file = "lupashin_update.pdf", overwrite = T, suppress_zoom_legend = T)
# Zhan
tic()
make_proteinorm_report(normList = norm_zhan$normList,
                       groups = norm_zhan$targets$group,
                       file = "zhan_update.pdf", overwrite = T)
toc()

# Rebello
make_proteinorm_report(normList = norm_reb$normList,
                       groups = norm_reb$targets$group,
                       file = "rebello_update.pdf",
                       overwrite = T)

# Make QC report ----------------------------------------------------------
write_qc_report(processed_data = norm_higgs,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                enrich = "protein",
                file = "higgs_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_ndu,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                enrich = "protein",
                file = "ndu_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_lupashin,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                enrich = "protein",
                file = "lupashin_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_zhan,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                enrich = "protein",
                file = "zhan_qc_update.pdf",
                overwrite = T)

write_qc_report(processed_data = norm_reb,
                chosen_norm_method = "vsn",
                grouping_column = "group",
                enrich = "protein",
                file = "rebello_qc_update.pdf",
                overwrite = T)



make_qc_report(normList = norm_ndu$normList, norm.meth = "vsn",
               groups = norm_ndu$targets$group,
               enrich = "protein", save = TRUE, file = "ndu_qc.pdf", overwrite = T)

make_qc_report(normList = norm_ndu$normList, norm.meth = "vsn",
               groups = norm_ndu$targets$group, legend = F,
               enrich = "protein", save = TRUE, file = "ndu_qc_no_legend.pdf", overwrite = T)


make_qc_report(normList = norm_lupashin$normList, norm.meth = "vsn",
               groups = norm_lupashin$targets$group,
               enrich = "protein", save = TRUE, file = "lupashin_qc.pdf", overwrite = T)


make_qc_report(normList = norm_zhan$normList, norm.meth = "vsn",
               groups = norm_zhan$targets$group,
               enrich = "protein", save = TRUE, file = "zhan_qc.pdf", overwrite = T)


make_qc_report(normList = norm_reb$normList, norm.meth = "vsn",
               groups = norm_reb$targets$group,
               enrich = "protein", save = TRUE, file = "rebello_qc.pdf", overwrite = T)

make_qc_report(normList = norm_reb$normList, norm.meth = "vsn",
               groups = norm_reb$targets$group,
               enrich = "protein", save = F)


# Make design -------------------------------------------------------------
des_higgs<- make_design(targets=norm_higgs$targets,
                        group_column = "group",
                        factor_columns = NULL,
                        paired_column = NULL)

des_ndu <- make_design(targets = norm_ndu$targets,
                       group_column = "group",
                       factor_columns = NULL,
                       paired_column =  NULL)

des_lupashin <- make_design(targets = norm_lupashin$targets,
                            group_column = "group",
                            factor_columns = NULL,
                            paired_column = NULL)

des_zhan <- make_design(targets = norm_zhan$targets,
                        group_column = "group",
                        factor_columns = NULL,
                        paired_column = NULL)

des_reb <- make_design(targets = norm_reb$targets,
                           group_column = "group",
                           factor_columns = NULL,
                           paired_column = NULL)


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

# Lupashin
contrasts_lupashin <- make_contrasts(file = "for_testing/Example Data/lupashin_030222/contrasts_bad.csv",
                                     design = des_lupashin$design)
contrasts_lupashin <- make_contrasts(file = "for_testing/Example Data/lupashin_030222/contrasts.csv",
                                     design = des_lupashin$design)


# Zhan
contrasts_zhan <- make_contrasts(file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/contrasts.txt",
                                 design = des_zhan$design)

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


fit_reb <- fit_limma_model(data = norm_reb$normList[["vsn"]],
                           design_obj = des_reb,
                           contrasts_obj = contrasts_rebello)

# Extract results ---------------------------------------------------------
results_lupashin <- extract_limma_DE_results(limma_fit = fit_lupashin)
results_ndu_brain <- extract_limma_DE_results(limma_fit = fit_ndu_brain)
results_ndu_intestine <- extract_limma_DE_results(limma_fit = fit_ndu_intestine)
results_ndu_kidney <- extract_limma_DE_results(limma_fit = fit_ndu_kidney)
results_zhan <- extract_limma_DE_results(limma_fit = fit_zhan)
results_reb <- extract_limma_DE_results(limma_fit = fit_reb)




# Write results -----------------------------------------------------------
write_limma_results(model_results = results_lupashin,
                    norm.method = "vsn",
                    annotation = ext_lupashin$annot,
                    ilab = "Lupashin_82928",
                    enrich = "protein",
                    overwrite = T)


write_limma_results(model_results = results_ndu_brain,
                    norm.method = "vsn",
                    annotation = ext_ndu$annot,
                    ilab = "ndu_brain_82928",
                    enrich = "protein")

write_limma_results(model_results = results_ndu_intestine,
                    norm.method = "vsn",
                    annotation = ext_ndu$annot,
                    ilab = "ndu_intestine_82928",
                    enrich = "protein")

write_limma_results(model_results = results_ndu_kidney,
                    norm.method = "vsn",
                    annotation = ext_ndu$annot,
                    ilab = "ndu_kidney_82928",
                    enrich = "protein")

write_limma_results(model_results = results_reb,
                    norm.method = "vsn",
                    annotation = ext_reb$annot,
                    ilab = "Rebello_82928",
                    enrich = "protein")


write_limma_results(model_results = results_zhan,
                    norm.method = "vsn",
                    annotation = ext_zhan$annot,
                    ilab = "zhan_982974",
                    enrich = "protein")




# testing report making ---------------------------------------------------
make_limma_reports(model_results = results_reb,
                   annotation = ext_reb$annot,
                   groups = norm_reb$targets$group,
                   output_dir = "output_rebello")

make_limma_reports(model_results = results_reb,
                   annotation = ext_reb$annot,
                   groups = norm_reb$targets$group,
                   output_dir = "output_rebello_wide",
                   width = 1500,
                   height = 1500)


make_limma_reports(model_results = results_lupashin,
                   annotation = ext_lupashin$annot,
                   groups = norm_lupashin$targets$group,
                   output_dir = "output_lupashin")


make_limma_reports(model_results = results_ndu_brain,
                   annotation = ext_ndu$annot,
                   groups = norm_ndu$targets$group,
                   output_dir = "output_ndu_brain")


make_limma_reports(model_results = results_zhan,
                   annotation = ext_zhan$annot,
                   groups = norm_zhan$targets$group,
                   output_dir = "output_zhan")


make_limma_reports(model_results = results_zhan,
                   annotation = ext_zhan$annot,
                   groups = norm_zhan$targets$group,
                   output_dir = "output_zhan_wide",
                   width = 2000)




# Testing phospho ---------------------------------------------------------
# extract data
thomas_tmt <- extract_data(file = "for_testing/Example Data/Thomas_03922_phos/proteinGroups.txt",
             pipe = "TMT",
             enrich = "protein")

thomas_phospho <- extract_data(file = "for_testing/Example Data/Thomas_03922_phos/Phospho (STY)Sites.txt",
                           pipe = "phosphoTMT",
                           enrich = "phospho")


tar_tmt <- make_targets(file = "for_testing/Example Data/Thomas_03922_phos/Thomas_032922_metafile_pro.csv",
                        sampleIDs = colnames(thomas_tmt$data),
                        pipe = "phosphoTMT",
                        enrich = "protein")

tar_phospho <- make_targets(file = "for_testing/Example Data/Thomas_03922_phos/Thomas_032922_metafile_phos.csv",
                        sampleIDs = colnames(thomas_tmt$data),
                        pipe = "phosphoTMT",
                        enrich = "phospho")
sub_protein <- subset_targets(tar_tmt, filter_column = "group", rm.vals = "Pool")
sub_phospho <- subset_targets(tar_phospho, filter_column = "group", rm.vals = "Pool")


norm_prot <- process_data(data = thomas_tmt$data, targets = sub_protein$targets,
                          min.reps = 3, min.grps = 2)

colnames(thomas_phospho$data) <- stringr::str_replace(colnames(thomas_phospho$data), ".phospho", ".lysate")
norm_phospho <- process_data(data = thomas_phospho$data, targets = sub_phospho$targets,
                          min.reps = 3, min.grps = 2)


make_proteinorm_report(normList = norm_prot$normList, groups = norm_prot$targets$group,
                       enrich = "protein", file = "thomas_protein.pdf", overwrite = T)
make_proteinorm_report(normList = norm_phospho$normList, groups = norm_phospho$targets$group,
                       enrich = "phospho", file = "thomas_phospho.pdf", overwrite = T)

make_qc_report(normList = norm_prot$normList, groups = norm_prot$targets$group,
               norm.method = "log2",
               enrich = "protein",
               file = "thomas_protein_qc.pdf",
               overwrite = T)

make_qc_report(normList = norm_phospho$normList, groups = norm_phospho$targets$group,
               norm.method = "log2",
               enrich = "phospho",
               file = "thomas_phospho_qc.pdf",
               overwrite = T)


des_prot <- make_design(targets = norm_prot$targets, group_column = "group")
des_phos <- make_design(targets = norm_phospho$targets, group_column = "group")

contrast_prot <- make_contrasts("for_testing/Example Data/Thomas_03922_phos/contrasts.csv",
                                design = des_prot$design)
contrast_phos <- make_contrasts("for_testing/Example Data/Thomas_03922_phos/contrasts.csv",
                                design = des_phos$design)


fit_prot <- fit_limma_model(data = norm_prot$normList[["vsn"]],
                            targets = des_prot$targets,
                            design = des_prot$design,
                            contrasts = contrast_prot$contrasts)

fit_phos <- fit_limma_model(data = norm_phospho$normList[["vsn"]],
                            targets = des_phos$targets,
                            design = des_phos$design,
                            contrasts = contrast_phos$contrasts)

results_prot <- extract_limma_DE_results(fit_prot)
results_phos <- extract_limma_DE_results(fit_phos)

write_limma_results(results_prot,
                    annotation = thomas_tmt$annot,
                    ilab = "thomas_09819",
                    norm.method = "vsn",
                    pipe = "TMT", enrich = "protein")

write_limma_results(results_phos,
                    annotation = thomas_phospho$annot,
                    ilab = "thomas_09819",
                    norm.method = "vsn",
                    pipe = "phosphoTMT", enrich = "phospho")

