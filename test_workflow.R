
library(pkgnet)
library(tictoc)
library(microbenchmark)
library(devtools)

CreatePackageReport("proteomicsDIA")
# source("bin/functions26_CURRENT_031022.r")

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


ext_kintler <- extract_data(file = "for_testing/Example Data/19_Kinter_120720_TMT_DIA_AG/Kinter_120720_DIA/Samples Report of Kinter_DIA_022521.csv",
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

target_ndu2 <- make_targets(#file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                           sampleIDs = colnames(ext_ndu$data),
                           pipe = "DIA",
                           enrich = "protein")

target_ndu3 <- make_targets(file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                           sampleIDs = c(colnames(ext_ndu$data), "sample_45"),
                           pipe = "DIA",
                           enrich = "protein") # Gave warning as expected

target_ndu4 <- make_targets(file = "for_testing/Example Data/NDu_030822_DIA/input_files/Du_030822_metafile_DIA.csv",
                            sampleIDs = colnames(ext_higgs$data),
                            pipe = "DIA",
                            enrich = "protein") # Gave error as expected,
                            # but need to improve the error message and figure this logic out

target_kintler <- make_targets(sampleIDs = colnames(ext_kintler$data),
                               pipe = "DIA",
                               enrich = "protein") # error because of no "Sample" in ID

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
sub_higgs <- subset_targets(targets=target_higgs, filter_column = "group", rm.vals = "Pool")

sub_ndu <- subset_targets(targets=target_ndu, filter_column = "group", rm.vals = "Pool")

sub_lupashin <- subset_targets(targets=target_lupashin, filter_column = "group", rm.vals = c("Pool", "input"))

sub_zhan <- subset_targets(targets = target_zhan, filter_column = "group", rm.vals = "pool")

sub_reb <- subset_targets(targets = target_reb, filter_column = "group", rm.vals = "Pool")


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
make_proteinorm_report(normList = norm_higgs$normList, groups = norm_higgs$targets$group, file = "higgs.pdf", overwrite = T)

# Ndu
make_proteinorm_report(normList = norm_ndu$normList, groups = norm_ndu$targets$group, file = "ndu.pdf", overwrite = T)

# Lupashin
make_proteinorm_report(normList=norm_lupashin$normList, groups = norm_lupashin$targets$group, file = "lupashin.pdf",
                       overwrite = T)
make_proteinorm_report(normList=norm_lupashin$normList, groups = norm_lupashin$targets$group, save = F)

# Zhan
tic()
make_proteinorm_report(normList = norm_zhan$normList, groups = norm_zhan$targets$group, file = "zhan.pdf")
toc()

make_proteinorm_report(normList = norm_reb$normList, groups = norm_reb$targets$group, file = "rebello.pdf", overwrite = T)



# Make QC report ----------------------------------------------------------
make_qc_report(normList = norm_higgs$normList, norm.method = "vsn",
               groups = norm_higgs$targets$group,
               batch = norm_higgs$targets$group,
               enrich = "protein", save = TRUE, file = "higgs_qc.pdf", overwrite = T)

make_qc_report(normList = norm_higgs$normList, norm.method = "vsn",
               groups = norm_higgs$targets$group,
               batch = norm_higgs$targets$group,
               enrich = "protein", save = FALSE)


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
contrasts_lupashin <- make_contrasts(file = "for_testing/Example Data/lupashin_030222/contrasts.csv",
                                     design = des_lupashin$design)
contrasts_lupashin <- make_contrasts(file = "for_testing/Example Data/lupashin_030222/contrasts_bad.csv",
                                     design = des_lupashin$design)

# Zhan
contrasts_zhan <- make_contrasts(file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/contrasts.txt",
                                 design = des_zhan$design)
# Rebello
contrasts_rebello <- make_contrasts(file = "for_testing/Example Data/rebello/contrasts.csv",
                                    design = des_reb$design)


test <- norm_reb$targets
test <- rbind(test, "test" = rep(NA, dim(test)[2]))

# Run the analysis --------------------------------------------------------
# Splitting up functionality
# First, fit the model

fit_lupashin <- fit_limma_model(data = norm_lupashin$normList[["vsn"]],
                                 targets = des_lupashin$targets,
                                 design = des_lupashin$design,
                                 contrasts = contrasts_lupashin$contrasts)

fit_ndu_brain <- fit_limma_model(data = norm_ndu$normList[["vsn"]],
                                targets = des_ndu$targets,
                                design = des_ndu$design,
                                contrasts = contrasts_ndu_brain$contrasts)
fit_ndu_intestine <- fit_limma_model(data = norm_ndu$normList[["vsn"]],
                                 targets = des_ndu$targets,
                                 design = des_ndu$design,
                                 contrasts = contrasts_ndu_intestine$contrasts)
fit_ndu_kidney <- fit_limma_model(data = norm_ndu$normList[["vsn"]],
                                 targets = des_ndu$targets,
                                 design = des_ndu$design,
                                 contrasts = contrasts_ndu_kidney$contrasts)

fit_zhan <- fit_limma_model(data = norm_zhan$normList[["vsn"]],
                            targets = des_zhan$targets,
                            design = des_zhan$design,
                            contrasts = contrasts_zhan$contrasts)


fit_reb <- fit_limma_model(data = norm_reb$normList[["vsn"]],
                           targets = des_reb$targets,
                           design = des_reb$design,
                           contrasts = contrasts_rebello$contrasts)



# Some testing of next steps ----------------------------------------------
z <- extract_limma_DE_results(limma_fit = fit_reb)


names(z)


efit <- fit_reb$eBayes_fit
adj.method <- "BH"
min.pval <- 0.05
min.lfc <- 1


## DECIDE TESTS
print(paste("extracting limma stat results for each comparison (statList)..."))
contrastNames <- colnames(efit$coefficients)
dt <- limma::decideTests(efit, adjust.method = adj.method, p.value = min.pval, lfc = min.lfc)
dtp <- limma::decideTests(efit, adjust.method = "none", p.value = min.pval, lfc = min.lfc)
sum.dt <- summary(dt)
sum.dtp <- summary(dtp)


stats <- limma::topTable(efit, coef = contrastNames[1], number = Inf, adjust.method = adj.method, sort = "none", p.value = 1,
                lfc = 0, confint = T)

cbind(dtp[, contrastNames[1]], dt[, contrastNames[1]])

## STAT RESULTS (statList)
statList <- list()
limmaStatColums <- c("logFC", "CI.L", "CI.R", "AveExpr", "t", "B", "P.Value", "adj.P.Val")
statList <- base::lapply(contrastNames, function(x) {
  stats <- limma::topTable(efit,
                           coef = x, number = Inf, adjust.method = adj.method,
                           sort.by = "none", p.value = 1, lfc = 0, confint = TRUE
  )
  df <- cbind(dtp[, x], dt[, x])
  colnames(df) <- c("sig.PVal", "sig.FDR")
  stats <- cbind(stats[, limmaStatColums], df[rownames(stats), ])
})
names(statList) <- contrastNames
print(paste("statList created. Success!!"))



test_statList <- statList[[4]]

dt[,names(statList)[4]][dt[,names(statList)[4]] > 1,]
dt[rownames(dt) %in% rownames(x$sig),]

x <- get_de(stats = test_statList, min.pval = min.pval, min.lfc = min.lfc, type = "p.adj", de.type = "limma")
x
x$sig






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
                       enrich = "protein", file = "thomas_protein.pdf")
make_proteinorm_report(normList = norm_phospho$normList, groups = norm_phospho$targets$group,
                       enrich = "phospho", file = "thomas_phospho.pdf")

make_qc_report(normList = norm_prot$normList, groups = norm_prot$targets$group,
               norm.method = "log2",
               enrich = "protein",
               file = "thomas_protein_qc.pdf")

make_qc_report(normList = norm_phospho$normList, groups = norm_phospho$targets$group,
               norm.method = "log2",
               enrich = "phospho",
               file = "thomas_phospho_qc.pdf")


des_prot <- make_design(targets = norm_prot$targets, group_column = "group")
des_phos <- make_design(targets = norm_phospho$targets, group_column = "group")

contrast_prot <- make_contrasts("for_testing/Example Data/Thomas_03922_phos/contrasts.csv",
                                design = des_prot$design)
contrast_phos <- make_contrasts("for_testing/Example Data/Thomas_03922_phos/contrasts.csv",
                                design = des_phos$design)
