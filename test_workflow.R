
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
                               enrich = "protein") # error because of now "Sample" in ID

target_lupashin <- make_targets(file = "for_testing/Example Data/lupashin_030222/Lupashin_030222_metafile_DIA.csv",
                                sampleIDs = colnames(ext_lupashin$data),
                                pipe = "DIA",
                                enrich = "protein") # Worked

target_zhan <- make_targets(file = "for_testing/Example Data/Zhan_DIA_217_samples/input_files/Zhan_111821_DIA_metadata.csv",
                            sampleIDs = colnames(ext_zhan$data),
                            pipe = "DIA",
                            enrich = "protein")

# Subset targets --------------------------------------------------------------
sub_higgs <- subset_targets(targets=target_higgs, filter_column = "group", rm.vals = "Pool")

sub_ndu <- subset_targets(targets=target_ndu, filter_column = "group", rm.vals = "Pool")

sub_lupashin <- subset_targets(targets=target_lupashin, filter_column = "group", rm.vals = c("Pool", "input"))

sub_zhan <- subset_targets(targets = target_zhan, filter_column = "group", rm.vals = "pool")




####### TO FIGURE OUT:
# I made new rlrNorm and giNorm functions. Did a lot of testing with them,
# and they seemed OK. though, I may have messed up the tolerance on my
# waldo compares
# In any case, I'm running into a weird bug:
# when I just use the decontaminated data as input, e.g.,
# with: ext_higg$data,
# the old and new giNorm functions are apparently giving me basically the same
# answer (down to ~15 decimal places):
waldo::compare(giNorm_old(logDat = logNorm(ext_higgs$data)),
               giNorm(dat = ext_higgs$data))
old <- giNorm_old(logDat = logNorm(ext_higgs$data))
new <- giNorm(dat = ext_higgs$data)
mean(new-old, na.rm = T) # basically 0.
# But, when I started running the proteionorm report function with this new function,
# I got very weird looking log2 ratios and the normalization metrics looked different
# So, I added the old normalization back into the process_data function and ran it:
norm_higgs <- process_data(data = ext_higgs$data,
                           targets = sub_higgs$targets,
                           min.reps = 5,
                           min.grps = 3)
mean(norm_higgs$normList$gi - norm_higgs$normList$gi2) # difference of 0.02!!
# So, somethign weird is going on here, and I don't know what
# In the old way, you longNorm and then take 2^ again to reverse ie:
# maybe numerical imprecision there? Should try doing some print statements
# in the bodies of the two functions, see if they are calculating the same
# col sums and getting the same median col sum.

# I modified the functions to output the filtered data that goes into the normalization
# steps:
old <- giNorm_old(logNorm(norm_higgs$filtered_data))
new <- giNorm(norm_higgs$filtered_data)
waldo::compare(old, new)
max(new - old) # biggest difference is to 15 decimals....
mean(new-old, na.rm = T) # again, the difference is 0 if I do it this way
# SO WHY ARE THEY DIFFERENT COMING OUT OF process_data?!?!?!?!?!
# I am fried and cannot think of a reason why? Something going wrong with
# row or column sorting????
# Add back in a ginorm_old version that just calls the normalyzer DE, not the
# underlying code, and see if tha tmakes a difference?
nrow(norm_higgs$filtered_data) < nrow(ext_higgs$data)

# Process data ------------------------------------------------------------
norm_higgs <- process_data(data = ext_higgs$data,
                           targets = sub_higgs$targets,
                           min.reps = 5,
                           min.grps = 3)
giNorm_old(2^logNorm(norm_higgs$filtered_data))
giNorm_old(logNorm(norm_higgs$filtered_data)) == giNorm(2^logNorm(norm_higgs$filtered_data))
mean(giNorm_old(logNorm(norm_higgs$filtered_data)) - giNorm(norm_higgs$filtered_data))

waldo::compare(2^logNorm(norm_higgs$filtered_data), as.matrix(norm_higgs$filtered_data))

mean(norm_higgs$normList$gi - norm_higgs$normList$gi2)

norm_ndu <- process_data(data = ext_ndu$data,
                         targets = sub_ndu$targets,
                         min.reps = 3,
                         min.grps = 1)
mean(norm_ndu$normList$gi - norm_ndu$normList$gi2, na.rm = T)

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


old_giNorm <- function(logDat) {
  giNormed <- NormalyzerDE::globalIntensityNormalization(as.matrix(logDat), noLogTransform = TRUE)
  colnames(giNormed) <- colnames(logDat)
  row.names(giNormed) <- rownames(logDat)
  return(as.matrix(giNormed))
}

old <- old_giNorm(logDat = logNorm(ext_higgs$data))
new <- giNorm(dat = ext_higgs$data)

colnames(dim(old))
colnames(dim(new))
unique(rownames(old) == rownames(new))
str(old)
str(new)

mean(new[old != new] - old[old != new], na.rm = T )

waldo::compare(old_giNorm(logDat = logNorm(ext_higgs$data)),
               giNorm(dat = ext_higgs$data))
sum(is.na(old))
sum(is.na(new))
typeof(new)
typeof(old)

class(new)
class(old)
giNorm

norm_higgs$normList$gi <- old

waldo::compare(norm_higgs$normList$gi, old)


# Normalization report ----------------------------------------------------
# Higgs
make_proteinorm_report(normList = norm_higgs$normList, groups = norm_higgs$targets$group, file = "higgs.pdf")
# Ndu
make_proteinorm_report(normList = norm_ndu$normList, groups = norm_ndu$targets$group, file = "ndu.pdf")

# Lupashin
make_proteinorm_report(normList=norm_lupashin$normList, groups = norm_lupashin$targets$group, file = "lupashin.pdf")
make_proteinorm_report(normList=norm_lupashin$normList, groups = norm_lupashin$targets$group, save = F)

# Zhan
tic()
make_proteinorm_report(normList = norm_zhan$normList, groups = norm_zhan$targets$group, file = "zhan.pdf")
toc()

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


