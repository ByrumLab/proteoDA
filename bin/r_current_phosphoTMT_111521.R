# create project directory in /home/uams/analysis/PI_DATE
#### in TERMINAL !!!!!!!!!!!!!!!!!

gsutil -m cp -r gs://uams-output/project_info_2021/Kim_102621_phos/Analysis/ . 

setwd("/home/uams/analysis/PI_DATE")

setwd("/home/uams/analysis/Kim_102621/Analysis")

##-------------------------
##   NEW FUNCTIONS
##-------------------------

source("/home/uams/uams-rscripts/current/functions25.R")


## PHOSPHO : BYRUM : (NB=4 | PX=10 | NS=36 | NP=4 | NG=10)
pro.file        = "/home/uams/analysis/Kim_102621/Analysis/txt/proteinGroups.txt"
pro.meta.file   = "/home/uams/analysis/Kim_102621/Analysis/Kim_102621_metafile_pro.csv"
contrast.file   = "/home/uams/analysis/Kim_102621/Analysis/contrasts.txt"



## EXTRACT DATA
ext <- extract_data(file=pro.file, pipe="phosphoTMT", enrich="protein")


## META DATA
meta <- make_targets(file=pro.meta.file, sampleIDs=colnames(ext$data), 
                     pipe="phosphoTMT", enrich="protein")


## REMOVE OUTLIER SAMPLES / GROUPS / BATCHES / FACTORS FROM ANALYSIS
table(meta$group)
sub <- subset_targets(targets=meta, factor="group",rm.vals="Pool")
table(sub$targets$group);table(sub$targets$batch)
sub <- subset_targets(targets=sub$targets, factor="batch",rm.vals=c(3,2,4))


## MANUALLY SUBSET METADATA TO REMOVE POOL SAMPLES, TMT BATCHES ETC.
# sub2<-meta[meta$group %in% "Pool"==FALSE,];dim(sub2)
# sub2<-sub2[sub2$batch %in% c(2,4)==FALSE,];dim(sub2)

## FILTER/NORMALIZE DATA
table(sub$targets$group)
norm <-process_data(data=ext$data, targets=sub$targets, min.reps=2,min.grps=1)

## PROTEINORM REPORT
table(norm$targets$batch)
make_norm_report(normList=norm$normList, groups=norm$targets$group,
                 batch=norm$targets$batch, sampleLabels=NULL, legend=TRUE,
                 enrich="protein", dir=NULL, save=TRUE)

## SELECT NORM METHOD
norm.methods ## available methods
norm.method = "vsn"

## MAKE QC REPORT
make_qc_report(normList=norm$normList, norm.method=norm.method, 
               groups=norm$targets$group, batch=NULL,
               sampleLabels=NULL, box.inset=-0.3, vio.inset=-0.3, pca.inset=-0.3,
               legend=FALSE,enrich="protein",dir=NULL,save=TRUE)

## DESIGN MATRIX
des <- make_design(targets=norm$targets, group="group", factors=NULL, paired=NULL)
des$designformula

## CONTRASTS
contrasts <- make_contrasts(file=contrast.file, design=des$design);contrasts

## LIMMA ANALYSIS
lim <- run_limma_analysis(data        = norm$normList[[norm.method]],
                           annot      = ext$annot[rownames(norm$normList[[norm.method]]),],
                           targets    = des$targets,
                           design     = des$design,
                           contrasts  = contrasts,
                           min.pval   = 0.055,
                           min.lfc    = 1,
                           adj.method = "BH",
                           paired     = FALSE,  ## TRUE if paired samples/mixed.effects model
                           pipe       = "phosphoTMT",
                           enrich     = "protein",
                           dir        = NULL, ## NULL creates 02_diff_exprssion directory
                           save       = TRUE,
                           ilab       = "Kim_102621"  ## PI_DATE
)





phos.file       = "/home/uams/analysis/Kim_102621/Analysis/txt/Phospho (STY)Sites.txt"
phos.meta.file  = "/home/uams/analysis/Kim_102621/Analysis/Kim_102621_metafile_phos.csv"
# combo.meta.file = "C:/Users/washamcharityl/Documents/PROJECTS/GOOGLE TESTING/TEST_PROJECT/Test_data/Byrum_022020_phos/sampleList3.combo.csv"
contrast.file   = "/home/uams/analysis/Kim_102621/Analysis/contrasts.txt"




## EXTRACT DATA
ext2 <- extract_data(file=phos.file, pipe="phosphoTMT", enrich="phospho")


## META DATA
meta2 <- make_targets(file=phos.meta.file, sampleIDs=colnames(ext2$data), 
                    pipe="phosphoTMT", enrich="phospho")


## REMOVE OUTLIER SAMPLES / GROUPS / BATCHES / FACTORS FROM ANALYSIS
table(meta2$group)
sub2 <- subset_targets(targets=meta2, factor="group",rm.vals="Pool")
table(sub2$targets$group);table(sub2$targets$batch)
sub2 <- subset_targets(targets=sub2$targets, factor="batch",rm.vals=c(3,2,4))


## MANUALLY SUBSET METADATA TO REMOVE POOL SAMPLES, TMT BATCHES ETC.
# sub2<-meta[meta$group %in% "Pool"==FALSE,];dim(sub2)
# sub2<-sub2[sub2$batch %in% c(2,4)==FALSE,];dim(sub2)

## FILTER/NORMALIZE DATA
table(sub2$targets$group)
norm2 <-process_data(data=ext2$data, targets=sub2$targets, min.reps=2,min.grps=1)

## PROTEINORM REPORT
make_norm_report(normList=norm2$normList, groups=norm2$targets$group,
               batch=norm2$targets$batch, sampleLabels=NULL, legend=TRUE,
               enrich="phospho", dir=NULL, save=TRUE)

## SELECT NORM METHOD
norm.methods ## available methods
norm.method2 = "vsn"

## MAKE QC REPORT
table(norm2$targets$batch)
make_qc_report(normList=norm2$normList, norm.method=norm.method2, 
               groups=norm2$targets$group, batch=NULL,
               sampleLabels=NULL, box.inset=-0.3, vio.inset=-0.3, pca.inset=-0.3,
               legend=FALSE,enrich="phospho",dir=NULL,save=TRUE)

## DESIGN MATRIX
des2 <- make_design(targets=norm2$targets, group="group", factors=NULL, paired=NULL)
des2$designformula

## CONTRASTS
contrasts2 <- make_contrasts(file=contrast.file, design=des2$design);contrasts2

## LIMMA ANALYSIS
lim2 <- run_limma_analysis(data      = norm2$normList[[norm.method2]],
                          annot      = ext2$annot[rownames(norm2$normList[[norm.method2]]),],
                          targets    = des2$targets,
                          design     = des2$design,
                          contrasts  = contrasts2,
                          min.pval   = 0.055,
                          min.lfc    = 1,
                          adj.method = "BH",
                          paired     = FALSE,  ## TRUE if paired samples/mixed.effects model
                          pipe       = "phosphoTMT",
                          enrich     = "phospho",
                          dir        = NULL, ## NULL creates 02_diff_exprssion directory
                          save       = TRUE,
                          ilab       = "Kim_102621"  ## PI_DATE
                          )




## FORMATTED LIMMA RESULTS
wb<-openxlsx::createWorkbook()
## add protein results worksheet
add_limma_results(wb=wb, statList=lim$statList, annot=lim$annot,data=lim$data, norm.method=norm.method, 
                  min.pval=lim$param$Value$min.pval, min.lfc=lim$param$Value$min.lfc, 
                  pipe=lim$param$Value$pipe, enrich=lim$param$Value$enrich)
## add phospho results worksheet
add_limma_results(wb=wb, statList=lim2$statList, annot=lim2$annot,data=lim2$data, norm.method=norm.method2, 
                  min.pval=lim2$param$Value$min.pval, min.lfc=lim2$param$Value$min.lfc, 
                  pipe=lim2$param$Value$pipe, enrich=lim2$param$Value$enrich)

filename<-paste0(lim$param$Value$ilab,"_Results.xlsx");filename
saveWorkbook(wb=wb,file=file.path("./protein_analysis",filename), overwrite=TRUE)


########## run in terminal and not console!!!!!!!!!!

# copy Scaffold file to the protein_analysis folder from the prot_fs
# sitting in /home/uams/analysis/PI_DATE
cp /home/uams/prot_fs/projects/Kaul_100621/Kaul_100621.sdia /home/uams/analysis/Kaul_100621/protein_analysis

mv Kim_102621.sf3 protein_analysis

# must be sitting in the directory with the "protein_analysis" folder visible 
# Ctrl+Alt+enter to send to terminal
python /home/uams/zip/createZip.py Kim 102621

# copy .zip file to uams-enduser bucket
gsutil -m cp Kim_102621.zip gs://uams-enduser

# copy everything to uams-output
# from working directory
gsutil -m cp -r /home/uams/analysis/Kim_102621/Analysis/ gs://uams-output/project_info_2021/Kim_102621_phos/Analysis 

