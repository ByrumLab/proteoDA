
# create project directory in /home/uams/analysis/PI_DATE
#### in TERMINAL !!!!!!!!!!!!!!!!!

gsutil -m cp -r gs://uams-output/project_info_2021/Gao_102021_TMT/Analysis/ . 


setwd("/home/uams/analysis/PI_DATE")

##-------------------------
##   NEW FUNCTIONS
##-------------------------

source("/home/uams/uams-rscripts/current/functions25.R")


## DIA : BALACHANDRAN : (NS=24 | NP=3 | NG=6)
pro.file  = "/home/uams/analysis/Kaul_100621/Samples Report of Kaul_100621.csv"
meta.file = "/home/uams/analysis/Kaul_100621/Kaul_100621_metafile.csv"
contrast.file = "/home/uams/analysis/Kaul_100621/contrasts.csv"



## EXTRACT DATA
ext <- extract_data(file=pro.file, pipe="DIA", enrich="protein")


## META DATA
meta <- make_targets(file=meta.file, sampleIDs=colnames(ext$data), 
                     pipe="DIA", enrich="protein")


## REMOVE OUTLIER SAMPLES / GROUPS / BATCHES / FACTORS FROM ANALYSIS
table(meta$group)
sub <- subset_targets(targets=meta, factor="group",rm.vals="Pool")
table(sub$targets$group)

## MANUALLY SUBSET METADATA TO REMOVE POOL SAMPLES, TMT BATCHES ETC.
# sub2<-meta[meta$group %in% "Pool"==FALSE,];dim(sub2)
# sub2<-sub2[sub2$batch %in% c(2,4)==FALSE,];dim(sub2)

## FILTER/NORMALIZE DATA
table(sub$targets$group)
norm <-process_data(data=ext$data, targets=sub$targets, min.reps=3,min.grps=1)

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

## select certain contrasts for a complex design!!!!!!!!!! run multiple different comparisons. 
############ Special change for multiple contrasts listed in same file
# contrasts_full <- read.csv("contrasts.csv", header=FALSE)
#    sink(file="contrast1.csv") 
#     cat(contrasts_full[1,1])
#    sink()
# 
#    sink(file="contrast2.csv")  # gender comparison
#     cat(contrasts_full[2,1])
#    sink()
# 
#    sink(file="contrast3.csv")  # group2 = separate groups by gender
#     cat(contrasts_full[3:6,1], sep ="\n")
#    sink()
# 
# ## CONTRASTS
# contrast.file = "/home/uams/analysis/Kaul_100621/contrast1.csv"
# contrasts <- make_contrasts(file=contrast.file, design=des$design);contrasts


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
                          pipe       = "DIA",
                          enrich     = "protein",
                          dir        = NULL, ## NULL creates 02_diff_exprssion directory
                          save       = TRUE,
                          ilab       = "Balachandran_061721"  ## PI_DATE
)



## FORMATTED LIMMA RESULTS
wb<-openxlsx::createWorkbook()
## add protein results worksheet
add_limma_results(wb=wb, statList=lim$statList, annot=lim$annot,data=lim$data, norm.method=norm.method, 
                  min.pval=lim$param$Value$min.pval, min.lfc=lim$param$Value$min.lfc, 
                  pipe=lim$param$Value$pipe, enrich=lim$param$Value$enrich)

filename<-paste0(lim$param$Value$ilab,"_Results.xlsx");filename
saveWorkbook(wb=wb,file=file.path("./protein_analysis",filename), overwrite=TRUE)

########## run in terminal and not console!!!!!!!!!!

# copy Scaffold file to the protein_analysis folder from the prot_fs
# sitting in /home/uams/analysis/PI_DATE
cp /home/uams/prot_fs/projects/Kaul_100621/Kaul_100621.sdia /home/uams/analysis/Kaul_100621/protein_analysis

# must be sitting in the directory with the "protein_analysis" folder visible 
# Ctrl+Alt+enter to send to terminal
python /home/uams/zip/createZip.py Kaul 100621

# copy .zip file to uams-enduser bucket
gsutil -m cp Kaul_100621.zip gs://uams-enduser

# copy everything to uams-output
# from working directory
gsutil -m cp -r /home/uams/analysis/Kaul_100621/ gs://uams-output/project_info_2021/Kaul_100621_DIA/Analysis 




