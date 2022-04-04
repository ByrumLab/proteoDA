
# create project directory in /home/uams/analysis/PI_DATE
#### in TERMINAL !!!!!!!!!!!!!!!!!

gsutil -m cp -r gs://uams-output/project_info_2021/Gao_102021_TMT/Analysis/ . 

setwd("/home/uams/analysis/PI_DATE")

setwd("/home/uams/analysis/Gao_102021/Analysis")


##-------------------------
##   NEW FUNCTIONS
##-------------------------

source("/home/uams/uams-rscripts/current/functions25.R")



## TMT : FEDOROV : (NB=1 | PX=10 | NS=10 | NP=0 | NG=2)
pro.file  = "/home/uams/analysis/Gao_102021/Analysis/txt/proteinGroups.txt"
meta.file = "/home/uams/analysis/Gao_102021/Analysis/Gao_102021_metafile.csv"
contrast.file = "/home/uams/analysis/Gao_102021/Analysis/contrasts.txt"



## EXTRACT DATA
ext <- extract_data(file=pro.file, pipe="TMT", enrich="protein")


## META DATA
meta <- make_targets(file=meta.file, sampleIDs=colnames(ext$data), 
                     pipe="TMT", enrich="protein")


## REMOVE OUTLIER SAMPLES / GROUPS / BATCHES / FACTORS FROM ANALYSIS
table(meta$group)
# remove pool
 sub <- subset_targets(targets=meta, factor="group",rm.vals="Pool")
 
 # subset by batch
 # table(sub$targets$group);table(sub$targets$batch)
 # sub <- subset_targets(targets=sub$targets, factor="batch",rm.vals=c(3,2,4))


## MANUALLY SUBSET METADATA TO REMOVE POOL SAMPLES, TMT BATCHES ETC.
# sub2<-meta[meta$group %in% "Pool"==FALSE,];dim(sub2)
# sub2<-sub2[sub2$batch %in% c(2,4)==FALSE,];dim(sub2)

## FILTER/NORMALIZE DATA
table(meta$group)
#norm <-process_data(data=ext$data, targets=meta, min.reps=2,min.grps=1)

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
               groups=norm$targets$group, batch=norm$targets$batch,
               sampleLabels=NULL, box.inset=-0.3, vio.inset=-0.3, pca.inset=-0.3,
               legend=FALSE,enrich="protein",dir=NULL,save=TRUE)

## DESIGN MATRIX
#des <- make_design(targets=norm$targets, group="group", factors=NULL, paired=NULL)
#des$designformula

# tmt batch

des <- make_design(targets=norm$targets, group="group", factors="batch", paired=NULL)
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
                          pipe       = "TMT",
                          enrich     = "protein",
                          dir        = NULL, ## NULL creates 02_diff_exprssion directory
                          save       = TRUE,
                          ilab       = "Gao_102021"  ## PI_DATE
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

mv Gao_102021.sf3 protein_analysis

# must be sitting in the directory with the "protein_analysis" folder visible 
# Ctrl+Alt+enter to send to terminal
python /home/uams/zip/createZip.py Gao 102021

# copy .zip file to uams-enduser bucket
gsutil -m cp Gao_102021.zip gs://uams-enduser

# copy everything to uams-output
# from working directory
gsutil -m cp -r /home/uams/analysis/Gao_102021/ gs://uams-output/project_info_2021/Gao_102021_TMT/Analysis 

