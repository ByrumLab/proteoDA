make_proteinorm_report <- function(normList, groups=NULL, batch=NULL, sampleLabels=NULL, legend=TRUE,
                                   enrich=c("protein","phospho"), dir=NULL, file=NULL, save=FALSE, keep.png=FALSE){

  enrich <- match.arg(enrich,choices=c("protein","phospho"), several.ok=FALSE)

  if(is.null(sampleLabels)){ sampleLabels<-colnames(normList[[1]]) };sampleLabels
  if(is.null(groups)){ groups <- rep("group",ncol(normList[[1]])) }
  groups<-make_factor(x=as.character(groups));groups
  if(!is.null(batch)){ batch <- make_factor(as.character(batch)) };batch
  if(save==FALSE){ keep.png<-FALSE; dir <- "." }

  ## CREATE QC OUTPUT DIRECTORY
  ## if use enrich type to create QC output directory
  if(save==TRUE){

    if(is.null(dir)){
      if(enrich=="protein"){ dir <- file.path("protein_analysis","01_quality_control") }
      if(enrich=="phospho"){ dir <- file.path("phospho_analysis","01_quality_control") }
    } ## DIR NULL
    if(!is.null(dir)){ if(!dir.exists(dir)){ dir.create(dir, recursive=TRUE) }};dir
    print(paste("QC output directory:", dir))

    ## FILENAME DEFINED
    if(!is.null(file)){
      if(file_ext(file) != "pdf"){
        stop("\nError! Invalid output file type...\nincorrect: file = '",file,"'",
             "\ncorrect:   file = 'proteiNorm_Report.pdf'") }
      pngdir=gsub(".pdf","",file)
      if(file.exists(file.path(dir,file))){
        file<-make_new_filename(x=file,dir=dir);file
        no<-sub(".pdf","",sub(".*_","",file));no
        pngdir=gsub(".pdf","",file)
      }
      if(!file.exists(file.path(dir,file))){ pngdir=gsub(".pdf","",file) }

    } ## FILE NOT NULL

    ## FILENAME NULL
    if(is.null(file)){
      file <-"proteiNorm_Report.pdf";no=""
      if(file.exists(file.path(dir,file))){
        file<-make_new_filename(x=file,dir=dir);file
        no<-sub(".pdf","",sub(".*_","",file));no
        pngdir=gsub(".pdf","",file)
      }
      if(!file.exists(file.path(dir,file))){ pngdir=gsub(".pdf","",file) }
    }  ## FILE NULL

    print(paste("QC output directory:", dir))
    print(paste("QC report file: ",file))
    print(paste("png directory: ",pngdir))

  } ## SAVE ==TRUE


  ## CREATE PLOTS
  nahm   <- plotNaHM(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,dir=dir,save=save)
  pcv    <- plotPCV(normList=normList,groups=groups,batch=batch,dir=dir,save=save)
  pmad   <- plotPMAD(normList=normList,groups=groups,batch=batch,dir=dir,save=save)
  pev    <- plotPEV(normList=normList, groups=groups,batch=batch,dir=dir,save=save)
  cor    <- plotCOR(normList=normList,groups=groups,batch=batch,dir=dir,save=save)
  lograt <- plotLogRatio(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,
                         zoom=FALSE,legend=TRUE,inset=0.02,dir=dir,save=save)
  plotLogRatio(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,
               zoom=TRUE,legend=TRUE,inset=0.02,dir=dir,save=save)
  nahm   <- plotNaHM(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,dir=dir,save=save)
  totint <- plotTotInten(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,dir=dir,save=save)
  if(save==FALSE){ draw(nahm$hm_batch); draw(nahm$hm_clust); draw(nahm$hm_group) }

  data2 <- list(pcv=pcv,pmad=pmad,pev=pev,cor=cor,lograt=lograt,nahm=nahm,totint=totint)
  if(save==TRUE){ data2<-list(pcv=pcv,pmad=pmad,pev=pev,cor=cor,lograt=lograt,
                              nahm=nahm,totint=totint,dir=dir,file=file) }


  ##  SAVE PROTEINORM REPORT PDF
  if(save==TRUE){

    pdf(file.path(dir,file), paper="USr", pagecentre=TRUE, pointsize=10,width=12,height=8)

    files <- c("PCVplot.png","PMADplot.png","PEVplot.png","CORplot.png","Log2RatioPlot.png",
               "Log2RatioPlot-zoom.png", "NaHMplot.png","NaHMplot_clust.png","NaHMplot_group.png",
               "NaHMplot_batch.png", "TotIntenPlot.png")
    pnglist<-paste0(paste0(file.path(dir),"/"),files);pnglist
    thePlots<-lapply(1:length(pnglist), function(i){grid::rasterGrob(png::readPNG(pnglist[i],native=F))})

    do.call(gridExtra::grid.arrange,c(thePlots[1:6], ncol=3))
    do.call(gridExtra::grid.arrange,c(thePlots[1], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[2], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[3], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[4], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[5], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[7], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[8], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[9], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[10], ncol=1))
    do.call(gridExtra::grid.arrange,c(thePlots[11], ncol=1))

    dev.off()

    ## REMOVE PNG FILES
    if(keep.png==FALSE){
      unlink(pnglist)
      print(file.path(dir,file))
      print("png files removed...")
    } ## remove png files.

    ## KEEP PNG FILES
    if(keep.png==TRUE){
      if(!dir.exists(file.path(dir,pngdir))){dir.create(file.path(dir,pngdir),recursive=TRUE)}
      lapply(files,function(x){
        file.copy(from=file.path(dir,x), to=file.path(dir,pngdir,x))
        file.remove(file.path(dir,x))
      })
      print(paste("png files moved to :",file.path(dir,pngdir)))
    } ## KEEP

  } ## SAVE == TRUE

  return(invisible(data2))


}
