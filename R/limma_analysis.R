run_limma_analysis <- function(data, annot, targets, design, contrasts, min.pval=0.055, min.lfc=1,
                               adj.method="BH", paired=FALSE, pipe="DIA", enrich="protein",
                               dir=NULL, save=TRUE, ilab="PI_DATE"){

  adj.method <- match.arg(arg=adj.method, choices=c("none","BH","BY","holm"),several.ok=FALSE);adj.method
  pipe   <- match.arg(arg=pipe,choices=c("DIA","TMT","phosphoTMT","LF"),several.ok=FALSE);pipe
  enrich <- match.arg(arg=enrich,choices=c("protein","phospho"), several.ok=FALSE);enrich

  param <- stats <- list()
  param[["min.lfc"]]    <- min.lfc
  param[["min.pval"]]   <- min.pval
  param[["adj.method"]] <- adj.method
  param[["paired"]]     <- paired
  param[["robust"]]     <- TRUE
  param[["pipe"]]       <- pipe
  param[["enrich"]]     <- enrich
  param[["ilab"]]       <- ilab


  ## MATCH TARGETS, DATA, ANNOTATION
  if(all(rownames(data) %in% rownames(annot))){ annot <- annot[rownames(data), ]
  } else { stop("Error! Row names of norm. data and row names of annotation
                    do not match.") }
  if(all(colnames(data) %in% rownames(targets))){
    data   <- data[, rownames(targets)]
    groups <- targets$group
  } else { stop("Error! Column names of norm. data and row names of targets
                    do not match.") }
  stopifnot(identical(colnames(data), rownames(targets)))
  stopifnot(identical(rownames(data), rownames(annot)))
  stopifnot(length(groups)==ncol(data))

  ## CREATE OUTPUT DIRECTORY
  if(save==TRUE){
    if(!is.null(dir)){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }

    if(is.null(dir)){
      if(enrich=="protein"){
        dir <- "./protein_analysis/02_diff_expression"
        if(!dir.exists(dir)){ dir.create(file.path(dir),recursive=TRUE) }
      }
      if(enrich=="phospho"){
        dir <- "./phospho_analysis/02_diff_expression"
        if(!dir.exists(dir)){ dir.create(file.path(dir), recursive=TRUE) }
      }
    }

  } ## SAVE == TRUE

  param[["dir"]] <- ifelse(save==TRUE, file.path(dir),"NULL")


  ##---------------------------
  ##  LIMMA EBAYES (NORMAL)
  ##---------------------------
  ## paired==FALSE is used for group comparisons,
  ## group comparisons correcting for e.g. batch or gender effect
  ## group comparisons for paired samples etc.
  if(paired==FALSE){

    ## fit limma models and perform DE
    fit     <- limma::lmFit(object=data, design=design)
    con.fit <- limma::contrasts.fit(fit=fit, contrasts=contrasts)
    efit    <- limma::eBayes(fit=con.fit, robust=TRUE)
    print(paste("limma DE analysis (paired==FALSE) complete. Success!!"))

    model<-list(efit=efit, con.fit=con.fit,fit=fit, design=design, contrasts=contrasts)

  } ## LIMMA EBAYES

  ##---------------------------------
  ##  LIMMA EBAYES (MIXED EFFECTS)
  ##---------------------------------
  ## paired==TRUE is used for group comparisons,
  ## when comparing within (paired samples) and across (groups)
  ## subjects
  if(paired==TRUE){

    message("performing paired analysis using mixed effects model ...")
    message("checking targets for 'paired' column ...")
    if("paired" %in% colnames(targets)==FALSE){
      stop("Error! Targets does not contain a column named 'paired'.
                 Correlation between sample pairs cannot be estimated.
                 To use the mixed effects model add a column named 'paired' (and set as factor)
                 to the targets file indicating paired sample info. or use make_design() function.
                 The design matrix should be created using a design formula that does
                 not include 'paired': e.g. ~0+group OR ~0+group+batch NOT ~0+group+paired")
    } ## FALSE
    if("paired" %in% colnames(targets)==TRUE){
      if(!is.factor(targets$paired)){
        if(is.character(targets$paired)){ targets$paired<-make_factor(targets$paired) }
        if(is.numeric(targets$paired)){ targets$paired<-make_factor(targets$paired, prefix="") }
      }} ## TRUE

    ## LIMMA EBAYES (MIXED EFFECTS)
    message("estimating correlation among sample pairs ...")
    corfit <- limma::duplicateCorrelation(object=data, design=design, block=targets$paired)
    cat("\n");message(paste("corfit = ", corfit$consensus.correlation)); cat("\n")
    if(corfit$consensus.correlation < 0.1){
      warning("Warning! The consensus correlation is either very small or has a negative value,
                    which may indicate little if any paired influence.")
    }
    fit     <- limma::lmFit(object=data, design=design, block=targets$paired,
                            correlation=corfit$consensus.correlation)
    con.fit <- limma::contrasts.fit(fit=fit, contrasts=contrasts)
    efit    <- limma::eBayes(fit=con.fit, robust=TRUE)
    print(paste("limma DE analysis (paired==TRUE) complete. Success!!"))

    model<-list(efit=efit, con.fit=con.fit,fit=fit, corfit=corfit,
                design=design, contrasts=contrasts)

  } ## MIXED EFFECTS


  ## DE stat results are extracted in 3 formats. statList is list object, where each
  ## item in the list is a data.frame of the stat results for a particular contrast.
  ## This list object is used to create/save DE plots and individual stat result files.
  ## comboStats = combined stat results in wide format. This data.frame is used to
  ## create results file for Big Query upload. This data.frame is used to create results file
  ## returned to the investigator.
  res <- extract_limma_results(efit=efit, annot=annot, data=data, min.pval=min.pval,
                               min.lfc=min.lfc, adj.method=adj.method, dir=dir, save=save,
                               enrich=enrich,ilab=ilab)

  ## SAVE DE PLOTS
  if(save==TRUE){ ## SAVE DE PLOTS
    print("saving limma plots ...")
    base::lapply(names(res$statList), function(x){

      ## VOLCANO PLOTS
      png(filename=file.path(dir,paste0(x,"_volcano_plot.png")), units="px",
          width=700,height=600, pointsize=15)
      volcanoPlot(stats=res$statList[[x]],comparison=x, min.pval=min.pval,
                  min.lfc=min.lfc, xlim=NULL,ylim=NULL,sig.type="p.adj",
                  top=NULL,labels=NULL,inset=-0.2,legend=TRUE)
      dev.off()
      png(filename=file.path(dir,paste0(x,"_volcano_plot_pvalue.png")), units="px",
          width=700,height=600, pointsize=15)
      volcanoPlot(stats=res$statList[[x]],comparison=x, min.pval=min.pval,
                  min.lfc=min.lfc, xlim=NULL,ylim=NULL,sig.type="pval",
                  top=NULL,labels=NULL,inset=-0.2,legend=TRUE)
      dev.off()

      ## MD PLOTS
      png(filename=file.path(dir,paste0(x,"_MD_plot.png")), units="px",
          width=700, height=600,pointsize=15)
      mdPlot(stats=res$statList[[x]], comparison=x, min.pval=min.pval,
             min.lfc=min.lfc, xlim=NULL, ylim=NULL, sig.type="p.adj",
             top=NULL, labels=NULL, inset=-0.2, legend=TRUE)
      dev.off()
      png(filename=file.path(dir,paste0(x,"_MD_plot_pvalue.png")), units="px",
          width=700, height=600,pointsize=15)
      mdPlot(stats=res$statList[[x]], comparison=x, min.pval=min.pval,
             min.lfc=min.lfc, xlim=NULL, ylim=NULL, sig.type="pval",
             top=NULL, labels=NULL, inset=-0.2, legend=TRUE)
      dev.off()

      ## P-VALUE HISTOGRAMS
      png(filename=file.path(dir,paste0(x,"_pvalue_histogram.png")), units="px",
          width=1400, height=600, pointsize=15)
      pvalueHistogram(stats=res$statList[[x]], comparison=x)
      dev.off()

      ## GLIMMA VOLCANO PLOTS
      glimmaVolcanoPlot(stats=res$statList[[x]], comparison=x, data=res$data,
                        res$annot, groups=groups, min.pval=min.pval,
                        min.lfc=min.lfc,sig.type="p.adj", pipe=pipe,
                        enrich=enrich, dir=dir)
      glimmaVolcanoPlot(stats=res$statList[[x]], comparison=x, data=res$data,
                        res$annot, groups=groups, min.pval=min.pval,
                        min.lfc=min.lfc,sig.type="pval", pipe=pipe,
                        enrich=enrich, dir=dir)

      ## GLIMMA MD PLOTS
      glimmaMDPlot(stats=res$statList[[x]], comparison=x, data=res$data,
                   annot=res$annot, groups=groups, min.pval=min.pval,
                   min.lfc=min.lfc, sig.type="p.adj", pipe=pipe,
                   enrich=enrich, dir=dir)
      glimmaMDPlot(stats=res$statList[[x]], comparison=x, data=res$data,
                   annot=res$annot, groups=groups, min.pval=min.pval,
                   min.lfc=min.lfc, sig.type="pval", pipe=pipe,
                   enrich=enrich, dir=dir)
    })

    print("All limma plots saved. Success!!")

  } ## SAVE DE PLOTS


  if(save==TRUE){

    ## SUMMARY OF DIFF EXPRESSION
    sumfile<-paste0("./",dir,"/summary.txt");sumfile
    sink(file=sumfile)
    cat(paste0("\n##",paste(rep("-",40),collapse="")))
    cat("\n##  Summary of Differential Expression")
    cat(paste0("\n##",paste(rep("-",40),collapse="")))
    cat("\n\n")
    cat(paste0("Significance Criteria: (|logFC| >= ",min.lfc, " & p-value <= ",min.pval,")"))
    cat("\n\n"); print(res$sum.dtp)
    cat(paste0("Significance Criteria: (|logFC| >= ",min.lfc, " & adj. p-value <= ",min.pval,")"))
    cat("\n\n"); print(res$sum.dt);cat("\n\n")
    sink()

  }## SAVE


  ## save summary DE results to log file
  if(!dir.exists("logs")){ dir.create("logs",recursive=TRUE) }
  sink(file="./logs/processing.log",append=TRUE)
  title = "LIMMA DE SUMMARY (P-VALUE)"
  cat(paste0("\n##",paste(rep("-",40),collapse="")))
  cat(paste0("\n##  ",title,"\n"))
  cat(paste0("##",paste(rep("-",40),collapse="")))
  cat("\n\n");print(res$sum.dtp); cat("\n\n")
  sink()

  sink(file="./logs/processing.log",append=TRUE)
  title = "LIMMA DE SUMMARY (ADJ. P-VALUE)"
  cat(paste0("\n##",paste(rep("-",40),collapse="")))
  cat(paste0("\n##  ",title,"\n"))
  cat(paste0("##",paste(rep("-",40),collapse="")))
  cat("\n\n");print(res$sum.dt); cat("\n\n")
  sink()

  ## save DE parameters to log file
  logs<-make_log(param=param, stats=stats, title="LIMMA DE ANALYSIS",save=TRUE)


  data2 <- list(statList=res$statList, comboStats=res$comboStats, comboStats.BQ=res$comboStats.BQ,
                de=res$de,dep=res$dep, dt=res$dt, dtp=res$dtp,sum.dt=res$sum.dt, sum.dtp=res$sum.dtp,
                contrastNames=res$contrastNames,model=model, targets=res$targets, groups=groups,
                data=res$data, annot=res$annot, param=logs$param,stats=logs$stats)

  return(data2)

} ## LIMMA



##----------------------------------------
##  EXTRACT LIMMA RESULTS (REQUIRED)
##----------------------------------------
extract_limma_results <- function(efit, data, annot, min.pval=0.055, min.lfc=1, adj.method="BH",
                                  dir=".", save=FALSE, enrich, ilab){


  ## DECIDE TESTS
  print(paste("extracting limma stat results for each comparison (statList)..."))
  contrastNames <- colnames(efit$coefficients);contrastNames
  dt  <- limma::decideTests(efit, adjust.method=adj.method, p.value=min.pval, lfc=min.lfc)
  dtp <- limma::decideTests(efit, adjust.method="none", p.value=min.pval, lfc=min.lfc)
  sum.dt  <- summary(dt);sum.dt
  sum.dtp <- summary(dtp);sum.dtp



  ## STAT RESULTS (statList)
  statList<-list()
  limmaStatColums <- c("logFC","CI.L","CI.R","AveExpr","t","B","P.Value","adj.P.Val")
  statList<-base::lapply(contrastNames,function(x){
    stats <- limma::topTable(efit, coef=x, number=Inf, adjust.method=adj.method,
                             sort.by="none", p.value=1, lfc=0, confint=TRUE)
    df<-cbind(dtp[,x],dt[,x]);colnames(df)<-c("sig.PVal","sig.FDR")
    stats<-cbind(stats[,limmaStatColums], df[rownames(stats),])
  })
  names(statList)<-contrastNames
  print(paste("statList created. Success!!"))

  ## COMBO STATS
  comboStats <- NULL
  for(x in names(statList)){
    tmp <- statList[[x]];head(tmp)
    colnames(tmp)<-paste(colnames(tmp),x,sep="_");colnames(tmp)
    if(!is.null(comboStats)){ comboStats <- cbind(comboStats, tmp[rownames(comboStats), ]) }
    if(is.null(comboStats)){ comboStats <- tmp }
  }

  ## SAVE LIMMA STAT RESULTS
  ## save stat results for individual contrasts as csv files.
  ## save combined stat results as a csv file in dir. also
  ## save BQ combined stat results (NAs /blanks replaced with zeros)
  ## as csv in project directory for Big Query upload.
  if(save==TRUE){ ## SAVE STAT RESULTS

    ## INDIVIDUAL STATS
    base::lapply(names(statList),function(x){
      stats <- statList[[x]]
      stats2 <- cbind(annot[rownames(stats),],data[rownames(stats),],stats);colnames(stats2)
      filename<-paste0(x,"_results.csv")
      utils::write.csv(stats2, file=file.path(dir,filename), row.names=FALSE)
    })
  }

  ## COMBINED STATS
  comboStats2<-cbind(annot[rownames(comboStats),],data[rownames(comboStats),],comboStats)
  if(save==TRUE){
    filename<-paste("combined_results.csv",sep="_");filename
    utils::write.csv(comboStats2,file=file.path(dir, filename),row.names=FALSE)
  }

  ## COMBINED STATS FOR BIG QUERY
  comboStats2[,][is.na(comboStats2)] <- 0
  comboStats2[,][comboStats2==""]    <- 0
  if(any(substr(colnames(comboStats2),start=1,stop=1)%in%c(0:9))==TRUE){
    colnames(comboStats2)<-paste0("X",colnames(comboStats2))
  }
  if(save==TRUE){
    filename <- paste(ilab, enrich, "results_BQ.csv",sep="_");filename
    if(file.exists(file.path(".",filename))){
      print("BQ file already exists. creating a new BQ filename...")
      filename <- make_new_filename(x=filename,dir=".");filename
      print(filename)
    }
    utils::write.csv(comboStats2,file=file.path(".",filename),row.names=FALSE)
    print("limma stat results for BQ saved. Success!!")
  }


  de <-lapply(names(statList), function(x){
    get_de(stats=statList[[x]], min.pval=min.pval, min.lfc=min.lfc, type="p.adj",de.type="limma")
  });names(de)<-names(statList)

  dep <-lapply(names(statList), function(x){
    get_de(stats=statList[[x]], min.pval=min.pval, min.lfc=min.lfc, type="pval",de.type="limma")
  });names(dep)<-names(statList)


  data2 <- list(statList=statList, comboStats=comboStats, comboStats.BQ=comboStats2, de=de,dep=dep,
                dt=dt, dtp=dtp, sum.dt=sum.dt, sum.dtp=sum.dtp, data=data, annot=annot,
                contrastNames=contrastNames)
  print("limma stat results extracted. Success!!")
  return(data2)


} ## GET RESULTS




##------------------
##    GET_DE
##------------------
## uses stats matrix to extract DE gene stuff including sig, up, dn stat matrices, ids for the sig lists
get_de <-function(stats, min.pval=0.055, min.lfc=1, type=c("p.adj","pval"), de.type=c("edger","limma")){

  if(de.type=="edger"){
    if(type=="p.adj"){cols=c("logFC","FDR")}
    if(type=="pval"){cols=c("logFC","PValue")}
  }
  if(de.type=="limma"){
    if(type=="p.adj"){cols=c("logFC","adj.P.Val")}
    if(type=="pval"){cols=c("logFC","P.Value")}
  }

  sig<-subset(stats, abs(stats[,cols[1]])>=min.lfc & stats[,cols[2]]<=min.pval)
  up<-subset(sig,sig[,cols[1]]>=min.lfc)
  dn<-subset(sig,sig[,cols[1]]<=min.lfc)
  percent <- round((nrow(sig)/nrow(stats))*100,2);percent

  info<-t(data.frame(list(no_up=nrow(up),no_dn=nrow(dn),no_sig=nrow(sig), pcnt_sig=paste0(percent,"%"),
                          no_genes=nrow(stats))));
  colnames(info)<-c("");print(info)


  ## TMM assumption that most genes are not DE messages
  if(percent>10){
    message(paste0("Warning: > 10% of the data is DE (",percent,"%  > 10%)"))
    message(paste0("Warning: the assumption that most genes/proteins/phospho are not DE may be violated (",
                   percent,"% is > 10%)"))
  }

  ids=list(sig=rownames(sig), up=rownames(up), dn=rownames(dn),all=rownames(stats))
  return(list(sig=sig, up=up, dn=dn, stats=stats, ids=ids, info=info,
              param=list(min.pval=min.pval, min.lfc=min.lfc)))


}
