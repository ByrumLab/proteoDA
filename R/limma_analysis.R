## ---------------------------
##  LIMMA EBAYES (NORMAL)
## ---------------------------
## paired==FALSE is used for group comparisons,
## group comparisons correcting for e.g. batch or gender effect
## group comparisons for paired samples etc.



## ---------------------------------
##  LIMMA EBAYES (MIXED EFFECTS)
## ---------------------------------
## paired==TRUE is used for group comparisons,
## when comparing within (paired samples) and across (groups)
## subjects

fit_limma_model <- function(data,
                            targets,
                            design,
                            contrasts,
                            paired = FALSE) {
  # Original fxn had a check that rownames in the data equaled rownames in the
  # annotation. Took that out for now, and removed annotation as an argument,
  # since we don't actually need the annotation to fit the model. But may want
  # to put that back in. We definitely need to check that at some point (maybe in
  # later functions?), but it may be worth doing it now to just not waste time running the
  # model on something that is incorrect.

  # On the other hand, depending on the various outputs, may be able to
  # reorder and match things back up if they ever are actually out of sync.

  # NEED TO ADD SOME SORT OF CHECK FOR THE CONTRASTS


  # Do some checks of the input data
  # All samples in data have targets
  if (!all(colnames(data) %in% rownames(targets))) {
    problemCols <- colnames(data)[colnames(data) %notin% rownames(targets)]
    cli::cli_abort(c("Not all column names in {.arg data} have a matching rowname in {.arg targets}",
                     "!" = "{cli::qty(length(problemCols))} Column{?s} without a match: {.val {problemCols}}"))
  }
  # All targets have corresponding data
  if (!all(rownames(targets) %in% colnames(data))) {
    problemTargets <- rownames(targets)[rownames(targets) %notin% colnames(data)]
    cli::cli_abort(c("Not all target rownames in {.arg targets} have a matching colname in {.arg data}",
                     "!" = "{cli::qty(length(problemTargets))} Row{?s} without a match: {.val {problemTargets}}"))
  }

  # Ensure data cols are in same order as rows in targets
  data <- data[, rownames(targets)]
  groups <- targets$group

  # Double-check that the order is all the same
  if (!identical(colnames(data), rownames(targets))) {
    cli::cli_abort(c("Colnames of {.arg data} and rownames of {.arg targets} still not the same after sorting",
                     "i" = "Not sure how this error can happen. Use the dubugger."))
  }

  cli::cli_rule()

  # On to model fitting
  # First step of fitting differs between paired and unpaired
  if (paired) {

    cli::cli_inform("Performing paired analysis with mixed effects model")
    # check for paired col
    if ("paired" %notin% colnames(targets)) {
      cli::cli_abort(c("{.arg targets} does not contain required column names {.val paired}",
                       "i" = "Add {.val paired} column, or using {.fun make_design}",
                       "i" = "Ensure that design formula does not include {.val paired}",
                       "i" = "e.g. ~0+group OR ~0+group+batch NOT ~0+group+paired"))
    }

    # Coerce paired col to factor, if it isn't
    if (!is.factor(targets$paired)) {
      targets$paired <- make_factor(targets$paired, prefix = "")
    }

    corfit <- limma::duplicateCorrelation(object = data, design = design, block = targets$paired)
    corfit_display <- round(corfit$consensus.correlation, 3)
    cli::cli_inform("Estimated inter-duplicate correlation = {.val {corfit_display}}")

    if (corfit$consensus.correlation < 0.1) {
      cli::cli_inform(cli::col_yellow("Estimated inter-duplicate correlation is low,
                                      which may indicate little or no paired influence"))
    }

    fit <- limma::lmFit(object = data, design = design,
                        block = targets$paired,
                        correlation = corfit$consensus.correlation)
  } else {
    cli::cli_inform("Performing standard un-paired model")
    fit <- limma::lmFit(object = data, design = design)
    }

  # contrasts fit and eBayes are the same across paired and not paired
  con.fit <- limma::contrasts.fit(fit = fit, contrasts = contrasts)
  efit <- limma::eBayes(fit = con.fit, robust = TRUE)


  # Set up return
  model <- list(efit = efit, con.fit = con.fit, fit = fit, design = design, contrasts = contrasts)
  if (paired) model[["corfit"]] <- corfit

  cli::cli_inform("limma DE analysis with {.arg paired} == {paired} complete")
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success!"))

  model
}




dummy_function <- function(data,
                           annot,
                           targets,
                           design,
                           contrasts,
                           min.pval = 0.055,
                           min.lfc = 1,
                           adj.method = "BH",
                           paired = FALSE,
                           pipe = "DIA",
                           enrich = "protein",
                           dir = NULL,
                           save = TRUE,
                           ilab = "PI_DATE") {






  enrich <- match.arg(arg = enrich, choices = c("protein", "phospho"), several.ok = FALSE)




  pipe <- match.arg(arg = pipe, choices = c("DIA", "TMT", "phosphoTMT", "LF"), several.ok = FALSE)

  adj.method <- match.arg(
    arg = adj.method,
    choices = c("none", "BH", "BY", "holm"), several.ok = FALSE
  )

  ## MATCH TARGETS, DATA, ANNOTATION
  # from original limma analysis
  if (all(rownames(data) %in% rownames(annot))) {
    annot <- annot[rownames(data), ]
  } else {
    stop("Error! Row names of norm. data and row names of annotation
                    do not match.")
  }

  stopifnot(identical(rownames(data), rownames(annot)))



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



  param <- stats <- list()
  param[["min.lfc"]]    <- min.lfc
  param[["min.pval"]]   <- min.pval
  param[["adj.method"]] <- adj.method
  param[["paired"]]     <- paired
  param[["robust"]]     <- TRUE
  param[["pipe"]]       <- pipe
  param[["enrich"]]     <- enrich
  param[["ilab"]]       <- ilab



  ## save DE parameters to log file
  logs<-make_log(param=param, stats=stats, title="LIMMA DE ANALYSIS",save=TRUE)


  data2 <- list(statList=res$statList, comboStats=res$comboStats, comboStats.BQ=res$comboStats.BQ,
                de=res$de,dep=res$dep, dt=res$dt, dtp=res$dtp,sum.dt=res$sum.dt, sum.dtp=res$sum.dtp,
                contrastNames=res$contrastNames,model=model, targets=res$targets, groups=groups,
                data=res$data, annot=res$annot, param=logs$param,stats=logs$stats)

  return(data2)

}



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
