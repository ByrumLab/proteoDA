make_design <- function(targets, group, factors=NULL, paired=NULL){

  param<-list()
  param[["group"]]<-group
  param[["factors"]]<-ifelse(is.null(factors),"NULL",paste(factors,collapse=", "))
  param[["paired"]]<-ifelse(is.null(paired),"NULL",paired)
  param<-t(data.frame(param))
  colnames(param)<-"make.design.parameters"
  param

  ## if input values are all column names in targets
  pass<-c(group,paired,factors) %in% colnames(targets);pass
  if(all(pass)==TRUE){

    ## subset targets to include input
    ## columns should be in this order: group,paired,followed by factors
    if(is.null(factors)&is.null(paired)){
      tar<-as.data.frame(targets[,c(group,paired,factors)],row.names=rownames(targets));tar
      colnames(tar)<-group
    } else{
      tar<-targets[,c(group,paired,factors)];head(tar)
    }
    ## check that each factor contains 2 or more levels
    levs <- apply(tar,2,FUN=function(x){length(unique(x))})>1;levs
    if(any(levs==FALSE)){
      stop(paste("Error! The following factors only have 1 level.",
                 "Factors must have 2 or more levels in order",
                 "to create design matrix. \nInvalid factors include: ",
                 paste(colnames(tar)[levs==FALSE],collapse=", ")))
    }

    ## make group a named list and create means model no intercept formula
    grp<-list(as.list.data.frame(tar[,group]));names(grp)<-group
    groupformula=paste0("~0+",names(grp));groupformula

    ## if paired factor supplied then design matrix and targets file for
    ## mixed effects model created. i.e. paired column in targets renamed
    ## "paired"
    if(!is.null(paired)){
      ## make paired a named list. if paired.type is paired then
      ## add to formula.
      p.col<-grep(paired, colnames(tar))
      stopifnot(colnames(tar)[p.col]==paired)
      colnames(tar)[p.col] <- "paired" ## rename column
      cat("\n")
      message("creating design matrix and targets file for paired sample design",
              "using limmas mixed effects model...");cat("\n")
      print(paste("input paired column name (",paired,") changed to 'paired' "))
      paired=NULL ## so paired is not included in design formula
      pairedformula<-names(paired)
    }
    ## no paired samples
    if(is.null(paired)){ pairedformula<-names(paired) }
    pairedformula

    ## named lists of factors. add to formula
    if(!is.null(factors)){
      ## single factor
      if(length(factors)==1){
        facs=list(as.list.data.frame(tar[,factors]))
        names(facs)=factors
      }
      ## multiple factors
      if(length(factors)>1){
        facs=as.list.data.frame(tar[,factors])
      }
    }
    ## no factors
    if(is.null(factors)){ facs<-NULL }
    additiveformula=names(facs);length(additiveformula)
    additiveformula


    ## make each a factor, if numeric or starts with number
    ## append column name to values to make it a character vector.
    for(i in base::seq_along(tar)){
      if(is.numeric(tar[,i]) | any(substr(tar[,i],1,1) %in% c(0:9))){
        tar[,i]<-make_factor(tar[,i], prefix=colnames(tar)[i])
      } else { tar[,i]<-make_factor(tar[,i]) }
    }


    ## DESIGN FORMULA
    if(length(groupformula)>0){designformula<-groupformula }
    if(length(pairedformula)>0){designformula<-paste(designformula,pairedformula,sep="+")}
    if(length(additiveformula)>0){
      additiveformula<-paste(additiveformula,collapse="+");additiveformula
      designformula<-paste(designformula,additiveformula,sep="+");designformula
    }
    designformula

    ## create design matrix
    design <- model.matrix(eval(parse(text=designformula,prompt="+")),data=tar); design

    ## change column names of design.
    desCols<-levels(tar[,group])
    for(x in c(paired,factors)){desCols<-c(desCols, levels(tar[,x])[-1])}
    print(colnames(design)); print(desCols)
    colnames(design)<-desCols;head(design)

  } else {
    ## if all input values are not column names in targets stop.
    if(any(pass==FALSE)){
      invalidCols <- c(group,paired,factors)[pass==FALSE];invalidCols
      stop("Error! Invalid input values: ", paste(invalidCols,collapse=", ") )
    }}


  print(head(tar));cat("\n"); print(head(design));cat("\n");
  print(tail(design));cat("\n"); print(designformula)
  cat("\n\n"); print("design matrix and targets created. Success!!")

  ## the targets file returned by this function should be used in the limma analysis
  data2 <- list(design=design, targets=tar, designformula=designformula)
  return(data2)


}
