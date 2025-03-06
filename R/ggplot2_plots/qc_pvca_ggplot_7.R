

# library(Biobase)
# library(pvca)



## y = filtered + normalized DGEList object (cpm filter/ TMM)
## factor_columns = column names of factors to test. factor columns to test
## factor columns must be columns listed in y$samples data.frame of DGEList object.
## if not listed then add the columns or if in separate metadata frame then add
## them to the y$samples metadata,
## threshold = numberic value of length 1L between 0 and 1
##
## example:
## qc_pvca_plot(y = proc$dgeLists$y.norm, factor_columns = c("treatment","batch"), threshold = 0.8)

qc_pvca_plot <- function(y, factor_columns, threshold = 0.8, title = NULL){

  if(class(y) != "DGEList"){stop("y is not a DGEList object.")}

## get normalized data
norm_lcpm <- as.matrix(edgeR::cpm(y=y, log=TRUE, prior.count=2, normalized.lib.sizes = T))

## make ExpressionSet object

## normalized data and metadata
assayData = as.matrix(norm_lcpm)
pData = y$samples[, factor_columns]
all.equal(rownames(pData), colnames(assayData))

## required nrow=ncol(pData)
varMetadata <- data.frame(labelDescription = c(factor_columns),
                          row.names = factor_columns)
varMetadata
#         labelDescription
# group     Disease/Status
# sex       Patient/Gender
# subtype  disease/Subtype

## make targets thing
phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=varMetadata)
Biobase::pData(phenoData)
#         group sex     subtype
# HDLEC-1 HDLEC   M     Control
# HDLEC-2 HDLEC   M     Control
# HDLEC-3 HDLEC   M     Control
# LMEC-4   LMEC   M Macrocystic
# LMEC-5   LMEC   M Macrocystic


## make expressionset object so you can run pvca
expSet <- Biobase::ExpressionSet(assayData=assayData, phenoData=phenoData)



## PVCA ANALYSIS

pvcaObj=pvca::pvcaBatchAssess(abatch = expSet,
                              batch.factors = factor_columns,
                              threshold = threshold)



## MAKE REGULAR BARPLOT

par(mar=c(8,6,2,2))
bp <- barplot(pvcaObj$dat, xlab = "", cex.lab=1,cex.axis = 1,
              ylab = "Weighted average proportion variance",
              ylim= c(0,1.1),col = c("dodgerblue3"), las=2,
              main=ifelse(is.null(title),"PVCA estimation bar chart", title))
axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 1, las=3)
values = pvcaObj$dat
new_values = round(values , 3)*100
text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 1)




## MAKE GGPLOT BARPLOT

input=data.frame(value=round(as.numeric(pvcaObj$dat),3),
                 component=pvcaObj$label,
                 dummy=c(rep("Factor",length(pvcaObj$dat))))
head(input)
#   value     component  dummy
# 1 0.006   sex:subtype Factor
# 2 0.042 group:subtype Factor
# 3 0.115     group:sex Factor
# 4 0.006       subtype Factor
# 5 0.045           sex Factor
# 6 0.565         group Factor

## GGPLOT BARPLOT PVCA ANALYSIS
# UofABioinformaticsHub/DataVisualisaton_BIS2016
#library(ggplot2)

title <- ifelse(is.null(title), "PCVA Estimation Bar Chart", title)
p <- ggplot2::ggplot(data = input, aes(x = component, y = value, title=" ", fill = dummy)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("dodgerblue3")) +
  labs(y = "Weighted average proportion variance", x="Effect", fill = "variable title") +
  theme_gray() +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 1,size=12),
        axis.text.y = element_text(size=12),
        plot.title = element_text(size=18)) +
  geom_bar(position = 'dodge', stat='identity') +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25)+
  theme(legend.position = "none") +
  ggtitle(title)

print(p)
p

}

# qc_pvca_plot(y= y, factor_columns = c("treatment","batch"), threshold = 0.9)











## FOR DAList from proteomics analyses
qc_pvca_plot2 <- function(DAList, factor_columns, threshold = 0.8, title = NULL){

  # suppressPackageStartupMessages(library(ggplot2))
  # suppressPackageStartupMessages(library(Biobase))


  DAList <- proteoDA:::validate_DAList(x = DAList)


   ## remove rows with missing values and sort by most variable proteins
   num_rows <- nrow(DAList$data)
   DAList$data <- DAList$data[!apply(is.na(DAList$data), 1, any), ]
   DAList$data <- DAList$data[order(proteoDA:::rowVars(as.matrix(DAList$data)), decreasing = TRUE), ]

  ## sort annotation to match filtered data
  DAList$annotation <- DAList$annotation[rownames(DAList$data), ]

  DAList <- proteoDA:::validate_DAList(x = DAList)
  cli::cli_inform("{.val {num_rows - nrow(DAList$data)}} rows with missing values were removed.")

  ## make ExpressionSet object

  ## normalized data and metadata
  assayData = as.matrix(DAList$data)
  pData = DAList$metadata[, factor_columns]
  all.equal(rownames(pData), colnames(assayData))

  ## required nrow=ncol(pData)
  varMetadata <- data.frame(labelDescription = c(factor_columns),
                            row.names = factor_columns)
  varMetadata
  #         labelDescription
  # group     Disease/Status
  # sex       Patient/Gender
  # subtype  disease/Subtype

  ## make targets thing
  # phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=varMetadata)
  phenoData <- Biobase::AnnotatedDataFrame(data = pData, varMetadata = varMetadata)
  # Biobase::pData(phenoData)
  #         group sex     subtype
  # HDLEC-1 HDLEC   M     Control
  # HDLEC-2 HDLEC   M     Control
  # HDLEC-3 HDLEC   M     Control
  # LMEC-4   LMEC   M Macrocystic
  # LMEC-5   LMEC   M Macrocystic


  ## make expressionset object so you can run pvca
  expSet <- Biobase::ExpressionSet(assayData=assayData, phenoData=phenoData)



  ## PVCA ANALYSIS

  pvcaObj=pvca::pvcaBatchAssess(abatch = expSet,
                                batch.factors = factor_columns,
                                threshold = threshold)



  ## MAKE REGULAR BARPLOT

  # par(mar=c(8,6,2,2))
  # bp <- barplot(pvcaObj$dat, xlab = "", cex.lab=1,cex.axis = 1,
  #               ylab = "Weighted average proportion variance",
  #               ylim= c(0,1.1),col = c("dodgerblue3"), las=2,
  #               main=ifelse(is.null(title),"PVCA estimation bar chart", title))
  # axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 1, las=3)
  # values = pvcaObj$dat
  # new_values = round(values , 3)*100
  # text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 1)




  ## MAKE GGPLOT BARPLOT

  input=data.frame(value=round(as.numeric(pvcaObj$dat),3),
                   component=pvcaObj$label,
                   dummy=c(rep("Factor",length(pvcaObj$dat))))
  head(input)
  #   value     component  dummy
  # 1 0.006   sex:subtype Factor
  # 2 0.042 group:subtype Factor
  # 3 0.115     group:sex Factor
  # 4 0.006       subtype Factor
  # 5 0.045           sex Factor
  # 6 0.565         group Factor

  ## GGPLOT BARPLOT PVCA ANALYSIS
  # UofABioinformaticsHub/DataVisualisaton_BIS2016
  # library(ggplot2)

  title <- ifelse(is.null(title), "PCVA Bar Chart", title)
  p <- ggplot2::ggplot(data = input, aes(x = component, y = value, title=" ", fill = dummy)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("dodgerblue3")) +
    labs(y = "Weighted average proportion variance", x="Effect", fill = "variable title",title = title) +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 1,size=12),
          axis.text.y = element_text(size=12),
          plot.title = element_text(size=18)) +
    geom_bar(position = 'dodge', stat='identity') +
    geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25)+
    theme(legend.position = "none") +
    ggtitle(title)

  print(p)
  list(p=p, DAList_NA =DAList)

}








## data = normalized data
## metadata = data.frame of the factors to investigate. only the factors you want to test can be included
qc_pvca_plot3 <- function(data, metadata, threshold = 0.8, title = NULL){

  if(!check_data(x = data)){
    cli::cli_abort("data is not a data.frame or matrix.")
  }


  if(!any(c("matrix","data.frame") %in% class(metadata))){
    stop("metadata is not a data frame of factors to test.")}

  if(!all(rownames(metadata)==colnames(data))){
    stop("column names of data and rownames of metadata do not match.")
  }
  factor_columns <- colnames(metadata)


  ## remove rows with missing values and sort by most variable proteins
  num_rows <- nrow(data)
  data <- data[!apply(is.na(data), 1, any), ]
  data <- data[order(proteoDA:::rowVars(as.matrix(data)), decreasing = TRUE), ]
  cli::cli_inform(c("{.val {num_rows-nrow(data)}} rows contain 1 or more missing
                    values and were removed. {.val {nrow(data)}} rows kept."))

  ## make ExpressionSet object

  ## normalized data and metadata
  assayData = as.matrix(data)
  pData = data.frame(metadata)
  all.equal(rownames(pData), colnames(assayData))

  ## required nrow=ncol(pData)
  varMetadata <- data.frame(labelDescription = c(factor_columns),
                            row.names = factor_columns)
  varMetadata
  #         labelDescription
  # group     Disease/Status
  # sex       Patient/Gender
  # subtype  disease/Subtype

  ## make targets thing
  # phenoData <- new("AnnotatedDataFrame", data=pData, varMetadata=varMetadata)
  phenoData <- Biobase::AnnotatedDataFrame(data = pData, varMetadata = varMetadata)
  # Biobase::pData(phenoData)
  #         group sex     subtype
  # HDLEC-1 HDLEC   M     Control
  # HDLEC-2 HDLEC   M     Control
  # HDLEC-3 HDLEC   M     Control
  # LMEC-4   LMEC   M Macrocystic
  # LMEC-5   LMEC   M Macrocystic


  ## make expressionset object so you can run pvca
  expSet <- Biobase::ExpressionSet(assayData=assayData, phenoData=phenoData)



  ## PVCA ANALYSIS

  pvcaObj=pvca::pvcaBatchAssess(abatch = expSet,
                                batch.factors = factor_columns,
                                threshold = threshold)



  ## MAKE REGULAR BARPLOT

  # par(mar=c(8,6,2,2))
  # bp <- barplot(pvcaObj$dat, xlab = "", cex.lab=1,cex.axis = 1,
  #               ylab = "Weighted average proportion variance",
  #               ylim= c(0,1.1),col = c("dodgerblue3"), las=2,
  #               main=ifelse(is.null(title),"PVCA estimation bar chart", title))
  # axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 1, las=3)
  # values = pvcaObj$dat
  # new_values = round(values , 3)*100
  # text(bp,pvcaObj$dat,labels = new_values, pos=3, cex = 1)




  ## MAKE GGPLOT BARPLOT

  input=data.frame(value=round(as.numeric(pvcaObj$dat), 3),
                   component=pvcaObj$label,
                   dummy=c(rep("Factor",length(pvcaObj$dat))))
  head(input)
  #   value     component  dummy
  # 1 0.006   sex:subtype Factor
  # 2 0.042 group:subtype Factor
  # 3 0.115     group:sex Factor
  # 4 0.006       subtype Factor
  # 5 0.045           sex Factor
  # 6 0.565         group Factor

  ## GGPLOT BARPLOT PVCA ANALYSIS
  # UofABioinformaticsHub/DataVisualisaton_BIS2016
  #library(ggplot2)

  title <- ifelse(is.null(title), "Principal Variant Component Analysis Estimation", title)
  p <- ggplot2::ggplot(data = input, aes(x = component, y = value, title=" ", fill = dummy)) +
    geom_bar(position = 'dodge', stat='identity') +
    geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
    scale_fill_manual(values = c("dodgerblue3")) +
    # labs(y = "Weighted average proportion variance", x="Effect", fill = "variable title") +
    labs(y = "Weighted average proportion variance", x="Effect", fill = "variable title") +
    theme_gray() +
    theme(axis.text.x = element_text(angle = 45,vjust = 0.9,hjust = 1,size=12),
          axis.text.y = element_text(size=12),
          plot.title = element_text(size=18),
          legend.position = "none") +
      ggtitle(label = title)


  plot(p)

  out <- list(p=p, data_na = data, metadata=metadata, threshold=threshold)
  out

}


