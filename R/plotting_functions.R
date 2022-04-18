
##------------------------------
##  PLOT RIGHT MARGIN
##------------------------------
right_margin <- function(x) {
  maxchar <- max(nchar(as.character(x))) + 10
  right <- maxchar/5
  return(right)
}


##------------------------------
##  PLOT LEFT MARGIN
##------------------------------
left_margin <- function(x) {
  maxchar <- max(nchar(as.character(x))) + 30
  left <- maxchar/5
  return(left)
}


##------------------------------
##  PLOT HEIGHT
##------------------------------
plot_height <- function(x){
  height <- (600/20)*length(x)
  height <- ifelse(height > 1000, 1000, height)
  return(height)
}


##------------------------------
##  [10] densityLog2Ratio
##------------------------------
densityLog2Ratio <- function(normList,
                             groups,
                             zoom = FALSE,
                             legend = TRUE,
                             inset = -0.2) {

  log2Ratio <- NULL
  groupList <- sort(unique(groups))

  for (method in names(normList)) {
    for (g1 in 1:(length(groupList) - 1)) {
      for (g2 in (g1 + 1):length(groupList)) {
        log2Ratio[[method]] <- c(log2Ratio[[method]],
                                 rowMeans(as.data.frame(normList[[method]][,groups == groupList[g1]]), na.rm = T) - rowMeans(as.data.frame(normList[[method]][, groups == groupList[g2]]), na.rm = T))
      }
    }
  }

  maxY=max(unlist(base::lapply(log2Ratio, FUN=function(x) max(density(x, na.rm=T)$y))));maxY
  minY=0

  if (!zoom) {
    minX <- 0.5 * min(unlist(min(density(log2Ratio[["vsn"]], na.rm = T)$x)))
    maxX <- 0.5 * max(unlist(max(density(log2Ratio[["vsn"]], na.rm = T)$x)))
  } else {
    minX <- -0.3
    maxX <- 0.3
    maxY <- maxY + (0.2*maxY)
    minY <- maxY - (0.5*maxY)
  }

  if (legend) {
    par(mar=c(5,5,4,3))
  } else {
    par(mar=c(5,5,3,4))
  }

  plot(NA, las = 1,
       xlim = c(minX, maxX),
       ylim = c(minY, maxY),
       xlab = "Log2 ratio",
       ylab = "Density",
       main = "Log2-ratio",
       cex.main = 1.5,
       cex.axis = 1.2,
       cex.lab = 1.3)
  abline(v = 0, lwd = 2, lty = 3, col = "grey")

  densityList <- list()
  for (method in names(normList)) {
    lines(density(log2Ratio[[method]], na.rm = T),
          col = binfcolors[which(names(log2Ratio) %in% method)],
          lwd = 3)

    densityList[[method]] <- density(log2Ratio[[method]], na.rm = T)
  }
  if (legend) {
    legend("topright",
           inset = c(inset, 0),
           names(log2Ratio),
           bty = "n",
           xpd = TRUE,
           box.col = "transparent",
           box.lwd = 0,
           bg = "transparent",
           border = "transparent",
           col = "transparent",
           pch = 22,
           pt.bg = binfcolors[1:length(log2Ratio)],
           pt.cex = 1.5,
           cex = 1,
           horiz = FALSE,
           ncol = 1)
  }

  return(invisible(densityList))

}



##------------------------------
##  [11] plotLogRatio
##------------------------------
plotLogRatio <- function(normList,
                         groups,
                         batch = NULL,
                         sampleLabels = NULL,
                         zoom = FALSE,
                         legend = TRUE,
                         inset = 0.02,
                         dir = ".",
                         save = FALSE) {

  groups <- make_factor(as.character(groups))

  if (is.null(batch)) {
    batch <- c(rep("1", ncol(normList[[1]])))
  }
  batch <- make_factor(as.character(batch), prefix = NULL)


  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(normList[[1]])
  }


  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive=TRUE)
    }
    filename <- file.path(dir, paste0("Log2RatioPlot", ifelse(zoom,"-zoom.png", ".png")))
    png(filename = filename,
        units = "px",
        width = 650,
        height = 650,
        pointsize = 15)
  }

  den <- densityLog2Ratio(normList, groups, zoom, legend, inset)
  if (save) dev.off()

  return(invisible(den))

}

### RESTYLED TO HERE

##------------------------------
##  [12] plotInten
##------------------------------
plotTotInten <- function(normList, groups, batch=NULL, sampleLabels=NULL,dir=".", save=FALSE){


  groups<-make_factor(groups);groups
  if(is.null(batch)){ batch <- c(rep("1",ncol(normList[[1]]))) };batch
  batch<-make_factor(as.character(batch),prefix=NULL)
  if(is.null(sampleLabels)){ sampleLabels <- colnames(normList[[1]]) };sampleLabels

  ## < 100 samples
  if(length(groups) < 100){
    width=round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02,0)
    height=800
    ncols=3
    par(oma=c(2, 1, 1, 1),mar=c(8,5,5,2))
  }

  ## >= 100 samples
  if(length(groups) >= 100){
    width=round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15,0);print(width)
    height=2400
    ncols=1
    par(oma=c(1,5,5,5),mar=c(8,2,2,2)) ## 100
  }


  if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
  if(save==TRUE){png(filename=file.path(dir,"TotIntenPlot.png"), units="px",
                     width=width, height=height, pointsize=15)}

  layout(matrix(1:9, ncol=ncols, byrow=TRUE))
  barList<-NULL
  for(i in names(normList)){
    barList[[i]] <- colSums(normList[[i]],na.rm = T)
    barplot(barList[[i]],
            # colSums(normList[[i]], na.rm=T),
            main="", las=2, yaxt="n",
            cex.main=1.5, cex.lab=1.2,
            col=colorGroup2(groups)[groups], names.arg=sampleLabels)
    title(main=i, font.main=1, cex.main=1.5, line=2)
    axis(side=2, cex.axis=1.2, las=2)
    if(i == "VSN") mtext(side=2, text="Total Intensity", line=6, cex=1.5)
    # abline(h=max(colSums(normList[[i]], na.rm=T)), lty=2)
  }
  names(barList)<-names(normList)
  if (save) dev.off()

  return(invisible(barList))

} ## PLOTTOTINTEN


##------------------------------
##  [13] heatmapClustered
##------------------------------
heatmapClustered <- function(missing, groups, batch=NULL, groupCol, batchCol, sampleLabels=NULL,legend=FALSE){

  if(is.null(batch)){ batch <- c(rep("1",ncol(missing))) }
  if(is.null(sampleLabels)){ sampleLabels <- colnames(missing) };sampleLabels

  ColAnn <- ComplexHeatmap::HeatmapAnnotation(Sample=groups, Batch=batch, col=list(Sample=groupCol, Batch=batchCol),
                              annotation_legend_param=list(Sample=list(title="Group", at=unique(groups),
                                                                       labels=paste("", unique(groups))),
                                                           Batch=list(title="Batch", at=unique(batch),
                                                                      labels=paste("Batch", unique(batch)))),
                              show_legend = ifelse(legend==FALSE,FALSE,TRUE)
  )
  hm_clust=ComplexHeatmap::Heatmap(missing+0, col=c("white", "black"), column_names_side="top",
                   column_title="Clustered", show_row_names=FALSE, show_column_names=TRUE,
                   name="Status", column_names_gp=grid::gpar(fontsize=7),
                   heatmap_legend_param=list(at=c(0, 1), labels=c("Missing", "Valid")),
                   show_heatmap_legend = ifelse(legend==FALSE, FALSE,TRUE),
                   top_annotation=ColAnn, column_labels=sampleLabels)
  # draw(hm_clust, heatmap_legend_side="right", annotation_legend_side="right")
  return(hm_clust)

} ## HEATMAPCLUSTERED


##------------------------------
##  [14] heatmapGroup
##------------------------------
heatmapGroup <- function(missing, groups, batch=NULL, groupCol, batchCol, sampleLabels=NULL, legend=FALSE){

  if(is.null(batch)){ batch <- c(rep("1",ncol(missing))) }
  if(is.null(sampleLabels)){ sampleLabels <- colnames(missing) };sampleLabels

  orderGroups=order(groups, batch)
  groups_groupSorted=groups[orderGroups]
  batch_groupSorted=batch[orderGroups]
  missing_groupSorted=missing[,orderGroups]

  ColAnn <- ComplexHeatmap::HeatmapAnnotation(Sample=groups_groupSorted, Batch=batch_groupSorted,
                              col=list(Sample=groupCol, Batch=batchCol),
                              annotation_legend_param=list(
                                Sample=list(title="Group", at=sort(unique(groups_groupSorted)),
                                            labels=paste("", sort(unique(groups_groupSorted)))),
                                Batch=list(title="Batch", at=sort(unique(batch_groupSorted)),
                                           labels=paste("Batch", sort(unique(batch_groupSorted))))
                              ),
                              show_legend=ifelse(legend==FALSE,FALSE,TRUE))
  hm_group=ComplexHeatmap::Heatmap(missing_groupSorted+0, col=c("white", "black"), column_names_side="top",
                   column_title="Sorted by Groups", show_row_names=FALSE,
                   show_column_names=TRUE, name="Status",
                   column_names_gp=grid::gpar(fontsize=7),
                   heatmap_legend_param=list(at=c(0, 1), labels=c("Missing", "Valids")),
                   show_heatmap_legend = ifelse(legend==FALSE, FALSE,TRUE),
                   top_annotation=ColAnn,
                   cluster_columns=FALSE, column_labels=sampleLabels[orderGroups])
  # draw(hm_group, heatmap_legend_side="right", annotation_legend_side="right")
  return(hm_group)

} ## HEATMAPGROUP


##------------------------------
##  [15] heatmapBatch
##------------------------------
heatmapBatch <- function(missing, groups, batch=NULL, groupCol, batchCol, sampleLabels=NULL, legend=FALSE){

  if(is.null(batch)){ batch <- c(rep("1",ncol(missing))) }
  if(is.null(sampleLabels)){ sampleLabels <- colnames(missing) };sampleLabels

  orderBatch=order(batch, groups)
  groups_batchSorted=groups[orderBatch]
  batch_batchSorted=batch[orderBatch]
  missing_batchSorted=missing[,orderBatch]
  ColAnn <- ComplexHeatmap::HeatmapAnnotation(Sample=groups_batchSorted, Batch=batch_batchSorted,
                              col=list(Sample=groupCol, Batch=batchCol),
                              annotation_legend_param=list(Sample=list(title="Group", at=sort(unique(groups_batchSorted)),
                                                                       labels=paste("", sort(unique(groups_batchSorted)))),
                                                           Batch=list(title="Batch", at=sort(unique(batch_batchSorted)),
                                                                      labels=paste("Batch", sort(unique(batch_batchSorted))))
                              ), show_legend=ifelse(legend==FALSE,FALSE,TRUE))
  hm_batch=ComplexHeatmap::Heatmap(missing_batchSorted+0, col=c("white", "black"), column_names_side="top",
                   column_title="Sorted By Batch", show_row_names=FALSE, show_column_names=TRUE,
                   name="Status", column_names_gp=grid::gpar(fontsize=7),
                   heatmap_legend_param=list(at=c(0, 1), labels=c("Missing", "Valid")),
                   show_heatmap_legend = ifelse(legend==FALSE,FALSE,TRUE),
                   top_annotation=ColAnn, cluster_columns=FALSE,
                   column_labels=sampleLabels[orderBatch])
  # draw(hm_batch, heatmap_legend_side="right", annotation_legend_side="right")
  return(hm_batch)

} ## HEATMAPBATCH


##------------------------------
##  [16] heatmapMissing
##------------------------------
heatmapMissing <- function(data, groups, batch=NULL, sampleLabels=NULL, showAllProtein=FALSE, legend=TRUE){

  if(is.null(batch)){ batch <- c(rep("1",ncol(data))) }
  if(is.null(sampleLabels)){ sampleLabels <- colnames(data) };sampleLabels

  missing=!is.na(data)
  if(!showAllProtein){
    complete=apply(missing, 1, all)
    completeNA=apply(!missing, 1, all)
    missing=missing[!complete & !completeNA,]
  }
  head(missing);dim(missing)

  batchCol = colorBatch(batch)
  groupCol = colorGroup2(groups)

  hm_clust = heatmapClustered(missing, groups, batch, groupCol, batchCol, sampleLabels, legend=FALSE)
  hm_group = heatmapGroup(missing, groups, batch, groupCol, batchCol, sampleLabels, legend=FALSE)
  hm_batch = heatmapBatch(missing, groups, batch, groupCol, batchCol, sampleLabels, legend=ifelse(legend==TRUE,TRUE,FALSE))

  ComplexHeatmap::draw(hm_clust+hm_group+hm_batch, heatmap_legend_side="right", annotation_legend_side="right",
       ht_gap=grid::unit(2, "cm"), column_title="Missing Values")

  return(invisible(missing))

} ## HEATMAPMISSING


##------------------------------
##  [17] plotNAHM
##------------------------------
plotNaHM <- function(normList, groups, batch=NULL, sampleLabels=NULL, dir=".", save=FALSE){

  if(is.null(batch)){ batch <- c(rep("1",ncol(normList[[1]]))) };batch
  if(is.null(sampleLabels)){ sampleLabels <- colnames(normList[[1]]) };sampleLabels
  if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)}

  batchCol = colorBatch(batch)
  groupCol = colorGroup2(groups)

  if(length(groups)<100){ width=round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02,0)
  width2=width/200;width=width/72 }
  if(length(groups)>=100){ width=round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15,0)
  width2=width/90; width=width/30 }


  if(save==TRUE){png(filename=file.path(dir, "NaHMplot.png"), units="in",
                     width=width, height=8, res=100,pointsize=8)}
  miss <- heatmapMissing(data=normList[["log2"]], groups=groups, batch=batch, sampleLabels=sampleLabels,legend=TRUE)
  if (save) dev.off()

  if(save==TRUE){png(filename=file.path(dir, "NaHMplot_clust.png"), units="in",
                     width=width2, height=8, res=100,pointsize=8)}
  hm_clust = heatmapClustered(missing=miss, groups, batch, groupCol, batchCol, sampleLabels, legend=TRUE)
  ComplexHeatmap::draw(hm_clust)
  if (save) dev.off()

  if(save==TRUE){png(filename=file.path(dir, "NaHMplot_group.png"),
                     units="in", width=width2, height=8, res=100,pointsize=8)}
  hm_group = heatmapGroup(missing=miss, groups, batch, groupCol, batchCol, sampleLabels, legend=T)
  ComplexHeatmap::draw(hm_group)
  if (save) dev.off()

  if(save==TRUE){png(filename=file.path(dir, "NaHMplot_batch.png"),
                     units="in", width=width2, height=8, res=100, pointsize=8)}
  hm_batch = heatmapBatch(missing=miss, groups, batch, groupCol, batchCol, sampleLabels, legend=T)
  ComplexHeatmap::draw(hm_batch)
  if (save) dev.off()

  data2 <- list(missing=miss, hm_clust=hm_clust,hm_group=hm_group,hm_batch=hm_batch)
  return(invisible(data2))


} ## PLOTNAHM




####################################
## QC REPORT FUNCTIONS BELOW????? ##
####################################

##------------------------------
##  [18] plotCorHM
##------------------------------
plotCorHM <- function(data, groups, batch=NULL, sampleLabels=NULL, dir=".",save=FALSE){

  groups <- make_factor(groups);groups
  if(is.null(batch)){ batch <- c(rep("1",ncol(data))) }
  batch <- make_factor(as.character(batch),prefix=NULL);batch
  if(is.null(sampleLabels)){ sampleLabels <- colnames(data) };sampleLabels


  if(length(groups) < 100){
    width=round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02,0)
    fontsize=12
  }

  if(length(groups) >= 100){
    width=round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15,0)
    width=round(width+width*0.3,0);print(width)
    fontsize=10
  }

  if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
  if(save==TRUE){png(filename=file.path(dir, "CorrHeatmap.png"), units="px",
                     width=width, height=width, pointsize=15)}
  cor_mat <- stats::cor(data, use="pairwise.complete.obs", method="pearson")

  ColAnn <- ComplexHeatmap::HeatmapAnnotation(Sample=groups,
                              col=list(Sample=colorGroup2(groups)),
                              annotation_legend_param=list(
                                Sample=list(title="Groups", at=levels(groups),
                                            labels=paste(levels(groups)))))

  RowAnn <- rowAnnotation(Batch=batch, col=list(Batch=colorBatch(batch)),
                          annotation_legend_param=list(
                            Batch=list(title="Batch", at=levels(batch),
                                       labels=paste("Batch", levels(batch)))))

  # col_fun = colorRamp2(c(-1,0, 1), colors=c("#47d604","white","#f54c57"), transparency = 0.7)
  hm_corr=ComplexHeatmap::Heatmap(cor_mat, name="Pearson correlation", border=TRUE,
                  col=circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7)),
                                           # colors=RColorBrewer::brewer.pal(8, "RdBu")[8:1],
                                           colors=RColorBrewer::brewer.pal(8, "Blues"),
                                           transparency=0.6),
                  # col=col_fun(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7))),
                  heatmap_legend_param=list(color_bar="continuous",
                                            legend_direction="horizontal",
                                            legend_width=grid::unit(5, "cm"),
                                            title_position="topcenter"),
                  column_names_gp=grid::gpar(fontsize=fontsize),
                  row_names_gp=grid::gpar(fontsize=fontsize),
                  top_annotation=ColAnn,
                  left_annotation=RowAnn,
                  column_labels=sampleLabels,
                  row_labels=sampleLabels,
                  ## function adds txt to cells of heatmap
                  cell_fun = function(j, i, x, y, width, height, fill){
                    grid.text(sprintf("%.1f",cor_mat[i, j]), x, y,gp=grid::gpar(fontsize=fontsize-2)) }
  )
  par(mar=c(10,10,10,10))
  ComplexHeatmap::draw(hm_corr, heatmap_legend_side="top",ht_gap=grid::unit(2, "cm"))
  if (save) dev.off()

  return(invisible(cor_mat))


} ## PLOT_CORHM


##------------------------------
##  [19] BOX PLOT
##------------------------------
plotBoxplot <- function(data, groups=NULL, sampleLabels=NULL, title=NULL, legend=TRUE, dir=".", save=FALSE){

  if(is.null(sampleLabels)){ sampleLabels <- as.character(colnames(data)) }
  if(is.null(groups)){ groups <- c(rep("group",ncol(data))) }
  groups <- make_factor(x=as.character(groups),prefix=NULL)

  ## remove rows with an NA
  data <- data[!apply(is.na(data), 1, any),];head(data);dim(data)

  ## plot margins
  x2<-left_margin(x=sampleLabels);x2
  x3<-right_margin(x=groups);x3

  width=round(0.0191*length(groups)^2 + 12.082*length(groups) + 671.75,0);print(width)

  if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
  if(save==TRUE){ png(filename=file.path(dir,"BoxPlot.png"), units="px",
                      width=width, height=750, pointsize=15)}

  op <- par(no.readonly = TRUE)
  par(mar = c(x2,6,3, x3), par(oma=c(0.5,0,0.5,x3/2+3)))
  box<-graphics::boxplot(data, col=colorGroup2(groups)[groups], names=sampleLabels,
                         notch=FALSE, horizontal=FALSE, outline=FALSE,
                         las=2, cex.axis=1, cex.labs=1, cex=1)
  mtext(side=2, text="Intensity", font=1, line=3, cex=1.2)
  title(main=ifelse(is.null(title), "", title), font.main=1, cex.main=1.3, line=1.1)
  if(legend==FALSE){ x3<-1};x3
  if(legend==TRUE){
    legend(par("usr")[2],par("usr")[4], bty="n",xpd=NA, legend=levels(groups), pch = 22, cex=1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg=colorGroup2(groups),
           col=colorGroup2(groups), pt.cex=1.2, pt.lwd=0, inset=0.02, horiz=F,ncol=1)
  }
  if (save) dev.off()

  return(invisible(box))


} ## BOXPLOT


##------------------------------
##  [20] VIOLIN PLOT
##------------------------------
plotViolin <- function(data, groups=NULL, sampleLabels=NULL, title=NULL, legend=TRUE, dir=".", save=FALSE){


  if(is.null(sampleLabels)){ sampleLabels <- as.character(colnames(data)) }
  if(is.null(groups)){ groups <- c(rep("group",ncol(data))) }
  groups <- make_factor(x=as.character(groups),prefix=NULL)

  ## remove rows with an NA
  data <- data[!apply(is.na(data), 1, any),];head(data);dim(data)


  ## convert data data.frame to a named list
  plotData <- list()
  for(k in 1:ncol(data)){plotData[[k]]<-data[,k]}
  names(plotData) <- colnames(data)

  ## plot margins
  x2<-left_margin(x=sampleLabels);x2
  x3<-right_margin(x=groups);x3


  width=round(0.0191*length(groups)^2 + 12.082*length(groups) + 671.75,0);print(width)

  if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
  if(save==TRUE){png(filename=file.path(dir,"ViolinPlot.png"), units="px",
                     width=width, height=750, pointsize=15)}


  op <- par(no.readonly = TRUE)
  par(mar = c(x2,6,3, x3), par(oma=c(0.5,0,0.5,x3/2+3)))
  vio<-vioplot::vioplot(plotData, col=colorGroup2(groups)[groups], names=sampleLabels,
                        main="", ylab="", xlab="",
                        las=2, font=1, cex.axis=1, cex.labs=1, cex=1)
  mtext(side=2, text="Density", line=3, cex=1.2)
  title(main=ifelse(is.null(title),"",title), font.main=1, cex.main=1.3, line=1.1)
  if(legend==FALSE){ x3<-1};x3
  if(legend==TRUE){
    legend(par("usr")[2],par("usr")[4], bty="n",xpd=NA, legend=levels(groups), pch = 22, cex=1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg=colorGroup2(groups),
           col=colorGroup2(groups), pt.cex=1.2, pt.lwd=0, inset=0.02, horiz=F,ncol=1)
  }
  if (save) dev.off()

  return(invisible(vio))


} ## VIOLIN



##------------------------------
##  [21] PCA PLOT
##------------------------------
## data=normalized intensities, or zscore intensities. if stdize=TRUE then it scales each row
## (each protein) to have a mean of zero and standard deviation = 1 if data matrix is already
## scaled then set stdize=FALSE, dims= vector of PC to plot.
plotPCA <- function(data, groups=NULL, sampleLabels=NULL, title=NULL, top=500, stdize=TRUE,
                    dims=c(1,2), cex.dot=2, xlim=NULL, ylim=NULL,legend=TRUE, dir=".", save=FALSE){


  if(is.null(sampleLabels)){ sampleLabels <- as.character(colnames(data)) }
  if(is.null(groups)){ groups <- c(rep("group",ncol(data))) }
  groups <- make_factor(x=as.character(groups),prefix=NULL)


  ## TOP VARIABLE PROTEINS
  ## remove rows with an NA and get top most variable proteins
  data <- data[!apply(is.na(data), 1, any),];head(data);dim(data)
  o=order(matrixStats::rowVars(as.matrix(data)), decreasing=TRUE)
  data=data[o,];head(data)
  top<-ifelse(nrow(data)>=top,top,nrow(data));print(top)
  data=data[1:top,];dim(data)
  ## center/scale rows (proteins) mean=0;stdev=1
  if(stdize==TRUE){ data <- t(scale(x=t(data),center=TRUE, scale=TRUE)) }

  ## PCA
  pca=stats::prcomp(t(data), scale=FALSE)
  pca$summary<- summary(pca)$importance
  ## % variance explained by each PC
  # eigs <- pca$sdev^2; eigs[1] / sum(eigs)
  # summary(pca)

  ## xlim
  max.char <- max(nchar(sampleLabels));max.char
  cex.char <- ifelse(max.char <= 10,1,0.8);cex.char

  if(is.null(xlim)){
    # max.char <- max(nchar(sampleLabels));max.char
    # cex.char <- ifelse(max.char <= 10,1,0.9);cex.char
    xmin <- abs(min(pca$x[,dims[1]])) + 0.10 * abs(min(pca$x[,dims[1]]));print(xmin)
    # xmax <- max(pca$x[,dims[1]]) + (0.05 * (chr * cex.char)) * max(pca$x[,dims[1]]);print(xmax)
    offset <- 0.6 * ((max.char/10) * cex.char);offset
    if(ncol(data) > 50){ offset <- 0.05;offset }
    xmax <- max(pca$x[,dims[1]]) + max(pca$x[,dims[1]])*offset;xmax
    xlim=c(-xmin,xmax);xlim
  }

  ## ylim
  if(is.null(ylim)){
    ymin <- abs(min(range(pca$x[,dims[2]]))) + 0.10 * abs(min(range(pca$x[,dims[2]])));print(ymin)
    ymax <- max(range(pca$x[,dims[2]]))+ 0.10*max(range(pca$x[,dims[2]]));print(ymax)
    ylim=c(-ymin,ymax);ylim
  }

  ## plot margins
  x2<-left_margin(x=sampleLabels);x2
  x3<-right_margin(x=groups);x3


  if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
  if(save==TRUE){png(filename=file.path(dir,"PCAplot.png"), units="px",
                     width=750, height=650, pointsize=15)}

  op <- par(no.readonly = TRUE)
  par(mar = c(x2,6,3, x3), par(oma=c(0.5,0,0.5,x3/2+3)))
  plot(x=pca$x[,dims[1]], y=pca$x[,dims[2]],
       pch=21, bg=colorGroup2(groups)[groups], col=colorGroup2(groups)[groups],
       cex=cex.dot, lwd = 1,  las=1, cex.axis=1.2, cex.lab=1.3,
       xlim=xlim, ylim=ylim,
       xlab=paste0("PC ",dims[1]," (", round(summary(pca)$importance["Proportion of Variance",dims[1]]*100, 2)," %)"),
       ylab=paste0("PC ",dims[2]," (", round(summary(pca)$importance["Proportion of Variance",dims[2]]*100, 2)," %)")
  )
  title(main=ifelse(is.null(title),"",title), font.main=1, cex.main=1.3, line=1.2)
  grid()
  if(ncol(data) <= 50){
    text(labels=sampleLabels, x=pca$x[,dims[1]], y=pca$x[,dims[2]],
         cex=cex.char, adj=ifelse(max.char < 3,-0.6, ifelse(max.char < 5, -0.4,-0.2)),
         col=colorGroup2(groups)[as.character(groups)])
  }
  ## mtext(paste("(glmpca = ", fam_list[j], " | ", plot_list[i],
  ## " | factor =  ", colnames(gpca$factors)[k],")", sep=""), font=3, cex=0.9, line=1)
  if(legend==FALSE){ x3<-1};x3
  if(legend==TRUE){
    legend(par("usr")[2],par("usr")[4], bty="n",xpd=NA, legend=levels(groups), pch = 22, cex=1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg=colorGroup2(groups),
           col=colorGroup2(groups), pt.cex=1.2, pt.lwd=0, inset=0.02, horiz=F,ncol=1)
  }
  if (save) dev.off()
  # par(op)

  data2 <- list(pca=pca, dat=data)
  return(invisible(data2))


} ## PCA PLOT


##------------------------------
##  [22] SAMPLE DENDROGRAM
##------------------------------
## normalized intensities input (cols=samples, rows=proteins.
## NAs removed, then top variable proteins identified.
## if stdize is true the rows of data are centered and scaled to 0 and 1 z-score,
## distance calc. then clustering.
plotDendrogram <- function(data, groups=NULL, sampleLabels=NULL, top=500, stdize=TRUE,
                           clust.metric="euclidean", clust.meth="complete",
                           cex.names=1, xlim=NULL, title=NULL, legend=TRUE,dir=".",save=FALSE){

  clust.metric <- match.arg(arg=clust.metric,
                            choices=c("pearson", "sqrt pearson", "spearman", "absolute pearson",
                                      "uncentered correlation", "weird", "cosine", "euclidean",
                                      "maximum", "manhattan", "canberra", "binary","minkowski"),
                            several.ok=FALSE)

  clust.meth <- match.arg(arg=clust.meth,
                          choices=c("ward.D","ward.D2","single","complete","average","mcquitty",
                                    "median","centroid"), several.ok=FALSE)

  if(is.null(sampleLabels)){ sampleLabels<-colnames(data) }
  if(is.null(groups)){ groups <- c(rep("group",ncol(data))) }
  groups <- make_factor(x=groups);groups


  ## remove rows with an NA,get top variable proteins
  data <- data[!apply(is.na(data), 1, any),];head(data);dim(data)
  o=order(matrixStats::rowVars(data), decreasing=TRUE)
  data=data[o,];head(data)
  top<-ifelse(nrow(data)>=top,top,nrow(data));print(top)
  data=data[1:top,]

  ## center and scale rows (protein) (mean=0, std=1)
  if(stdize==TRUE){ data <- t(scale(t(data))) }


  # if(is.null(xlim)){xlim<-c(0,100)}
  x2<-left_margin(x=sampleLabels);x2
  x3<-right_margin(x=groups);x3


  width=round(0.0191*length(groups)^2 + 12.082*length(groups) + 671.75,0);print(width)


  if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
  if(save==TRUE){png(filename=file.path(dir,"Dendrogram.png"), units="px",
                     width=width, height=650, pointsize=15)}


  op <- par(no.readonly = TRUE)
  par(xpd = T, mar = par()$mar + c(3,0,2,5))
  d <- ClassDiscovery::distanceMatrix(dataset=data, metric=clust.metric)
  hc <- stats::hclust(d, method=clust.meth)
  ClassDiscovery::plotColoredClusters(hc, labs=sampleLabels, cols=colorGroup2(groups)[groups],
                                      lwd=1.5, las=2, cex.axis=1.2, xlab="",ylab="", font=1, cex=cex.names,
                                      line=-0.6)
  title(main=ifelse(is.null(title),"",title), font.main=1, line=2.5,cex=1.5)
  ## mtext(paste0("(n = ", nrow(data), " | ",metric," | ",method,")"), font=3, cex=0.8, line=2.7)
  if(legend==TRUE){
    legend(par("usr")[2],par("usr")[4], bty="n",xpd=NA, legend=levels(groups), pch = 22, cex=1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg=colorGroup2(groups),
           col=colorGroup2(groups), pt.cex=1.2, pt.lwd=0, inset=0.02, horiz=F,ncol=1)
  }
  par(mar=c(5, 4, 4, 2) + 0.1)

  if (save) dev.off()

  data2<- list(hc=hc, d=d, dat=data, clust.meth=clust.meth, clust.metric=clust.metric, stdize=stdize, top=top)

  return(invisible(data2))


} ## SAMPLE DENDROGRAM



##------------------------------
##  [23] CORRELATION SCATTER
##------------------------------
plotCorrScatter <- function(data, method=c("pearson","spearman","kendall"),
                            alpha=0.05, pch=".", dir=".", save=FALSE){

  data <- data[!apply(is.na(data), 1, any),];head(data);dim(data)

  if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
  if(save==TRUE){png(filename=file.path(dir,"corrScatter.png"), units="px",
                     width=1000, height=1000, pointsize=10)}

  cor_plot <- psych::pairs.panels(x        = data,
                                  smooth   = TRUE,          ## If TRUE, draws loess smooths
                                  scale    = FALSE,         ## If TRUE, scales the correlation text font
                                  density  = TRUE,          ## If TRUE, adds density plots and histograms
                                  ellipses = TRUE,          ## If TRUE, draws ellipses
                                  method   = method,        ## correlation method (also "spearman" or "kendall")
                                  pch      = ".",           ## pch symbol
                                  lm       = FALSE,         ## If TRUE, plots linear fit rather than the LOESS (smoothed) fit
                                  cor      = TRUE,          ## If TRUE, reports correlations
                                  jiggle   = FALSE,         ## If TRUE, data points are jittered
                                  factor   = 2,             ## Jittering factor
                                  hist.col = "dodgerblue1", ## Histograms color
                                  stars    = TRUE,          ## If TRUE, adds significance level with stars
                                  ci       = TRUE,          ## If TRUE, adds confidence intervals
                                  alpha    = alpha)         ## alpha level for confidence regions

  return(invisible(cor_plot))


}



colorBatch <- function(batch){

  batchCol <- unlist(ifelse(length(unique(batch)) == 1, grDevices::rainbow(1, start=0.5),
                            list(grDevices::rainbow(length(unique(batch))))))
  names(batchCol) <- unique(batch)

  return(batchCol)

} ## COLOR BATCH


##-----------------------------
##  [07] GROUP PLOT COLORS 1
##-----------------------------
colorGroup <- function(group){

  if(length(unique(group)) < 9){
    groupCol <- yarrr::piratepal(palette="basel")[1:length(unique(group))]
    names(groupCol) <- unique(group)
  } else {
    if(length(unique(group)) >= 9){
      groupCol <- grDevices::rainbow(length(unique(group)))
      names(groupCol) <- unique(group)
    }}

  return(groupCol)

} ## COLOR GROUP 1


##-----------------------------
##  [07] GROUP PLOT COLORS 2
##-----------------------------
colorGroup2 <- function(group){

  blueberry <- "#1F5EDC";  cherry <- "#EE0010"
  apple     <- "#32CD32";  barbie <- "#FF1493"
  fanta     <- "#FF7F00";  grape  <- "#A342FC"
  ocean     <- "#00C8FF";  mtndew <- "#ADFF2F"
  gold      <- "#FFE100";  orchid <- "#E36EF6"
  aceblue   <- "#009ACE";  poop   <- "#996633"

  binfcolors <-c(blueberry,cherry,apple,barbie,fanta,grape,ocean,mtndew,gold,orchid,aceblue,poop)
  names(binfcolors)<-c("blueberry","cherry","apple","barbie","fanta",
                       "grape","ocean","mtndew","gold","orchid","aceblue","poop")
  ## plots color palette
  colors <- unikn::newpal(col=binfcolors, names=names(colors))
  colors <- unikn::usecol(pal=binfcolors, n="all", use_names=TRUE)
  # unikn::seecol(pal=binfcolors, title="color palette: (binfcolors)")

  ## group=6-12
  if(length(unique(group)) > 5){
    groupCol <- binfcolors[1:length(unique(group))]
    names(groupCol) <- unique(group)
    # unikn::seecol(binfcolors)
  }
  if(length(unique(group)) > 12){
    groupCol <- c(grDevices::rainbow(length(unique(group)))); groupCol
    names(groupCol) <- unique(group)
    # unikn::seecol(pal=grDevices::rainbow(20))
  }
  if(length(unique(group))<=5){
    groupCol<-binfcolors[1:length(unique(group))]
    names(groupCol) <- unique(group)
    # unikn::seecol(groupCol)
  }
  if(length(unique(group))==3){
    groupCol<-c(binfcolors[c("blueberry","barbie","apple")])
    names(groupCol) <- unique(group)
    # unikn::seecol(groupCol)
  }
  if(length(unique(group))==2){
    groupCol<-binfcolors[c("aceblue","cherry")]
    names(groupCol) <- unique(group)
    # unikn::seecol(pal=groupCol)
  }
  return(groupCol)
} ## COLOR GROUP 2


##--------------------
##  [08] MERGE PDFS
##--------------------
## pdflist=c("path/to/file1.pdf", "path/to/file2.pdf", "path/to/file3.pdf")
## newpdf="path/to/name_of_new_combined_pdf_file.pdf
## delete.files==FALSE removes individual pdf files
merge_pdf_files <- function(pdflist, newpdf, delete.files=FALSE){

  ## determine if pdf files exist. logical T/F
  pass <- file.exists(pdflist);pass

  if(newpdf %in% pdflist==TRUE){
    stop("Error! newpdf cannot be included in pdflist.\n",
         "\nincorrect: newpdf='./dir/file3.pdf', pdflist=c('./dir/file1.pdf','./dir/file3.pdf')",
         "\n  correct: newpdf='./dir/file4.pdf,  pdflist=c('./dir/file1.pdf','./dir/file3.pdf')")
  }

  ## if all pdfs exist combine files into a single pdf
  ## remove individual pdf files if delete.files is FALSE
  if(all(pass==TRUE)){
    ## combine pdfs
    pdftools::pdf_combine(input=pdflist, output = newpdf)
    pdftools::pdf_length(newpdf)
    print(paste(newpdf, " created. Success!!"))

    ## remove individual pdfs
    if(delete.files==TRUE){
      file.remove(pdflist, recursive=TRUE)
      print("individual pdf files removed...")
    }
  }

  ## stop if any pdf files do not exist
  if(any(pass==FALSE)){
    stop(paste("Error! Some files in pdflist do not exist. Invalid files include: ",
               paste(pdflist[pass==FALSE],collapse="\n")))
  }


} ## MERGE PDF FILES


##---------------------------
##  [09] MERGE PNGS TO PDF
##---------------------------
## pnglist=c("./path/to/file.png","./path/to/file2.png")
## newpdf="./path/to/name.of.new.pdf.filename.pdf"
## delete.files=TRUE; removes list of png files.
merge_png_pdf <- function(pnglist, newpdf, delete.files=FALSE) {

  grDevices::pdf(newpdf)
  n <- length(pnglist)
  for(i in 1:n){
    pngfile <- pnglist[i]
    pngraster <- png::readPNG(pngfile)
    grid::grid.raster(pngraster, width=grid::unit(0.9,"npc"), height= grid::unit(0.8,"npc"))
    if(i < n){ plot.new() }
  }
  dev.off()
  if(delete.files){ unlink(pnglist) }

}
