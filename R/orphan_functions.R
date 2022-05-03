
##------------------------------
##  [01] findDensityCutoff
##------------------------------

## NOT SURE THIS FUNCTION GETS USED ANYWHERE?

## Finds upper limit of the histogram so each box contains at least cutoffPercent of the data
findDensityCutoff <- function(longdat, cutoffPercent = 0.001) {
  densityCutoff <- 0
  newDensityCutoff <- max(longdat, na.rm = TRUE)
  totalNum <- length(longdat)
  while(newDensityCutoff != densityCutoff) {
    densityCutoff <- newDensityCutoff
    breakSeq <- seq(from = max(min(longdat, na.rm = TRUE), 0),
                    to = densityCutoff,
                    by = (densityCutoff - max(min(longdat, na.rm=TRUE), 0)) / 30)
    freqs <- graphics::hist(longdat[longdat < densityCutoff], breaks = breakSeq, plot = FALSE)
    newDensityCutoff <- freqs$breaks[which((freqs$counts / totalNum) < cutoffPercent)[1] + 1]
  }
  return(densityCutoff)

}


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


## function adds txt to cells of heatmap
cell_fun = function(j, i, x, y, width, height, fill) {
  grid::grid.text(sprintf("%.1f", cor_mat[i, j]), x, y, gp = grid::gpar(fontsize = fontsize - 2))
}


# Get colors for groups
#
# Used to get colors for the groups in our missing value heatmaps.
#
# @param group A vector of group names.
#
# @return A vector of colors for each unique group
# @export
#
# @examples
# # No examples yet
#
# colorGroup <- function(group){
#
#   if(length(unique(group)) < 9){
#     groupCol <- yarrr::piratepal(palette="basel")[1:length(unique(group))]
#     names(groupCol) <- unique(group)
#   } else {
#     if(length(unique(group)) >= 9){
#       groupCol <- grDevices::rainbow(length(unique(group)))
#       names(groupCol) <- unique(group)
#     }}
#
#   return(groupCol)
# }
