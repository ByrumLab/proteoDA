
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
    freqs <- hist(longdat[longdat < densityCutoff], breaks = breakSeq, plot = FALSE)
    newDensityCutoff <- freqs$breaks[which((freqs$counts / totalNum) < cutoffPercent)[1] + 1]
  }
  return(densityCutoff)

}

