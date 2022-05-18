

# ALLOW GROUPS TO VARY HERE.

make_limma_de_report <- function(data, annot, targets, design, contrasts, min.pval = 0.055, min.lfc = 1,
                                 adj.method = "BH", paired = FALSE, pipe = "DIA", enrich = "protein",
                                 dir = NULL, save = TRUE, ilab = "PI_DATE") {


  # I think data here is just the raw data that got passed into the original function: an element from the normList
  param[["dir"]] <- ifelse(save == TRUE, file.path(dir), "NULL")

  ## CREATE OUTPUT DIRECTORY
  if (save == TRUE) {
    if (!is.null(dir)) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
    }

    if (is.null(dir)) {
      if (enrich == "protein") {
        dir <- "./protein_analysis/02_diff_expression"
        if (!dir.exists(dir)) {
          dir.create(file.path(dir), recursive = TRUE)
        }
      }
      if (enrich == "phospho") {
        dir <- "./phospho_analysis/02_diff_expression"
        if (!dir.exists(dir)) {
          dir.create(file.path(dir), recursive = TRUE)
        }
      }
    }
  } ## SAVE == TRUE


  ## SAVE DE PLOTS
  if (save == TRUE) { ## SAVE DE PLOTS
    print("saving limma plots ...")
    base::lapply(names(res$statList), function(x) {

      ## VOLCANO PLOTS
      grDevices::png(
        filename = file.path(dir, paste0(x, "_volcano_plot.png")), units = "px",
        width = 700, height = 600, pointsize = 15
      )
      volcanoPlot(
        stats = res$statList[[x]], comparison = x, min.pval = min.pval,
        min.lfc = min.lfc, xlim = NULL, ylim = NULL, sig.type = "p.adj",
        top = NULL, labels = NULL, inset = -0.2, legend = TRUE
      )
      grDevices::dev.off()
      grDevices::png(
        filename = file.path(dir, paste0(x, "_volcano_plot_pvalue.png")), units = "px",
        width = 700, height = 600, pointsize = 15
      )
      volcanoPlot(
        stats = res$statList[[x]], comparison = x, min.pval = min.pval,
        min.lfc = min.lfc, xlim = NULL, ylim = NULL, sig.type = "pval",
        top = NULL, labels = NULL, inset = -0.2, legend = TRUE
      )
      grDevices::dev.off()

      ## MD PLOTS
      grDevices::png(
        filename = file.path(dir, paste0(x, "_MD_plot.png")), units = "px",
        width = 700, height = 600, pointsize = 15
      )
      mdPlot(
        stats = res$statList[[x]], comparison = x, min.pval = min.pval,
        min.lfc = min.lfc, xlim = NULL, ylim = NULL, sig.type = "p.adj",
        top = NULL, labels = NULL, inset = -0.2, legend = TRUE
      )
      grDevices::dev.off()
      grDevices::png(
        filename = file.path(dir, paste0(x, "_MD_plot_pvalue.png")), units = "px",
        width = 700, height = 600, pointsize = 15
      )
      mdPlot(
        stats = res$statList[[x]], comparison = x, min.pval = min.pval,
        min.lfc = min.lfc, xlim = NULL, ylim = NULL, sig.type = "pval",
        top = NULL, labels = NULL, inset = -0.2, legend = TRUE
      )
      grDevices::dev.off()

      ## P-VALUE HISTOGRAMS
      grDevices::png(
        filename = file.path(dir, paste0(x, "_pvalue_histogram.png")), units = "px",
        width = 1400, height = 600, pointsize = 15
      )
      pvalueHistogram(stats = res$statList[[x]], comparison = x)
      grDevices::dev.off()

      ## GLIMMA VOLCANO PLOTS
      glimmaVolcanoPlot(
        stats = res$statList[[x]], comparison = x, data = res$data,
        res$annot, groups = groups, min.pval = min.pval,
        min.lfc = min.lfc, sig.type = "p.adj", pipe = pipe,
        enrich = enrich, dir = dir
      )
      glimmaVolcanoPlot(
        stats = res$statList[[x]], comparison = x, data = res$data,
        res$annot, groups = groups, min.pval = min.pval,
        min.lfc = min.lfc, sig.type = "pval", pipe = pipe,
        enrich = enrich, dir = dir
      )

      ## GLIMMA MD PLOTS
      glimmaMDPlot(
        stats = res$statList[[x]], comparison = x, data = res$data,
        annot = res$annot, groups = groups, min.pval = min.pval,
        min.lfc = min.lfc, sig.type = "p.adj", pipe = pipe,
        enrich = enrich, dir = dir
      )
      glimmaMDPlot(
        stats = res$statList[[x]], comparison = x, data = res$data,
        annot = res$annot, groups = groups, min.pval = min.pval,
        min.lfc = min.lfc, sig.type = "pval", pipe = pipe,
        enrich = enrich, dir = dir
      )
    })

    print("All limma plots saved. Success!!")
  }


  invisible(NULL)
}
















####################################
##    LIMMA DE PLOT FUNCTIONS     ##
####################################
{ ## DE PLOT FUNCTIONS


  ## ------------------------
  ##  [01]  VOLCANO PLOT
  ## ------------------------
  volcanoPlot <- function(stats, comparison = NULL, min.pval = 0.055, min.lfc = 1,
                          xlim = NULL, ylim = NULL, sig.type = c("pval", "p.adj"),
                          top = NULL, labels = NULL, inset = -0.3, legend = FALSE) {
    stats <- data.frame(stats)


    ## row names for DIA,TMT,phosphoTMT protein data are a combination
    ## of UniprotID,GN, and protein group id. phospho data = uniprotid, GN,
    ## and phospho site (e.g. S37). This info is parsed into a matrix.
    ncols <- length(unlist(strsplit(rownames(stats)[1], split = "_")))

    # mat  <- matrix(unlist(strsplit(rownames(stats),"_")), ncol=ncols, byrow=TRUE);head(mat)
    # rownames(mat)=rownames(stats)

    ## rownames(stats)[1489] ## "D3ZHA7_RGD1560334_predicted_1943"
    mat <- matrix(nrow = nrow(stats), ncol = 3)
    mat[, 1] <- gsub("_.*", "", rownames(stats)) ## before first _
    mat[, 2] <- sub("_[^_]+$", "", sub("^[^_]*_", "", rownames(stats))) ## after first _ & before last _
    mat[, 3] <- sub(".*_", "", rownames(stats)) ## after last _
    rownames(mat) <- rownames(stats)


    ## annotation info extracted from rownames of protein data
    ## Gene_name will be used to label significant points if top >0
    if (length(grep("S", mat[, 3])) == 0) {
      colnames(mat) <- c("UniprotID", "Gene_name", "id")
      mat <- data.frame(mat)
      mat$labels <- mat$Gene_name
      stats <- cbind(stats, mat[rownames(stats), ])
    }

    ## annotation info extracted from rownames of phospho data
    ## Gene_name(phosphosite) will be used to label sig. points if top >0
    if (length(grep("S", mat[, 3])) > 0) {
      colnames(mat) <- c("UniprotID", "Gene_name", "phosAAPosition")
      mat <- data.frame(mat)
      mat$labels <- paste0(mat$Gene_name, "(", mat$phosAAPosition, ")")
      stats <- cbind(stats, mat[rownames(stats), ])
    }

    ## column named 'p' added to stats data.frame based on sig.type selected. texted for y-axis label also defined
    if (sig.type == "p.adj") {
      stats$p <- stats$adj.P.Val
      ylab <- expression(paste("-", log[10], " (adj. p-value)", sep = ""))
      filename <- paste(comparison, "volcano_plot.png", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }

    if (sig.type == "pval") {
      stats$p <- stats$P.Value
      ylab <- expression(paste("-", log[10], " (p-value)", sep = ""))
      filename <- paste(comparison, "volcano_plot_pvalue.png", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }


    ## determine x/y limits for the plot.
    if (length(xlim) != 2) {
      xlim <- c(-1, 1) * (max(abs(stats$logFC), na.rm = TRUE) * 1.1)
    }

    if (length(ylim) != 2) {
      ylim <- c(0, max(-log10(stats$p) * 1.1, na.rm = TRUE))
    }



    ## global par plot parameters
    graphics::par(
      font.axis = 1, lwd = 1, font.main = 2, cex.main = 1.5, fg = "black",
      col.axis = "black", cex.axis = 1.2, cex.lab = 1.3, font.axis = 1, xaxs = "i", yaxs = "i"
    )
    if (legend == TRUE) {
      graphics::par(mar = c(6, 7, 4, 6))
    }
    if (legend == FALSE) {
      graphics::par(mar = c(6, 7, 4, 3))
    }

    ## create blank plot area
    plot(x = 1, las = 1, type = "n", main = "", xlab = "", ylab = "", xlim = xlim, ylim = ylim)
    graphics::title(line = 1, main = comparison, cex.main = 1.3)
    graphics::mtext(side = 1, line = 3, cex = 1.3, text = expression(paste(log[2], " (fold-change)", sep = "")))
    graphics::mtext(side = 2, line = 4.5, cex = 1.3, text = ylab)
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = "gray94", lwd = 0)
    graphics::grid(col = "white", lty = "dotted", lwd = 1.2)
    graphics::box(which = "plot", lwd = 4, lty = 1, col = "gray94")

    ## subset based on significance cutoffs
    up <- subset(stats, p <= min.pval & logFC >= min.lfc)
    down <- subset(stats, p <= min.pval & logFC <= -min.lfc)

    ## add points
    with(stats, graphics::points(logFC, -log10(p),
      pch = 21, bg = scales::alpha(colour = "grey60", alpha = 0.5), ## black
      col = scales::alpha(colour = "grey50", alpha = 0.5), lwd = 0.5, cex = 1.1
    ))
    with(up, graphics::points(logFC, -log10(p),
      pch = 21, bg = scales::alpha(colour = "#FF0000", alpha = 1), ## red
      col = scales::alpha(colour = "#990000", alpha = 1), lwd = 0.5, cex = 1.1
    ))
    with(down, graphics::points(logFC, -log10(p),
      pch = 21, bg = "#1E90FF", col = "#003399", lwd = 0.5, ## blue
      cex = 1
    ))

    graphics::abline(v = c(-min.lfc, min.lfc), col = "dimgray", lty = 2, lwd = 0.8)
    graphics::abline(h = -log10(min.pval), col = "dimgray", lty = 2, lwd = 0.8)


    ## ADD LABELS

    ## user  can supply a vector of labels to add to the points on the volcano plot.
    ## for protein data the list can be Gene_name, UniprotID, UniprotID_Gene_name, proKey,
    ## for phospho data phosKey, Gene_name or Gene_name(site).
    ## labels added to corresponding points in orange.
    if (!is.null(labels)) {
      if (length(labels) > 0) {

        ## get row numbers in stats corresponding to the labels list
        ## remove any labels that were not identified and subset
        ## stats. divide subdata into logFC>0 and logFC<0 so labels can be
        ## positioned on the R and L side of the point
        n <- keep <- c()
        for (i in 1:length(labels)) {
          tmp <- grep(paste0(labels[i], "_"), rownames(stats))

          tmp2 <- ifelse(length(tmp) == 0, FALSE, TRUE)
          keep <- c(keep, tmp2)
          table(keep)
          n <- c(n, tmp)
        }

        ## ADD USER INPUT LABELS TO POINTS
        # labs <- labels#[keep==T];labs
        subdata <- stats[n, ]

        ## ADD USER INPUT LABELS TO POS logFC > 0
        ## subset data for positive logFC values, add labels
        subup <- subset(subdata, subdata$logFC > 0)
        if (nrow(subup) > 0) {
          graphics::points(
            x = subup$logFC, y = -log10(subup$p), pch = 21, cex = 1.2,
            bg = "orange", col = "darkorange3", lwd = 0.5,
            graphics::text(
              x = subup$logFC, y = -log10(subup$p), col = "orange",
              labels = subup$labels, pos = 4, cex = 0.6, font = 2
            )
          )
        }

        ## ADD USER INPUT LABELS TO NEG logFC < 0
        subdn <- subset(subdata, subdata$logFC < 0)
        if (nrow(subdn > 0)) {
          graphics::points(
            x = subdn$logFC, y = -log10(subdn$p), pch = 21, cex = 1.2,
            bg = "orange", col = "darkorange3", lwd = 0.5,
            graphics::text(
              x = subdn$logFC, y = -log10(subdn$p), col = "orange",
              labels = subdn$labels, pos = 2, cex = 0.6, font = 2
            )
          )
        }
      } ## LEN >0
    } ## LABELS


    ## TOP

    ## the user can use 'top' to label the most sig. up and down regulated
    ## features (proteins=gene name, phospho=gene name(site).
    ## e.g. top=5, will label the top 5 up and down reg.
    if (!is.null(top)) { ## TOP
      if (length(top) == 1 & top > 0) {
        ## sig. up-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC Z-> A
        sig <- up[with(up, order(logFC, order(p, decreasing = FALSE), decreasing = TRUE)), ]
        top2 <- ifelse(top > nrow(sig), nrow(sig), top)
        sig <- sig[1:top2, ]
        graphics::points(
          x = sig$logFC, y = -log10(sig$p), pch = 21, cex = 1.2,
          bg = "green", col = "forestgreen", lwd = 0.5,
          graphics::text(
            x = sig$logFC, y = -log10(sig$p), col = "green3",
            labels = sig$labels, pos = 4, cex = 0.6, font = 2
          )
        )

        ## sig. down-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC A->Z
        sig <- down[with(down, order(logFC, order(p, decreasing = FALSE), decreasing = FALSE)), ]
        top2 <- ifelse(top > nrow(sig), nrow(sig), top)
        sig <- sig[1:top2, ]
        graphics::points(
          x = sig$logFC, y = -log10(sig$p), pch = 21, cex = 1.2,
          bg = "green", col = "forestgreen", lwd = 0.5,
          graphics::text(
            x = sig$logFC, y = -log10(sig$p), col = "green3",
            labels = sig$labels, pos = 2, cex = 0.6, font = 2
          )
        )
      }
    } ## TOP


    if (legend == TRUE) {
      legend("topright", c("up", "not sig.", "down"),
        box.col = "transparent",
        box.lwd = 0, bg = "transparent", border = "black", pch = 22,
        pt.bg = c("#FF0000", "grey50", "#1E90FF"), pt.cex = 1.4,
        col = c("#FF0000", "grey50", "#1E90FF"), horiz = FALSE,
        inset = c(inset, 0), bty = "n", xpd = TRUE, ncol = 1, cex = 0.9
      )
    }
  } ## VOLCANO


  ## ------------------
  ##  [02]  MD PLOT
  ## ------------------
  mdPlot <- function(stats, comparison = NULL, min.pval = 0.055, min.lfc = 1,
                     xlim = NULL, ylim = NULL, sig.type = c("pval", "p.adj"),
                     top = NULL, labels = NULL, inset = -0.3, legend = FALSE) {
    stats <- data.frame(stats)

    ## row names for DIA,TMT,phosphoTMT protein data are a combination
    ## of UniprotID,GN, and protein group id. phospho data = uniprotid, GN,
    ## and phospho site (e.g. S37). This info is parsed into a matrix.
    ncols <- length(unlist(strsplit(rownames(stats)[1], split = "_")))
    # mat  <- matrix(unlist(strsplit(rownames(stats),"_")), ncol=ncols, byrow=TRUE);head(mat)
    # rownames(mat)=rownames(stats)

    ## rownames(stats)[1489] ## "D3ZHA7_RGD1560334_predicted_1943"
    mat <- matrix(nrow = nrow(stats), ncol = 3)
    mat[, 1] <- gsub("_.*", "", rownames(stats)) ## before first _
    mat[, 2] <- sub("_[^_]+$", "", sub("^[^_]*_", "", rownames(stats))) ## after first _ & before last _
    mat[, 3] <- sub(".*_", "", rownames(stats)) ## after last _
    rownames(mat) <- rownames(stats)


    ## annotation info extracted from rownames of protein data
    ## Gene_name will be used to label significant points if top >0
    if (length(grep("S", mat[, 3])) == 0) {
      colnames(mat) <- c("UniprotID", "Gene_name", "id")
      mat <- data.frame(mat)
      mat$labels <- mat$Gene_name
      stats <- cbind(stats, mat[rownames(stats), ])
    }

    ## annotation info extracted from rownames of phospho data
    ## Gene_name(phosphosite) will be used to label sig. points if top >0
    if (length(grep("S", mat[, 3])) > 0) {
      colnames(mat) <- c("UniprotID", "Gene_name", "phosAAPosition")
      mat <- data.frame(mat)
      mat$labels <- paste0(mat$Gene_name, "(", mat$phosAAPosition, ")")
      stats <- cbind(stats, mat[rownames(stats), ])
    }

    ## column named 'p' added to stats data.frame based on sig.type selected. texted for y-axis label also defined
    if (sig.type == "p.adj") {
      stats$p <- stats$adj.P.Val
      main.title <- paste(comparison, "(adj. p-value)")
      filename <- paste(comparison, "MD_plot.png", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }

    if (sig.type == "pval") {
      stats$p <- stats$P.Value
      main.title <- paste(comparison, "(p-value)")
      filename <- paste(comparison, "MD_plot_pvalue.png", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }


    ## determine x/y limits for the plot.
    if (length(xlim) != 2) {
      xlim <- c(min(stats$AveExpr, na.rm = TRUE), max(stats$AveExpr, na.rm = TRUE) * 1.1)
    }
    if (length(ylim) != 2) {
      ylim <- c(-1, 1) * max(abs(stats$logFC), na.rm = TRUE) * 1.1
    }

    ## global par plot parameters
    graphics::par(
      font.axis = 1, lwd = 1, font.main = 2, cex.main = 1.5, fg = "black",
      col.axis = "black", cex.axis = 1.2, cex.lab = 1.3, font.axis = 1, xaxs = "i", yaxs = "i"
    )
    if (legend == TRUE) {
      graphics::par(mar = c(6, 7, 4, 6))
    }
    if (legend == FALSE) {
      graphics::par(mar = c(6, 7, 4, 3))
    }

    ## create blank plot area
    plot(x = 1, las = 1, type = "n", main = "", xlab = "", ylab = "", xlim = xlim, ylim = ylim)
    graphics::title(line = 1, main = main.title, cex.main = 1.3)
    graphics::mtext(side = 1, line = 3, cex = 1.3, text = "avg. expression")
    graphics::mtext(side = 2, line = 3.5, cex = 1.3, text = expression(paste(log[2], " (fold-change)", sep = "")))

    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2], graphics::par("usr")[4], col = "gray94")
    graphics::grid(col = "white", lty = "dotted", lwd = 1.2)
    graphics::box(which = "plot", lwd = 4, lty = 1, col = "gray94")

    ## subset based on significance cutoffs
    up <- subset(stats, p <= min.pval & logFC >= min.lfc)
    down <- subset(stats, p <= min.pval & logFC <= -min.lfc)

    ## add points
    with(stats, graphics::points(AveExpr, logFC,
      pch = 21, bg = scales::alpha(colour = "grey60", alpha = 0.5), ## black
      col = scales::alpha(colour = "grey50", alpha = 0.5), lwd = 0.5, cex = 1.1
    ))
    with(up, graphics::points(AveExpr, logFC,
      pch = 21, bg = scales::alpha(colour = "#FF0000", alpha = 1), ## red
      col = scales::alpha(colour = "#990000", alpha = 1), lwd = 0.5, cex = 1.1
    ))
    with(down, graphics::points(AveExpr, logFC,
      pch = 21, bg = "#1E90FF", col = "#003399", lwd = 0.5, ## blue
      cex = 1
    ))
    graphics::abline(h = c(-min.lfc, min.lfc), col = "dimgray", lty = 2, lwd = 0.8)


    ## ADD LABELS

    ## user  can supply a vector of labels to add to the points on the volcano plot.
    ## for protein data the list can be Gene_name, UniprotID, UniprotID_Gene_name, proKey,
    ## for phospho data phosKey, Gene_name or Gene_name(site).
    ## labels added to corresponding points in orange.
    if (!is.null(labels)) {
      if (length(labels) > 0) {

        ## get row numbers in stats corresponding to the labels list
        ## remove any labels that were not identified and subset
        ## stats. divide subdata into logFC>0 and logFC<0 so labels can be
        ## positioned on the R and L side of the point
        n <- keep <- c()
        for (i in 1:length(labels)) {
          tmp <- grep(paste0(labels[i], "_"), rownames(stats))
          tmp2 <- ifelse(length(tmp) == 0, FALSE, TRUE)
          keep <- c(keep, tmp2)
          table(keep)
          n <- c(n, tmp)
        }

        ## ADD USER INPUT LABELS TO POINTS
        # labs <- labels#[keep==T];labs
        subdata <- stats[n, ]

        ## ADD USER INPUT LABELS TO POS logFC > 0
        ## subset data for positive logFC values, add labels
        subup <- subset(subdata, subdata$logFC > 0)
        if (nrow(subup) > 0) {
          graphics::points(
            x = subup$AveExpr, y = subup$logFC, pch = 21, cex = 1.2,
            bg = "orange", col = "darkorange3", lwd = 0.5,
            graphics::text(
              x = subup$AveExpr, y = subup$logFC, col = "orange",
              labels = subup$labels, pos = 4, cex = 0.6, font = 2
            )
          )
        }

        ## ADD USER INPUT LABELS TO NEG logFC < 0
        subdn <- subset(subdata, subdata$logFC < 0)
        if (nrow(subdn > 0)) {
          graphics::points(
            x = subdn$AveExpr, y = subdn$logFC, pch = 21, cex = 1.2,
            bg = "orange", col = "darkorange3", lwd = 0.5,
            graphics::text(
              x = subdn$AveExpr, y = subdn$logFC, col = "orange",
              labels = subdn$labels, pos = 2, cex = 0.6, font = 2
            )
          )
        }
      } ## LEN >0
    } ## LABELS


    ## TOP

    ## the user can use 'top' to label the most sig. up and down regulated
    ## features (proteins=gene name, phospho=gene name(site).
    ## e.g. top=5, will label the top 5 up and down reg.
    if (!is.null(top)) { ## TOP
      if (length(top) == 1 & top > 0) {
        ## sig. up-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC Z-> A
        sig <- up[with(up, order(logFC, order(p, decreasing = FALSE), decreasing = TRUE)), ]
        top2 <- ifelse(top > nrow(sig), nrow(sig), top)
        sig <- sig[1:top2, ]
        graphics::points(
          x = sig$AveExpr, y = sig$logFC, pch = 21, cex = 1.2,
          bg = "green", col = "forestgreen", lwd = 0.5,
          graphics::text(
            x = sig$AveExpr, y = sig$logFC, col = "green3",
            labels = sig$labels, pos = 4, cex = 0.6, font = 2
          )
        )

        ## sig. down-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC A->Z
        sig <- down[with(down, order(logFC, order(p, decreasing = FALSE), decreasing = FALSE)), ]
        top2 <- ifelse(top > nrow(sig), nrow(sig), top)
        sig <- sig[1:top2, ]
        graphics::points(
          x = sig$AveExpr, y = sig$logFC, pch = 21, cex = 1.2,
          bg = "green", col = "forestgreen", lwd = 0.5,
          graphics::text(
            x = sig$AveExpr, y = sig$logFC, col = "green3",
            labels = sig$labels, pos = 4, cex = 0.6, font = 2
          )
        )
      }
    } ## TOP


    if (legend == TRUE) {
      legend("topright", c("up", "not sig.", "down"),
        box.col = "transparent",
        box.lwd = 0, bg = "transparent", border = "black", pch = 22,
        pt.bg = c("#FF0000", "grey50", "#1E90FF"), pt.cex = 1.4,
        col = c("#FF0000", "grey50", "#1E90FF"), horiz = FALSE,
        inset = c(inset, 0), bty = "n", xpd = TRUE, ncol = 1, cex = 0.9
      )
    }
  } ## MD PLOT


  ## -------------------------------
  ##  [03]  GLIMMA VOLCANO PLOT
  ## -------------------------------
  glimmaVolcanoPlot <- function(stats, comparison, data, annot, groups, min.pval = 0.055, min.lfc = 1,
                                sig.type = c("p.adj", "pval"), pipe = "DIA",
                                enrich = "protein", dir = ".") {
    pipe <- match.arg(arg = pipe, choices = c("DIA", "TMT", "phosphoTMT", "LF"), several.ok = FALSE)
    enrich <- match.arg(arg = enrich, choices = c("protein", "phospho"), several.ok = FALSE)
    # if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)}

    stats <- data.frame(stats)

    ## column named 'p' added to stats data.frame based on sig.type selected.
    ## text for y-axis label also defined
    if (sig.type == "p.adj") {
      stats$p <- stats$adj.P.Val
      ylab <- expression(paste("-", log[10], " (adj. p-value)", sep = ""))
      filename <- paste(comparison, "volcano_plot", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }
    if (sig.type == "pval") {
      stats$p <- stats$P.Value
      ylab <- expression(paste("-", log[10], " (p-value)", sep = ""))
      filename <- paste(comparison, "volcano_plot_pvalue", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }

    ## remove rows with adj.PValue==NA or P.Value==NA
    stats <- stats[!is.na(stats$p), ]

    ## sort data, annotation to match order in stats
    ## replace NAs in data with zeros (visualization purposes)
    ## make sure the number of groups entered matches the number of sample columns in data
    if (all(rownames(stats) %in% rownames(data))) {
      data <- data[rownames(stats), ]
      data[is.na(data)] <- 0
    } else {
      stop("Error! data and stats do not match. check rownames.")
    }
    if (all(rownames(stats) %in% rownames(annot))) {
      annot <- annot[rownames(stats), ]
    } else {
      stop("Error! annot and stats do not match. check rownames.")
    }
    if (all(rownames(annot) == rownames(data)) == FALSE) {
      stop("Error! Rownames of data and annotation do not match.")
    }
    if (length(groups) != ncol(data)) {
      stop("Error! The number of groups does not equal the number of data columns.")
    }


    ## add status column to stats (determines colors of significant dots in plot)
    stats$status <- rep(0, nrow(stats))
    stats$status[(stats$p <= min.pval & stats$logFC >= min.lfc) == TRUE] <- 1
    stats$status[(stats$p <= min.pval & stats$logFC <= -min.lfc) == TRUE] <- -1


    ## DETERMINE ANNOTATION COLUMNS OF INTEREST
    if (pipe == "DIA" & enrich == "protein") {
      glimmaAnnotationColums <- c("id", "UniprotID", "Gene_name", "Description")
    }
    if ((pipe == "TMT" | pipe == "LF") & enrich == "protein") {
      glimmaAnnotationColums <- c("id", "UniprotID", "Gene_name", "Description")
    }
    if (pipe == "phosphoTMT" & enrich == "protein") {
      glimmaAnnotationColums <- c("id", "UniprotID", "Gene_name", "Description")
    }
    if (pipe == "phosphoTMT" & enrich == "phospho") {
      glimmaAnnotationColums <- c("id", "proGroupID", "UniprotID", "Gene_name", "Description", "phosAAPosition")
    }
    # if(any(glimmaAnnotationColums %in% colnames(annot))==FALSE){
    #    stop("Error! glimmaAnnotationColums are not all present in annot.")}
    if (all(glimmaAnnotationColums %in% colnames(annot))) {
      glimmaAnnot <- annot[, glimmaAnnotationColums]
    } else {
      stop("Error! glimmaAnnotationColums are not all present in annotation.")
    }

    Glimma::glXYPlot(
      x = stats$logFC,
      y = -log(stats$p, 10),
      counts = data[rownames(stats), ],
      groups = groups,
      samples = colnames(data),
      status = stats$status,
      transform = FALSE,
      anno = glimmaAnnot[rownames(stats), ],
      display.columns = colnames(glimmaAnnot),
      xlab = "logFC",
      # ylab      = paste0(" ","-","log10 (",ifelse(sig.type=="p.adj","adj. p-value","p-value"),")"),
      ylab = paste0("negLog10(", ifelse(sig.type == "p.adj", "adj.P", "P"), ")"),
      # ylab      = paste("neg log10",sig.type),
      side.main = "Gene_name",
      side.xlab = "Group",
      side.ylab = "Normalized Intensity",
      side.log = FALSE,
      sample.cols = colorGroup2(groups)[groups],
      cols = c("#00bfff", "#858585", "#ff3030"),
      jitter = 10,
      path = file.path(dir),
      folder = "Volcano-Plots",
      html = filename,
      main = paste0(comparison, " (", ifelse(sig.type == "p.adj", "adj. p-value", "p-value"), ")"),
      launch = FALSE
    )


    fixfile <- file.path(dir, "Volcano-Plots/js", paste0(filename, ".js"))
    tmptxt <- readLines(fixfile)
    tmptxt <- gsub("logCPM", "Intensity", tmptxt)
    writeLines(as.character(tmptxt), file(fixfile), sep = "\n")
  } ## GLIMMA VOLCANO PLOT


  ## --------------------------
  ##  [04]  GLIMMA MD PLOT
  ## --------------------------
  glimmaMDPlot <- function(stats, comparison, data, annot, groups, min.pval = 0.05, min.lfc = 1,
                           sig.type = c("p.adj", "pval"), pipe = "DIA",
                           enrich = "protein", dir = ".") {
    pipe <- match.arg(arg = pipe, choices = c("DIA", "TMT", "phosphoTMT", "LF"), several.ok = FALSE)
    enrich <- match.arg(arg = enrich, choices = c("protein", "phospho"), several.ok = FALSE)
    # if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)}

    stats <- data.frame(stats)

    ## column named 'p' added to stats data.frame based on sig.type selected.
    ## text for y-axis label also defined
    if (sig.type == "p.adj") {
      stats$p <- stats$adj.P.Val
      filename <- paste(comparison, "MD_plot", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }
    if (sig.type == "pval") {
      stats$p <- stats$P.Value
      filename <- paste(comparison, "MD_plot_pvalue", sep = ifelse(is.null(comparison) | comparison == "", "", "_"))
    }

    ## remove rows with adj.PValue==NA or P.Value==NA
    stats <- stats[!is.na(stats$p), ]


    ## sort data, annotation to match order in stats
    ## replace NAs in data with zeros (visualization purposes)
    ## make sure the number of groups entered matches the number of sample columns in data
    if (all(rownames(stats) %in% rownames(data))) {
      data <- data[rownames(stats), ]
      data[is.na(data)] <- 0
    } else {
      stop("Error! data and stats do not match. check rownames.")
    }
    if (all(rownames(stats) %in% rownames(annot))) {
      annot <- annot[rownames(stats), ]
    } else {
      stop("Error! annot and stats do not match. check rownames.")
    }
    if (all(rownames(annot) == rownames(data)) == FALSE) {
      stop("Error! Rownames of data and annotation do not match.")
    }
    if (length(groups) != ncol(data)) {
      stop("Error! The number of groups does not equal the number of data columns.")
    }


    ## add status column to stats (determines colors of significant dots in plot)
    stats$status <- rep(0, nrow(stats))
    stats$status[(stats$p <= min.pval & stats$logFC >= min.lfc) == TRUE] <- 1
    stats$status[(stats$p <= min.pval & stats$logFC <= -min.lfc) == TRUE] <- -1


    ## DETERMINE ANNOTATION COLUMNS OF INTEREST
    if (pipe == "DIA" & enrich == "protein") {
      glimmaAnnotationColums <- c("id", "UniprotID", "Gene_name", "Description")
    }
    if ((pipe == "TMT" | pipe == "LF") & enrich == "protein") {
      glimmaAnnotationColums <- c("id", "UniprotID", "Gene_name", "Description")
    }
    if (pipe == "phosphoTMT" & enrich == "protein") {
      glimmaAnnotationColums <- c("id", "UniprotID", "Gene_name", "Description")
    }
    if (pipe == "phosphoTMT" & enrich == "phospho") {
      glimmaAnnotationColums <- c("id", "proGroupID", "UniprotID", "Gene_name", "Description", "phosAAPosition")
    }
    if (all(glimmaAnnotationColums %in% colnames(annot))) {
      glimmaAnnot <- annot[, glimmaAnnotationColums]
    } else {
      stop("Error! glimmaAnnotationColums are not all present in annotation.")
    }


    Glimma::glXYPlot(
      x = stats$AveExpr,
      y = stats$logFC,
      counts = data[rownames(stats), ],
      groups = groups,
      samples = colnames(data),
      status = stats$status,
      transform = FALSE,
      anno = glimmaAnnot[rownames(stats), ],
      display.columns = colnames(glimmaAnnot),
      xlab = "AveExpr",
      ylab = "logFC",
      # ylab      = paste("neg log10",sig.type),
      side.main = "Gene_name",
      side.xlab = "Group",
      side.ylab = "Normalized Intensity",
      side.log = FALSE,
      sample.cols = colorGroup2(groups)[groups],
      cols = c("#00bfff", "#858585", "#ff3030"),
      jitter = 10,
      path = file.path(dir),
      folder = "MD-Plots",
      html = filename,
      main = paste0(comparison, " (", ifelse(sig.type == "p.adj", "adj. p-value", "p-value"), ")"),
      launch = FALSE
    )

    fixfile <- file.path(dir, "MD-Plots/js", paste0(filename, ".js"))
    tmptxt <- readLines(fixfile)
    tmptxt <- gsub("logCPM", "Intensity", tmptxt)
    writeLines(as.character(tmptxt), file(fixfile), sep = "\n")
  } ## GLIMMA MD PLOT



  ## --------------------------
  ##  [05] pvalueHistogram
  ## --------------------------
  pvalueHistogram <- function(stats, comparison = NULL) {
    stats <- data.frame(stats)

    graphics::layout(matrix(1:2, ncol = 2, byrow = TRUE))
    ##   P-VALUE HISTOGRAM
    ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6164648/#B1-high-throughput-07-00023
    ## inspect histogram of unadjusted P.Values for the differential expression results.
    ## Figure shows histogram of p-values from gene-by-gene statistical tests.
    graphics::par(mar = c(6, 7, 5, 5))
    graphics::hist(stats$P.Value,
      col = "dodgerblue", breaks = 50,
      main = "", xlab = "", ylab = "",
      xlim = c(0, 1), las = 1, cex.axis = 1.2
    )
    graphics::title(line = 1, main = comparison, cex.main = 1.3)
    graphics::mtext(side = 1, line = 3, cex = 1.3, text = "p-value")
    graphics::mtext(side = 2, line = 4, cex = 1.3, text = "Frequency")

    ##   ADJUSTED P-VALUE HISTOGRAM
    graphics::par(mar = c(6, 7, 5, 5))
    graphics::hist(stats$adj.P.Val,
      col = "dodgerblue", breaks = 50,
      main = "", xlab = "", ylab = "",
      xlim = c(0, 1), las = 1, cex.axis = 1.2
    )
    graphics::title(line = 1, main = comparison, cex.main = 1.3)
    graphics::mtext(side = 1, line = 3, cex = 1.3, text = "adj. p-value")
    graphics::mtext(side = 2, line = 4, cex = 1.3, text = "Frequency")
  } ## PVALUE HISTOGRAM
} ## DE PLOT FUNCTIONS
