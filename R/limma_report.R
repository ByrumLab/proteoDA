make_limma_spreadsheets <-


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




# Some results p.csv saving taken out of extract_limma_results.

dummy_function_2 <- function(){

## COMBO STATS
comboStats <- NULL
for (x in names(statList)) {
  tmp <- statList[[x]]
  colnames(tmp) <- paste(colnames(tmp), x, sep = "_")
  if (!is.null(comboStats)) {
    comboStats <- cbind(comboStats, tmp[rownames(comboStats), ])
  }
  if (is.null(comboStats)) {
    comboStats <- tmp
  }
}

## SAVE LIMMA STAT RESULTS
## save stat results for individual contrasts as csv files.
## save combined stat results as a csv file in dir. also
## save BQ combined stat results (NAs /blanks replaced with zeros)
## as csv in project directory for Big Query upload.
if (save) { ## SAVE STAT RESULTS

  ## INDIVIDUAL STATS
  base::lapply(names(statList), function(x) {
    stats <- statList[[x]]
    stats2 <- cbind(annot[rownames(stats), ], data[rownames(stats), ], stats)
    filename <- paste0(x, "_results.csv")
    utils::write.csv(stats2, file = file.path(dir, filename), row.names = FALSE)
  })
}

## COMBINED STATS
comboStats2 <- cbind(annot[rownames(comboStats), ], data[rownames(comboStats), ], comboStats)
if (save) {
  filename <- paste("combined_results.csv", sep = "_")
  utils::write.csv(comboStats2, file = file.path(dir, filename), row.names = FALSE)
}

## COMBINED STATS FOR BIG QUERY
comboStats2[, ][is.na(comboStats2)] <- 0
comboStats2[, ][comboStats2 == ""] <- 0
if (any(substr(colnames(comboStats2), start = 1, stop = 1) %in% c(0:9)) == TRUE) {
  colnames(comboStats2) <- paste0("X", colnames(comboStats2))
}
if (save) {
  filename <- paste(ilab, enrich, "results_BQ.csv", sep = "_")
  if (file.exists(file.path(".", filename))) {
    print("BQ file already exists. creating a new BQ filename...")
    filename <- make_new_filename(x = filename, dir = ".")
  }
  utils::write.csv(comboStats2, file = file.path(".", filename), row.names = FALSE)
  print("limma stat results for BQ saved. Success!!")
}

}



