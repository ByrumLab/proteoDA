make_limma_de_report <- function(data, annot, targets, design, contrasts, min.pval=0.055, min.lfc=1,
                                          adj.method="BH", paired=FALSE, pipe="DIA", enrich="protein",
                                          dir=NULL, save=TRUE, ilab="PI_DATE") {



  param[["dir"]] <- ifelse(save==TRUE, file.path(dir),"NULL")

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

  }


  invisible(NULL)
}

