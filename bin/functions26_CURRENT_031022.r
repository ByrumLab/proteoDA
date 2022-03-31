


#####################################
##       REQUIRED PACKAGES         ##
##################################### 
{ ## LIBRARIES
   
   ##------------------------------
   ##  LOAD REQUIRED PACKAGES
   ##------------------------------
   ## get_pkgs() function: install and load multiple R packages.
   ## checks to see if packages are installed. Install them if they are not, 
   ## then load them into the R session.
   if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
   
   pkgs <- c("circlize",         ## heatmaps
             "ClassDiscovery",   ## sample dendrograms
             "biomaRt",          ## ensembl annotation info from biomart
             # "clusterProfiler",## ORA and GSEA
             "ComplexHeatmap",   ## heatmaps
             # "dendextend",     ## dendrograms
             # "DOSE",           ## clusterProfiler plots
             "dplyr",            ## used to rank normalization methods 
             # "enrichplot",     ## plots from clusterProfiler
             # "edgeR",          ## coolmap() other functions
             "genefilter",       ## proteiNorm plots
             # "gplots",         ## plots, heatmaps.2
             "Glimma",           ## interactive volcano/MD plots
             "grid",             ## plots
             "gridExtra",        ## arranging plot grids
             "impute",           ## imputation methods
             "imputeLCMD",       ## imputation methods
             # "leapR",          ## data integration
             "limma",            ## DE analysis moderated t-test
             "matrixStats",      ##
             # "mixOmics",       ## multiomics data integration
             # "mogsa",          ## multiomics pathway
             # "msigdbr",        ## MSigDB package genesets
             "NormalyzerDE",     ## normalization methods
             "openxlsx",         ## excel input/output files
             "pcaMethods",       ## PCA plots
             # "pdftools",       ## combine pdfs into a single file (QC_Results.pdf)
             "png",              ## used to create QC plots 
             # "phosR",          ## phospho site analysis
             "psych",            ## correlation scatter plot matrix
             "RColorBrewer",     ## plot colors
             "scales",           ## alpha transparency level of plot colors
             "statmod",          ## required by other packages
             "stringr",          ## used to parse Fasta.header
             "tidyverse",        ## miscellaneous
             "tools",            ## file extensions, basenames, file import
             "unikn",            ## p-value histogram
             "utils",            ## read/write csv,txt,tsv files
             "vioplot",          ## violin plots
             "yarrr"             ## plot colors 
   )
   
   get_pkgs <- function(pkg){
      new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
      if (length(new.pkg)){
         # install.packages(new.pkg, dependencies = TRUE)
         BiocManager::install(new.pkg, dependencies=TRUE, force=TRUE)
         sapply(pkg, require, character.only = TRUE)
      }
   }
   get_pkgs(pkgs)
   
   
} ## LIBRARIES


#####################################
##       GENERAL PARAMETERS        ##
#####################################
{ ## GENERAL
   
   op <- par(no.readonly = TRUE) ## restores par() plot parameters
   `%notin%` <- Negate(`%in%`)   ## opposite of %in% (not in list = TRUE)
   
} ## GENERAL


#####################################
##          VARIOUS LISTS          ##
#####################################
{ ## LISTS
   
   ## TMT/phosphoTMT/peptide: annotation/contamination columns in MaxQuant peptide.txt file
   peptideAnnotationColums <- c("id", "Protein.group.IDs", "Leading.razor.protein", "Gene.names")
   peptideContamColums     <- c("Reverse", "Potential.contaminant")
   
   ## TMT/protein: annotation/contamination columns in MaxQuant proteinGroups.txt file
   tmtAnnotationColums <- c("id", "Majority.protein.IDs", "Fasta.headers", "Score")
   tmtContamColums     <- c("Reverse", "Potential.contaminant", "Only.identified.by.site")
   
   ## phosphoTMT/protein: annotation/contamination columns in MaxQuant proteinGroups.txt file
   proteinAnnotationColums <- c("id", "Majority.protein.IDs", "Fasta.headers", "Score", 
                                "Phospho..STY..site.IDs")
   proteinContamColums     <- c("Reverse", "Potential.contaminant", "Only.identified.by.site")
   
   ## phosphoTMT/phosphopeptide: annotation/contamination/site probability columns in 
   ## MaxQuant Phospho (STY)Sites.txt file   
   phosphoAnnotationColums <- c("Proteins", "Positions.within.proteins", "Leading.proteins",
                                "Protein", "Protein.names", "Gene.names", "Fasta.headers",
                                "Localization.prob", "Score.diff", "PEP","Score", 
                                "Amino.acid", "Sequence.window","Phospho..STY..Probabilities",
                                "Charge", "id", "Protein.group.IDs", "Positions", "Position",
                                "Peptide.IDs", "Mod..peptide.IDs", "Evidence.IDs")
   phosphoContamColums     <- c("Reverse", "Potential.contaminant")
   phosphoLocalProbColum   <- c("Localization.prob")
   
   ## DIA: protein annotation/contaminantion columns from Scaffold DIA Sample Report csv file
   diaAnnotationColums <- c("id", "Protein.Name","Accession.Number","Molecular.Weight",
                            "Protein.Group.Score", "Identified.Peptide.Count", 
                            "Exclusivity","Quantitative.Value")
   diaContamColums     <- c("Protein.Name") ## Group by / DECOY_
   
   ## METADATA COLUMNS
   tmtMetaColums <- c("sample", "batch", "tag", "group")
   diaMetaColums <- c("sample", "number","group")
   lfMetaColums  <- c("sample","group")
   targetColums  <- c("sampleIDs", "sample", "group")
   
   ## named lists of normalization and imputation methods
   norm.methods        <- c("log2", "median", "mean", "vsn","quantile", "cycloess", "rlr", "gi")
   names(norm.methods) <- c("Log2", "Median", "Mean", "VSN","Quantile", "Cyclic Loess", "RLR", "Global Intensity")
   impute.methods          <- c("none", "knn", "qrilc", "mindet", "minprob", "min", "zero")
   names(impute.methods)   <- c("None","KNN", "QRILC", "MinDet", "MinProb", "Min", "Zero")
   
   ## list of reporter ions (tags) for various TMT kits
   TMT6plex  <- c("126","127","128","129","130","131");names(TMT6plex) <-c(1:6);TMT6plex
   TMT10plex <- c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N");names(TMT10plex) <-c(1:10);TMT10plex
   TMT11plex <- c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C");names(TMT11plex) <- c(1:11);TMT11plex
   TMT16plex <- c("126","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N");names(TMT16plex)<-c(1:16);TMT16plex
   
   ## column names of limma stat output
   limmaStatColums <- c("logFC","CI.L","CI.R","AveExpr","t","B","P.Value","adj.P.Val")
   
   
   ## create a color palette for plots called binfcolors
   ## color names courtesy of my nieces and nephews (ages 6, 9, and 12)
   binfcolors <- c("#1F5EDC","#EE0010","#32CD32","#FF1493","#FF7F00",
                   "#A342FC","#00C8FF","#ADFF2F","#FFE100","#E36EF6","#009ACE","#996633")
   names(binfcolors) <- c("blueberry","cherry","apple","barbie","fanta",
                          "grape","ocean","mtndew","gold","orchid","aceblue","poop")
   ## unikn::seecol(pal=binfcolors,title="color palette: (binfcolors)")
   
   lightbinfcolors <- c("#c8d8f7","#ffdadd","#d0f3d0","#ffb1db","#ffe1c4",
                        "#dbb6fe","#c4f2ff","#e3ffb8","#fff8c4","#f1b8fb",
                        "#baeeff","#eddcca")
   names(lightbinfcolors) <- c("lightblueberry","lightcherry","lightapple","lightbarbie","lightfanta",
                               "lightgrape","lightocean","lightmtndew","lightgold","lightorchid",
                               "lightaceblue","lightpoop")
   ## unikn::seecol(pal=lightbinfcolors,title="color palette: (lightbinfcolors)")
   
   
} ## LISTS


#####################################
##        HELPER FUNCTIONS         ##
#####################################
{ ## HELPER FUNCTIONS
   
   ##-----------------
   ##  [01] ADD EXCEL
   ##-----------------
   ## FUNCTIONS
   ## writes data to excel worksheet 
   ## wb=existing wb object or if NULL then a new wb is created using createWorkbook()
   ## sheet="mysheetname" or sheet=1, if NULL then Sheet1,Sheet2 etc.
   ## x=data to write, rownames are included
   ## open=logical to open the wb after data is written (this is a temp file) T wb is opened F=wb is not
   ## wb object is returned if wb is NULL then use mywb <- add_excel() to store results 
   ## in mywb rather than wb which would overwrite existing wb in GlobalEnv
   ## USE: wb<-createWorkbook()
   ##      add_excel(wb=wb,sheet="mysheet", x=mydf, open=F)
   ##      wb2<-add_excel(wb=NULL,sheet=NULL,x=mydf,open=F) 
   ##      saveWorkbook(wb2, file="workbookname.xlsx",overwrite=TRUE,returnValue=TRUE)
   add_excel <- function(x, wb=NULL, sheet=NULL, rowNames=FALSE, open=FALSE){
      library(openxlsx)
      colStyle <- openxlsx::createStyle(fgFill="#c0c0c0", halign="center", valign="center", 
                                        textDecoration="bold", fontColour="black", border="TopBottomLeftRight",
                                        fontSize=11, numFmt="TEXT")
      if(is.null(wb)){ wb <- createWorkbook() }
      if(is.null(sheet)){ sheet=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheet
      ifelse(nchar(sheet) > 30, substring(sheet,1,30), sheet)
      
      addWorksheet(wb,sheetName=sheet, tabColour="dodgerblue2")
      writeData(wb, sheet=sheet, x=x, borders="columns", rowNames=rowNames,
                headerStyle=colStyle, keepNA=TRUE, na.string="NA", withFilter=TRUE)
      setRowHeights(wb, sheet, rows=1, heights=20)
      print(paste("data was added to ", sheet))
      activeSheet(wb) <- sheet
      
      if(open==TRUE){ openXL(wb) }
      return(wb)
      
   } ## ADD EXCEL
   
   
   ##---------------------
   ##  [02] MAKE FACTOR
   ##---------------------
   ## x = vector of values that will be converted to a factor 
   ## prefix = character string; text added to values in x if x is numeric
   ## e.g. prefix="snap" x=c(1,1,2,2,3) == c(snap1,snap1,snap1,snap2,snap2,snap3) 
   ## if x is numeric vector and prefix=NULL then default text "A" added
   ## e.g. prefix=NULL x=c(1,1,2,2,3) == c(A1,A1,A1,A2,A2,A3) 
   ## if x is already a character vector prefix is ignored.
   ## e.g. prefix="crackle" x=c(B1,B1,B1,B2,B2,B3) == c(B1,B1,B1,B2,B2,B3)
   ## e.g. prefix=NULL x=c(B1,B1,B1,B2,B2,B3) == c(B1,B1,B1,B2,B2,B3)
   ## character vector is then converted to a factor, factor levels are 
   ## ordered according to order of unique values in x rather than alphabetically.
   ## e.g. c(C5,C5,C5, C1,C1,C1, C3,C3,C3) levels == c(C5,C1,C3) NOT c(C1,C3,C5)
   ## e.g. c(E,E,E,E, A,A,A, C,C,C) levels == c(E,A,C) NOT c(A,C,E)
   make_factor <- function(x, prefix=NULL){
      if(is.numeric(x)){
         prefix <- ifelse(is.null(prefix), "X", prefix)
         x <- paste(prefix, x, sep="")
      }
      x <- factor(x, levels=ordered(unique(x)))
      return(x)
   } ## MAKE FACTOR
   
   
   ##---------------------------------------
   ##  [03] MAKE ALL CONTRASTS (OPTIONAL)
   ##---------------------------------------
   # Adapted from Gordon Smyth https://support.bioconductor.org/p/9228/
   ## creates contrast matrix for all unique pairwise group combinations
   make_all_contrasts <- function(design, groups){
      
      designCols <- colnames(design)
      groups <- make_factor(groups)
      groupLevels <- levels(groups)
      
      n <- length(groupLevels);n
      stopifnot(identical(groupLevels, designCols[1:n]))
      contrasts <- matrix(0, length(designCols), choose(n, 2))
      rownames(contrasts) <- designCols    
      colnames(contrasts) <- 1:choose(n,2)
      
      k = 0
      for (i in 1:(n-1)){
         for (j in (i+1):n){
            k = k + 1
            contrasts[j, k] = 1
            contrasts[i, k] = -1
            colnames(contrasts)[k] = paste0(groupLevels[j], "_vs_", groupLevels[i])
            print(contrasts)
         }
      }        
      return(contrasts)
      
   } ## MAKE ALL CONTRASTS
   
   
   ##--------------------------------------------
   ##  [04] EXTRACT CONTRAST GROUPS (OPTIONAL)
   ##--------------------------------------------
   ## function called by make_contrasts() function.
   ## input is a character vector of contrasts with format g2_vs_g1=g2-g1
   ## returns a named list of the groups involved in each contrast
   ## for semi-complex contrasts such as g2_vs_g1=g2-g1 or diff.g3=(g3-g2)-(g2-g1)
   ## output format: 
   ## $g2_vs_g1
   ## [1] g2 g1
   extract_contrast_groups <-function(contrast.vec){
      
      contrast.vec <- gsub(" ", "", as.character(contrast.vec)) ## remove blank spaces
      mat   <- matrix(unlist(strsplit(contrast.vec,"=")), ncol=2, byrow=TRUE)
      names <- mat[,1]
      
      vec <- mat[,2]
      vec <- gsub(" ","",as.character(vec)) ## removes blank spaces
      vec <- gsub("\\/+[[:digit:]]","", as.character(vec))  ## removes division by one digit number (/2)
      vec <- gsub("\\-+[[:digit:]]", "", as.character(vec)) ## removes subtraction by one digit number (-2)
      vec <- gsub("\\++[[:digit:]]", "", as.character(vec)) ## removes addition of one digit number (+2)
      vec <- gsub("\\*+[[:digit:]]", "", as.character(vec)) ## removes multiply by one digit number (*2)
      vec <- gsub("\\+", "-", as.character(vec)) ## replaces + with - (+ -> -)
      vec <- gsub("\\(|\\)", "", as.character(vec)) ## removes parentheses ()
      vec <- gsub("\\/", "-", as.character(vec)) ## replace division sign / with -
      vec <- gsub(" ", "", as.character(vec)) ## removes blank spaces
      vec2 <- strsplit(vec, "-")
      names(vec2) <- names
      
      vec3 <- NULL
      for(i in 1:length(vec2)){ vec3[[i]] <- unique(vec2[[i]]) }
      names(vec3) <- names
      
      return(vec3)
      
   } ## CONTRAST GROUPS
   
   
   ##-----------------------------
   ##  [05] COPY TO ENVIRONMENT
   ##-----------------------------
   ## from = env1, to=env2, names=c(character vector of variable names
   ## to copy from env1 to the new env (env2))
   ## env2<-new.env()
   ## copy2env(from=.GlobalEnv, to=env2, names=c("targets, "y.org","norm_counts","ilab"))
   ## ls(env2)
   ## mylist<-as.list.environment(env2) ## converts environment ot a list object 
   copy2env <- function(from, to, names=ls(from, all.names=TRUE)){
      
      mapply(assign, names, mget(names, from), list(to), 
             SIMPLIFY = FALSE, USE.NAMES = FALSE)
      invisible(NULL)
      
   } ## COPY2ENV
   
   
   ##---------------------------
   ##  [06] BATCH PLOT COLORS 
   ##---------------------------
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
         grid::grid.raster(pngraster, width=unit(0.9,"npc"), height= unit(0.8,"npc"))
         if(i < n){ plot.new() }
      }
      dev.off()
      if(delete.files){ unlink(pnglist) }
      
   } ## MERGE PNG PDF
   

   ##---------------------------
   ##  [10] MAKE NEW FILENAME
   ##---------------------------
   ## when saving a new file if the file already exists in a directory. 
   ## create a new filename by adding _## to end.
   make_new_filename <- function(x,dir="."){
      kk=0
      # library(xfun,quietly=TRUE)
      filelist<-c(list.files(path=dir));filelist
      ext <- tools::file_ext(x);ext
      xx=x;xx
      success=FALSE
      while(success==FALSE){
         kk=kk+1
         # xx=sub("_[^_]+$", "", xx);xx# run the function again
         xx <- paste0(gsub(paste0(".",ext),"",xx),"_",kk,".",ext);xx
         success=(file.path(dir,xx)%in%file.path(dir,filelist)==FALSE);success
         if(success==FALSE){
            xx=sub("_[^_]+$", "", xx);xx# run the function again
            xx=x
         };xx
      };xx
      if(success==TRUE){return(xx)}
   } ## MAKE NEW FILE NAME
   
   
   ##-------------------------
   ##  [11] TABLES DT/KABLE
   ##-------------------------
   ## DT TABLE ABOVE BUT AS A FUNCTION
   dt_table <- function(x,cap_text="") {
      # suppressPackageStartupMessages(library(DT))
      DT::datatable(x, class="cell-border stripe", rownames=FALSE, filter="none", 
                    editable=FALSE, extensions="Buttons", 
                    options=list(initComplete=JS("function(settings, json) {","$('body').css({'font-family': 'Calibri'});", "}"), 
                                 dom="Bfrtip", buttons=c("copy", "csv", "excel", "pdf", "print"),
                                 options=list(pageLength = nrow(x), scrollX = TRUE),
                                 caption=htmltools::tags$caption(style = "caption-side: bottom; text-align: left;","",
                                                                 htmltools::em(paste(cap_text," \n", sep = "")))))
      
   } ## DT_TABLE
   
   
   ##---------------------
   ##  [12] DT TABLES 2
   ##---------------------
   dt_table2 <- function(x){
      # suppressPackageStartupMessages(library(DT))
      datatable(x, rownames=FALSE, filter="none", 
                options=list(pageLength=nrow(x), scrollX=TRUE,
                             initComplete = JS("function(settings, json) {","$('body').css({'font-family': 'Arial'});", "}")),
                caption=htmltools::tags$caption(style='caption-side: top; text-align: left;', '', htmltools::em('')))
      
   } ## DT_TABLE 2
   
   
   ##-------------------
   ##  [13] KAB_TABLE
   ##-------------------
   kab_table <- function(x, cap_text=""){
      # suppressPackageStartupMessages(library(kableExtra))      
      kableExtra::kable_styling(kable(x, row.names=FALSE, full_width=TRUE, 
                                      position="center", caption=paste(cap_text,"\n\n"),
                                      format="html"), 
                                c("striped", "bordered", "responsive"), 
                                html_font = "Arial", font_size=14)
      
   } ## KABLE_TABLE
   
   
} ## HELPER FUNCTIONS


#####################################
##    MAIN PROCESSING FUNCTIONS    ##
#####################################
{ ## MAIN FUNCTIONS
   
   
   ##--------------------------------
   ##  [01] IMPORT DATA (REQUIRED)
   ##--------------------------------
   import_data <-function(file, pipe, enrich){
      
      ## check that file exists in specified location 
      ## and is csv/tsv/txt
      file <-gsub("\\./","",file);file
      file.dir <- dirname(file);file.dir
      file <- file.path(file.dir, match.arg(basename(file), list.files(file.dir),
                                            several.ok=FALSE));file
      filext <- tools::file_ext(file);filext
      filext <- match.arg(filext, c("csv","txt","tsv"),several.ok=FALSE)
      if(filext=="txt" | filext=="tsv"){ sep="\t" }
      if(filext=="csv"){ sep=","}
      
      ## vector of required column names based on analysis pipeline 
      ## and type of enrichment. these columns are processed by other
      ## functions that extract annotation, remove contaminants, etc.
      ## values for each variable are defined according to R syntax
      ## rules. Values for each variable are listed at the top of the
      ## functions.R file.
      if(all(pipe=="DIA" & enrich=="protein")){ 
         reqCols<-c(diaAnnotationColums,diaContamColums) 
      }
      if(all(pipe=="TMT" & enrich=="protein")){ 
         reqCols<-c(tmtAnnotationColums, tmtContamColums) 
      }
      if(all(pipe=="phosphoTMT" & enrich=="protein")){ 
         reqCols <- c(proteinAnnotationColums, proteinContamColums) 
      }
      if(all(pipe=="phosphoTMT" & enrich=="phospho")){ 
         reqCols <- c(phosphoAnnotationColums, phosphoContamColums, phosphoLocalProbColum) 
      }
      if(all(pipe=="LF" & enrich=="protein")){ 
         reqCols<-c(tmtAnnotationColums, tmtContamColums) 
      }
      
      ## IMPORT DIA DATA
      if(pipe=="DIA"){
         
         ## sample report generated by scaffold DIA is imported as data.frame. 
         ## Sample Reports generated by Scaffold DIA contain meta data in the 
         ## first 10 lines, and text in the last row to mark the end of the file. 
         ## Thus, when reading in the file the first 10 lines are skipped and the 
         ## last row is removed. column name of column 1 (X. on import) is changed to id
         ## column 1 must be changed to "id" to allow extraction and other processing.
         ## Note: column names of the input file are defined according to R syntax rules.
         ## special characters, spaces, are converted to periods.
         data <- utils::read.csv(file=file, sep=sep, stringsAsFactors=FALSE, 
                                 header=TRUE, check.names=TRUE, skip=10)
         if(data[nrow(data),1]=="END OF FILE"){ data <- data[-nrow(data),] };data[nrow(data),1]
         
         ## name of column 1 is changed to 'id'
         if(colnames(data)[1]=="X."){ colnames(data)[1] <- "id" }
         data <- data.frame(data)
         no_input <- nrow(data)
         
         ## checks to make sure the input file contains the required annotation/contamination
         ## columns. This ensures that the file imported correctly and helps to identify
         ## sampleIDs columns and extract targets info. later.
         if(all(reqCols %in% colnames(data))==TRUE){
            
            print(paste(pipe, enrich, "file imported. Success!!"))
            print(paste("the file contains",no_input, enrich, "entries"));cat("\n")
            
            data2<-list(data=data, no_input=no_input, reqCols=reqCols)
            return(data2)
            
         } else { stop(paste0("Error! The input file may be missing one or more",
                              "required annotation/contamination \ncolumns or ",
                              "id not import properly. \n\nFor DIA experiments ",
                              "the input sample report file should include the following", 
                              "column names:\n", paste(reqCols, collapse="', '"), "\n", 
                              " The annotation/sample data should also begin in row 11 ",
                              " of the report."))
         }
      } ## DIA IMPORT
      
      
      ## IMPORT TMT/phosp DATA
      if(any(pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF")){
         
         ## imports file as data.frame. Note: column names of the input file
         ## are defined according to R syntax rules. special characters, spaces,
         ## are converted to periods.
         data <- utils::read.csv(file=file, sep=sep, stringsAsFactors=FALSE, 
                                 header=TRUE, check.names=TRUE)
         no_input <- nrow(data)
         
         ## checks that file contains required annotation/contaminant columns, 
         ## the id column of MaxQuant output files is a unique key, if the values
         ## are non-numeric or duplicated the file may contain some rows with trash 
         ## information that needs to be removed. if this is the case a msg is returned
         ## alerting the user that the file needs manual inspection. A data.frame 
         ## of the input data is returned if it contains all the required columns
         ## and the id column contains numeric unique values (0-n)
         if(all(reqCols %in% colnames(data))==TRUE){
            if(all(!duplicated(data$id) & is.integer(data$id))){
               
               print(paste(pipe, enrich, "file imported. Success!!"))
               print(paste("the file contains",no_input, enrich, "entries"));cat("\n")
               
               data2 <- list(data=data, no_input=no_input)
               return(data2)
               
            } else { 
               stop(message(paste("Error! 'id' column values. The id' column of",
                                  "maxQuant output files should contain unique",
                                  "integer values (e.g. 0-n). Inspect the file manually.",
                                  "e.g. sort id column a -> z. Remove junk information",
                                  "typically at the end of the file.")))
            }
         } else {
            stop(message(paste("Error! the input file is missing one or more required", 
                               "annotation/contamination columns. For this pipeline",
                               "the file should include the following column names: ",
                               paste(reqCols, collapse=", "))))
         }
         
      } ## TMT IMPORT
      
      
   } ## IMPORT DATA
   
   
   ##----------------------------------------
   ##  [02] REMOVE CONTAMINANTS (REQUIRED)
   ##----------------------------------------
   remove_contaminants <- function(data, pipe, enrich){
      
      ## vector of required contaminant column names based on 
      ## analysis pipeline and type of enrichment. these columns 
      ## are processed by other functions that remove contaminants, etc.
      ## values for each variable are defined according to R syntax
      ## rules. Values for each variable are listed at the top of the
      ## functions.R file.
      if(pipe=="TMT" & enrich=="protein"){ contamColums <- tmtContamColums }
      if(pipe=="phosphoTMT" & enrich=="protein"){ contamColums <- proteinContamColums }
      if(pipe=="phosphoTMT" & enrich=="phospho"){ contamColums <- phosphoContamColums }
      if(pipe=="DIA" & enrich=="protein"){ contamColums <- diaContamColums }
      if(pipe=="LF" & enrich=="protein"){ contamColums <- tmtContamColums }
      
      if(pipe=="DIA"){
         
         no_rows <- nrow(data);no_rows
         data <- data[!grepl("DECOY", data[, contamColums]),];nrow(data)
         data <- data[!grepl("Group of", data[, contamColums]),];nrow(data)
         
         no_contam  <- no_rows - nrow(data)
         no_qfilter <- nrow(data)
         
         print(paste(no_contam, "contaminants removed. Success!!"))
         print(paste(no_qfilter, pipe, enrich, "entries retained ..."))
         
         data2 <- list(data=data, no_contam=no_contam, no_qfilter=no_qfilter)
         return(data2)
         
      } ## DIA CONTAM
      
      if(pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF"){
         
         no_rows <- nrow(data)
         for(x in contamColums){
            data[,x][is.na(data[,x])] <- ""
            remove<-data[,x]=="+"
            data <- data[remove==FALSE, ]
         }
         
         no_contam  <- no_rows - nrow(data)
         no_qfilter <- nrow(data)
         
         print(paste(no_contam, "contaminants removed. Success!!"))
         print(paste(no_qfilter, pipe, enrich, "entries retained ..."))
         
         data2 <- list(data=data, no_contam=no_contam, no_qfilter=no_qfilter)
         return(data2)
         
      } ## TMT/phosphoTMT/LF
      
   } ## REMOVE CONTAM
   
   
   ##---------------------------------------------
   ##  [03] LOCAL PROBABILITY FILTER (REQUIRED)
   ##---------------------------------------------
   local_prob_filter <- function(data, min.prob, pipe, enrich){
      
      no_rows  <- nrow(data)
      min.prob <- ifelse(is.numeric(min.prob),min.prob,0.75)
      remove   <- data[, phosphoLocalProbColum] < min.prob
      data     <- data[remove==FALSE, ]
      
      no_localprob_removed <- no_rows - nrow(data)
      no_localprob_kept    <- nrow(data)
      
      print(paste(no_localprob_removed, pipe, "entries with localization",
                  "probabilities < ", min.prob, "removed. Success!!"))
      print(paste(no_localprob_kept, pipe, "entries kept for further analysis..."))
      
      data2 <- list(data=data, min.prob=min.prob, no_localprob_removed=no_localprob_removed, 
                    no_localprob_kept=no_localprob_kept)
      return(data2)
      
      
   } ## LOCAL PROB FILTER
   
   
   ##---------------------------------------
   ##  [04] EXTRACT SAMPLE IDS (REQUIRED)
   ##---------------------------------------
   ## colNames = input file column names to search. For phospho data if sampleIDs
   ## are input by user they should be Reporter.intensity.corrected.X. without
   ## the ___1 class extension. This allows the function to obtain all intensity
   ## data for single, double, and tripple phosphopeptides.
   extract_sampleIDs <- function(colNames, sampleIDs=NULL, pipe, enrich){
      
      ## the first 11 columns of DIA sample reports contain
      ## 2 col. of junk & 9 annotation info. sample data should
      ## begin in column 12, with pool samples used to create the 
      ## library in the first 3 columns (12-14) followed by the 
      ## experimental samples. 
      ## NOTE: the 2 columns of the file are removed when the file
      ## is imported. Thus samples will start in column 10.
      if(pipe=="DIA" & enrich=="protein"){
         
         if(!is.null(sampleIDs)){
            
            IDs<-sampleIDs
            keep<-colNames%in%IDs;keep;table(keep)
            sampleIDs<-colNames[keep==TRUE]
            stopifnot(length(sampleIDs)==length(IDs));sampleIDs
            if(identical(sampleIDs, character(0))){
               stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
            }
         } ## NOT NULL
         
         if(is.null(sampleIDs)){
            ## columns containing sampleIDs should end with .mzML extension.
            keep <- grep("mzML", colNames);keep
            sampleIDs <- colNames[keep];sampleIDs
            
            ## if any annotation/contaminant column names are still included in list
            ## sampleIDs remove them.
            remove    <- sampleIDs %in% unique(diaAnnotationColums,diaContamColums);remove
            sampleIDs <- sampleIDs[remove==FALSE];sampleIDs
            
         } ## NULL
         
         
         data2 <- list(sampleIDs=sampleIDs)
         return(data2)
         
      } ## DIA
      
      ## label free looking for iBAQ.XXX columns. user should input ids
      if(pipe=="LF" & enrich=="protein"){
         
         if(!is.null(sampleIDs)){
            
            IDs<-sampleIDs
            pattern="NULL"
            keep<-colNames%in%IDs;keep;table(keep)
            sampleIDs<-colNames[keep==TRUE]
            stopifnot(length(sampleIDs)==length(IDs));sampleIDs
            if(identical(sampleIDs, character(0))){
               stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
            }
         } ## NOT NULL
         
         if(is.null(sampleIDs)){
            
            pattern="iBAQ."
            sampleIDs <- grep(pattern, colNames, value=TRUE)
            sampleIDs <- sampleIDs[-grep("iBAQ.peptides",sampleIDs)];sampleIDs
            
            if(identical(sampleIDs, character(0))){
               stop(paste0("\nsampleIDs could not be identified based on pattern='",pattern,
                           "' sampleIDs==character(0)"))
            }
            
         } ## NULL
         
         data2<-list(sampleIDs=sampleIDs, pattern=pattern)
         return(data2)
         
      } ## LF
      
      ## TMT/phosphoTMT
      if(pipe=="TMT" | pipe=="phosphoTMT"){
         
         if(enrich=="protein"){
            
            if(!is.null(sampleIDs)){
               
               pattern="NULL"               
               IDs<-sampleIDs
               
               keep<-colNames%in%IDs;keep;table(keep)
               sampleIDs<-colNames[keep==TRUE]
               stopifnot(length(sampleIDs)==length(IDs));sampleIDs
               if(identical(sampleIDs, character(0))){
                  stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
               }
            } ## NOT NULL
            
            
            ## MaxQuant TMT sample data columns of interest. X is the index 
            ## of the TMT tag starting with 0.
            ## https://adinasarapu.github.io/posts/2020/01/blog-post-tmt/
            ## proteinGroups.txt		   "Reporter intensity corrected X <experiment name(s)>"
            ## peptides.txt			   "Reporter intensity corrected X <experiment name(s)>"
            ## Phospho (STY)Sites.txt	"Reporter intensity corrected X <experiment name(s)>"
            
            ## get the corrected reporter intensity column names for each protein
            ## lysate sample 
            ## "Reporter.intensity.corrected.X.<experiment name(s)> _Lysate"
            
            if(is.null(sampleIDs)){
               pattern <- "Reporter\\.intensity\\.corrected\\.[[:digit:]]+.*[L/l]ysate*"
               sampleIDs <- grep(pattern, colNames, value=TRUE)
               
               ## if pattern doesn't match any column names then reduce pattern to more 
               ## generic format. 
               ## "Reporter.intensity.corrected.X" 
               if(identical(sampleIDs, character(0))){
                  pattern = "Reporter\\.intensity\\.corrected\\.[[:digit:]]+"
                  sampleIDs = grep(pattern, colNames, value=TRUE);sampleIDs
                  cat("sampleIDs identified using generic pattern.")
               }
               
            }## NULL
            
            data2 <- list(sampleIDs=sampleIDs, pattern=pattern)
            return(data2)
            
         } ## PROTEIN 
         
         if(enrich=="phospho"){
            
            if(!is.null(sampleIDs)){
               
               pattern="NULL"               
               IDs<-sampleIDs
               keep<-colNames%in%IDs;keep;table(keep)
               
               classIDs  <- colNames[keep==TRUE]
               stopifnot(length(classIDs)==length(IDs));classIDs
               classNums <- gsub("^.*_{3}","",classIDs) ## keep everything after ____ (class number)
               sampleIDs <- gsub("_{3}.*","",classIDs)  ## keep everything before ____ (should match targets)
               
               if(identical(sampleIDs, character(0))){
                  stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
               }
            } ## NOT NULL
            
            
            ## MaxQuant TMT sample data columns of interest. X is the index 
            ## of the 'tag', starting with 0. ___XX indicates peptides with 
            ## single, double, triple phospho sites) 
            ## better to analyze ___1, ___2, ___3, than summarized site 
            ## intensities because in a biological system the same protein 
            ## may have distinct functions when it is differentially phosphorylated.
            ## https://adinasarapu.github.io/posts/2020/01/blog-post-tmt/
            ## proteinGroups.txt		   "Reporter intensity corrected X <experiment name(s)>"
            ## peptides.txt			   "Reporter intensity corrected X <experiment name(s)>"
            ## Phospho (STY)Sites.txt	"Reporter intensity corrected X <experiment name(s)>"
            ## Phospho (STY)Sites.txt	"Reporter intensity corrected X <experiment name(s)>___XX" 
            
            ## get column names for the corrected reporter intensities for each 
            ## phospho class 
            ## "Reporter.intensity.corrected.X.<experiment name(s)> _Phospho___1"
            if(is.null(sampleIDs)){
               
               pattern   <- "Reporter\\.intensity\\.corrected\\.[[:digit:]]+.*[P/p]hos*"
               classIDs  <- grep(pattern, colNames, value = TRUE)
               classNums <- gsub("^.*_{3}","",classIDs) ## keep everything after ____ (class number)
               sampleIDs <- gsub("_{3}.*","",classIDs)  ## keep everything before ____ (should match targets)
               
               ## if pattern doesn't match any column names then reduce pattern 
               ## to more generic format.
               ## "Reporter.intensity.corrected.X" 
               if(identical(classIDs, character(0))){
                  pattern   <- "Reporter\\.intensity\\.corrected\\.*[[:digit:]].*\\_{3}[[:digit:]]$"
                  classIDs  <- grep(pattern, colNames, value = TRUE);classIDs
                  classNums <- gsub("^.*_{3}","", classIDs) ## keep everything after ____
                  sampleIDs <- unique(gsub("_{3}.*","",classIDs)) ## keep everything before ____
                  cat("classIDs identified using generic pattern.")
               }
               
               
            } ## NULL
            
            
            data2 <- list(classIDs=classIDs, sampleIDs=sampleIDs, classNums=classNums, 
                          pattern=pattern)
            return(data2)
          
         } ## PHOSPHO
         
      } ## TMT/PHOSPHOTMT
      
      
   } ## EXTRACT SAMPLEIDS
   
   
   ##-----------------------------------------
   ##  [05] EXTRACT PROTEIN DATA (REQUIRED)
   ##-----------------------------------------
   extract_protein_data <- function(data, sampleIDs=NULL, pipe, enrich){
      
      ## extract protein data and annotation from DIA experiment
      if(all(pipe=="DIA" & enrich=="protein")){ ## DIA/PROTEIN
         
         ## get sampleIDs defined by user
         if(!is.null(sampleIDs)){ 
            IDs <- extract_sampleIDs(colNames=colnames(data),sampleIDs=sampleIDs, pipe=pipe, enrich=enrich)
            pattern   <- IDs$pattern
            sampleIDs <- IDs$sampleIDs
         }
         
         ## if sampleIDs are not defined by user then use all column names
         ## of the input data in the search.
         if(is.null(sampleIDs)){ 
            IDs <- extract_sampleIDs(colNames=colnames(data), sampleIDs=sampleIDs, pipe=pipe, enrich=enrich)
            pattern   <- IDs$pattern
            sampleIDs <- IDs$sampleIDs
         }
         
         ## extract gene name, gene symbol, and uniprot id info. from the Fasta.header
         ## append to qfilterData add unique protein ID (proKey = uniprot_GN_id 
         data$Accession.Number <- gsub("(.+?)(\\ .*)","\\1",stringr::str_extract(data$Protein.Name,"(?<=\\|)[^\\|]+(?=\\ )"))
         data$UniprotID    <- stringr::str_extract(data$Protein.Name, "(?<=\\|)[^\\|]+(?=\\|)")
         data$Gene_name    <- stringr::str_extract(data$Protein.Name, "(?<=GN\\=)[^\\|]+(?= PE\\=)")
         data$Description  <- stringr::str_extract(data$Protein.Name, "(?<= )[^\\|]+(?= OS\\=)")
         
         if(pipe=="DIA"){ annotColums <- diaAnnotationColums }
         annotColums <- c(annotColums, "UniprotID", "Gene_name","Description")
         rownames(data) <- paste(data$UniprotID, data$Gene_name, data$id, sep="_")
         
         ## extract annotation columns and protein intensity data 
         rawAnnot <- data[, annotColums]
         rawData  <- data[, sampleIDs]
         
         ## intensity data is not a numeric data.frame and requires additional processing
         ## numbers are character strings separated by commas (4,000,000), some cells
         ## contain text 'Missing Value' if the protein does not have detectable expression. 
         ## replace 'Missing Value' with zeros, remove comma's, convert data in each column
         ## to numeric values
         rawData[,][rawData[,]=="Missing Value"] <- 0
         replaceCommas <- function(x){ x<-as.numeric(gsub("\\,", "", x)) }
         for(i in 1:ncol(rawData)){ rawData[,i] <- replaceCommas(rawData[,i]) }
         
         no_samples   <- ncol(rawData)
         no_extracted <- nrow(rawData)
         
         print(paste("intensity data for", no_extracted, 
                     pipe, enrich, "entries and ",no_samples, "samples extracted.",
                     "Success!!"));cat("\n")
         
         data2 <- list(rawData=rawData, rawAnnot=rawAnnot, sampleIDs=sampleIDs,
                       annotColums=annotColums, no_samples=no_samples, 
                       no_extracted=no_extracted, pipe=pipe, enrich=enrich, pattern=pattern)
         return(data2)
         
         
      } ## DIA/PROTEIN            
      
      
      ## after filtering contaminants append additional protein annotation info.,
      ## extract protein data, parse into annotation and sample protein intensities.
      if(all((pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF") & enrich=="protein")){ ## TMT/PROTEIN
         
         ## get sampleIDs if defined by user.
         if(!is.null(sampleIDs)){ 
            IDs <- extract_sampleIDs(colNames=colnames(data),sampleIDs=sampleIDs, pipe=pipe, enrich=enrich)
            pattern   <- IDs$pattern
            sampleIDs <- IDs$sampleIDs
         }
         
         ## if sampleIDs are not defined by user then use all column names
         ## of the input data in the search.
         if(is.null(sampleIDs)){ 
            IDs <- extract_sampleIDs(colNames=colnames(data), sampleIDs=sampleIDs,pipe=pipe, enrich=enrich)
            pattern   <- IDs$pattern
            sampleIDs <- IDs$sampleIDs
         }
         
         if(pipe=="TMT" | pipe=="LF"){ annotColums <- tmtAnnotationColums }
         if(pipe=="phosphoTMT"){ annotColums <- proteinAnnotationColums }
         
         ## extract gene name, gene symbol, and uniprot id info. from the Fasta.header
         ## append to qfilterData add unique protein ID (proKey = uniprot_GN_id 
         data$UniprotID   <- stringr::str_extract(data$Fasta.headers, "(?<=\\|)[^\\|]+(?=\\|)")
         data$Gene_name   <- stringr::str_extract(data$Fasta.headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)")
         data$Description <- stringr::str_extract(data$Fasta.headers, "(?<= )[^\\|]+(?= OS\\=)")
         annotColums <- c(annotColums, "UniprotID", "Gene_name","Description")
         
         ## set row names as a combination of uniprotid_gene.name_id (protein key)
         rownames(data) <- paste(data$UniprotID, data$Gene_name, data$id, sep="_")
         
         ## extract annotation columns and corrected protein intensity data 
         ## intensity corrected columns) for the protein lysate samples.
         rawAnnot <- data[, annotColums]
         rawData  <- data[, sampleIDs]
         
         no_samples   <- ncol(rawData) 
         no_extracted <- nrow(rawData)
         
         print(paste("corrected reporter intensity data for", no_extracted, 
                     pipe, enrich, "entries and ",no_samples, "samples extracted.",
                     "Success!!"));cat("\n")
         
         data2 <- list(rawData=rawData, rawAnnot=rawAnnot, sampleIDs=sampleIDs,
                       annotColums=annotColums, no_samples=no_samples, 
                       no_extracted=no_extracted, pipe=pipe, enrich=enrich, pattern=pattern)
         return(data2)
         
      } ## TMT/PROTEIN
      
   } ## EXTRACT PROTEIN
   
   
   ##-----------------------------------------
   ##  [06] EXTRACT PHOSPHO DATA (REQUIRED)
   ##-----------------------------------------
   extract_phospho_data <- function(data, sampleIDs=NULL, pipe, enrich){
      
      if(pipe=="phosphoTMT" & enrich=="phospho"){
         
         annotColums <- phosphoAnnotationColums
         
         ## after filtering contaminants & low prob sites parse the Fasta.header 
         ## column to obtain gene name (Description), gene symbol (Gene_name), 
         ## UniprotID, and Flanking aa sequence (7aa-site-7aa), extract first id 
         ## from protein.group.IDs, first a.a. position within protein from 
         ## Positions.within.proteins, and create extract gene name,
         ## gene symbol, uniprot id info. from phospho(STY)sites Fasta_header and 
         ## Flanking a.a.sequence around phospho sites from phospho(Group)STY)sites
         ## Sequence_window and add to the qfilterData data.frame                            
         data$UniprotID   <- stringr::str_extract(data$Fasta.headers, "(?<=\\|)[^\\|]+(?=\\|)")
         data$Gene_name   <- stringr::str_extract(data$Fasta.headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)")
         data$Description <- stringr::str_extract(data$Fasta.headers, "(?<= )[^\\|]+(?= OS\\=)")
         data$Flanking    <- gsub("\\;.*$", "", data$Sequence.window) %>% str_sub(9,23) %>% paste0("-p")
         
         ## get protein group ID (id col. proteinGroups) and site position for 
         ## the first protein entry i.e. corresponds to Majority.protein.ID
         ## to get the first entry remove everything after the first semicolon.
         data$proGroupID     <- sub(";.*", "", data$Protein.group.IDs)
         data$phosPosition   <- sub(";.*", "", data$Positions.within.proteins) 
         data$phosAAPosition <- paste(data$Amino.acid, data$phosPosition,sep="")
         
         annotColums <- c(annotColums, "proGroupID", "UniprotID", "Gene_name", 
                          "Description", "Flanking", "phosPosition", "phosAAPosition");annotColums 
         rownames(data) <- paste(data$UniprotID, data$Gene_name, data$phosAAPosition,sep="_")
         
         ## GET SAMPLE IDS
         ## use user defined sampleIDs to extract the corrected intensity data for 
         ## single, double, and tripple phosphopeptides ___1, ___2, ____3.
         if(!is.null(sampleIDs)){ 
            table(sampleIDs %in% colnames(data))
            IDs <- extract_sampleIDs(colNames=sampleIDs, pipe=pipe, enrich=enrich)
         }
         
         ## use extract_sampleIDs() function to extract the corrected intensity data for 
         ## single, double, and tripple phosphopeptides ___1, ___2, ____3.
         if(is.null(sampleIDs)){ IDs <- extract_sampleIDs(colNames=colnames(data),
                                                          pipe=pipe, enrich=enrich) }
         classIDs  <- IDs$classIDs; classIDs
         classNums <- IDs$classNums;classNums
         pattern   <- IDs$pattern
         
         ## extract annotation and corrected intensity data 
         ## for single, double, tripple phosphopeptides
         rawAnnot <- data[, annotColums]; dim(rawAnnot)
         rawData  <- data[, classIDs]; dim(rawData)
         
         no_samples   <- ncol(rawData)
         no_extracted <- nrow(rawData)
         
         print(paste("intensity data for", no_extracted,pipe, enrich,
                     "entries and", no_samples, "samples extracted. Success!!"))
         cat("\n")
         
         ## EXTRACT CLASS DATA
         classList <- vector("list", max(IDs$classNums))
         names(classList) <- c(paste("class_",1:max(IDs$classNums),sep=""))
         classData <- NULL; classAnnot <- NULL
         
         for(i in 1:length(classList)){ ## INDIV CLASSES
            ## get class == i
            keep      <- IDs$classNums==i; table(keep)
            classIDs  <- IDs$classIDs[keep==TRUE]   ## sampleID___1
            sampleIDs <- IDs$sampleIDs[keep==TRUE]  ## sampleID
            classNums <- IDs$classNums[keep]        ## 1
            stopifnot(all(classNums==i))
            
            ## INDIVIDUAL CLASS DATA
            classData <- rawData[, classIDs]  ##__1
            stopifnot(colnames(classData)==classIDs)
            colnames(classData) <- sampleIDs  ## rename minus __1
            
            ## INDIVIDUAL CLASS ANNOTATION
            classAnnot <- rawAnnot[, annotColums]
            classAnnot$Class <- rep(i, nrow(classAnnot))  ## Class
            annotColums2 <- c(annotColums, "Class")
            stopifnot(rownames(classData)==rownames(classAnnot))
            
            no_class_samples <- ncol(classData) 
            no_class_phospho <- nrow(classData)
            
            print(paste("classData",i,":","intensity data for",no_class_phospho, 
                        "class", i, pipe,enrich, "entries and", no_class_samples,
                        "samples extracted. Success!!",sep=" "))
            
            ## id list from extract_sampleIDs() for class i
            idList <- list(classIDs=classIDs, sampleIDs=sampleIDs, classNums=classNums, 
                           annotColums=annotColums2, pattern=pattern)
            classList[[i]] <- list(rawData=classData, rawAnnot=classAnnot[rownames(classData), ],
                                   idList=idList, no_class_samples=no_class_samples, 
                                   no_class_phospho=no_class_phospho)
         }## INDIV CLASSES
         
         data2 <- list(classList=classList)
         return(data2)
         
      }
      
   } ## EXTRACT PHOSPHO
   
   
   ##---------------------------------
   ##  [07] EXTRACT DATA (REQUIRED)
   ##---------------------------------
   extract_data <- function(file=NULL, sampleIDs=NULL, min.prob=0.75,
                            pipe=c("DIA","TMT","phosphoTMT","LF"),
                            enrich=c("protein","phospho")){
      
      ## if no file is input by user the open file dialog box is 
      ## called to allow the user to navigate to and select the file.
      if(is.null(file)){ file=file.choose() }
      
      ## input parameters are added to param variable
      param <- stats <- list()
      param[["file"]]      <- file
      param[["pipe"]]      <- pipe
      param[["enrich"]]    <- enrich
      param[["sampleIDs"]] <- ifelse(is.null(sampleIDs),"NULL", paste(sampleIDs,collapse=", "))
      param[["min.prob"]]  <- ifelse(enrich=="protein", "NULL",min.prob)
      
      
      ## DIA PROTEIN
      if(pipe=="DIA" & enrich=="protein"){
         
         ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
         # if(!dir.exists("protein_analysis/01_quality_control")){ 
         #    dir.create("protein_analysis/01_quality_control", recursive=TRUE)
         #    print("quality control directory created...")
         # }
         
         ## IMPORT DATA
         ## the first 10 lines are skipped and last line removed. column 1 (.X)
         ## is renamed "id"
         DIADATA <- import_data(file=file, pipe=pipe, enrich=enrich)
         diaData <- DIADATA$data
         reqCols <- DIADATA$reqCols
         
         param[["file"]] <- DIADATA$file
         stats[["total_input_rows"]] <- DIADATA$no_input
         
         ## REMOVE IRRELEVANT COLUMNS 
         ## irrelevant columns, columns 2-3 (visible,star) are removed from 
         ## the sample report
         if(colnames(diaData)[1] != "id"){ colnames(diaData)[1] <-"id" }
         remove  <- colnames(diaData) %in% c("Visible", "Star");remove
         diaData <- diaData[, remove==FALSE]
         
         
         ## REPLACE MISSING VALUES WITH ZERO
         ## intensity data is not a numeric data.frame and requires 
         ## additional processing numbers are character strings separated
         ## by commas (4,000,000), some cells contain text 'Missing Value' 
         ## if the protein does not have detectable expression. 
         ## replace 'Missing Value' with zeros, remove comma's, convert 
         ## data in each column to numeric values
         diaData[,][diaData[,]=="Missing Value"] <- "0"
         replaceCommas  <- function(x){ x <- as.numeric(gsub("\\,", "", x)) }
         replaceNumbers <- colnames(diaData) %in% reqCols
         replaceColums  <- colnames(diaData)[replaceNumbers==FALSE];replaceColums
         for(i in replaceColums){ diaData[,i] <- replaceCommas(diaData[,i]) }
         
         
         ## SAVE BIG QUERY INPUT FILE
         ## save a copy of samples report for upload into Big Query. 
         ## This file cannot contain NA or blank values, column names 
         ## cannot begin with a number, replace all NA and blank cells
         ## with zeros. If column names begin with a number add X to 
         ## beginning of column name. change file name to ilab_Samples_Report_BQ.csv
         bqData <- data.frame(diaData)
         bqData[,][is.na(bqData[,])] <- 0
         bqData[,][bqData[,]==""]    <- 0
         
         ## if column names begin with a number append X to the start of each name.
         testColumNumber <-substr(colnames(bqData),1,1);testColumNumber
         if(length(grep("[[:digit:]]",testColumNumber)) > 0){
            colnames(bqData) <- paste0("X",colnames(bqData))
         }
         bn   <- basename(file);bn
         bn   <- gsub("Samples Report of ","",bn);bn
         ilab <- gsub(paste0(".",tools::file_ext(bn)),"",bn);ilab
         ## save DIA Big Query Input File
         filename <- paste0(ilab,"_Samples_Report_BQ.csv")
         utils::write.csv(bqData, file.path(".",filename), row.names=FALSE)
         print(paste("Big Querry Samples Report (",filename,
                     ") saved to project directory. Success!!"))
         param[["ilab"]] <- ilab
         
         
         ## REMOVE CONTAMINANTS
         qfilterData <- remove_contaminants(data=diaData, pipe=pipe, enrich=enrich)
         stats[["no_contam_removed"]] <- qfilterData$no_contam
         
         ## EXTRACT PROTEIN DATA
         extProt <- extract_protein_data(data=qfilterData$data, sampleIDs=sampleIDs,
                                         pipe=pipe, enrich=enrich)
         
         stats[["total_input_samples"]] <- extProt$no_samples
         stats[["no_extracted_rows"]]   <- extProt$no_extracted
         
         ## SAVE PARAM/STATS TO LOG FILE
         logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)
         
         ## return a list of data.frame containing the extracted quality 
         ## filtered intensity data, corresponding extracted annotation, 
         ## stats, and input parameters
         data2 <- list(data=extProt$rawData, annot=extProt$rawAnnot, 
                       param=logs$param, stats=logs$stats[c(3,1,2,4),], ilab=ilab)
         return(data2)
         
         
      } ## DIA PROTEIN
      
      
      ## TMT PROTEIN
      if((pipe=="TMT" | pipe=="LF") & enrich=="protein"){
         
         ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
         # if(!dir.exists("protein_analysis/01_quality_control")){ 
         #    dir.create("protein_analysis/01_quality_control", recursive=TRUE) 
         #    print("quality control directory created...")
         # }
         
         ## IMPORT DATA
         maxQuant <- import_data(file=file, pipe=pipe, enrich=enrich)
         stats[["total_input_rows"]] <- maxQuant$no_input
         
         ## REMOVE CONTAMINANTS
         qfilterData <- remove_contaminants(data=maxQuant$data, pipe=pipe, enrich=enrich)
         stats[["no_contam_removed"]] <- qfilterData$no_contam
         
         
         ## EXTRACT PROTEIN DATA AND ANNOTATION
         extProt <- extract_protein_data(data=qfilterData$data, sampleIDs=sampleIDs, 
                                         pipe=pipe, enrich=enrich)
         
         stats[["total_input_samples"]] <- extProt$no_samples
         stats[["no_extracted_rows"]]   <- extProt$no_extracted
         param[["pattern"]] <- extProt$pattern
         
         ## SAVE PARAM/STATS TO LOG FILE
         logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)
         
         ## return a list of data.frame containing the extracted quality 
         ## filtered intensity data, corresponding extracted annotation, 
         ## stats, and input parameters
         data2 <- list(data=extProt$rawData, annot=extProt$rawAnnot, 
                       param=logs$param,stats=logs$stats)
         return(data2)
         
      } ## TMT PROTEIN
      
      
      ## phosphoTMT PROTEIN
      if(pipe=="phosphoTMT" & enrich=="protein"){
         
         ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
         # if(!dir.exists("protein_analysis/01_quality_control")){ 
         #    dir.create("protein_analysis/01_quality_control", recursive=TRUE) 
         #    print("quality control directory created...")
         # }
         
         ## IMPORT DATA
         maxQuant <- import_data(file=file, pipe=pipe, enrich=enrich)
         stats[["total_input_rows"]] <- maxQuant$no_input
         
         
         ## REMOVE CONTAMINANTS
         qfilterData <- remove_contaminants(data=maxQuant$data, pipe=pipe, enrich=enrich)
         stats[["no_contam_removed"]] <- qfilterData$no_contam
         
         
         ## EXTRACT PROTEIN DATA
         extProt <- extract_protein_data(data=qfilterData$data, sampleIDs=sampleIDs, 
                                         pipe=pipe, enrich=enrich)
         
         stats[["total_input_samples"]] <- extProt$no_samples
         stats[["no_extracted_rows"]] <- extProt$no_extracted
         param[["pattern"]] <- extProt$pattern
         
         ## SAVE PARAM/STATS TO LOG FILE
         logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)
         
         
         ## return a list of data.frame containing the extracted quality 
         ## filtered intensity data, corresponding extracted annotation, 
         ## stats, and input parameters
         data2 <- list(data=extProt$rawData, annot=extProt$rawAnnot, 
                       param=logs$param, stats=logs$stats[c(3,1,2,4),])
         return(data2)
         
      } ## phosphoTMT PROTEIN
      
      
      ## phosphoTMT PHOSPHO
      if(pipe=="phosphoTMT" & enrich=="phospho"){
         
         ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
         # if(!dir.exists("phospho_analysis/01_quality_control")){ 
         #    dir.create("phospho_analysis/01_quality_control", recursive=TRUE) 
         #    print("quality control directory created...")
         # }
         
         ## IMPORT DATA
         maxQuant <- import_data(file=file, pipe=pipe, enrich=enrich)
         stats[["total_input_rows"]] <- maxQuant$no_input
         
         ## REMOVE CONTAMINANTS
         qfilterData <- remove_contaminants(data=maxQuant$data, pipe=pipe, enrich=enrich)
         stats[["no_contam_removed"]] <- qfilterData$no_contam
         
         ## LOCAL PROBABILITY FILTER
         qfilterData <- local_prob_filter(data=qfilterData$data, min.prob=min.prob,
                                          pipe=pipe, enrich=enrich)
         param[["min.prob"]] <- qfilterData$min.prob
         stats[["no_localprob_removed"]] <- qfilterData$no_localprob_removed
         
         ## EXTRACT PHOSPHO DATA
         extPhos <- extract_phospho_data(data=qfilterData$data, sampleIDs=sampleIDs, 
                                         pipe=pipe, enrich=enrich)
         stats[["total_input_samples"]] <- extPhos$classList$class_1$no_class_samples
         stats[["no_extracted_rows"]] <- extPhos$classList$class_1$no_class_phospho
         param[["pattern"]] <- extPhos$classList$class_1$idList$pattern
         
         ## SAVE PARAM/STATS TO LOG FILE
         logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)
         
         
         ## return a list of data.frame containing the extracted quality 
         ## filtered intensity data for class 1 phosphopeptides, 
         ## corresponding extracted annotation, stats, and input parameters
         ## classList corresponds to data/annotation for each phosphopeptide
         ## class e.g. class___1, class___2, class___3 
         
         data2 <- list(data=extPhos$classList$class_1$rawData,   ## class 1 data
                       annot=extPhos$classList$class_1$rawAnnot, ## class 1 annot
                       classList=extPhos$classList, ## class_1$data/annot, class_2$data/annot etc.
                       sampleIDs=extPhos$classList$class_1$idList$sampleIDs, 
                       pattern=extPhos$classList$class_1$idList$pattern,
                       param=logs$param,stats=logs$stats[c(4,1,2,3,5), ])
         return(data2)
         
         
      } ## phosphoTMT PHOSPHO
      
      
   } ## EXTRACT DATA
   
   
   
   ##------------------------------------
   ##  [08] IMPORT METADATA (OPTIONAL)
   ##------------------------------------
   import_meta <-function(file, sampleIDs, pipe){
      
      ## check that file exists in specified location. 
      ## if file exits and is csv/tsv/txt
      file <-gsub("\\./","",file)
      file.dir <- dirname(file);file.dir
      file <- file.path(file.dir, match.arg(basename(file), list.files(file.dir),
                                            several.ok=FALSE))
      filext <- tools::file_ext(file);filext
      filext <- match.arg(filext, c("csv","txt","tsv"),several.ok=FALSE)
      if(filext=="txt" | filext=="tsv"){ sep="\t" }
      if(filext=="csv"){ sep=","}
      
      ## the metadata file should contain the following required columns
      if(pipe=="phosphoTMT" | pipe=="TMT"){ reqCols<- tmtMetaColums }
      if(pipe=="DIA"){ reqCols <- diaMetaColums }
      if(pipe=="LF"){ reqCols <- lfMetaColums }
      
      ## import file
      data <- utils::read.csv(file=file, sep=sep, stringsAsFactors=FALSE, 
                              header=TRUE, check.names=TRUE)
      
      if(all(reqCols %in% colnames(data))==TRUE){
         print("metadata file imported. Success!!")
         return(data)
      }
      
      if(all(reqCols %in% colnames(data))==FALSE){
         stop(message(paste("Error! For",pipe,"experiments the sample metadata file",
                            "should include the following required columns:",
                            paste(reqCols,collapse=", "), ". The following columns are",
                            "missing from the metadata file: ",
                            paste(reqCols[reqCols %in% colnames(data)==FALSE], 
                                  collapse=", "))))
      }
      
   } ## IMPORT META DATA
   
   
   
   ##-------------------------------------
   ##  [09] GET BATCHES (REQUIRED) (NEW)
   ##-------------------------------------
   get_tmt_batches <- function(sampleIDs, pipe, enrich){
      
      ## remove Reporter.intensity.corrected.XX text from sampleIDs
      pattern <- gsub("Reporter\\.intensity\\.corrected\\.[[:digit:]]+","",sampleIDs)
      
      ## determine if remaining experiment info. text contains a number
      ## if a number is present we assume this refers to different batches
      hasNumber <- grep("[[:digit:]]",pattern);hasNumber
      
      ## if experiment info. part of sampleID contains a number
      ## then isolate number by removing periods, underscores, and letters
      ## from the text
      if(all(length(hasNumber)==length(sampleIDs))==TRUE){
         
         b <- gsub("\\.","",pattern);b   ## removes period
         b <- gsub("_","",b);b           ## removes underscore
         b <- gsub("[A-Z]","",b);b       ## removes capital letters
         b <- gsub("[a-z]","",b);b       ## removes lowercase letters
         
         ## check that b is now a number and free of other
         ## characters and symbols by checking against a vector of 1:100
         if(all(b %in% 1:100)){ b <- as.numeric(b);b } else { 
            stop(paste("b is either not a number or is a number > 100. :(",
                       "unique batch values include:",paste(unique(b),collapse=", "))) }
      }
      
      ## if experiment text does not contain a number then use unique
      ## values of the text to assign batches assuming that the first
      ## unique value corresponds to experiment 1
      if(all(length(hasNumber)==length(sampleIDs))==FALSE){
         b=pattern;b
         for(i in 1:length(unique(b))){ b[b %in% unique(b)[i]] <- i };b
         if(all(b %in% 1:100)){ b <- as.numeric(b);b } else {
            stop(paste("b is either not a number or is a number > 100. :(",
                       "unique batch values include:",paste(unique(b),collapse=", "))) }
      }
      cat("\n");print(b);cat("\n");print(unique(b))
      
      return(b)
      
      
   } ## GET TMT BATCHES
   
   
   ##------------------------------------------
   ##  [10] GET DIA SAMPLE NUMBER (REQUIRED)
   ##------------------------------------------
   ## sample number extracted from sampleIDs is matched 
   ## to sample number in meta file. pool samples in metadata 
   ## number column should be labeled P1, P2, P3 samples in metadata 
   ## number column should be labeled 1,2,3,4, etc.
   get_dia_sample_number <- function(sampleIDs){
      
      ## scaffold DIA sample reports should contain 
      ## Pool_1/Sample_1 as part of the sampleIDs
      ## e.g. "Balachandran_061721_Pool_3.mzML", 
      ## e.g. "Balachandran_061721_Sample_01_CS1.mzML"
      pattern1 <- "[P/p]ool_[[:digit:]]"   ## POOL PATTERN (Pool_1)
      pattern2 <- "[S/s]ample_[[:digit:]]" ## SAMPLE PATTERN (Sample_1)
      pat1     <- grep(pattern1,sampleIDs)
      pat2     <- grep(pattern2,sampleIDs)
      pattern  <- c(pat1, pat2);pattern
      
      if(all(length(pattern)==length(sampleIDs))==TRUE){
         
         ## pool samples
         number1 <- gsub("\\.","_", sampleIDs[pat1]);number1
         number1 <- gsub(".*[P/p]ool_","",number1)  ## removes everything b4 Pool_
         number1 <- gsub("_.*","",number1);number1  ## removes everything after number _mzML
         number1 <- paste0("P",as.numeric(number1)) ## adds capital P to number == P1, P2, P3
         
         ## samples
         number2 <- gsub("\\.","_",sampleIDs[pat2])
         number2 <- gsub(".*[S/s]ample_","",number2) ## removes everything b4 Sample_
         number2 <- gsub("_.*","",number2);number2   ## removes everthing after number _mzML
         number2 <- as.numeric(number2);number2      ## if 01,02 == 1,2
         number  <- c(number1, number2);number       ## P1,P2,P3, 1,2,3,..., n
         
         print("sample numbers extracted from DIA sampleIDs. Success!!")
         return(number)
      }
      
      if(all(length(pattern)!=length(sampleIDs))==FALSE){
         stop(message(paste("Error! Sample numbers could not be extracted", 
                            "from the input sampleIDs. DIA sampleIDs should",
                            "contain the following patterns:", pattern1,
                            "(pool samples e.g. P/pool_1/P/pool_01) & ",
                            pattern2,"(samples e.g. Sample_1/Sample_01 (not S1)).")))
      } 
      
   } ## GET SAMPLE NUMBER
   
   
   ##---------------------------------
   ##  [11] MAKE TARGETS (OPTIONAL)
   ##---------------------------------
   make_targets <- function(file=NULL, sampleIDs, pipe=c("DIA","TMT","phosphoTMT"),
                            enrich=c("protein","phospho")){
      
      ## METADATA FILE INPUT BY USER
      if(!is.null(file)){
         
         ## DIA METADATA FILE
         ## sample number extracted from sampleIDs is matched to sample number 
         ## in meta file. pool samples in metadata number column should be labeled
         ## P1, P2, P3. samples in metadata number column should be labeled 1,2,3,4,
         ## if combining the two is successful a file named targets.csv is saved to
         ## the project directory.
         if(pipe=="DIA"){ ## DIA
            
            ## import metadata file. use sample number as metadata key
            metadata <- import_meta(file=file, sampleIDs=sampleIDs, pipe=pipe)
            
            ## extract sample number from sampleIDs and use as data key
            number <- get_dia_sample_number(sampleIDs=sampleIDs);number
            
            ## if the extracted sample numbers and metadata sample numbers match.
            ## create basic targets info. and combine with metadata. save targets
            ## to csv file and return combined metadata data.frame.
            if(all(number %in% metadata$number)){
               ## create basic targets info.
               targets <- data.frame(sampleIDs=sampleIDs, number=number, 
                                     pipe=rep(pipe,length(sampleIDs)),
                                     enrichment=rep(enrich,length(sampleIDs)))
               rownames(targets) <- targets$sampleIDs;head(targets)
               
               ## use sample number info. to sort targets info. so that it is in the same order
               ## as the meta data file. combine targets and metadata into a single data.frame
               m <- match(metadata$number, number);m
               stopifnot(number[m]==metadata$number)
               targets <- targets[m, ];head(targets)
               remove  <- colnames(metadata)%in%colnames(targets);remove
               targets <- cbind(targets, metadata[,remove==FALSE]);head(targets)
               
               ## save targets as a csv file in the main project directory.
               filename <- "targets.csv"
               utils::write.csv(targets, file=filename, row.names=FALSE)
               print(paste(pipe, "targets file: ", filename, "saved to project directory. Success!!"))
               print("metadata and targets info. combined. Success!!")
               return(targets)
               
            } ## COMBINED DIA
            
            ## sample numbers do not match. return basic extracted targets info. and 
            ## imported metadata as separate data.frames. integrate manually.
            ## NOTE: make sure the sample numbers for the pool samples in the imported
            ## metadata file are P1, P2, P3 and samples are 1,2,3,4 ...
            if(all(number %in% metadata$number)==FALSE){
               message("Error! The extracted sample numbers do not match the metadata file.\n", 
                       "The following sample numbers are not present in the metadata file:", 
                       paste(number[number %in% metadata$number==F],collapse=", "),"\n")
               message("Basic targets info. (e.g. sampleIDs/number) extracted from the input\n",
                       " sampleIDs and the imported metadata could not be combined.\n",
                       "Returning metadata and basic targets info. as separate data.frames.\n",
                       "Integration failed. Combine manually :(")
               
               ## create basic targets info.
               targets <- data.frame(sampleIDs=sampleIDs, number=number, 
                                     pipe=rep(pipe,length(sampleIDs)),
                                     enrichment=rep(enrich,length(sampleIDs)))
               rownames(targets) <- targets$sampleIDs;head(targets)
               
               data2 <- list(targets=targets, metadata=metadata)
               return(data2)
               
            } ## NOT COMBINED DIA
            
            
            # } ## PATTERN 1
            
         } ## DIA
         
         
         ## TMT/phosphoTMT METADATA FILE
         ## batch number and channel indexes are extracted from sampleIDs 
         ## this combined key is used to match the combined batch/tag/channel info. 
         ## in the metadata file.remember the batch column refers to TMT batch number
         ## in the metadata file.
         if(pipe=="TMT" | pipe=="phosphoTMT"){ ## TMT/phosphoTMT
            
            ## import meta data file. set row names as batch.channel key
            metadata <- import_meta(file=file, sampleIDs=sampleIDs, pipe=pipe);metadata
            rownames(metadata) <- paste(metadata$batch, metadata$tag,sep=".");metadata
            
            ## list of TMT tags and TMT batches in metadata file
            unq_tags    <- unique(metadata$tag);unq_tags   ## list of unique TMT tags
            unq_bats    <- unique(metadata$batch);unq_bats ## list of unique TMT batches
            no_unq_tags <- length(unq_tags);no_unq_tags    ## no of unique TMT tags
            no_unq_bats <- length(unq_bats);no_unq_bats    ## no of unique batches
            ## stop if each batch does not contain same no samples
            stopifnot(length(sampleIDs)==no_unq_tags * no_unq_bats) 
            
            ## use no. unique tags to determine type of TMT kit used in experiment.
            if(all(no_unq_tags==6  & all(unq_tags %in% TMT6plex))){ plex  <- TMT6plex }
            if(all(no_unq_tags==10 & all(unq_tags %in% TMT10plex))){ plex <- TMT10plex }
            if(all(no_unq_tags==11 & all(unq_tags %in% TMT11plex))){ plex <- TMT11plex }
            if(all(no_unq_tags==16 & all(unq_tags %in% TMT16plex))){ plex <- TMT16plex }
            if(no_unq_tags %in% c(6,10,11,16)==FALSE){
               return(message("Error! TMTplex could not be determined."))
            }; plex
            
            ## use TMT tags to add channel index number to metadata
            ## create named vector 
            ## combine each TMT batch with the known TMT tag list 
            ## e.g. (1.126, 1.130C, 2.126,2.130C) like metadata rownames)
            ## add channel column to metadata 
            metadata$channel <- rep(NA,nrow(metadata))
            for(i in base::seq_along(plex)){
               m <- metadata$tag %in% plex[i]
               metadata$channel[m==TRUE] <-names(plex)[i]
            }
            if(any(is.na(metadata$channel))==TRUE){
               stop(paste("Error! adding channel index to metadata failed.",
                          "NAs in channel column remain:",
                          paste(metadata$channel,collapse=", ")))
            }
            
            
            
            ##--------------------------
            ##    PROCESS SAMPLEIDS
            ##--------------------------
            print("processing sampleIDs ...")
            
            ## extract channel index number from sampleIDs
            channel <- gsub("Reporter\\.intensity\\.corrected\\.","",sampleIDs); channel
            channel <- as.numeric(sub("\\..*","",channel));channel
            
            
            ## GET BATCH NUMBERS
            b <- get_tmt_batches(sampleIDs=sampleIDs, pipe=pipe, enrich=enrich);b
            print(paste0("batch info. for ", length(unique(b)), " TMT batches created. Success!!"))
            
            ## unique TMT tags and TMT batches extracted from sampleIDs
            unq_idx     <- unique(channel);unq_idx    ## list of unique TMT tags
            unq_b       <- unique(b);unq_b            ## list of unique TMT batches
            no_unq_idx  <- length(unq_idx);no_unq_idx ## no of unique TMT tags
            no_unq_b    <- length(unq_b);no_unq_b     ## no of unique batches
            ## stop if each batch does not contain same no samples
            stopifnot(length(sampleIDs)==no_unq_idx * no_unq_b) 
            
            ## use no. unique tags to determine type of TMT kit used in experiment.
            if(all(no_unq_idx==6  & all(unq_idx %in% names(TMT6plex)))){  plex2  <- TMT6plex }
            if(all(no_unq_idx==10 & all(unq_idx %in% names(TMT10plex)))){ plex2 <- TMT10plex }
            if(all(no_unq_idx==11 & all(unq_idx %in% names(TMT11plex)))){ plex2 <- TMT11plex }
            if(all(no_unq_idx==16 & all(unq_idx %in% names(TMT16plex)))){ plex2 <- TMT16plex }
            if(no_unq_idx %in% c(6,10,11,16)==FALSE){
               return(stop(paste("Error! TMTplex tags could not be determined from channel",
                                 "index extracted from sampleIDs. :("))) }
            plex2
            
            ## channel indexes from sampleIDs used to obtain TMT tag info. 
            tag <- plex2[names(plex2)[channel]];tag
            
            ## basic targets info. created setting rownames as batch.tag (1.131C)
            ## https://stackoverflow.com/questions/23534066/cbind-warnings-row-names-were-found-from-a-short-variable-and-have-been-discar
            targets <- data.frame(sampleIDs=sampleIDs, batch=b, tag=tag, channel=channel, 
                                  pipe=rep(pipe,length(sampleIDs)),
                                  enrichment=rep(enrich, length(sampleIDs)), row.names=NULL)
            rownames(targets) <- paste(targets$batch,targets$tag,sep=".");head(targets)
            
            
            ## TMT COMBINED
            if(all(rownames(targets) %in% rownames(metadata))){
               ## merge metadata with targets information removing duplicate columns in targets
               targets <- targets[rownames(metadata), ]
               remove  <- colnames(targets) %in% colnames(metadata);table(remove)
               targets2<-as.data.frame(targets[,remove==FALSE])
               colnames(targets2)<-colnames(targets)[remove==FALSE]; head(targets2)
               targets <- cbind(metadata,targets2);head(targets)
               stopifnot(targets$channel==metadata$channel)
               rownames(targets) <- targets$sampleIDs;head(targets)
               
               ## save targets file to project directory
               if(pipe=="phosphoTMT" & enrich=="protein"){ filename <- "targets.pro.csv" }
               if(pipe=="phosphoTMT" & enrich=="phospho"){ filename <- "targets.phos.csv" }
               if(pipe=="TMT" & enrich=="protein"){ filename <- "targets.csv" }
               utils::write.csv(targets, file=filename, row.names=FALSE)
               print(paste("targets file: ", filename, "saved to project directory. Success!!"))
               
               ## return the combined targets info. and imported metadata
               cat("\n");print("targets metadata created. Success!!")
               return(targets)
               
            } ## TMT COMBINED
            
            
            ## TMT NOT COMBINED
            ## batch.channel index info. do not match. return basic extracted targets info.
            ## and imported metadata as separate data.frames. integrate manually.
            if(all(rownames(targets) %in% rownames(metadata))==FALSE){
               message("Error! The extracted batch.channel info. does not match the metadata file.\n", 
                       "The following batch.channel info. does not match the metadata file:\n", 
                       paste(rownames(targets)[rownames(targets) %in% rownames(metadata)==F],collapse="\n "),"\n")
               message("Basic targets info. (e.g. sampleIDs/batch/tag/channel)\n",
                       "and the imported metadata could not be combined.\n",
                       "Returning metadata and basic targets info. as separate data.frames.\n",
                       "Integration failed. Combine the two data.frames manually. :(")
               
               data2 <- list(targets=targets, metadata=metadata)
               return(data2)
               
               
            } ## TMT NOT COMBINED
            
         } ## TMT/PHOSPHOTMT
         
         
      } ## NOT NULL
      
      
      ## NO METADATA FILE
      if(is.null(file)){
         
         ## DIA (NO META DATA FILE)
         if(pipe=="DIA"){ ## DIA NULL
            ## extract sample number from sampleIDs and use as data key
            number  <- get_dia_sample_number(sampleIDs=sampleIDs);number
            targets <-data.frame(sampleIDs=sampleIDs, number=number, 
                                 pipe=rep(pipe,length(sampleIDs)),
                                 enrichment=rep(enrich,length(sampleIDs)), row.names=NULL)
            rownames(targets)<-targets$sampleIDs;head(targets)
            print(paste(enrich, "targets template for", pipe, "experiment created. Success!!",
                        "Import metadata and combine."))
            
            return(targets)
            
         } ## DIA NULL
         
         
         ## TMT/phosphoTMT (NO METADATA FILE)
         if(pipe=="TMT" | pipe=="phosphoTMT"){ ## TMT/PHOSPHO NULL
            
            ## extract channel index number from sampleIDs
            channel <- gsub("Reporter\\.intensity\\.corrected\\.","",sampleIDs); channel
            channel <- sub("\\..*","",channel);channel
            channel<-as.numeric(channel);channel
            ## extract batch number from sampleIDs
            b <- get_tmt_batches(sampleIDs=sampleIDs, pipe=pipe, enrich=enrich);b
            
            ## list of TMT tags and TMT batches
            unq_idx     <- unique(channel);unq_idx   ## list of unique TMT tags
            unq_b       <- unique(b);unq_b           ## list of unique TMT batches
            no_unq_idx  <- length(unq_idx);no_unq_idx## no of unique TMT tags
            no_unq_b    <- length(unq_b);no_unq_b    ## no of unique batches
            ## stop if each batch does not contain same no samples
            stopifnot(length(sampleIDs)==no_unq_idx * no_unq_b) 
            
            ## use no. unique tags to determine type of TMT kit used in experiment.
            if(all(no_unq_idx==6  & all(unq_idx %in% names(TMT6plex)))){ plex2  <- TMT6plex }
            if(all(no_unq_idx==10 & all(unq_idx %in% names(TMT10plex)))){ plex2 <- TMT10plex }
            if(all(no_unq_idx==11 & all(unq_idx %in% names(TMT11plex)))){ plex2 <- TMT11plex }
            if(all(no_unq_idx==16 & all(unq_idx %in% names(TMT16plex)))){ plex2 <- TMT16plex }
            plex2
            ## channel indexes from sampleIDs used to obtain TMT tag info. 
            if(!is.null(plex2)){ tag <- plex2[names(plex2)[channel]];tag }
            if(is.null(plex2)){
               tag=rep("", length(sampleIDs))
               message(paste("Warning! TMTplex tags could not be determined from channel",
                             "index values extracted from sampleIDs. :(")) 
            }
            
            ## create tagets basic info. data.frame.
            ## https://stackoverflow.com/questions/23534066/cbind-warnings-row-names-were-found-from-a-short-variable-and-have-been-discar
            targets <- data.frame(sampleIDs=sampleIDs, batch=b, tag=tag, channel=channel, 
                                  pipe=rep(pipe,length(sampleIDs)),
                                  enrichment=rep(enrich, length(sampleIDs)), row.names=NULL);head(targets)
            rownames(targets) <- targets$sampleIDs;targets;dim(targets)
            print(paste(enrich, "basic targets info. for", pipe, "experiment created. Success!!",
                        "Import metadata and combine the two data.frames."))
            
            return(targets)
            
         } ## TMT/PHOSPHO NULL
         
         
      } ## NULL
      
      
   } ## MAKE TARGETS
   
   
   ##---------------------------------
   ##  [12] GET TARGET COLUMNS (OPTIONAL)
   ##---------------------------------
   #x="sex", targets=sub$targets, required=FALSE
   #x="group", targets=sub$targets, required=TRUE
   get_target_cols <- function(targets, x, required=TRUE){
      
      ## get sample id column
      if(required==TRUE){
         if(any(is.na(x) || is.null(x) || x=="" || x==" ")){
            stop(paste("x is a required meta data column and is either NULL, NA, or blank.", 
                       " Please enter the column name or column number \nin the targets.file",
                       " identifying the sample.ids/experimental groups. ", 
                       "\n\nAvailable column names include: \n\n  ", 
                       paste(colnames(targets),collapse="'\n '") ,sep=""));cat("\n")
         }}
      
      ## x is not defined and is defined as NULL
      if(any(is.na(x) || is.null(x) || x=="" || x==" ")){ x=NULL }
      
      x2=NULL;x3<-NA
      ## x is column name
      len<-ifelse(required==TRUE, 1, length(x));len
      several.ok=ifelse(required==TRUE, FALSE,TRUE);several.ok
      if(all(any(is.character(x)) & length(x)==len & any(x%in%colnames(targets)))){
         x2 <- match.arg(x, colnames(targets), several.ok=several.ok);x2
         k<-x%in%x2;k
         if(FALSE %in% k==TRUE){ x3<-as.numeric(x[k==FALSE])};x3 } 
      
      if(all(is.null(x2) & all(is.na(x3)))){ x3<-as.numeric(x) 
      } else { 
         if(is.null(x2) & any(!is.na(x3))){ x3<-as.numeric(x3) }}
      
      ## if x is column number
      len<-ifelse(required==TRUE, 1, length(x3));len
      if(len==0){x3<-NA};x3
      several.ok=ifelse(required==TRUE, FALSE,TRUE);several.ok
      if(any(!is.na(x3))){
         if(all(any(is.numeric(x3) & ttutils::isInteger(x3)) & length(x3)==len)){
            if(any(x3 %in% c(1:ncol(targets)))){
               keep<-x3 %in% c(1:ncol(targets))
               x3 <-x3[keep==TRUE];x3
               x3 <-colnames(targets)[x3]
               stopifnot(x3 %in% colnames(targets))
            } else {
               if(all(x3 %in% c(1:ncol(targets)))==FALSE){ x3 <- NULL  }}
         }} else {
            if(is.na(x3)){ x3<-NULL }}
      
      x4=c(x2,x3);x4
      
      if(required==TRUE & is.null(x4)){
         stop(paste("\n\n x is a required meta data column and is NULL i.e. is either not a valid column name 
                 \n or column number or is a vector with multiple values. Please enter one column name or 
                 \n one column number (multiple values are not allowed) in the targets.file identifying 
                 \n the sample.ids/experimental groups/paired samples etc. \n\n Available column names include:
                 \n\n  ", paste(colnames(targets),collapse="'\n '") ,sep=""));cat("\n\n")
      }
      
      return(x4)
      
      
   } ## GET TARGET COLUMNS
   
   
   
   ##--------------------------------------------------
   ##  [12] MAKE TARGETS 2 (FROM RNASEQ) (OPTIONAL)
   ##--------------------------------------------------
   ## MAKE TARGETS FILE
   # targets=sub$targets, sample.ids="sampleIDs",samples="sample", group="subtype",factors=c("group","sex")
   make_targets2 <- function(targets, sample.ids, samples=NULL, group, factors=NULL){
      
      sampleIDs <- get_target_cols(targets, x=sample.ids,required=TRUE)
      samples    <- get_target_cols(targets, x=samples, required=FALSE)
      group      <- get_target_cols(targets, x=group, required=TRUE)
      factors    <- get_target_cols(targets, x=factors, required=FALSE)
      
      my_targets<-NULL
      
      ## DEFINE SAMPLE IDS
      my_targets <- data.frame(sampleIDs=as.character(targets[,sample.ids]),row.names=targets[,sample.ids])
      # print(head(my_targets));cat("\n\n")
      
      ## DEFINE SAMPLE NAMES
      if(!is.null(samples)){
         my_targets<-data.frame(cbind(my_targets, sample=targets[,samples]), stringsAsFactors=FALSE)
      }
      if(is.null(samples)){
         my_targets<-data.frame(cbind(my_targets, sample=targets[,sample.ids]), stringsAsFactors=FALSE) }
      # print(head(my_targets));cat("\n\n")
      
      ## [5.1] DEFINE GROUP OF INTEREST
      my_targets <- data.frame(cbind(my_targets, group=targets[,group]), stringsAsFactors=FALSE)
      # stopifnot(my_targets$group==targets[,group])
      # print(head(my_targets));cat("\n\n")
      
      my_targets$group <- factor(my_targets$group, levels=ordered(unique(as.character(my_targets$group))))
      # cat("\n\n");print(head(my_targets));cat("\n\n")
      # print(table(my_targets$group));cat("\n")
      # print(levels(my_targets$group));cat("\n\n")
      
      
      ## [5.2] DEFINE FACTORS OF INTEREST
      if(!is.null(factors)){
         for(i in 1:length(factors)){
            my_targets <- data.frame(cbind(my_targets, factor=targets[,factors[i]]), stringsAsFactors=FALSE)
            stopifnot(my_targets$factor==targets[,factors[i]])
            
            my_targets$factor <- factor(my_targets$factor, levels=ordered(unique(as.character(my_targets$factor))))
            # print(paste("factors: ", factors[i]))
            # print(head(my_targets$factor));cat("\n")
            colnames(my_targets)[ncol(my_targets)] <- factors[i]
            
         }
      }
      
      colnames(my_targets)<-make.unique(colnames(my_targets))
      # print(head(my_targets));cat("\n\n")
      # print(table(my_targets[,ncol(my_targets)]));cat("\n\n")
      # print(levels(my_targets[,ncol(my_targets)]));cat("\n\n")
      
      return(my_targets)
      
   } ## MAKE TARGETS
   
   
   
   ##------------------------
   ##  [11] MAKE LOG FILES
   ##------------------------
   make_log <- function(param, stats, title="", save=TRUE){
      
      param<-data.frame(t(as.data.frame(t(param))))
      colnames(param)<-"";param
      
      stats<-data.frame(t(as.data.frame(t(stats))))
      colnames(stats)<-"";stats
      
      if(save==TRUE){
         
         if(!dir.exists("logs")){dir.create("logs",recursive=TRUE)}
         
         if(all(length(param) > 0 & nrow(param)>0)){
            
            if(!file.exists(file.path("logs","parameters.log"))){ 
               sink(file="./logs/parameters.log",append=FALSE)
            } else { sink(file="./logs/parameters.log",append=TRUE) }
            
            # sink(file="./param.txt",append=TRUE)
            cat(paste0("\n##",paste(rep("-",40),collapse="")))
            cat(paste0("\n##  ",title,"\n"))
            cat(paste0("##",paste(rep("-",40),collapse="")))
            print(param); cat("\n\n")
            sink()
         }
         colnames(param)<-"Value";param
         
         
         if(all(length(stats) > 0 & nrow(stats) > 0)){
            
            if(!file.exists(file.path("logs","processing.log"))){ 
               sink(file="./logs/processing.log",append=FALSE)
            } else { sink(file="./logs/processing.log",append=TRUE) }
            
            # sink(file="./log.txt",append=TRUE)
            cat(paste0("\n##",paste(rep("-",40),collapse="")))
            cat(paste0("\n##  ",title,"\n"))
            cat(paste0("##",paste(rep("-",40),collapse="")))
            print(stats); cat("\n\n")
            sink()
         }
         colnames(stats)<-"Value";stats
      } ## SAVE
      
      data2<-list(param=param, stats=stats)
      return(data2)
      
   } ## MAKE LOG
   
   
   ##-----------------------------------------
   ##  [11] MAKE LOG FILES 2 (FROM RNASEQ)
   ##-----------------------------------------
   ## MAKE LOG FILES FROM RNA-SEQ FUNCTIONS
   make_log2 <- function(param, stats, title="", save=TRUE){
      
      # param<-data.frame(parameter=names(param),value=t(as.data.frame(t(param))));param
      # names<-names(param)
      param<-data.frame(value=t(as.data.frame(t(param))));param
      param$step<-rep(title,nrow(param));param
      colnames(param)<-NULL;param
      
      stats<-data.frame(value=t(as.data.frame(t(stats))));stats
      stats$step<-rep(title,nrow(stats));stats
      colnames(stats)<-NULL;stats
      
      if(save==TRUE){
         
         
         if(!dir.exists("logs")){dir.create("logs",recursive=TRUE)}
         
         if(length(param) > 0){
            
            if(!file.exists(file.path("logs","parameters.log"))){ 
               sink(file="./logs/parameters.log",append=FALSE)
               cat(c("parameter","value","step"))
            } else { sink(file="./logs/parameters.log",append=TRUE) }
            print(param);# cat("\n\n")
            sink()
         }
         # colnames(param)<-c("parameter","value","step")#"Value";param
         
         if(length(stats) > 0){
            
            if(!file.exists(file.path("logs","processing.log"))){
               sink(file="./logs/processing.log",append=FALSE) 
               cat(c("parameter","value","step"))
            } else { sink(file="./logs/processing.log",append=TRUE) }
            
            print(stats)#; cat("\n\n")
            sink()
         }
         
      } ## SAVE
      
      data2<-list(param=param, stats=stats)
      return(data2)
      
   } ## MAKE LOG
   
   
   
   ##--------------------------------
   ##  [12] SUBSET DATA (OPTIONAL) 
   ##--------------------------------
   ## factor = column name in targets to filter on (1 column name)
   ## rm.vals=character vector of factor levels used to remove sample rows 
   subset_targets <- function(targets, factor, rm.vals){
      
      lst <- c(); rm <- 0; 
      param <- stats <- list()
      
      if(any(is.null(factor) | length(factor)!=1)){
         stop(paste("Error! factor is either NULL or contains multiple values",
                    "filtering cannot be performed.","'factor' should refer to",
                    " 1 column name in targets."))
      }
      
      no_input <- nrow(targets)
      stats[[factor]] <- lst
      stats[["total_input_samples"]] <- no_input
      
      ## REMOVE ROWS
      if(!is.null(factor)){
         for(x in rm.vals){
            no_rows<-nrow(targets);no_rows
            if(any(targets[,factor] %in% x)){
               remove <- targets[,factor] %in% x
               targets <- targets[remove==FALSE,]
               lst <- c(lst,x)
               rm_rows <- no_rows-nrow(targets)
               stats[[paste0(factor," = ",x)]] <- rm_rows
               rm<-rm + rm_rows
               message(paste0("[",factor," == ",x,"]", " : ", rm_rows," samples removed", 
                              " : ", nrow(targets), " samples kept"))
            }}
         
      } ## RM ROWS
      
      stopifnot(no_input == rm + nrow(targets))
      param[["factor"]] <-factor
      param[["rm.vals"]] <- paste(lst,collapse=", ")
      stats[["total_samples_removed"]] <- rm
      stats[["total_samples_kept"]] <- nrow(targets)

      title<-paste0("SUBSET TARGETS STATS (",factor,")");title
      logs <-make_log(param=param, stats=stats,title=title,save=TRUE)
      
      
      message("metadata for ",rm," of the ",no_input," samples were removed from the targets file. Success!!\n")
      data2 <- list(targets=targets, stats=logs$stats,param=logs$param)
      return(data2)   
      
      
   } ## SUBSET
   
   
   
   ##----------------------------------------------------------
   ##  [13] FILTER MULTI-BATCH TMT (COMBAT) (OPTIONAL)
   ##----------------------------------------------------------
   ## removes rows that are not detected in any TMT batch
   ## removes rows where all samples in a TMT batch have NAs
   ## (i.e. not detected in one of the TMT batches)
   ## filtering NAs function to allow comBat correction
   ## (i.e. comBat requires that at least 1 sample/condition in each batch has a value)
   ## this function removes rows that are not detected in one of the TMT batches.
   ## returned data contains the proteins that were detected in ALL TMT batches 
   ## in at least one of the samples.
   ## NOTE: batch is vector and can be any blocking factor e.g batch, group or paste(group+batch),etc.
   filter_multi_batch_tmt <- function(data, batch){
      
      rm<-keep<-idList <- stats<-param<-NULL
      param[["batch"]] <- paste(names(table(batch)),"=",table(batch))
      stats[["no_input_rows"]] <- nrow(data)
      
      ## replace zeros with NAs
      data[,][data[,] == 0] <- NA
      
      ## logical matrix   T=NA/F=value
      temp=is.na(data);table(temp)
      
      ## remove rows with all NAs
      allNA=apply(temp,1,all);head(allNA) ## logical row=all NAs ==TRUE
      rm<-rownames(temp)[allNA==TRUE];head(rm);length(rm)
      pass.data=temp[allNA==FALSE,];dim(pass.data) ## rows with all NAs removed
      idList[["missing_in_all_batches"]] <- c(rm)
      
      stats[["no_na_rows_all"]] <- length(rm)
      
      ## remove rows where all samples in a group = NA
      for(x in unique(batch)){
         
         grp.data=pass.data[, batch %in% x];head(grp.data) ## subset data for a particular group
         grpNA=apply(grp.data,1,all);table(grpNA) ## logical identifying rows with all NAs (TRUE)
         rm.grp=rownames(pass.data)[grpNA==TRUE];length(rm.grp)
         
         idList[[paste0("missing_in_batch_",x)]] <- c(rm.grp)
         
         if(is.null(rm)){rm<-rm.grp}
         if(!is.null(rm)){rm<-c(rm,rm.grp)}
         pass.data<-pass.data[grpNA==FALSE, ];dim(pass.data)
         
         stats[[paste0("no_na_rows_",x)]] <- length(rm.grp)
         
      }
      head(rm);length(rm) ## row names that were removed
      keep <- rownames(pass.data);length(keep) ## rownames that were kept
      idList[["present_in_all_batches"]] <- c(keep)
      
      stats[["no_rows_removed"]] <- length(rm)
      stats[["no_rows_kept"]] <- length(keep)
      
      logs<-make_log(param=param, stats=stats, title="FILTER MULTI-BATCH TMT", save=TRUE)
      
      stopifnot(nrow(data)==(length(rm)+length(keep)))
      stopifnot(rownames(data) %in% c(rm,keep))
      data2<- list(data=data[keep, ], dataNA=data[rm, ], ids=idList, param=logs$param, stats=logs$stats)
      
      ## filtering processing stats
      print(paste("A total of",length(rm), "entries were removed from the data set leaving",
                  length(keep),"entries for further analysis. Success!!"))
      
      
      return(data2)
      
      
   } ## FILTER NAS (NULTIBATCH)
   
   
   
   ##------------------------------------------------
   ##  [13] FILTER (X REPS OF Y GROUPS) (OPTIONAL) UPDATED NEW
   ##------------------------------------------------
   ## group= column name in targets to guide filtering criteria.
   ## targets=data.frame of the sample info to extract must include 'sample' column
   ## data columns will be matched/subsetted to targets info based on target rownames.
   ## data = dataframe of raw intensities (no normalization/transformations)
   ## min,grps based on group column in targets
   filter_data <- function(data, targets, group="group", min.reps=3, min.grps=1){
      
      param <- stats <- list()
      
      stats[["total_input_samples"]] <- ncol(data)
      
      ## reorder/subset data matrix to match the row order in targets
      if(any(rownames(targets)%in%colnames(data))==FALSE){
         stop("Erorr! some row names of targets do not match any column names in data.")
      }
      data <- data[, rownames(targets)]
      stopifnot(rownames(targets)==colnames(data))
      
      stats[["no_removed_samples"]] <- stats[["total_input_samples"]] - ncol(data)
      stats[["no_filtered_samples"]] <- ncol(data)
      stats[["total_input_rows"]]=nrow(data)
      
      print(paste(stats[["no_removed_samples"]], "samples were removed from the data matrix..."))
      
      ## change data column names and targets row names
      ## to sample name i.e. sample column in targets
      if("sample" %in% colnames(targets)){
         stopifnot(rownames(targets)==colnames(data))
         print("column names of data and row names of targets converted to targets sample names...");cat("\n")
         rownames(targets) <- targets$sample
         colnames(data)    <- rownames(targets)
         stopifnot(colnames(data)==rownames(targets))
      }
      
      ## check that the group name provided is a valid column name in targets.
      ## extract group column as character vector and make a factor. 
      if(group%in%colnames(targets)==FALSE){
         stop("Error! group should be a column name in targets (e.g. group='group')") }
      groups <- make_factor(as.character(targets[, group]))
      
      param[["group"]] <- group
      param[["groups"]] <- paste(paste0(names(table(groups)),"=",table(groups)),collapse=", ")
      
      nreps <- table(groups); nreps              ##  number of samples in each group
      ngrps <- length(unique(groups));ngrps      ##  number of groups
      # print(group.names);print(nreps);print(ngrps)
      
      
      
      ## if min.reps is NULL set min.reps to the smallest group size
      ## if min.grps is NULL set the min.grps to 1
      if(is.null(min.reps)){ min.reps <- min(nreps);min.reps }
      if(is.null(min.grps)){ min.grps <- 1;  min.grps }
      
      ## if input no min.reps exceeds the max no. reps / group
      ## min.reps value is lowered to the smallest samle group.
      repsCutoff <- min.reps <= nreps;repsCutoff
      if(all(repsCutoff)==FALSE){
         warning(paste0("The min.reps threshold min.reps = ", min.reps, " exceeds the max. ",
                        "number of replicates per group: ", paste(names(nreps),nreps,sep="=",collapse=", "), 
                        ". \nThe min.reps ",
                        "threshold was lowered to equal the smallest sample group.\n",
                        "min.reps = ", min(nreps)))
         min.reps <- min(nreps); min.reps
         repsCutoff <- min.reps <= nreps;repsCutoff
      }
      
      ## repsCutoff is used to determine how many groups meet the min.reps criteria.
      ## if min.grps parameter is greater than grpsCutoff, then the threshold is lowered
      ## to one (i.e. )
      grpsCutoff <- sum(as.numeric(repsCutoff));grpsCutoff
      if(min.grps > grpsCutoff){
         message(paste0("Based on the min.reps parameter ",min.reps,". The min.grps threshold : ",
                        min.grps," exceeds what is allowed by the dataset: ", grpsCutoff,".\n",
                        "The min.grps threshold has been lowered to one ", 
                        "group. \nmin.grps = ", 1))
         min.grps <- 1
      }
      
      param[["min.reps"]] <- min.reps
      param[["min.grps"]] <- min.grps
      print(paste0("extracting entries with intensity > 0 in at least ", min.reps,
                   " of the samples in ", min.grps," or more groups..."))
      
      ## FILTERING 
      ## zeros are replaced with NA    
      tmpData <- data; dim(tmpData)
      tmpData[,][tmpData[,] == 0] <- NA; head(tmpData)
      
      ## calc. no samples in each group with intensities > 0 
      noSamplesPerGroup <- NULL
      for(x in unique(groups)){
         keep<-groups %in% x;table(keep)
         grpData<- data.frame(tmpData[,keep]);head(grpData)
         ## no samples in each group with intensities > 0 
         noSamplesPerGroup <- cbind(noSamplesPerGroup, apply(grpData, 1, FUN = function(x){sum(!is.na(x))}))
      }
      noSamplesPerGroup<-data.frame(noSamplesPerGroup)
      colnames(noSamplesPerGroup) <- unique(groups)
      rownames(noSamplesPerGroup) <- rownames(grpData)
      head(noSamplesPerGroup);dim(noSamplesPerGroup)
      
      ## min X samples in at least Y groups
      aboveCutoff <- apply(noSamplesPerGroup, 1 , FUN=function(x){ sum(x >= min.reps) >= min.grps });table(aboveCutoff)
      noSamplesPerGroup$aboveCutoff <-aboveCutoff;head(noSamplesPerGroup)
      filterData <- tmpData[aboveCutoff == TRUE, ]; dim(filterData)
      removeData <- tmpData[aboveCutoff == FALSE, ];dim(removeData)
      
      stats[["no_removed_rows"]]=nrow(removeData)
      stats[["no_filtered_rows"]]=nrow(filterData)
      
      ## filtering processing stats
      print(paste("A total of",nrow(removeData), "entries were removed from the data set leaving",
                  nrow(filterData),"entries for further analysis. Success!!"))
      
      logs<-make_log(param,stats,title="FILTERING X in Y",save=TRUE)
      
      data2 <- list(data=filterData, targets=targets, noSamplesPerGroup=noSamplesPerGroup, 
                    rm.data=removeData, param=logs$param, stats=logs$stats)
      return(data2)
      
      
   } ## FILTER DATA
   
   
   
   ##-----------------------------------
   ##  [14] NORMALIZE DATA (REQUIRED) (UPDATED NEW)
   ##-----------------------------------
   normalize_data <- function(data, targets){
      
      if(ncol(data) > nrow(targets)){
         warning("Warning! data < targets. the data set will be subsetted to match", 
                 "the samples listed in the targets file prior to normalization...")
      }
      if(any(rownames(targets)%in%colnames(data))==FALSE){
         stop("some row names in targets do not match the column names in data...")
      }
      data <- data[,rownames(targets)]
      stopifnot(colnames(data)==rownames(targets))
      
      ## change target row names to sample name
      ## change column names of data to target row names
      ## i.e. sample name
      if("sample"%in%colnames(targets)){
         print("column names of data and row names of targets converted to targets sample names...");cat("\n")
         rownames(targets) <- targets$sample
         colnames(data)<-rownames(targets)
         stopifnot(colnames(data)==rownames(targets))
      }
      
      ## create named list object to hold output dataframes of normalized data.
      ## names= normalization methods
      normList <- vector("list", length(norm.methods))
      names(normList) <- norm.methods
      
      ## apply normalization methods using functions listed above.  
      ## NOTE: most of the normalizations use log2(intensity) as input except
      ## VSN which normalizes using raw intensities with no log2 transformation.
      normList[["log2"]]     <- logNorm(dat=data)
      normList[["median"]]   <- medianNorm(logDat=normList[["log2"]])
      normList[["mean"]]     <- meanNorm(logDat=normList[["log2"]])
      normList[["vsn"]]      <- vsnNorm(dat=data)
      normList[["quantile"]] <- quantNorm(logDat=normList[["log2"]])
      normList[["cycloess"]] <- cycLoessNorm(logDat=normList[["log2"]])
      normList[["rlr"]]      <- rlrNorm(logDat=normList[["log2"]])
      normList[["gi"]]       <- giNorm(logDat=normList[["log2"]])
      
      data2 <- list(normList=normList, targets=targets)
      print(paste("Data normalized... Success!!"))
      return(data2)
      
   } ## NORM DATA
   
   
   ##---------------------------------
   ##  [15] PROCESS DATA (OPTIONAL) (UPDATED NEW)
   ##---------------------------------
   process_data <- function(data, targets, group="group", min.reps=3, min.grps=1){
      
      ## conditional filter applied to the data
      filt <- filter_data(data=data, targets=targets, group=group, 
                          min.reps=min.reps, min.grps=min.grps)
      
      ## apply 8 norm. methods to the filtered data set.
      ## output returned includes list object containing df for each norm. tech.
      ## filtered unnormalized dataset (data) is also returned 
      norm <- normalize_data(data=filt$data, targets=filt$targets)
      
      param<-filt$param;colnames(param)<-NULL
      stats<-filt$stats;colnames(stats)<-NULL
      logs<-make_log(param=as.list(unlist(param)), stats=as.list(unlist(stats)), 
                     title="DATA PROCESSING",save=TRUE)
      
      print("data processing complete... returning filtered and normalized intensity data...")
      data2<-list(normList=norm$normList, targets=norm$targets, filt=filt,
                  param=logs$param,stats=logs$stats)
      
      return(data2)
      
      
   } ## PROCESS DATA
   
   
   ##--------------------------------
   ##  RANK NORMALIZATION METHODS    
   ##--------------------------------
   rank_norm_methods <- function(normList, groups, method, sum.meths){
      
      ## ss = list object of summary stats
      ## mat = matrix of summary stats
      ## colms = column names of summary stats
      ## 
      # groups<-make_factor(targets$group);groups
      
      # sink(file="rank_norm.txt")
      ## named list object with pooled values for each group
      if(method=="PCV"){  qcdata=base::lapply(normList, function(x) PCV(x, groups=groups));decreasing=F;neg=1 }
      if(method=="PMAD"){ qcdata=base::lapply(normList, function(x) PMAD(x, groups=groups));decreasing=F;neg=1 }
      if(method=="PEV"){  qcdata=base::lapply(normList, function(x) PEV(x, groups));decreasing=F;neg=1 }
      if(method=="COR"){  qcdata=base::lapply(normList, function(x) unlist(COR(x, groups)));decreasing=T;neg=-1 }
      
      
      ## list containing data.frame of summary stats (min, 1stQu, median, mean, 3rdQu, IQR)
      ## for each norm. method 1 QC method.
      ss.lst <- base::lapply(qcdata, FUN=function(x){ 
         tmp <- IQR(x); names(tmp) <-"IQR"
         sum <- c(summary(x),tmp)
         names(sum) <- c("min","q1","med","mean","q3","max", "iqr")
         sum
      });ss.lst
      
      ## list object of summary stat metrics calc for 
      ## each norm method for a particular QC method i.e. PCV
      ## is converted to data.frame
      ss.mat <- matrix(nrow=length(names(ss.lst)), ncol=length(ss.lst[[1]]));ss.mat
      colnames(ss.mat) <- names(ss.lst[[1]]);ss.mat
      rownames(ss.mat) <- names(ss.lst);ss.mat
      for(x in names(ss.lst)){ ss.mat[x,]<-ss.lst[[x]] }
      ss.mat <- data.frame(ss.mat);ss.mat
      
      ## summary stat column names (named vector)
      ## subsetted to match sum.meths input in select_norm_method()
      colms <- colnames(ss.mat)[colnames(ss.mat) %in% sum.meths==T];colms
      # rw <- names(ss.mat);rw
      
      ## list object of ranks for each summary stat method
      ## sort summary stat matrix A -> Z based on one of the sum.meths 
      ## e.g. sorted A->Z based on median PCV values in col. 3.
      ## e.g. sort Z->A based on median COR values. norm. method with 
      ## lowest value is in row 1. a new tmp col. with 1-8 rank is added
      ## to end. 1 = best 8 = worst; duplicate values receive the same rank.
      rnk.lst <- base::lapply(c(1:length(colms)),function(x){
         o <- order(ss.mat[,colms[x]],decreasing=decreasing);o
         ss.mat<-ss.mat[o,];ss.mat
         ss.mat$rnk <- dplyr::dense_rank(neg*ss.mat[,colms[x]]);ss.mat ## 1:nrow(ss.mat)
         ss.mat <- ss.mat[, c(colms[x], "rnk")];ss.mat
         # colnames(ss.mat)[ncol(ss.mat)] <- "rnk";ss.mat
      });rnk.lst
      names(rnk.lst) <- colms;rnk.lst
      # [[COR]]
      # $min
      #                min rnk
      # log2     0.7811268   4
      # median   0.7818913   3
      # mean     0.7818913   3
      # 
      # $q1
      #                 q1 rnk
      # log2     0.8419107   2
      # median   0.8412799   4
      
      
      rnks <- data.frame(rnk.lst);rnks
      rnkcolms <- colnames(rnks)[grep(".rnk",colnames(rnks))];rnkcolms
      stcolms <- colnames(rnks)[colnames(rnks) %in% paste0(sum.meths,".",sum.meths)];stcolms
      if(length(sum.meths)==1){
         # df2 <- data.frame(rnks[,stcolms], row.names=rownames(rnks));df2
         rnks <- data.frame(rnks[, rnkcolms], row.names=rownames(rnks));rnks
         colnames(rnks) <- rnkcolms;rnks 
      }
      if(length(sum.meths) > 1){ 
         # df2<-rnks[,stcolms];df2
         # colnames(df2)<-paste0(method,".",gsub("\\..*","",colnames(df2)));df2
         rnks<-rnks[,rnkcolms];rnks
      }
      # sink()
      sink(file="./ranks.txt")
      print("ss.lst");print(ss.lst);cat("\n")
      print("ss.mat:");  print(ss.mat);cat("\n")
      # print("df2:");  print(df2);cat("\n")
      print("rnk.lst");print(rnk.lst);cat("\n")
      print("rnks:");   print(rnks);cat("\n")
      # sink()
      
      return(rnks)
      
   } ## RANK NORM
   
   
   ##---------------------------------------
   ##   SELECT TOP NORMALIZATION METHOD
   ##---------------------------------------
   select_norm_method <- function(normList, targets, method=c("PCV","PMAD","PEV","COR"), 
                                  sum.meths=c("min","q1","med","mean","q3","max","iqr")){
      
      sum.meths <- match.arg(arg=sum.meths, several.ok=TRUE);sum.meths
      method    <- match.arg(arg=method, several.ok=TRUE);method
      
      # sink(file="./select_norm.txt")
      param <- list()
      param[["method"]] <- paste(method,collapse=", ")
      param[["sum.meths"]] <- paste(sum.meths,collapse=", ")
      
      groups <- targets$group;groups
      
      sink(file="./ranks.txt", append=TRUE)
      print("dd")
      dd <- base::lapply(method, function(x){
         rank_norm_methods(normList=normList,groups,method=x, sum.meths=sum.meths)
      })
      names(dd)<-method;dd
      
      #####################################
      tst<-dd;tst
      neg <- ifelse(x=="COR",-1,1);neg
      for(x in names(tst)){
         
         tst[[x]]["sum"] <-rowSums(tst[[x]], na.rm=F);tst
         tst[[x]]["rank"] <- dplyr::dense_rank(tst[[x]][,"sum"]);tst
      };tst
      
      
      mt<-matrix(nrow=8,ncol=length(names(tst)))
      rownames(mt)<-rownames(tst[[1]]);mt
      colnames(mt)<-names(tst);mt
      
      ## 
      for(x in colnames(mt)){ mt[,x]<-tst[[x]][,"rank"] }
      mt<-data.frame(mt);mt
      mt$cmb.sum<-rowSums(mt,na.rm=F);mt
      mt$cmb.rank<-dplyr::dense_rank(neg*mt$cmb.sum);mt
      mt<-data.frame(mt);mt
      
      sink()
      
      # proteiNormPlots(normList,targets,save=F)
      
      ##################################
      # $PCV
      #          med.rnk iqr.rnk
      # log2           3       3
      # median         8       6
      # mean           7       7
      # vsn            1       8
      # quantile       4       4
      # cycloess       6       1
      # rlr            5       5
      # gi             2       2
      # 
      # $PMAD
      #          med.rnk iqr.rnk
      # log2           5       2
      # median         8       6
      # mean           7       7
      # vsn            1       8
      # quantile       6       3
      # cycloess       2       4
      # rlr            3       5
      # gi             4       1
      
      # ## SINGLE SUMMARY METHOD TESTED
      # if(length(sum.meths)==1){ 
      #    dd <- data.frame(dd) 
      #    colnames(dd) <- paste0(method,".",sum.meths);dd
      # }
      # 
      # ## MULTIPLE SUMMARY METHODS TESTED
      # if(length(sum.meths) > 1){
      #    dd <- data.frame(dd);dd
      #    tmp2 <- c() 
      #    for(x in sum.meths){ 
      #       tmp <- paste0(method,".",x)
      #       tmp2<- c(tmp2,tmp)
      #    };tmp2
      #    if(all(colnames(dd)%in%tmp2 & tmp2%in%colnames(dd))){ dd <- dd[,tmp2];head(dd) }
      # }
      # head(dd)
      # 
      # 
      # ## SINGLE QC METHOD TESTED
      # dat <- dat2 <- data.frame(row.names=rownames(dd));dat
      # if(length(method)==1){
      #    for(i in 1:length(sum.meths)){
      #       
      #       dd2 <- as.data.frame(dd[,grep(sum.meths[i],colnames(dd))]);dd2
      #       colnames(dd2)<-colnames(dd)[grep(sum.meths[i],colnames(dd))];dd2
      #       rownames(dd2)<-rownames(dd);dd2
      #       dd2$SUM <- rowSums(dd2, na.rm=F);dd2
      #       o <- order(dd2$SUM,decreasing=F);o
      #       dd2 <- dd2[o,];dd2
      #       dd2$RANK <- dplyr::dense_rank(dd2[,"SUM"]);dd2#1:nrow(dd2);dd2
      #       keep <- colnames(dd2)%in%c("SUM","RANK")
      #       colnames(dd2)[keep==TRUE] <- paste0(colnames(dd2)[keep==TRUE],".",sum.meths[i])
      #       colnames(dd2)
      #       dd2 <- dd2[rownames(dd),];dd2
      #       dat <- cbind(dat,dd2[,ncol(dd2)]);dat
      #       colnames(dat)[i] <- colnames(dd2)[ncol(dd2)];dat
      #    };dat
      #    
      # } ## SINGLE QC METHOD
      # 
      # ## MULTIPLE QC METHODS TESTED
      # if(length(method) > 1){
      # for(i in 1:length(sum.meths)){
      #    
      #    dd2 <- dd[,grep(sum.meths[i],colnames(dd))];dd2
      #    dd2$SUM <- rowSums(dd, na.rm=F);dd2
      #    # o <- order(dd2$SUM,decreasing=F);o
      #    # dd2 <- dd2[o,];dd2
      #    dd2$RANK <- dplyr::dense_rank(dd2[,"SUM"]);dd2#1:nrow(dd2);dd2
      #    keep <- colnames(dd2)%in%c("SUM","RANK")
      #    colnames(dd2)[keep==TRUE] <- paste0(colnames(dd2)[keep==TRUE],".",sum.meths[i])
      #    colnames(dd2)
      #    dd2 <- dd2[rownames(dat),];dd2
      #    
      #    dat2<-cbind(dat2,dd2);dat2 ## all c(PEV.med.rnk,COR.med.rnk,SUM.med,RANK.med)
      #    
      #    dat <- cbind(dat,dd2[,ncol(dd2)]);dat ## RANK.med,RNK.iqr only
      #    colnames(dat)[i] <- colnames(dd2)[ncol(dd2)];dat
      # }
      # };dat
      # 
      # dat$SUM.final <- rowSums(dat, na.rm=T);dat
      # o <- order(dat$SUM.final,decreasing=F);o
      # dat <- dat[o,];dat
      # dat$RANK.final <- dplyr::dense_rank(dat[,"SUM.final"]);dat ## c(1:nrow(dat))
      
      
      ## SELECT TOP RANKED NORMALIZATION METHOD
      
      ## SINGLE RANKED 1
      ## one norm. method is ranked number 1
      topNorm <- rownames(dat)[dat$RANK.final==1];topNorm
      if(length(topNorm)==1){
         topNorm <- rownames(dat)[1];topNorm
         print("one normalization method was ranked number 1 ...")
         print(paste("top ranked normalization method:", topNorm))
      }
      
      ## MULTIPLE RANKED 1
      ## multiple norm. methods ranked number 1. cycloess and vsn
      ## given priority over other methods. when cycloess and vsn
      ## are both ranked number 1 cycloess is given priorty over
      ## vsn.
      if(length(topNorm) > 1){
         print("multiple normalization methods were ranked number 1.")
         print(paste("candidates include: ", paste(topNorm,collapse=", ")))
         print("Note: cycloess and vsn are given priority over other methods...")
         ## if cycloess and vsn are in the list rated no. 1, select cycloess
         if(all("cycloess"%in% topNorm==TRUE & "vsn" %in% topNorm==TRUE)){
            print("cycloess selected over vsn")
            topNorm <- "cycloess"
         }
         if(all("cycloess"%in% topNorm==TRUE & "vsn" %in% topNorm==FALSE)){
            topNorm<-"cycloess"
         }
         if(all("cycloess"%in% topNorm==FALSE & "vsn" %in% topNorm==TRUE)){
            topNorm<-"vsn"
         }
         cat("\n");print(paste("top ranked normalization method:", topNorm))
      }
      
      param[["top_norm_methods"]] <- paste(rownames(data)[1:3], collapse=", ")
      param <- t(data.frame(param))
      colnames(param) <- "select.norm.parameters"
      sink(file="./param.txt",append=TRUE)
      print(param);cat("\n")
      sink()
      
      top3Norm <-rownames(dat)[dat$RANK.final%in%c(1,2,3)==TRUE];top3Norm
      names(top3Norm)<-dat$RANK.final[dat$RANK.final%in%c(1,2,3)==TRUE];top3Norm
      cat("\n");print("top 3 normalization methods: ")
      print(dat[top3Norm,c("SUM.final","RANK.final")])
      
      data2 <- list(dat=dat, topNorm=topNorm, dat2=dat2, dd=dd, param=param)
      return(data2)
      
   } ## SELECT NORM
   
   
   
   ##----------------------------
   ##  IMPUTE DATA (OPTIONAL)
   ##----------------------------
   tmp <- expand.grid(norm.methods, impute.methods);tmp  ## DO NOT REMOVE REQUIRED
   combo.methods <-paste(tmp[,1], tmp[,2],sep="_");combo.methods ## DO NOT REMOVE REQUIRED
   impute_data <- function(normList, targets, norm.meth=c("log2","median","mean","vsn","quantile",
                                                           "cycloess","rlr","gi"),
                           impute.meth=c("none","knn","qrilc","mindet","minprob","min","zero")){ 
      
      ##  CHECK INPUT
      norm.meth   <- match.arg(arg=norm.meth, 
                               choices=c(unique(names(normList),c("log2","median","mean","vsn","quantile",
                                                    "cycloess","rlr","gi"))), several.ok=TRUE)
      
      impute.meth <- match.arg(argimpute.meth, 
                               choices=c("none","knn","qrilc","mindet","minprob","min","zero"),
                               several.ok=TRUE)
      
      ## IMPUTE DATA
      impList <- list()
      for(j in 1:length(norm.meth)){
         
         normData <- normList[[norm.meth[j]]]
         normData <- normData[,rownames(targets)]
         stopifnot(colnames(normData)==rownames(targets))
         
         ## get imputed normData         
         for(i in 1:length(impute.meth)){
            
            if(impute.meth[i] == "none"){ 
               impData <- normData
               nm <- paste(norm.meth[j], impute.meth[i],sep="_") ##combo name
               impList[[nm]] <- impData
            } 
            
            if(impute.meth[i] == "knn"){ 
               impData <- impute::impute.knn(normData, rowmax = 0.9)$data
               nm <- paste(norm.meth[j], impute.meth[i],sep="_") ##combo name
               impList[[nm]] <- impData
            } 
            
            if(impute.meth[i] == "qrilc"){
               impData <- imputeLCMD::impute.QRILC(dataSet.mvs = normData)[[1]]
               nm <- paste(norm.meth[j], impute.meth[i],sep="_")
               impList[[nm]] <- impData
            }
            
            if(impute.meth[i] == "mindet"){
               impData <- imputeLCMD::impute.MinDet(dataSet.mvs = normData)
               nm <- paste(norm.meth[j], impute.meth[i],sep="_")
               impList[[nm]] <- impData
            }
            
            if(impute.meth[i] == "minprob"){
               impData <- imputeLCMD::impute.MinProb(dataSet.mvs = normData)
               nm <- paste(norm.meth[j], impute.meth[i],sep="_")
               impList[[nm]] <- impData
            }
            
            if(impute.meth[i] == "min"){
               normData[is.na(normData)] = min(normData, na.rm = TRUE)
               impData <- normData
               nm <- paste(norm.meth[j], impute.meth[i],sep="_")
               impList[[nm]] <- impData
            }
            
            if(impute.meth[i] == "zero"){
               normData[is.na(normData)] <- 0
               impData <- normData
               nm <- paste(norm.meth[j], impute.meth[i],sep="_")
               impList[[nm]] <- impData
            }
            ## "bpca" = pcaMethods::bpca(Matrix = normData)
            
         }  ## I LOOP
      } ## J LOOP
      
      stopifnot(rownames(impData)==rownames(normData))
      stopifnot(colnames(impData)==colnames(normData))
      
      return(list(impList=impList))
      
   } ## IMPUTE DATA
   
   
   
   
} ## MAIN FUNCTIONS


#####################################
##      PROTEINORM FUNCTIONS       ##
#####################################
{ ## PROTEINORM NORMALIZATION FUNCTIONS
   
   logNorm <- function(dat) {
      logInt <- log2(dat)
      #logInt <- replace(is.infinite(logInt), NA)
      logInt[is.infinite(as.matrix(logInt))] <- NA
      return(as.matrix(logInt))
   }
   
   medianNorm <- function(logDat) {
      # Find medians of each sample
      # Divide by median
      # Multiply by mean of medians
      sampleMed <- apply(logDat, 2, median, na.rm=TRUE)
      meanMed <- mean(sampleMed, na.rm=TRUE)
      out <- t(t(logDat) / sampleMed)
      out <- out * meanMed
      return(as.matrix(out))
   }
   
   meanNorm <- function(logDat) {
      # Find means of each sample
      # Divide by mean
      # Multiply by mean of means
      sampleMean <- apply(logDat, 2, mean, na.rm=TRUE)
      meanMean <- mean(sampleMean, na.rm=TRUE)
      out <- t(t(logDat) / sampleMean)
      out <- out * meanMean
      return(as.matrix(out))
   }
   
   vsnNorm <- function(dat) {
      vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
      colnames(vsnNormed) <- colnames(dat)
      row.names(vsnNormed) <- rownames(dat)
      return(as.matrix(vsnNormed))
   }
   
   quantNorm <- function(logDat) {
      quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy=FALSE)
      colnames(quantNormed) <- colnames(logDat)
      row.names(quantNormed) <- rownames(logDat)
      return(as.matrix(quantNormed))
   }
   
   cycLoessNorm <- function(logDat) {
      cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method="fast")
      colnames(cycLoessNormed) <- colnames(logDat)
      row.names(cycLoessNormed) <- rownames(logDat)
      return(as.matrix(cycLoessNormed))
   }
   
   rlrNorm <- function(logDat) {
      rlrNormed <- NormalyzerDE::performGlobalRLRNormalization(as.matrix(logDat), noLogTransform=TRUE)
      colnames(rlrNormed) <- colnames(logDat)
      row.names(rlrNormed) <- rownames(logDat)
      return(as.matrix(rlrNormed))
   }
   
   giNorm <- function(logDat) {
      giNormed <- NormalyzerDE::globalIntensityNormalization(as.matrix(logDat), noLogTransform=TRUE)
      colnames(giNormed) <- colnames(logDat)
      row.names(giNormed) <- rownames(logDat)
      return(as.matrix(giNormed))
   }
    
   
} ## PROTEINORM NORM FUNCTIONS



#####################################
##    PROTEINORM PLOT FUNCTIONS    ##
#####################################
{ ## PROTEINORM QC PLOT FUNCTIONS
   
   
   ##------------------------------
   ##  PLOT RIGHT MARGIN
   ##------------------------------
   right_margin <- function(x){
      maxchar <- max(nchar(as.character(x)))+10;#print(maxchar)
      right=maxchar/5;right
      return(right)
   }
   
   
   ##------------------------------
   ##  PLOT LEFT MARGIN
   ##------------------------------
   left_margin <- function(x){
      maxchar <- max(nchar(as.character(x)))+30;#print(maxchar)
      left=maxchar/5;left
      return(left)
   }
   
   
   ##------------------------------
   ##  PLOT HEIGHT
   ##------------------------------
   plot_height <- function(x){
      len<-length(x);len ## no samples
      height=(600/20)*length(x);#print(height)
      height<-ifelse(height>1000,1000,height)
      return(height)
      
   }
   
   
   ##------------------------------
   ##  [01] findDensityCutoff
   ##------------------------------
   ## Finds upper limit of the histogram so each box contains at least cutoffPercent of the data
   findDensityCutoff <- function(longdat, cutoffPercent=0.001) {
      densityCutoff <- 0
      newDensityCutoff <- max(longdat, na.rm=TRUE)
      totalNum <- length(longdat)
      while(newDensityCutoff != densityCutoff) {
         densityCutoff <- newDensityCutoff
         breakSeq <- seq(max(min(longdat, na.rm=TRUE), 0), densityCutoff,
                         by=(densityCutoff - max(min(longdat, na.rm=TRUE), 0)) / 30)
         freqs <- hist(longdat[longdat < densityCutoff], breaks=breakSeq, plot=FALSE)
         newDensityCutoff <- freqs$breaks[which((freqs$counts / totalNum) < cutoffPercent)[1]+1]
      }
      return(densityCutoff)
      
   } ## FINDDENSITYCUTOFF
   
   
   ##------------------------------
   ##  [02] PCV 
   ##------------------------------
   PCV <- function(data, groups){
      PCV=NULL
      for(group in unique(groups)){
         tempData=as.matrix(data[,groups %in% group])
         CVs=genefilter::rowSds(tempData, na.rm=FALSE) /
            rowMeans(tempData, na.rm=FALSE)
         PCV[group]=mean(CVs, na.rm=T)
      }
      return(PCV)
      
   } ## PCV
   
   
   ##------------------------------
   ##  [03] PMAD 
   ##------------------------------
   PMAD <- function(data, groups){
      PMAD=NULL
      for(group in unique(groups)){
         tempData=as.matrix(data[, groups %in% group])
         MAD=matrixStats::rowMads(tempData, na.rm=FALSE)
         PMAD[group]=mean(MAD, na.rm=T)
      }
      return(PMAD)
      
   } ## PMAD
   
   
   ##------------------------------
   ##  [04] PEV 
   ##------------------------------
   PEV <- function(data, groups){
      PEV=NULL
      for(group in unique(groups)){
         tempData=as.matrix(data[,groups %in% group])
         
         rowNonNACnt=rowSums(!is.na(tempData)) - 1
         EV=rowNonNACnt * matrixStats::rowVars(tempData, na.rm=FALSE)
         PEV[group]=sum(EV, na.rm=TRUE)/sum(rowNonNACnt, na.rm=TRUE)
      }
      return(PEV)
      
   } ## PEV
   
   
   ##------------------------------
   ##  [05] COR 
   ##------------------------------
   COR <- function(data, groups){
      COR=NULL
      for(group in unique(groups)){
         corGroup=NULL
         tempData=as.matrix(data[,groups %in% group])
         if(ncol(tempData)==1){ corGroup <- 1}
         if(ncol(tempData) > 1){
         corVals=stats::cor(tempData, use="pairwise.complete.obs", method="pearson")
         for (index in seq_len(ncol(corVals) - 1)) {
            corGroup <- c(corGroup, corVals[index, -(seq_len(index)), drop="FALSE"])
         }}
         COR[[group]]=corGroup
      }
      return(COR)
   } ## COR
   
   
   ##------------------------------
   ##  [06] plotPCV 
   ##------------------------------
   plotPCV <- function(normList, groups, batch=NULL, dir=".", save=FALSE){
      
      groups<-make_factor(as.character(groups));groups
      if(is.null(batch)){ batch <- c(rep("1",ncol(normList[[1]]))) }
      batch<-make_factor(as.character(batch),prefix=NULL)
      
      if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
      if(save==TRUE){png(filename=file.path(dir,"PCVplot.png"), units="px", 
                         width=650, height=650, pointsize=15)}
      
      par(mar=c(8,6,4,3))
      plotData=base::lapply(normList, function(x) PCV(x, groups))
      boxplot(plotData, main="PCV", las=2, col=binfcolors[1:length(normList)], 
              boxlwd=1, yaxt="n", xaxt="n", cex.main=1.5)
      axis(2,cex.axis=1.2, las=2)
      axis(side=1, at=base::seq_along(names(normList)), labels=names(normList), cex.axis=1.3, las=2)
      mtext(side=2, text="Pooled Coefficient of Variation", line=4.5, cex=ifelse(save==TRUE, 1.2,0.9))
      points(rep(base::seq_along(normList), each=length(plotData[[1]])), unlist(plotData), pch="*", cex=1)
      if(save==TRUE){dev.off()}
      
      return(invisible(plotData))
      
      
   } ## PLOTPCV
   
   
   ##------------------------------
   ##  [07] plotPMAD 
   ##------------------------------
   plotPMAD <- function(normList, groups, batch=NULL, dir=".", save=FALSE){
      
      if(is.null(batch)){ batch <- c(rep("1",ncol(normList[[1]]))) }
      
      if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
      if(save==TRUE){png(filename=file.path(dir, "PMADplot.png"), units="px", 
                         width=650, height=650, pointsize=15)}
      
      par(mar=c(8,6,4,3))
      plotData=base::lapply(normList, function(x) PMAD(x, groups))
      boxplot(plotData, main="PMAD", las=2, col=binfcolors[1:length(normList)], boxlwd=1,
              yaxt="n", xaxt="n", cex.main=1.5)
      axis(2,cex.axis=1.2, las=2)
      axis(side=1, at=base::seq_along(names(normList)), labels=names(normList), cex.axis=1.3, las=2)
      mtext(side=2, text="Median Absolute Deviation", line=4.5, cex=ifelse(save==TRUE, 1.2,0.9))
      points(rep(base::seq_along(normList), each=length(plotData[[1]])), unlist(plotData), pch="*", cex=1)
      
      if(save==TRUE){dev.off()}
      
      return(invisible(plotData))
      
      
   } ## PLOTPMAD
   
   
   ##------------------------------
   ##  [08] plotPEV 
   ##------------------------------
   plotPEV <- function(normList, groups, batch=NULL, dir=".", save=FALSE){
      
      groups<-make_factor(as.character(groups));groups
      if(is.null(batch)){ batch <- c(rep("1",ncol(normList[[1]]))) }
      batch<-make_factor(as.character(batch),prefix=NULL)
      
      if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
      if(save==TRUE){png(filename=file.path(dir,"PEVplot.png"), units="px", 
                         width=650, height=650, pointsize=15)}
      par(mar=c(8,6,4,3))
      plotData=base::lapply(normList, function(x) PEV(x, groups))
      boxplot(plotData, main="PEV", las=2, col=binfcolors[1:length(normList)], boxlwd=1, 
              yaxt="n", xaxt="n", cex.main=1.5)
      axis(2,cex.axis=1.2, las=2)
      axis(side=1, at=base::seq_along(names(normList)), labels=names(normList), cex.axis=1.3, las=2)
      mtext(side=2, text="Pooled Estimate of Variance", line=4.5, cex=ifelse(save==TRUE, 1.2,0.9))
      points(rep(base::seq_along(normList), each=length(plotData[[1]])), unlist(plotData), pch="*", cex=1)
      if(save==TRUE){dev.off()}
      
      return(invisible(plotData))
      
   } ## PLOTPEV
   
   
   ##-----------------------------
   ##  [09] plotCOR 
   ##-----------------------------
   plotCOR <- function(normList, groups, batch=NULL, dir=".", save=FALSE){
      
      groups<-make_factor(as.character(groups));groups
      if(is.null(batch)){ batch <- c(rep("1",ncol(normList[[1]]))) }
      batch<-make_factor(as.character(batch),prefix=NULL)
      
      if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
      if(save==TRUE){png(filename=file.path(dir, "CORplot.png"), units="px", 
                         width=650, height=650, pointsize=15)}
      
      par(mar=c(8,6,4,3))
      plotData=base::lapply(normList, function(x) unlist(COR(x, groups)))
      boxplot(plotData, main="Cor", las=2, col=binfcolors[1:length(normList)],
              boxlwd=1, yaxt="n", xaxt="n", cex.main=1.5)
      axis(side=2,cex.axis=1.2, las=2)
      axis(side=1, at=base::seq_along(names(normList)), labels=names(normList), cex.axis=1.3, las=2)
      mtext(side=2, text="Intragroup Correlation", line=4.5, cex=ifelse(save==TRUE,1.2, 0.9))
      points(rep(base::seq_along(normList), each=length(plotData[[1]])), unlist(plotData), pch="*", cex=1)
      if(save==TRUE){dev.off()}
      
      return(invisible(plotData))
      
      
   } ## PLOTCOR
   
   
   ##------------------------------
   ##  [10] densityLog2Ratio 
   ##------------------------------
   densityLog2Ratio <- function(normList, groups, zoom=FALSE, legend=TRUE, inset=-0.2){

      log2Ratio=NULL
      groupList=sort(unique(groups));groupList
      for(method in names(normList)){
         for(g1 in 1:(length(groupList)-1)){
            for(g2 in (g1+1):length(groupList)){
               log2Ratio[[method]]=c(log2Ratio[[method]],
                                     rowMeans(as.data.frame(normList[[method]][,groups==groupList[g1]]), na.rm=T) - 
                                        rowMeans(as.data.frame(normList[[method]][,groups == groupList[g2]]), na.rm=T))
            }}}
      
      ## minX=min(unlist(base::lapply(log2Ratio, FUN=function(x) min(density(x, na.rm=T)$x))))
      ## maxX=max(unlist(base::lapply(log2Ratio, FUN=function(x) max(density(x, na.rm=T)$x))))
      maxY=max(unlist(base::lapply(log2Ratio, FUN=function(x) max(density(x, na.rm=T)$y))));maxY
      minY=0
      
      if(zoom==FALSE){
         minX=0.5 * min(unlist(min(density(log2Ratio[["vsn"]], na.rm=T)$x)))
         maxX=0.5 * max(unlist(max(density(log2Ratio[["vsn"]], na.rm=T)$x)))
         maxY=maxY
         minY=minY
      }
      if(zoom==TRUE){
         minX=-0.3
         maxX=0.3
         maxY=maxY + (0.2*maxY)
         minY=maxY - (0.5*maxY)
      }
      
      if(legend==TRUE){par(mar=c(5,5,4,3))}
      if(legend==FALSE){par(mar=c(5,5,3,4))}
      plot(NA, las=1, xlim=c(minX, maxX), ylim=c(minY,maxY), xlab="Log2 ratio", ylab="Density",
           main="Log2-ratio", cex.main=1.5, cex.axis=1.2, cex.lab=1.3)
      abline(v=0, lwd=2,lty=3, col="grey")
      
      densityList<-list()
      for(method in names(normList)){
         lines(density(log2Ratio[[method]], na.rm=T), 
               ## col=grDevices::rainbow(length(log2Ratio))[which(names(log2Ratio) %in% method)],
               col=binfcolors[which(names(log2Ratio) %in% method)], 
               ## lty=ifelse(which(names(log2Ratio) %in% method) %% 2 == 0, 2, 1),
               lwd=3)
         
         densityList[[method]] <- density(log2Ratio[[method]],na.rm=T)
         
         # den=density(log2Ratio[[method]],na.rm=T)
         # keep=den$y>=max(den$y)/2 ## points above half-height ==TRUE
         # half.y<-den$y[keep==TRUE] ## yvalues forming above half height
         # half.x<-den$x[keep==TRUE] ## x-values forming above half height
         
         # lf<-grep(min(half.x),half.x);lf
         # lower<-c(half.x[lf],half.y[lf]);lower ## left
         # rt <-grep(max(half.x),half.x);rt    
         # upper<-c(half.x[rt],half.y[rt]);upper ## right
         # tp <- grep(max(half.y),half.y)
         # center=c(half.x[tp],half.y[tp]);center
         # width=upper[1]-lower[1] ## distance across x at half height (width at half height)
         # print(paste0("meth = ",method,"      width = ",round(width,2),
         #              "      max.y = ",round(center[2],2),"      max.x = ",round(center[1],2)))
      }
      if(legend==TRUE){
         legend("topright", inset=c(inset, 0), names(log2Ratio), bty="n", xpd=TRUE,
                box.col="transparent", box.lwd=0, bg="transparent", border="transparent", 
                col="transparent", pch=22, pt.bg=binfcolors[1:length(log2Ratio)], 
                pt.cex=1.5, cex=1, horiz=FALSE, ncol=1)
      }
      
      return(invisible(densityList))
      
   } ## DENSITYLOG2RATIO
   
   
   ##------------------------------
   ##  [11] plotLogRatio 
   ##------------------------------
   plotLogRatio <- function(normList, groups, batch=NULL, sampleLabels=NULL, 
                            zoom=FALSE, legend=TRUE, inset=0.02, dir=".", save=FALSE){
      
      groups<-make_factor(as.character(groups));groups
      if(is.null(batch)){ batch <- c(rep("1",ncol(normList[[1]]))) };batch
      batch<-make_factor(as.character(batch),prefix=NULL)
      if(is.null(sampleLabels)){ sampleLabels <- colnames(normList[[1]]) };sampleLabels
      
      if(save==TRUE){ if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)} }
      if(save==TRUE){
         filename <- file.path(dir, paste0("Log2RatioPlot",ifelse(zoom==T,"-zoom.png",".png")));filename
         png(filename=filename, units="px", width=650, height=650,pointsize=15) 
      }
      den <- densityLog2Ratio(normList, groups, zoom, legend, inset)
      if(save==TRUE){dev.off()}
      
      return(invisible(den))
      
   } ## PLOTLOGRATIO
   
   
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
         width=round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02,0);print(width)
         height=800
         ncols=3;print(ncols)
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
      if(save==T){dev.off()}
      
      return(invisible(barList))
      
   } ## PLOTTOTINTEN
   
   
   ##------------------------------
   ##  [13] heatmapClustered 
   ##------------------------------
   heatmapClustered <- function(missing, groups, batch=NULL, groupCol, batchCol, sampleLabels=NULL,legend=FALSE){
      
      if(is.null(batch)){ batch <- c(rep("1",ncol(missing))) }
      if(is.null(sampleLabels)){ sampleLabels <- colnames(missing) };sampleLabels
      
      ColAnn <- HeatmapAnnotation(Sample=groups, Batch=batch, col=list(Sample=groupCol, Batch=batchCol), 
                                  annotation_legend_param=list(Sample=list(title="Group", at=unique(groups), 
                                                                           labels=paste("", unique(groups))),
                                                               Batch=list(title="Batch", at=unique(batch), 
                                                                          labels=paste("Batch", unique(batch)))),
                                  show_legend = ifelse(legend==FALSE,FALSE,TRUE)
      )
      hm_clust=Heatmap(missing+0, col=c("white", "black"), column_names_side="top",  
                       column_title="Clustered", show_row_names=FALSE, show_column_names=TRUE, 
                       name="Status", column_names_gp=gpar(fontsize=7), 
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
      
      ColAnn <- HeatmapAnnotation(Sample=groups_groupSorted, Batch=batch_groupSorted, 
                                  col=list(Sample=groupCol, Batch=batchCol), 
                                  annotation_legend_param=list(
                                     Sample=list(title="Group", at=sort(unique(groups_groupSorted)),
                                                 labels=paste("", sort(unique(groups_groupSorted)))),
                                     Batch=list(title="Batch", at=sort(unique(batch_groupSorted)),
                                                labels=paste("Batch", sort(unique(batch_groupSorted))))
                                  ),
                                  show_legend=ifelse(legend==FALSE,FALSE,TRUE))
      hm_group=Heatmap(missing_groupSorted+0, col=c("white", "black"), column_names_side="top",  
                       column_title="Sorted by Groups", show_row_names=FALSE, 
                       show_column_names=TRUE, name="Status", 
                       column_names_gp=gpar(fontsize=7), 
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
      ColAnn <- HeatmapAnnotation(Sample=groups_batchSorted, Batch=batch_batchSorted, 
                                  col=list(Sample=groupCol, Batch=batchCol),
                                  annotation_legend_param=list(Sample=list(title="Group", at=sort(unique(groups_batchSorted)),
                                                                           labels=paste("", sort(unique(groups_batchSorted)))),
                                                               Batch=list(title="Batch", at=sort(unique(batch_batchSorted)),
                                                                          labels=paste("Batch", sort(unique(batch_batchSorted))))
                                  ), show_legend=ifelse(legend==FALSE,FALSE,TRUE))
      hm_batch=Heatmap(missing_batchSorted+0, col=c("white", "black"), column_names_side="top",  
                       column_title="Sorted By Batch", show_row_names=FALSE, show_column_names=TRUE, 
                       name="Status", column_names_gp=gpar(fontsize=7), 
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
      
      draw(hm_clust+hm_group+hm_batch, heatmap_legend_side="right", annotation_legend_side="right",
           ht_gap=unit(2, "cm"), column_title="Missing Values")
      
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
      if(save==TRUE){dev.off()}
      
      if(save==TRUE){png(filename=file.path(dir, "NaHMplot_clust.png"), units="in", 
                         width=width2, height=8, res=100,pointsize=8)}
      hm_clust = heatmapClustered(missing=miss, groups, batch, groupCol, batchCol, sampleLabels, legend=TRUE)
      draw(hm_clust)
      dev.off()
      
      if(save==TRUE){png(filename=file.path(dir, "NaHMplot_group.png"), 
                         units="in", width=width2, height=8, res=100,pointsize=8)}
      hm_group = heatmapGroup(missing=miss, groups, batch, groupCol, batchCol, sampleLabels, legend=T)
      draw(hm_group)
      dev.off()
      
      if(save==TRUE){png(filename=file.path(dir, "NaHMplot_batch.png"), 
                         units="in", width=width2, height=8, res=100, pointsize=8)}
      hm_batch = heatmapBatch(missing=miss, groups, batch, groupCol, batchCol, sampleLabels, legend=T)
      draw(hm_batch)
      dev.off()

      data2 <- list(missing=miss, hm_clust=hm_clust,hm_group=hm_group,hm_batch=hm_batch)
      return(invisible(data2))
      
      
   } ## PLOTNAHM
   
   
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
      
      ColAnn <- HeatmapAnnotation(Sample=groups, 
                                  col=list(Sample=colorGroup2(groups)),
                                  annotation_legend_param=list(
                                     Sample=list(title="Groups", at=levels(groups),
                                                 labels=paste(levels(groups)))))
      
      RowAnn <- rowAnnotation(Batch=batch, col=list(Batch=colorBatch(batch)),
                              annotation_legend_param=list(
                                 Batch=list(title="Batch", at=levels(batch),
                                            labels=paste("Batch", levels(batch)))))
      
      # col_fun = colorRamp2(c(-1,0, 1), colors=c("#47d604","white","#f54c57"), transparency = 0.7)
      hm_corr=Heatmap(cor_mat, name="Pearson correlation", border=TRUE,
                      col=circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7)),
                                               # colors=RColorBrewer::brewer.pal(8, "RdBu")[8:1],
                                               colors=RColorBrewer::brewer.pal(8, "Blues"),
                                               transparency=0.6),
                      # col=col_fun(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7))),
                      heatmap_legend_param=list(color_bar="continuous", 
                                                legend_direction="horizontal", 
                                                legend_width=unit(5, "cm"), 
                                                title_position="topcenter"),
                      column_names_gp=gpar(fontsize=fontsize), 
                      row_names_gp=gpar(fontsize=fontsize), 
                      top_annotation=ColAnn, 
                      left_annotation=RowAnn,
                      column_labels=sampleLabels, 
                      row_labels=sampleLabels,
                      ## function adds txt to cells of heatmap
                      cell_fun = function(j, i, x, y, width, height, fill){
                         grid.text(sprintf("%.1f",cor_mat[i, j]), x, y,gp=gpar(fontsize=fontsize-2)) } 
      )
      par(mar=c(10,10,10,10))
      draw(hm_corr, heatmap_legend_side="top",ht_gap=unit(2, "cm"))
      if(save==TRUE){dev.off()}
      
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
      if(save==TRUE){dev.off()}
      
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
      if(save==TRUE){dev.off()}
      
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
      if(save==TRUE){ dev.off() }
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
      
      if(save==TRUE){ dev.off() }
      
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
      
      
   } ## CORR SCATTER
   
   
   ##------------------------------
   ##  [24] PROTEINORM REPORT 
   ##------------------------------
   make_proteinorm_report <- function(normList, groups=NULL, batch=NULL, sampleLabels=NULL, legend=TRUE,
                                      enrich=c("protein","phospho"), dir=NULL, file=NULL, save=FALSE, keep.png=FALSE){
      
      enrich <- match.arg(enrich,choices=c("protein","phospho"), several.ok=FALSE)
      
      if(is.null(sampleLabels)){ sampleLabels<-colnames(normList[[1]]) };sampleLabels
      if(is.null(groups)){ groups <- rep("group",ncol(normList[[1]])) }
      groups<-make_factor(x=as.character(groups));groups
      if(!is.null(batch)){ batch <- make_factor(as.character(batch)) };batch
      if(save==FALSE){ keep.png<-FALSE; dir <- "." }
      
      ## CREATE QC OUTPUT DIRECTORY
      ## if use enrich type to create QC output directory
      if(save==TRUE){
         
         if(is.null(dir)){
            if(enrich=="protein"){ dir <- file.path("protein_analysis","01_quality_control") }
            if(enrich=="phospho"){ dir <- file.path("phospho_analysis","01_quality_control") }
         } ## DIR NULL
         if(!is.null(dir)){ if(!dir.exists(dir)){ dir.create(dir, recursive=TRUE) }};dir
         print(paste("QC output directory:", dir))
         
         ## FILENAME DEFINED
         if(!is.null(file)){
            if(file_ext(file) != "pdf"){ 
               stop("\nError! Invalid output file type...\nincorrect: file = '",file,"'",
                    "\ncorrect:   file = 'proteiNorm_Report.pdf'") }
            pngdir=gsub(".pdf","",file)
            if(file.exists(file.path(dir,file))){
               file<-make_new_filename(x=file,dir=dir);file
               no<-sub(".pdf","",sub(".*_","",file));no
               pngdir=gsub(".pdf","",file)
            }
            if(!file.exists(file.path(dir,file))){ pngdir=gsub(".pdf","",file) }
            
         } ## FILE NOT NULL
         
         ## FILENAME NULL
         if(is.null(file)){ 
            file <-"proteiNorm_Report.pdf";no=""
            if(file.exists(file.path(dir,file))){
               file<-make_new_filename(x=file,dir=dir);file
               no<-sub(".pdf","",sub(".*_","",file));no
               pngdir=gsub(".pdf","",file)
            }
            if(!file.exists(file.path(dir,file))){ pngdir=gsub(".pdf","",file) }
         }  ## FILE NULL
         
         print(paste("QC output directory:", dir))
         print(paste("QC report file: ",file))
         print(paste("png directory: ",pngdir))
         
      } ## SAVE ==TRUE
      
      
      ## CREATE PLOTS
      nahm   <- plotNaHM(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,dir=dir,save=save)
      pcv    <- plotPCV(normList=normList,groups=groups,batch=batch,dir=dir,save=save)
      pmad   <- plotPMAD(normList=normList,groups=groups,batch=batch,dir=dir,save=save)
      pev    <- plotPEV(normList=normList, groups=groups,batch=batch,dir=dir,save=save)
      cor    <- plotCOR(normList=normList,groups=groups,batch=batch,dir=dir,save=save)
      lograt <- plotLogRatio(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,
                             zoom=FALSE,legend=TRUE,inset=0.02,dir=dir,save=save)
      plotLogRatio(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,
                   zoom=TRUE,legend=TRUE,inset=0.02,dir=dir,save=save)
      # nahm   <- plotNaHM(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,dir=dir,save=save))
      totint <- plotTotInten(normList=normList,groups=groups,batch=batch,sampleLabels=sampleLabels,dir=dir,save=save)      
      if(save==FALSE){ draw(nahm$hm_batch); draw(nahm$hm_clust); draw(nahm$hm_group) }

      data2 <- list(pcv=pcv,pmad=pmad,pev=pev,cor=cor,lograt=lograt,nahm=nahm,totint=totint)
      if(save==TRUE){ data2<-list(pcv=pcv,pmad=pmad,pev=pev,cor=cor,lograt=lograt,
                                  nahm=nahm,totint=totint,dir=dir,file=file) }
      
      
      ##  SAVE PROTEINORM REPORT PDF 
      if(save==TRUE){
         
         pdf(file.path(dir,file), paper="USr", pagecentre=TRUE, pointsize=10,width=12,height=8)
         
         files <- c("PCVplot.png","PMADplot.png","PEVplot.png","CORplot.png","Log2RatioPlot.png",
                    "Log2RatioPlot-zoom.png", "NaHMplot.png","NaHMplot_clust.png","NaHMplot_group.png",
                    "NaHMplot_batch.png", "TotIntenPlot.png")
         pnglist<-paste0(paste0(file.path(dir),"/"),files);pnglist
         thePlots<-lapply(1:length(pnglist), function(i){grid::rasterGrob(readPNG(pnglist[i],native=F))})
         
         do.call(gridExtra::grid.arrange,c(thePlots[1:6], ncol=3))
         do.call(gridExtra::grid.arrange,c(thePlots[1], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[2], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[3], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[4], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[5], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[7], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[8], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[9], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[10], ncol=1))
         do.call(gridExtra::grid.arrange,c(thePlots[11], ncol=1))
         
         dev.off()
         
         ## REMOVE PNG FILES
         if(keep.png==FALSE){
            unlink(pnglist)
            print(file.path(dir,file))
            print("png files removed...")
         } ## remove png files.
         
         ## KEEP PNG FILES
         if(keep.png==TRUE){
            if(!dir.exists(file.path(dir,pngdir))){dir.create(file.path(dir,pngdir),recursive=TRUE)}
            lapply(files,function(x){
               file.copy(from=file.path(dir,x), to=file.path(dir,pngdir,x))
               file.remove(file.path(dir,x))
            })
            print(paste("png files moved to :",file.path(dir,pngdir)))
         } ## KEEP
         
      } ## SAVE == TRUE
      
      return(invisible(data2))
      
      
   } ## PROTEIN NORM REPORT
   
   
   ##------------------------------
   ##  [25] QC REPORT  NEWEST WORKING
   ##------------------------------
   ## normList = either a named list object (list of df norm. intensities) or a data.frame or matrix
   ## containing normalized intensities
   ## if a list object then norm.meth is used to extract the corresponding data.frame. 
   ## file = pdf file name
   make_qc_report <- function(normList, norm.meth="cycloess", groups=NULL, batch=NULL, sampleLabels=NULL, 
                              stdize=TRUE, top=500, dims=c(1,2), cex.dot=2, clust.metric="euclidean",
                              clust.meth="complete",cex.names=1, xlim=NULL, ylim=NULL,
                              legend=TRUE, enrich=c("protein","phospho"), dir=NULL, file=NULL, save=FALSE, 
                              keep.png=FALSE){
      
      param<-stats<-list()
      
      norm.meth <- match.arg(arg=norm.meth, 
                             choices=unique(c(names(normList),"log2","median","mean","vsn","quantile",
                                              "cycloess","rlr","gi")), several.ok=FALSE);norm.meth
      
      clust.metric <- match.arg(arg=clust.metric, 
                                choices=c("pearson", "sqrt pearson", "spearman", "absolute pearson", 
                                          "uncentered correlation", "weird", "cosine", "euclidean", 
                                          "maximum", "manhattan", "canberra", "binary","minkowski"), 
                                several.ok=FALSE);clust.metric   
      
      clust.meth <- match.arg(arg=clust.meth, 
                              choices=c("ward.D","ward.D2","single","complete","average","mcquitty",
                                        "median","centroid"), several.ok=FALSE);clust.meth
      
      enrich <- match.arg(enrich,choices=c("protein","phospho"), several.ok=FALSE);enrich 
      
      
      ## NORMALIZED INTENSITY DATA
      ## list object or dataframe/matrix of norm. intensities
      if(any(class(normList) == "list")){ data<-normList[[norm.meth]] }
      if(any(class(normList) %in% "data.frame")){ data=as.matrix(normList) }
      if(any(class(normList) %in% "matrix")){ data=as.matrix(normList) };head(data);dim(data)
      
      top <- ifelse(top > nrow(data), nrow(data),top);top
      
      if(is.null(sampleLabels)){ sampleLabels<-colnames(data) };sampleLabels
      if(is.null(groups)){ groups <- rep("group",ncol(data)) }
      groups<-make_factor(x=as.character(groups));groups
      if(!is.null(batch)){ batch <- make_factor(as.character(batch)) };batch
      if(save==FALSE){ keep.png <- FALSE }
      
      
      ## CREATE OUTPUT DIRECTORIES, QC REPORT FILE NAME
      if(save==TRUE){
         
         ## CREATE QC OUTPUT DIRECTORY
         ## if use enrich type to create QC output directory
         if(is.null(dir)){
            if(enrich=="protein"){ dir <- file.path("protein_analysis","01_quality_control") }
            if(enrich=="phospho"){ dir <- file.path("phospho_analysis","01_quality_control") }
         }
         if(!is.null(dir)){ if(!dir.exists(dir)){ dir.create(dir, recursive=TRUE) }};dir
         
         
         ## FILENAME NOT NULL
         if(!is.null(file)){
            if(file_ext(file) != "pdf"){ 
               stop("\nError! Invalid output file type...\nincorrect: file = '",file,"'",
                    "\ncorrect:   file = 'QC_Report.pdf'") }
            pngdir=gsub(".pdf","",file)
            if(file.exists(file.path(dir,file))){
               file<-make_new_filename(x=file,dir=dir);file
               no<-sub(".pdf","",sub(".*_","",file));no
               pngdir=gsub(".pdf","",file)
            }
            if(!file.exists(file.path(dir,file))){ pngdir=gsub(".pdf","",file) }
            
         } ## FILE NOT NULL
         
         ## FILENAME NULL
         if(is.null(file)){ 
            file <-"QC_Report.pdf";no=""
            if(file.exists(file.path(dir,file))){
               file<-make_new_filename(x=file,dir=dir);file
               no<-sub(".pdf","",sub(".*_","",file));no
               pngdir=gsub(".pdf","",file)
            }
            if(!file.exists(file.path(dir,file))){ pngdir=gsub(".pdf","",file) }
         }  ## FILE NULL
         
         print(paste("QC output directory:", dir))
         print(paste("QC report file: ",file))
         print(paste("png directory: ",pngdir))
         
      } ## SAVE == TRUE
      
      
      ##-----------------
      ##  CREATE PLOTS
      ##-----------------
      ## BOX PLOTS
      box <- plotBoxplot(data=data,groups=groups, sampleLabels=sampleLabels,title="Box Plot",legend=legend,dir=dir,save=save)
      if(!is.null(batch)){
         if(save==TRUE){ 
            width=round(0.0191*length(groups)^2 + 12.082*length(groups) + 671.75,0);print(width)
            png(filename=file.path(dir,"BoxPlot2.png"), units="px", width=width, height=750, pointsize=15) }
            box2 <- plotBoxplot(data=data,groups=batch,sampleLabels=sampleLabels,title="Box Plot",legend=legend,dir=dir,save=FALSE)
         if(save==TRUE){ dev.off() }
      }
      
      ## VIOLIN PLOTS
      vio <- plotViolin(data,groups=groups,sampleLabels, title="Violin Plot", legend, dir, save)
      if(!is.null(batch)){
         width=round(0.0191*length(groups)^2 + 12.082*length(groups) + 671.75,0);print(width)
         if(save==TRUE){ png(filename=file.path(dir,"violinPlot2.png"), units="px", width=width, height=750, pointsize=15) }
         vio2 <- plotViolin(data=data,groups=batch,sampleLabels=sampleLabels, title="Violin Plot", legend=legend, dir=dir, save=FALSE)
         if(save==TRUE){ dev.off() }
      }
      
      ## PCA PLOTS
      pca <- plotPCA(data=data,groups=groups,sampleLabels=sampleLabels,title="PCA Plot", top=top, stdize=stdize, dims=dims,
                     cex.dot=cex.dot,xlim=xlim,ylim=ylim,legend=legend,dir=dir,save=save)
      if(!is.null(batch)){
         if(save==TRUE){ png(filename=file.path(dir,"PCAplot2.png"), units="px",width=750, height=650, pointsize=15) }
         pca2 <- plotPCA(data=data,groups=batch,sampleLabels=sampleLabels,title="PCA Plot", top=top,stdize=stdize, dims=dims,
                         cex.dot=cex.dot,xlim=xlim,ylim=xlim,legend=legend,dir=dir,save=FALSE)
         if(save==TRUE){ dev.off() }
      }
      
      ## CLUSTER DENDROGRAMS
      dendro <- plotDendrogram(data=data,groups=groups,sampleLabels=sampleLabels,top=top,stdize=stdize,
                               clust.metric=clust.metric,clust.meth=clust.meth,
                               cex.names=1,xlim=NULL,title="Cluster Dendrogram",legend=legend,dir=dir,save=save)
      if(!is.null(batch)){
         if(save==TRUE){ 
            width=round(0.0191*length(groups)^2 + 12.082*length(groups) + 671.75,0);print(width)
            png(filename=file.path(dir,"Dendrogram2.png"), units="px", width=width, height=650, pointsize=15) }
               dendro2 <- plotDendrogram(data=data,groups=batch,sampleLabels=sampleLabels,top=top,stdize=stdize,
                                         clust.metric=clust.metric,clust.meth=clust.meth,
                                         cex.names=1,xlim=NULL,title="Cluster Dendrogram",legend=legend,dir=dir,save=FALSE)
         if(save==TRUE){ dev.off() }
      }
      
      ## SAMPLE CORRELATION HEATMAP (PEARSON)
      corhm <- plotCorHM(data=data,groups=groups,batch=batch,sampleLabels=sampleLabels,dir=dir,save=save)
      
      
      ## PLOTS LIST OBJECT
      if(is.null(batch)){  plotList <- list(box=box, vio=vio, pca=pca, dendro=dendro, corhm=corhm,norm.meth=norm.meth,
                                            dir=ifelse(save==TRUE,dir,NULL), file=ifelse(save==TRUE,file,NULL)) }
      if(!is.null(batch)){ plotList <- list(box=box,box2=box2,vio=vio,vio2=vio2,pca=pca,pca2=pca2,
                                            dendro=dendro,dendro2=dendro2,corhm=corhm,norm.meth=norm.meth,
                                            dir=ifelse(save==TRUE,dir,NULL),file=ifelse(save==TRUE,file,NULL)) }
      # if(save==TRUE){
      # if(is.null(batch)){  plotList <- list(box=box, vio=vio, pca=pca, dendro=dendro, corhm=corhm,norm.meth=norm.meth,
      #                                       dir=dir,file=file) }
      # if(!is.null(batch)){ plotList <- list(box=box,box2=box2,vio=vio,vio2=vio2,pca=pca,pca2=pca2,
      #                                       dendro=dendro,dendro2=dendro2,corhm=corhm,norm.meth=norm.meth,
      #                                       dir=dir,file=file) }
      # } ## SAVE==TRUE
      
      
      ##-----------------------------
      ##  SAVE QC REPORT PDF FILE
      ##-----------------------------
      if(save==TRUE){
         
         ##  MAKE PDF FILE
         pdf(file.path(dir,file), paper="USr", pagecentre=TRUE, pointsize=15,width=12,height=8)
         
         ## QC_REPORT.PDF (NO BATCH)
         if(is.null(batch)){
            
            files <- c("BoxPlot.png","ViolinPlot.png","PCAplot.png","Dendrogram.png","CorrHeatmap.png")
            pnglist<-paste0(paste0(file.path(dir),"/"),files);pnglist
            thePlots<-lapply(1:length(pnglist), function(i){grid::rasterGrob(readPNG(pnglist[i],native=F))})
            
            if(ncol(data)<=50){
               do.call(gridExtra::grid.arrange,c(thePlots[1:4], ncol=2))
            }
            do.call(gridExtra::grid.arrange,c(thePlots[1], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[2], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[3], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[4], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[5], ncol=1))
            
         } ## BATCH IS NULL
         

         ## QC_REPORT.PDF (BATCH INCLUDED)
         if(!is.null(batch)){
            
            ## QC_REPORT.PDF
            files <-c("BoxPlot.png","BoxPlot2.png", "ViolinPlot.png","ViolinPlot2.png", 
                      "PCAplot.png","PCAplot2.png", "Dendrogram.png","Dendrogram2.png", "CorrHeatmap.png")
            pnglist<-paste0(paste0(file.path(dir),"/"),files);pnglist
            thePlots<-lapply(1:length(pnglist), function(i){grid::rasterGrob(readPNG(pnglist[i],native=F))})
            
            ## QC REPORT PDF FILE 
            if(ncol(data) <= 50){
            do.call(gridExtra::grid.arrange,c(thePlots[1:4], ncol=2))
            do.call(gridExtra::grid.arrange,c(thePlots[5:8], ncol=2))
            }
            do.call(gridExtra::grid.arrange,c(thePlots[1], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[2], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[3], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[4], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[5], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[6], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[7], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[8], ncol=1))
            do.call(gridExtra::grid.arrange,c(thePlots[9], ncol=1))
            
         } ## BATCH IS NOT NULL
         
         dev.off()
         
         ## REMOVE PNG FILES
         if(keep.png==FALSE){
            unlink(pnglist)
            print(file.path(dir,file))
            print("png files removed...")
         } ## remove png files.
         
         ## KEEP PNG FILES
         if(keep.png==TRUE){
            if(!dir.exists(file.path(dir,pngdir))){dir.create(file.path(dir,pngdir),recursive=TRUE)}
            lapply(files,function(x){
               file.copy(from=file.path(dir,x), to=file.path(dir,pngdir,x))
               file.remove(file.path(dir,x))
            })
            print(paste("png files moved to :",file.path(dir,pngdir)))
         } ## KEEP
         
      } ## SAVE == TRUE
      
      
      param[["norm.meth"]] <- norm.meth
      param[["batch"]] <- ifelse(is.null(batch),"NULL",paste(unique(batch),collapse=", "))
      param[["stdize"]] <-stdize
      param[["top"]] <- top
      param[["dims"]] <- paste(dims,collapse=", ")
      param[["clust.metric"]]<-clust.metric
      param[["clust.meth"]]<-clust.meth
      if(save==TRUE){ param[["dir"]]<-dir }
      if(save==TRUE){ param[["file"]]<-file }
      if(keep.png==TRUE){ param[["png.dir"]] <- pngdir }
      param[["enrich"]]<-enrich
      
      stats[["no_samples"]]  <- ncol(data)
      stats[["no_groups"]]   <- length(unique(groups))
      stats[["no_batches"]]  <- ifelse(is.null(batch),0, length(unique(batch)))
      stats[["tot_no_rows"]] <- nrow(data)
      stats[["top_no_rows"]] <- top
      
      logs<-make_log(param,stats, title="QC REPORT",save=TRUE);logs
      
      data2 <- list(plots=plotList, param=logs$param, stats=logs$stats)
      
      # merge_png_pdf(pnglist = pnglist, newpdf = file.path(dir,"plots2.pdf"), delete.files=TRUE)
      # merge_pdf_files(pdflist=c(file.path(dir,"plots1.pdf"),file.path(dir,"plots2.pdf")), 
      #                 newpdf=file.path(dir,"QC_Results.pdf"),delete.files=TRUE)
      
      print("QC report files created. Success!!")
      
      return(invisible(data2))

      
   } ## MAKE QC PLOTS
   
   
   
} ## PROTEINORM PLOT FUNCTIONS



#####################################
##       LIMMA DE FUNCTIONS        ##
#####################################
{ ## LIMMA DE FUNCTIONS
   
   ##-----------------------------------------
   ##  MAKE DESIGN MATRIX (REQUIRED) (NEW)
   ##-----------------------------------------
   ## paired column name in paired argument if limma's mixed effects model 
   ## will be used for DE analysis. (for complicated paired sample designs)
   ## ~0+group; dupCorr(block=targets$paired)
   ## for regular paired sample design include paired column name
   ## in factors argument ~0+group+paired
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
      
      
   } ## MAKE DESIGN (NEW)
   
   
   
   ##----------------------------------
   ##   MAKE CONTRASTS (REQUIRED)
   ##----------------------------------
   make_contrasts <- function(file=NULL, design){
      
      param<-stats<-list()
      
      if(is.null(file)){file=file.choose()}
      ## check that file exists in specified location. if file exits and is csv/tsv/txt
      file <-gsub("\\./","",file)
      file.dir <- dirname(file);file.dir
      file <- file.path(file.dir, match.arg(basename(file), list.files(file.dir),several.ok=F))
      filext <- tools::file_ext(file);filext
      filext <- match.arg(filext, c("csv","txt","tsv"),several.ok=F)
      if(filext=="txt" | filext=="tsv"){ sep="\t" }
      if(filext=="csv"){ sep=","}
      
      ## imports contrast file
      contrast.vec <- utils::read.csv(file=file, sep=sep, stringsAsFactors=F, header=F)
      contrast.vec <- contrast.vec[,1]
      contrast.vec <-gsub(" ", "", contrast.vec)
      
      ## extract groups included in each contrast
      contrast.grps <-extract_contrast_groups(contrast.vec=contrast.vec)
      
      ## if the groups defined in the contrast file do no match the design 
      ## then return the imported contrast.vec. 
      if(all(unique(unlist(contrast.grps))%in%colnames(design))==FALSE){
         message("the contrast file imported successfully, however, some of the groups ",
                 "included in the contrasts do not match the design matrix.",
                 "The contrast matrix will need to be defined manually using the makeContrasts() function.")
         message("contrasts <- makeContrasts(contrasts=contrast.vec, levels=design)\n",
                 "colnames(contrasts) <- gsub('=.*','',colnames(contrasts))")
         message("Returning the imported contrast.vec for additional processing ...")
         return(contrast.vec)
      }
      ## define contrasts
      contrasts <- makeContrasts(contrasts=contrast.vec, levels=design)
      colnames(contrasts) <- gsub("=.*","",colnames(contrasts))
      
      param[["file"]]<-file
      stats[["no_contrasts"]] <-ncol(contrasts)
      title="CONTRASTS"
      logs<-make_log(param=param, stats=stats, title=title, save=TRUE)
      
      data2=list(contrasts=contrasts, contrast.vec=contrast.vec,param=logs$param, stats=logs$stats)
      return(data2)
      
   } ## CONTRASTS
   
   
   
   
   
   ##------------------------------------
   ##  MAKE ALL CONTRASTS (OPTIONAL)
   ##------------------------------------
   # Adapted from Gordon Smyth https://support.bioconductor.org/p/9228/
   ## creates contrast matrix for all unique pairwise group combinations
   make_all_contrasts <- function(design, targets){
      
      designCols <- colnames(design)
      targets$group <- factor(targets$group, levels=ordered(unique(targets$group))); targets$group
      groupLevels <- levels(targets$group); groupLevels
      
      n = length(groupLevels);n
      stopifnot(identical(groupLevels, designCols[1:n]))
      contrasts <- matrix(0, length(designCols), choose(n, 2));contrasts
      rownames(contrasts) = designCols    
      colnames(contrasts) = 1:choose(n,2); contrasts
      k = 0
      for (i in 1:(n-1)){
         for (j in (i+1):n){
            k = k + 1
            contrasts[j, k] = 1
            contrasts[i, k] = -1
            colnames(contrasts)[k] = paste0(groupLevels[j], "_vs_", groupLevels[i])
            print(contrasts)
         }
      }        
      return(contrasts)
   }
   
   
   ##------------------------------------
   ##   LIMMA DE ANALYSIS (REQUIRED) (NEW)
   ##------------------------------------
   ## data=norm$normList[["vsn"]];annot=ext$annot[rownames(norm$normList[["vsn"]]), ];
   ## targets=des$targets; design=des$design; contrasts=contrasts; min.pval=0.055;
   ## min.lfc=1;adj.method="BH"; 
   ## mixed.effects=FALSE;pipe="DIA";enrich="protein"
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
      
      
   } ## GET_DE
   
   
   
   
   ##---------------------
   ##  ADD HYPERLINKS
   ##---------------------
   ## dataframe, url.col = column name with ids you want to make hyperlinks for
   ## url = url used for the ids, e.g. uniprot site, ncbi, mirbase etc.
   ## the ids column is replaced with the hyperlink formula format for excel
   make_excel_hyperlinks <- function(data, url.col, url){
      
      ids<-data[,url.col]
      tmp<-is.na(ids);table(tmp)
      url2<-paste0(url, ids)
      ids2=paste0("HYPERLINK(\"", url2,"\", \"", ids, "\")")
      ids2[tmp==TRUE]<-NA
      data[,url.col]<-ids2
      class(data[,url.col])<-"formula"
      
      return(data)
   }
   
   
   
   ##-----------------------
   ##  ADD LIMMA RESULTS
   ##-----------------------
   add_limma_results <- function(wb, sheetName=NULL, statList, annot, data, norm.method, 
                                 min.pval, min.lfc, pipe, enrich){
      
      ## merge column colors
      binfcolors <-c("#1F5EDC","#EE0010","#32CD32","#FF1493","#FF7F00",
                     "#A342FC","#00C8FF","#ADFF2F","#FFE100","#E36EF6","#009ACE","#996633")
      names(binfcolors)<-c("blueberry","cherry","apple","barbie","fanta",
                           "grape","ocean","mtndew","gold","orchid","aceblue","poop")
      
      ## column heading colors
      lightbinfcolors<-c("#c8d8f7","#ffdadd","#d0f3d0","#ffb1db","#ffe1c4",
                         "#dbb6fe","#c4f2ff","#e3ffb8","#fff8c4","#f1b8fb",
                         "#baeeff","#eddcca")
      names(lightbinfcolors)<-c("lightblueberry","lightcherry","lightapple","lightbarbie","lightfanta",
                                "lightgrape","lightocean","lightmtndew","lightgold","lightorchid",
                                "lightaceblue","lightpoop")
      
      
      normName <- names(norm.methods)[grep(norm.method,norm.methods)];normName
      
      if(pipe=="DIA" & enrich=="protein"){
         
         annotCols=c("id", "Protein.Name","Accession.Number","Molecular.Weight", "Protein.Group.Score",
                     "Identified.Peptide.Count","Exclusivity", 
                     "UniprotID", "Gene_name","Description")
         newNames=c("id", "Protein Name","Accession Number","Molecular Weight", "Protein Group Score",
                    "Identified Peptide Count","Exclusivity", 
                    "UniProt ID", "Gene Name","Description")
         annot.title <- "Protein Annotation"
         data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Intensities");data.title 
         
         if(is.null(sheetName)){ sheetName <- "Protein Results" } 
         ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
         if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
         print(sheetName)
      }
      
      if(pipe=="TMT" & enrich=="protein"){
         annotCols= c("id", "Fasta.headers","Majority.protein.IDs", "Score", "UniprotID",       
                      "Gene_name","Description")
         newNames= c("id", "Fasta headers","Majority protein IDs", "Score", "UniProt ID",       
                     "Gene Name","Description")
         annot.title <- "Protein Annotation"
         data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Exclusive MS1 Intensities");data.title
         # sheetName   <- "Protein Results"
         
         if(is.null(sheetName)){ sheetName <- "Protein Results" } 
         ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
         if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
         print(sheetName)
         
      }
      
      
      if(pipe=="LF" & enrich=="protein"){
         annotCols= c("id", "Fasta.headers","Majority.protein.IDs", "Score", "UniprotID",       
                      "Gene_name","Description")
         newNames= c("id", "Fasta headers","Majority protein IDs", "Score", "UniProt ID",       
                     "Gene Name","Description")
         annot.title <- "Protein Annotation"
         data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Intensities");data.title
         # sheetName   <- "Protein Results"
         
         if(is.null(sheetName)){ sheetName <- "Protein Results" } 
         ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
         if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
         print(sheetName)
      }
      
      
      if(pipe=="phosphoTMT" & enrich=="protein"){
         
         annotCols= c("id", "Fasta.headers", "Majority.protein.IDs", "Phospho..STY..site.IDs", "Score", 
                      "UniprotID","Gene_name","Description")
         newNames= c("id", "Fasta headers", "Majority protein IDs", "Phospho (STY) site IDs", "Score", 
                     "UniProt ID","Gene Name","Description")
         annot.title <- "Protein Annotation"
         data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Exclusive MS1 Intensities");data.title
         # sheetName   <- "Protein Results"
         
         if(is.null(sheetName)){ sheetName <- "Protein Results" } 
         ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
         if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
         print(sheetName)
         
      }
      
      if(pipe=="phosphoTMT" & enrich=="phospho"){
         
         annotCols <- c("id","Fasta.headers", "proGroupID","UniprotID","Gene_name",
                        "Description", "PEP","Score", "Localization.prob", "Phospho..STY..Probabilities",
                        "Flanking","phosAAPosition","Class")
         newNames <- c("id","Fasta headers", "Protein Group ID","UniProt ID","Gene Name", 
                       "Description", "PEP","Score", "Localization prob", "Phospho (STY) Probabilities",
                       "Flanking","Site","Class")
         
         annot.title <- "Protein/Phospho Annotation"
         data.title  <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Exclusive MS1 Intensities");data.title
         # sheetName   <- "Phospho Results"
         
         if(is.null(sheetName)){ sheetName <- "Phospho Results" } 
         ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
         if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
         print(sheetName)
         
      }
      
      
      ## get subset of annotation columns and rename them
      annot <- annot[, annotCols];head(annot)
      colnames(annot) <- newNames;head(annot)
      annot <- make_excel_hyperlinks(data=annot, url.col="UniProt ID", url="https://www.uniprot.org/uniprot/")
      
      
      ## add worksheet to workbook
      openxlsx::addWorksheet(wb=wb, sheetName=sheetName, tabColour="#FFE100")
      openxlsx::setRowHeights(wb=wb, sheet=sheetName, rows = 3, heights = 30);## increase row height for the header
      openxlsx::setRowHeights(wb=wb, sheet=sheetName, rows = 4, heights = 20);## increase row height for the header
      
      
      ##--------------
      ##  ANNOTATION
      ##--------------
      ## merge annotation columns
      annot.start=1
      annot.end=ncol(annot);annot.end
      
      openxlsx::mergeCells(wb=wb,sheet=sheetName, cols=annot.start:annot.end, rows=3)
      openxlsx::writeData(wb=wb,sheet=sheetName, x=annot.title, startCol=1, startRow=3, colNames=FALSE)
      mergeStyle <- openxlsx::createStyle(fgFill="#636262", halign="center", valign="center", textDecoration="bold",
                                          fontColour="white", fontSize=14, numFmt="TEXT")
      openxlsx::addStyle(wb=wb,sheet=sheetName, style=mergeStyle, rows=3, cols=1, stack=TRUE)
      
      
      ## write annotation data to worksheet
      colStyle <- openxlsx::createStyle(fgFill = "#d8d8d8", halign="left", valign="center", textDecoration="bold",
                                        fontColour="black", fontSize=11, numFmt="TEXT")
      openxlsx::addStyle(wb=wb,sheet=sheetName, style=colStyle, rows=4, cols=annot.start:annot.end, stack=TRUE)
      openxlsx::writeData(wb=wb, sheet=sheetName, x=annot,
                          startCol = annot.start,
                          startRow = 4,
                          colNames = TRUE,
                          rowNames = FALSE,
                          borders="columns",
                          keepNA=TRUE, na.string="NA", 
                          sep="\t")
      
      ## add hyperlink style (blue font, underlined)
      hyperlinkStyle <- openxlsx::createStyle(fontColour="#0000FF",halign="left",valign="center",textDecoration="underline")
      cols <- grep("UniProt ID",colnames(annot));cols
      openxlsx::addStyle(wb=wb, sheet=sheetName, style=hyperlinkStyle, cols=cols, rows=5:(nrow(annot)+4),
                         gridExpand=TRUE,stack=TRUE)
      
      
      
      ##--------------
      ##  DATA
      ##--------------
      ## start/stop column positions for data
      data.start=annot.end+1;data.start
      data.end=annot.end+ncol(data);data.end
      
      ## data merged column 
      openxlsx::mergeCells(wb=wb,sheet=sheetName, cols=data.start:data.end, rows=3)
      openxlsx::writeData(wb=wb,sheet=sheetName, x=data.title, startCol=data.start, startRow=3, colNames=FALSE)
      mergeStyle <- openxlsx::createStyle(fgFill = "#516285", halign="center", valign="center", textDecoration="bold",
                                          fontColour="white", fontSize=14, numFmt="TEXT")
      openxlsx::addStyle(wb=wb,sheet=sheetName, style=mergeStyle, rows=3, cols=data.start, stack=TRUE)
      
      ## write data to worksheet
      colStyle <- openxlsx::createStyle(fgFill = "#cdd4e1", halign="left", valign="center", textDecoration="bold",
                                        fontColour="black", fontSize=11, numFmt="TEXT")
      openxlsx::addStyle(wb=wb,sheet=sheetName, style=colStyle, rows=4, cols=data.start:data.end, stack=TRUE)
      openxlsx::writeData(wb=wb, sheet=sheetName, x=data,
                          startCol = data.start,
                          startRow = 4,
                          colNames = TRUE,
                          rowNames = FALSE,
                          borders="columns",
                          keepNA=TRUE, na.string="NA", 
                          sep="\t")
      
      
      ##--------------
      ##  STATS
      ##--------------
      ## colors for stat columns
      if(length(names(statList)) > length(binfcolors)){
         numColors<- ceiling(length(names(statList))/length(binfcolors))+1
         colors2 <- rep(binfcolors,numColors)
         lightcolors2 <- rep(lightbinfcolors,numColors)
      } else {
         if(length(names(statList))<=length(binfcolors)){
            colors2<-binfcolors
            lightcolors2<-lightbinfcolors
         }
      }
      stopifnot(length(colors2) > length(names(statList)))
      
      
      ## write stat data to worksheet
      stat.start=NULL;stat.end=NULL
      for(i in base::seq_along(names(statList))){
         stats <- statList[[i]]
         comparison <- names(statList)[i];comparison
         mycolor<-colors2[i];mycolor
         mylightcolor<-lightcolors2[i];mylightcolor
         
         ## stats start
         if(is.null(stat.start)){ 
            stat.start=data.end+1;stat.start 
         } else {
            if(!is.null(stat.start)){ 
               stat.start=stat.end+1;stat.start }
         }
         
         ## stats end
         if(is.null(stat.end)){ 
            stat.end=data.end+ncol(stats);stat.end 
         } else {
            if(!is.null(stat.end)){ 
               stat.end=stat.end+ncol(stats);stat.end }
         }
         
         
         openxlsx::mergeCells(wb=wb,sheet=sheetName, cols=stat.start:stat.end, rows=3)
         openxlsx::writeData(wb=wb,sheet=sheetName, x=comparison, startCol=stat.start, startRow=3, colNames=FALSE)
         mergeStyle <- openxlsx::createStyle(fgFill = mycolor, halign="center", valign="center", textDecoration="bold",
                                             fontColour="white", fontSize=14, numFmt="TEXT")
         openxlsx::addStyle(wb=wb,sheet=sheetName, style=mergeStyle, rows=3, cols=stat.start, stack=TRUE)
         
         
         ## write stats to worksheet
         colStyle <- openxlsx::createStyle(fgFill = mylightcolor, halign="left", valign="center", textDecoration="bold",
                                           fontColour="black", fontSize=11, numFmt="TEXT")
         openxlsx::addStyle(wb=wb,sheet=sheetName, style=colStyle, rows=4, cols=stat.start:stat.end, stack=TRUE)
         openxlsx::writeData(wb=wb, sheet=sheetName, x=stats,
                             startCol = stat.start,
                             startRow = 4,
                             colNames = TRUE,
                             rowNames = FALSE,
                             borders="columns",
                             keepNA=TRUE, na.string="NA", 
                             sep="\t")
         
         
         posStyle  <- openxlsx::createStyle(fontColour="#990000", bgFill="#FFC7CE")
         negStyle  <- openxlsx::createStyle(fontColour="#006100", bgFill="#C6EFCE")
         normStyle <- openxlsx::createStyle(fontColour="#000000", bgFill="#FFFFFF")
         
         fc.col    <- grep("logFC", colnames(stats))+stat.start-1;fc.col
         fc.rule1  <- paste0(">=",min.lfc);fc.rule1 
         fc.rule2  <- paste0("<=",-min.lfc);fc.rule2
         fdr.col   <- grep("adj.P.Val", colnames(stats))+stat.start-1; fdr.col
         fdr.rule  <- paste0("<=",min.pval);fdr.rule
         pval.col  <- grep("P.Value", colnames(stats))+stat.start-1;pval.col
         pval.rule <- paste0("<=",min.pval);pval.rule
         
         openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fc.col, rows=5:(nrow(stats)+4), 
                                         type="expression", rule=fc.rule1, style=posStyle)
         openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fc.col, rows=5:(nrow(stats)+4), 
                                         type="expression", rule=fc.rule2, style=negStyle)
         openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fc.col, rows=5:(nrow(stats)+4), 
                                         type="contains", rule="NA", style=normStyle)
         openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4), 
                                         type="expression", rule=fdr.rule, style=posStyle)
         openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4), 
                                         type="contains", rule="NA", style=normStyle)
         openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4), 
                                         type="expression", rule=pval.rule, style=posStyle)
         openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4), 
                                         type="contains", rule="NA", style=normStyle)
         
         
      } ## FOR LOOP
      
      
      ## add global styles
      headerStyle <- openxlsx::createStyle(halign = "center", valign="center", #fontColour = "#000000",  
                                           border = "TopBottomLeftRight", borderColour = c("black","black","black","black"),
                                           textDecoration = "bold")
      openxlsx::addStyle(wb=wb, sheet=sheetName, style=headerStyle, cols=1:stat.end, rows=3, gridExpand=TRUE, stack=TRUE)
      openxlsx::addStyle(wb=wb, sheet=sheetName, style=headerStyle, cols=1:stat.end, rows=4, gridExpand=TRUE, stack=TRUE)
      openxlsx::addFilter(wb=wb, sheet=sheetName, cols=1:stat.end, rows=4)
      
      
   } ## ADD LIMMA RESULTS
   
   
   
} ## LIMMA DE FUNCTIONS



####################################
##    LIMMA DE PLOT FUNCTIONS     ##
####################################
{ ## DE PLOT FUNCTIONS
   
   
   ##------------------------
   ##  [01]  VOLCANO PLOT
   ##------------------------
   volcanoPlot <- function(stats, comparison=NULL, min.pval=0.055, min.lfc=1, 
                           xlim=NULL, ylim=NULL, sig.type=c("pval","p.adj"), 
                           top=NULL, labels=NULL, inset=-0.3, legend=FALSE){
      
      stats<-data.frame(stats);head(stats)
      
      ## row names for DIA,TMT,phosphoTMT protein data are a combination 
      ## of UniprotID,GN, and protein group id. phospho data = uniprotid, GN, 
      ## and phospho site (e.g. S37). This info is parsed into a matrix. 
      ncols=length(unlist(strsplit(rownames(stats)[1],split="_")));ncols
      # mat  <- matrix(unlist(strsplit(rownames(stats),"_")), ncol=ncols, byrow=TRUE);head(mat)
      # rownames(mat)=rownames(stats)
      
      ## rownames(stats)[1489] ## "D3ZHA7_RGD1560334_predicted_1943"
      mat<-matrix(nrow=nrow(stats), ncol=3)
      mat[,1] <-gsub("_.*","",rownames(stats)) ## before first _
      mat[,2] <- sub("_[^_]+$","", sub("^[^_]*_", "", rownames(stats))) ## after first _ & before last _
      mat[,3] <- sub(".*_","",rownames(stats)) ## after last _
      rownames(mat)<-rownames(stats)
      
      
      ## annotation info extracted from rownames of protein data
      ## Gene_name will be used to label significant points if top >0
      if(length(grep("S",mat[,3]))==0){ 
         colnames(mat)<-c("UniprotID","Gene_name","id")
         mat=data.frame(mat);head(mat)
         mat$labels<-mat$Gene_name
         stats<-cbind(stats,mat[rownames(stats),])
      } 
      
      ## annotation info extracted from rownames of phospho data
      ## Gene_name(phosphosite) will be used to label sig. points if top >0
      if(length(grep("S",mat[,3]))>0){
         colnames(mat)<-c("UniprotID","Gene_name","phosAAPosition")
         mat<-data.frame(mat)
         mat$labels <- paste0(mat$Gene_name,"(",mat$phosAAPosition,")")
         head(mat)
         stats<-cbind(stats,mat[rownames(stats),]);head(stats)
      } 
      head(stats)
      
      ## column named 'p' added to stats data.frame based on sig.type selected. texted for y-axis label also defined
      if(sig.type=="p.adj"){
         stats$p <- stats$adj.P.Val
         ylab=expression(paste("-",log[10]," (adj. p-value)",sep=""))
         filename<-paste(comparison,"volcano_plot.png",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      
      if(sig.type=="pval"){
         stats$p<-stats$P.Value
         ylab=expression(paste("-",log[10]," (p-value)",sep=""))
         filename<-paste(comparison,"volcano_plot_pvalue.png",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      
      
      ## determine x/y limits for the plot. 
      if(length(xlim) !=2){ xlim <- c(-1,1)*(max(abs(stats$logFC), na.rm=TRUE)*1.1) }; xlim
      if(length(ylim) !=2){ ylim <- c(0, max(-log10(stats$p)*1.1, na.rm=TRUE)) };ylim
      
      
      ## global par plot parameters
      par(font.axis=1, lwd=1, font.main=2, cex.main=1.5, fg="black",
          col.axis="black", cex.axis=1.2, cex.lab=1.3, font.axis=1, xaxs="i", yaxs="i")
      if(legend==TRUE){par(mar=c(6,7,4,6))}
      if(legend==FALSE){par(mar=c(6,7,4,3))}
      
      ## create blank plot area
      plot(x=1, las=1, type="n", main="",xlab="",ylab="", xlim=xlim,ylim=ylim)
      title(line=1, main=comparison, cex.main=1.3)
      mtext(side=1, line=3, cex=1.3, text=expression(paste(log[2]," (fold-change)",sep="")))
      mtext(side=2, line=4.5,cex=1.3, text=ylab)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray94",lwd=0)
      grid(col="white", lty="dotted", lwd=1.2)
      box(which="plot", lwd=4, lty=1, col="gray94")
      
      ## subset based on significance cutoffs
      up   <- subset(stats, p <= min.pval & logFC >= min.lfc)
      down <- subset(stats, p <= min.pval & logFC <= -min.lfc)
      
      ## add points
      with(stats, points(logFC, -log10(p), pch=21, bg=alpha(colour="grey60",alpha=0.5), ## black
                         col=alpha(colour="grey50",alpha=0.5), lwd=0.5, cex=1.1)) 
      with(up, points(logFC, -log10(p), pch=21, bg=alpha(colour="#FF0000",alpha=1),     ## red
                      col=alpha(colour="#990000",alpha=1), lwd=0.5, cex=1.1)) 
      with(down, points(logFC, -log10(p), pch=21, bg="#1E90FF", col="#003399", lwd=0.5, ## blue
                        cex=1)) 
      
      abline(v= c(-min.lfc,min.lfc), col="dimgray", lty=2, lwd=0.8)
      abline(h= -log10(min.pval), col="dimgray", lty=2, lwd=0.8)
      
      
      ## ADD LABELS
      
      ## user  can supply a vector of labels to add to the points on the volcano plot.
      ## for protein data the list can be Gene_name, UniprotID, UniprotID_Gene_name, proKey,
      ## for phospho data phosKey, Gene_name or Gene_name(site). 
      ## labels added to corresponding points in orange.
      if(!is.null(labels)){
         if(length(labels)>0){
            
            ## get row numbers in stats corresponding to the labels list
            ## remove any labels that were not identified and subset
            ## stats. divide subdata into logFC>0 and logFC<0 so labels can be
            ## positioned on the R and L side of the point
            n<-keep<-c()
            for(i in 1:length(labels)){
               tmp=grep(paste0(labels[i],"_"),rownames(stats));tmp
               tmp2<-ifelse(length(tmp)==0,FALSE,TRUE)
               keep<-c(keep,tmp2);table(keep)
               n <- c(n,tmp);n
            };n
            
            ## ADD USER INPUT LABELS TO POINTS 
            # labs <- labels#[keep==T];labs
            subdata <- stats[n,];dim(subdata)
            
            ## ADD USER INPUT LABELS TO POS logFC > 0
            ## subset data for positive logFC values, add labels 
            subup<- subset(subdata, subdata$logFC>0);dim(subup)
            if(nrow(subup)>0){
               points(x=subup$logFC,y=-log10(subup$p),pch=21, cex=1.2,
                      bg="orange", col="darkorange3",lwd=0.5,
                      text(x=subup$logFC, y=-log10(subup$p), col="orange",
                           labels=subup$labels, pos=4,cex=0.6,font=2))
            }
            
            ## ADD USER INPUT LABELS TO NEG logFC < 0
            subdn<- subset(subdata, subdata$logFC<0);dim(subdn)
            if(nrow(subdn>0)){
               points(x=subdn$logFC,y=-log10(subdn$p),pch=21, cex=1.2,
                      bg="orange", col="darkorange3",lwd=0.5,
                      text(x=subdn$logFC, y=-log10(subdn$p), col="orange",
                           labels=subdn$labels, pos=2,cex=0.6,font=2))
            }
         } ## LEN >0 
      } ## LABELS
      
      
      ## TOP
      
      ## the user can use 'top' to label the most sig. up and down regulated
      ## features (proteins=gene name, phospho=gene name(site).
      ## e.g. top=5, will label the top 5 up and down reg.
      if(!is.null(top)){ ## TOP
         if(length(top)==1 & top>0){
            ## sig. up-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC Z-> A
            sig <- up[with(up,order(logFC,order(p,decreasing=FALSE),decreasing=TRUE)),]
            top2<-ifelse(top > nrow(sig), nrow(sig),top);top2
            sig<- sig[1:top2, ]
            points(x=sig$logFC,y=-log10(sig$p),pch=21, cex=1.2,
                   bg="green", col="forestgreen",lwd=0.5,
                   text(x=sig$logFC, y=-log10(sig$p), col="green3",
                        labels=sig$labels, pos=4,cex=0.6,font=2))
            
            ## sig. down-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC A->Z
            sig <- down[with(down,order(logFC,order(p,decreasing=FALSE),decreasing=FALSE)),]
            top2<-ifelse(top > nrow(sig), nrow(sig),top);top2
            sig<- sig[1:top2, ]
            points(x=sig$logFC,y=-log10(sig$p),pch=21, cex=1.2,
                   bg="green", col="forestgreen",lwd=0.5,
                   text(x=sig$logFC, y=-log10(sig$p), col="green3",
                        labels=sig$labels, pos=2,cex=0.6,font=2))   
            
         }
      } ## TOP
      
      
      if(legend==TRUE){
         legend("topright", c("up", "not sig.", "down"), box.col="transparent",
                box.lwd=0, bg="transparent", border="black", pch=22,
                pt.bg=c("#FF0000","grey50","#1E90FF"), pt.cex=1.4,
                col=c("#FF0000","grey50","#1E90FF"), horiz=FALSE,
                inset=c(inset, 0), bty="n", xpd=TRUE, ncol=1, cex=0.9)
      }
      
      
   } ## VOLCANO
   
   
   ##------------------
   ##  [02]  MD PLOT
   ##------------------
   mdPlot <- function(stats, comparison=NULL, min.pval=0.055, min.lfc=1, 
                      xlim=NULL, ylim=NULL, sig.type=c("pval","p.adj"), 
                      top=NULL, labels=NULL, inset=-0.3, legend=FALSE){
      
      stats<-data.frame(stats);head(stats)
      
      ## row names for DIA,TMT,phosphoTMT protein data are a combination 
      ## of UniprotID,GN, and protein group id. phospho data = uniprotid, GN, 
      ## and phospho site (e.g. S37). This info is parsed into a matrix. 
      ncols=length(unlist(strsplit(rownames(stats)[1],split="_")));ncols
      # mat  <- matrix(unlist(strsplit(rownames(stats),"_")), ncol=ncols, byrow=TRUE);head(mat)
      # rownames(mat)=rownames(stats)
      
      ## rownames(stats)[1489] ## "D3ZHA7_RGD1560334_predicted_1943"
      mat<-matrix(nrow=nrow(stats), ncol=3)
      mat[,1] <-gsub("_.*","",rownames(stats)) ## before first _
      mat[,2] <- sub("_[^_]+$","", sub("^[^_]*_", "", rownames(stats))) ## after first _ & before last _
      mat[,3] <- sub(".*_","",rownames(stats)) ## after last _
      rownames(mat)<-rownames(stats)
      
      
      ## annotation info extracted from rownames of protein data
      ## Gene_name will be used to label significant points if top >0
      if(length(grep("S",mat[,3]))==0){ 
         colnames(mat)<-c("UniprotID","Gene_name","id")
         mat=data.frame(mat);head(mat)
         mat$labels<-mat$Gene_name
         stats<-cbind(stats,mat[rownames(stats),])
      } 
      
      ## annotation info extracted from rownames of phospho data
      ## Gene_name(phosphosite) will be used to label sig. points if top >0
      if(length(grep("S",mat[,3]))>0){
         colnames(mat)<-c("UniprotID","Gene_name","phosAAPosition")
         mat<-data.frame(mat)
         mat$labels <- paste0(mat$Gene_name,"(",mat$phosAAPosition,")")
         head(mat)
         stats<-cbind(stats,mat[rownames(stats),]);head(stats)
      } 
      
      ## column named 'p' added to stats data.frame based on sig.type selected. texted for y-axis label also defined
      if(sig.type=="p.adj"){
         stats$p <- stats$adj.P.Val
         main.title = paste(comparison, "(adj. p-value)")
         filename<-paste(comparison,"MD_plot.png",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      
      if(sig.type=="pval"){
         stats$p<-stats$P.Value
         main.title = paste(comparison, "(p-value)")
         filename<-paste(comparison,"MD_plot_pvalue.png",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      
      
      ## determine x/y limits for the plot. 
      if(length(xlim) !=2){ xlim <- c(min(stats$AveExpr, na.rm=TRUE),max(stats$AveExpr, na.rm=TRUE)*1.1) }; xlim
      if(length(ylim) !=2){ ylim <- c(-1,1)*max(abs(stats$logFC), na.rm=TRUE)*1.1 };ylim
      
      
      ## global par plot parameters
      par(font.axis=1, lwd=1, font.main=2, cex.main=1.5, fg="black",
          col.axis="black", cex.axis=1.2, cex.lab=1.3, font.axis=1, xaxs="i", yaxs="i")
      if(legend==TRUE){par(mar=c(6,7,4,6))}
      if(legend==FALSE){par(mar=c(6,7,4,3))}
      
      ## create blank plot area
      plot(x=1, las=1, type="n", main="",xlab="",ylab="", xlim=xlim, ylim=ylim)
      title(line=1, main=main.title, cex.main=1.3)
      mtext(side=1, line=3,cex=1.3, text="avg. expression")
      mtext(side=2, line=3.5, cex=1.3, text=expression(paste(log[2]," (fold-change)",sep="")))
      
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray94")
      grid(col="white", lty="dotted", lwd=1.2)
      box(which="plot", lwd=4, lty=1, col="gray94")
      
      ## subset based on significance cutoffs
      up   <- subset(stats, p <= min.pval & logFC >= min.lfc)
      down <- subset(stats, p <= min.pval & logFC <= -min.lfc)
      
      ## add points
      with(stats, points(AveExpr, logFC, pch=21, bg=alpha(colour="grey60",alpha=0.5), ## black
                         col=alpha(colour="grey50",alpha=0.5), lwd=0.5, cex=1.1)) 
      with(up, points(AveExpr,logFC, pch=21, bg=alpha(colour="#FF0000",alpha=1),     ## red
                      col=alpha(colour="#990000",alpha=1), lwd=0.5, cex=1.1)) 
      with(down, points(AveExpr,logFC, pch=21, bg="#1E90FF", col="#003399", lwd=0.5, ## blue
                        cex=1)) 
      abline(h= c(-min.lfc,min.lfc), col="dimgray", lty=2, lwd=0.8)
      
      
      ## ADD LABELS
      
      ## user  can supply a vector of labels to add to the points on the volcano plot.
      ## for protein data the list can be Gene_name, UniprotID, UniprotID_Gene_name, proKey,
      ## for phospho data phosKey, Gene_name or Gene_name(site). 
      ## labels added to corresponding points in orange.
      if(!is.null(labels)){
         if(length(labels)>0){
            
            ## get row numbers in stats corresponding to the labels list
            ## remove any labels that were not identified and subset
            ## stats. divide subdata into logFC>0 and logFC<0 so labels can be
            ## positioned on the R and L side of the point
            n<-keep<-c()
            for(i in 1:length(labels)){
               tmp=grep(paste0(labels[i],"_"),rownames(stats));tmp
               tmp2<-ifelse(length(tmp)==0,FALSE,TRUE)
               keep<-c(keep,tmp2);table(keep)
               n <- c(n,tmp);n
            };n
            
            ## ADD USER INPUT LABELS TO POINTS 
            # labs <- labels#[keep==T];labs
            subdata <- stats[n,];dim(subdata)
            
            ## ADD USER INPUT LABELS TO POS logFC > 0
            ## subset data for positive logFC values, add labels 
            subup<- subset(subdata, subdata$logFC>0);dim(subup)
            if(nrow(subup)>0){
               points(x=subup$AveExpr,y=subup$logFC,pch=21, cex=1.2,
                      bg="orange", col="darkorange3",lwd=0.5,
                      text(x=subup$AveExpr, y=subup$logFC, col="orange",
                           labels=subup$labels, pos=4,cex=0.6,font=2))
            }
            
            ## ADD USER INPUT LABELS TO NEG logFC < 0
            subdn<- subset(subdata, subdata$logFC<0);dim(subdn)
            if(nrow(subdn>0)){
               points(x=subdn$AveExpr,y=subdn$logFC,pch=21, cex=1.2,
                      bg="orange", col="darkorange3",lwd=0.5,
                      text(x=subdn$AveExpr, y=subdn$logFC, col="orange",
                           labels=subdn$labels, pos=2,cex=0.6,font=2))
            }
         } ## LEN >0 
      } ## LABELS
      
      
      ## TOP
      
      ## the user can use 'top' to label the most sig. up and down regulated
      ## features (proteins=gene name, phospho=gene name(site).
      ## e.g. top=5, will label the top 5 up and down reg.
      if(!is.null(top)){ ## TOP
         if(length(top)==1 & top>0){
            ## sig. up-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC Z-> A
            sig <- up[with(up,order(logFC,order(p,decreasing=FALSE),decreasing=TRUE)),]
            top2<-ifelse(top > nrow(sig), nrow(sig),top);top2
            sig<- sig[1:top2, ]
            points(x=sig$AveExpr,y=sig$logFC,pch=21, cex=1.2,
                   bg="green", col="forestgreen",lwd=0.5,
                   text(x=sig$AveExpr, y=sig$logFC, col="green3",
                        labels=sig$labels, pos=4,cex=0.6,font=2))
            
            ## sig. down-reg proteins/sites are sorted by pval (A->Z), then sorted by logFC A->Z
            sig <- down[with(down,order(logFC,order(p,decreasing=FALSE),decreasing=FALSE)),]
            top2<-ifelse(top > nrow(sig), nrow(sig),top);top2
            sig<- sig[1:top2, ]
            points(x=sig$AveExpr,y=sig$logFC,pch=21, cex=1.2,
                   bg="green", col="forestgreen",lwd=0.5,
                   text(x=sig$AveExpr, y=sig$logFC, col="green3",
                        labels=sig$labels, pos=4,cex=0.6,font=2))   
            
         }
      } ## TOP
      
      
      if(legend==TRUE){
         legend("topright", c("up", "not sig.", "down"), box.col="transparent",
                box.lwd=0, bg="transparent", border="black", pch=22,
                pt.bg=c("#FF0000","grey50","#1E90FF"), pt.cex=1.4,
                col=c("#FF0000","grey50","#1E90FF"), horiz=FALSE,
                inset=c(inset, 0), bty="n", xpd=TRUE, ncol=1, cex=0.9)
      }
      
      
   } ## MD PLOT
   
   
   ##-------------------------------
   ##  [03]  GLIMMA VOLCANO PLOT
   ##-------------------------------
   glimmaVolcanoPlot <- function(stats, comparison, data, annot, groups, min.pval=0.055, min.lfc=1, 
                                 sig.type=c("p.adj","pval"), pipe="DIA", 
                                 enrich="protein", dir="."){
      
      pipe   <- match.arg(arg=pipe, choices=c("DIA", "TMT","phosphoTMT","LF"),several.ok=FALSE);pipe
      enrich <- match.arg(arg=enrich, choices=c("protein","phospho"),several.ok=FALSE);enrich
      # if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)}
      
      stats <- data.frame(stats);dim(stats)
      
      ## column named 'p' added to stats data.frame based on sig.type selected. 
      ## text for y-axis label also defined
      if(sig.type=="p.adj"){
         stats$p <- stats$adj.P.Val
         ylab=expression(paste("-",log[10]," (adj. p-value)",sep=""))
         filename<-paste(comparison,"volcano_plot",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      if(sig.type=="pval"){
         stats$p<-stats$P.Value
         ylab=expression(paste("-",log[10]," (p-value)",sep=""))
         filename<-paste(comparison,"volcano_plot_pvalue",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      
      ## remove rows with adj.PValue==NA or P.Value==NA
      stats <- stats[!is.na(stats$p),];dim(stats)
      
      
      ## sort data, annotation to match order in stats
      ## replace NAs in data with zeros (visualization purposes)
      ## make sure the number of groups entered matches the number of sample columns in data
      if(all(rownames(stats) %in% rownames(data))){ 
         data <- data[rownames(stats), ]
         data[is.na(data)] <- 0
      } else { stop("Error! data and stats do not match. check rownames.") }
      if(all(rownames(stats) %in% rownames(annot))){ 
         annot <- annot[rownames(stats), ]
      } else { stop("Error! annot and stats do not match. check rownames.") }
      if(all(rownames(annot)==rownames(data))==FALSE){ 
         stop("Error! Rownames of data and annotation do not match.")}
      if(length(groups) != ncol(data)){stop("Error! The number of groups does not equal the number of data columns.")}
      
      
      ## add status column to stats (determines colors of significant dots in plot)
      stats$status <- rep(0,nrow(stats))
      stats$status[(stats$p <= min.pval & stats$logFC >= min.lfc)==TRUE] <- 1
      stats$status[(stats$p <= min.pval & stats$logFC <= -min.lfc)==TRUE] <- -1
      
      
      ## DETERMINE ANNOTATION COLUMNS OF INTEREST
      if(pipe=="DIA" & enrich=="protein"){ 
         glimmaAnnotationColums <- c("id", "UniprotID","Gene_name","Description") }
      if((pipe=="TMT" | pipe=="LF") & enrich=="protein"){ 
         glimmaAnnotationColums <- c("id", "UniprotID","Gene_name","Description") }
      if(pipe=="phosphoTMT" & enrich=="protein"){ 
         glimmaAnnotationColums <- c("id","UniprotID","Gene_name","Description") }
      if(pipe=="phosphoTMT" & enrich=="phospho"){ 
         glimmaAnnotationColums <- c("id", "proGroupID", "UniprotID", "Gene_name", "Description","phosAAPosition") }
      # if(any(glimmaAnnotationColums %in% colnames(annot))==FALSE){ 
      #    stop("Error! glimmaAnnotationColums are not all present in annot.")}
      if(all(glimmaAnnotationColums %in% colnames(annot))){ 
         glimmaAnnot <- annot[, glimmaAnnotationColums]
      } else { stop("Error! glimmaAnnotationColums are not all present in annotation.")}
      
      glXYPlot(x         = stats$logFC,
               y         = -log(stats$p,10),
               counts    = data[rownames(stats), ],
               groups    = groups,
               samples   = colnames(data),
               status    = stats$status, 
               transform = FALSE, 
               anno      = glimmaAnnot[rownames(stats), ],
               display.columns = colnames(glimmaAnnot),
               xlab      = "logFC",
               # ylab      = paste0(" ","-","log10 (",ifelse(sig.type=="p.adj","adj. p-value","p-value"),")"),
               ylab = paste0("negLog10(",ifelse(sig.type=="p.adj","adj.P","P"),")"),
               # ylab      = paste("neg log10",sig.type),
               side.main = "Gene_name",
               side.xlab = "Group",
               side.ylab = "Normalized Intensity",
               side.log  = FALSE,
               sample.cols = colorGroup2(groups)[groups],
               cols      = c("#00bfff", "#858585", "#ff3030"),
               jitter    = 10,
               path      = file.path(dir),
               folder    = "Volcano-Plots",
               html      = filename,
               main      = paste0(comparison," (",ifelse(sig.type=="p.adj","adj. p-value","p-value"),")"),
               launch    = FALSE)
      
      
      fixfile <- file.path(dir,"Volcano-Plots/js",paste0(filename,".js"));print(fixfile)
      tmptxt <- readLines(fixfile)
      tmptxt <- gsub("logCPM","Intensity",tmptxt)
      writeLines(as.character(tmptxt), file(fixfile), sep="\n")
      
      
   } ## GLIMMA VOLCANO PLOT
   
   
   ##--------------------------
   ##  [04]  GLIMMA MD PLOT
   ##--------------------------
   glimmaMDPlot <- function(stats, comparison, data, annot, groups, min.pval=0.05, min.lfc=1, 
                            sig.type=c("p.adj","pval"), pipe="DIA", 
                            enrich="protein", dir="."){
      
      pipe   <- match.arg(arg=pipe, choices=c("DIA", "TMT","phosphoTMT","LF"),several.ok=FALSE);pipe
      enrich <- match.arg(arg=enrich, choices=c("protein","phospho"),several.ok=FALSE);enrich
      # if(!dir.exists(dir)){dir.create(dir,recursive=TRUE)}
      
      stats <- data.frame(stats);dim(stats)
      
      ## column named 'p' added to stats data.frame based on sig.type selected. 
      ## text for y-axis label also defined
      if(sig.type=="p.adj"){
         stats$p <- stats$adj.P.Val
         filename<-paste(comparison,"MD_plot",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      if(sig.type=="pval"){
         stats$p<-stats$P.Value
         filename<-paste(comparison,"MD_plot_pvalue",sep=ifelse(is.null(comparison) | comparison=="","","_"))
      }
      
      ## remove rows with adj.PValue==NA or P.Value==NA
      stats <- stats[!is.na(stats$p),];dim(stats)
      
      
      ## sort data, annotation to match order in stats
      ## replace NAs in data with zeros (visualization purposes)
      ## make sure the number of groups entered matches the number of sample columns in data
      if(all(rownames(stats) %in% rownames(data))){ 
         data <- data[rownames(stats), ]
         data[is.na(data)] <- 0
      } else { stop("Error! data and stats do not match. check rownames.") }
      if(all(rownames(stats) %in% rownames(annot))){ 
         annot <- annot[rownames(stats), ]
      } else { stop("Error! annot and stats do not match. check rownames.") }
      if(all(rownames(annot)==rownames(data))==FALSE){
         stop("Error! Rownames of data and annotation do not match.")}
      if(length(groups) != ncol(data)){
         stop("Error! The number of groups does not equal the number of data columns.")}
      
      
      ## add status column to stats (determines colors of significant dots in plot)
      stats$status <- rep(0,nrow(stats))
      stats$status[(stats$p <= min.pval & stats$logFC >= min.lfc)==TRUE] <- 1
      stats$status[(stats$p <= min.pval & stats$logFC <= -min.lfc)==TRUE] <- -1
      
      
      ## DETERMINE ANNOTATION COLUMNS OF INTEREST
      if(pipe=="DIA" & enrich=="protein"){ 
         glimmaAnnotationColums <- c("id", "UniprotID","Gene_name","Description") }
      if((pipe=="TMT" | pipe=="LF") & enrich=="protein"){ 
         glimmaAnnotationColums <- c("id", "UniprotID","Gene_name","Description") }
      if(pipe=="phosphoTMT" & enrich=="protein"){ 
         glimmaAnnotationColums <- c("id","UniprotID","Gene_name","Description") }
      if(pipe=="phosphoTMT" & enrich=="phospho"){ 
         glimmaAnnotationColums <- c("id", "proGroupID", "UniprotID", "Gene_name", "Description","phosAAPosition") }
      if(all(glimmaAnnotationColums %in% colnames(annot))){ 
         glimmaAnnot <- annot[, glimmaAnnotationColums]
      } else { stop("Error! glimmaAnnotationColums are not all present in annotation.")}
      
      
      glXYPlot(x         = stats$AveExpr,
               y         = stats$logFC,
               counts    = data[rownames(stats), ],
               groups    = groups,
               samples   = colnames(data),
               status    = stats$status, 
               transform = FALSE, 
               anno      = glimmaAnnot[rownames(stats), ],
               display.columns = colnames(glimmaAnnot),
               xlab      = "AveExpr",
               ylab      = "logFC",
               # ylab      = paste("neg log10",sig.type),
               side.main = "Gene_name",
               side.xlab = "Group",
               side.ylab = "Normalized Intensity",
               side.log  = FALSE,
               sample.cols = colorGroup2(groups)[groups],
               cols      = c("#00bfff", "#858585", "#ff3030"),
               jitter    = 10,
               path      = file.path(dir),
               folder    = "MD-Plots",
               html      = filename,
               main      = paste0(comparison," (",ifelse(sig.type=="p.adj","adj. p-value","p-value"),")"),
               launch    = FALSE)
      
      fixfile <- file.path(dir,"MD-Plots/js",paste0(filename,".js"));print(fixfile)
      tmptxt <- readLines(fixfile)
      tmptxt <- gsub("logCPM","Intensity",tmptxt)
      writeLines(as.character(tmptxt), file(fixfile), sep="\n")

      
   } ## GLIMMA MD PLOT
   
   
   
   ##--------------------------
   ##  [05] pvalueHistogram     
   ##--------------------------
   pvalueHistogram <- function(stats, comparison=NULL){
      
      stats <- data.frame(stats);dim(stats)
      
      layout(matrix(1:2, ncol=2, byrow=TRUE))
      ##   P-VALUE HISTOGRAM
      ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6164648/#B1-high-throughput-07-00023
      ## inspect histogram of unadjusted P.Values for the differential expression results.  
      ## Figure shows histogram of p-values from gene-by-gene statistical tests.
      par(mar=c(6,7,5,5))
      graphics::hist(stats$P.Value, col="dodgerblue", breaks=50, 
                     main="", xlab="", ylab="", 
                     xlim=c(0,1), las=1,cex.axis=1.2)
      title(line=1, main=comparison, cex.main=1.3)
      mtext(side=1, line=3, cex=1.3, text="p-value")
      mtext(side=2, line=4,cex=1.3, text="Frequency")
      
      ##   ADJUSTED P-VALUE HISTOGRAM
      par(mar=c(6,7,5,5))
      graphics::hist(stats$adj.P.Val, col="dodgerblue", breaks=50, 
                     main="", xlab="", ylab="", 
                     xlim=c(0,1), las=1,cex.axis=1.2)
      title(line=1, main=comparison, cex.main=1.3)
      mtext(side=1, line=3, cex=1.3, text="adj. p-value")
      mtext(side=2, line=4,cex=1.3, text="Frequency")
      
      
   } ## PVALUE HISTOGRAM
   
   
} ## DE PLOT FUNCTIONS











