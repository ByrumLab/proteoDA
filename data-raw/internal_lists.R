######################################
## Lists of column names and colors ##
## that are internally used by      ##
## other functions                  ##
######################################

##Not wholly convinced this is the best way to do it, but its a start.

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


usethis::use_data(x = binfcolors,
                  diaAnnotationColums, diaContamColums, diaMetaColums,
                  impute.methods,
                  lfMetaColums,
                  lightbinfcolors,
                  limmaStatColums,
                  norm.methods,
                  peptideAnnotationColums, peptideContamColums,
                  phosphoAnnotationColums, phosphoContamColums, phosphoLocalProbColum,
                  proteinAnnotationColums, proteinContamColums,
                  targetColums,
                  TMT10plex, TMT11plex, TMT16plex, TMT6plex,
                  tmtAnnotationColums, tmtContamColums, tmtMetaColums,
                  internal = T, overwrite = T)
