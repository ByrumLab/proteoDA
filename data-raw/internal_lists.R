######################################
## Lists of column names and colors ##
## that are internally used by      ##
## other functions                  ##
######################################


## DIA: protein annotation/contaminantion columns from Scaffold DIA Sample Report csv file
diaAnnotationColums <- c("id", "Protein.Name","Accession.Number","Molecular.Weight",
                         "Protein.Group.Score", "Identified.Peptide.Count",
                         "Exclusivity","Quantitative.Value")
diaContamColums     <- c("Protein.Name") ## Group by / DECOY_

## METADATA COLUMNS
diaMetaColums <- c("sample", "number","group")

## named lists of normalization and imputation methods
norm.methods        <- c("log2", "median", "mean", "vsn","quantile", "cycloess", "rlr", "gi")
names(norm.methods) <- c("Log2", "Median", "Mean", "VSN","Quantile", "Cyclic Loess", "RLR", "Global Intensity")


## create a color palette for plots called binfcolors
## color names courtesy of my nieces and nephews (ages 6, 9, and 12)
binfcolors <- c("#1F5EDC","#EE0010","#32CD32","#FF1493","#FF7F00",
                "#A342FC","#00C8FF","#ADFF2F","#FFE100","#E36EF6","#009ACE","#996633")
names(binfcolors) <- c("blueberry","cherry","apple","barbie","fanta",
                       "grape","ocean","mtndew","gold","orchid","aceblue","poop")


lightbinfcolors <- c("#c8d8f7","#ffdadd","#d0f3d0","#ffb1db","#ffe1c4",
                     "#dbb6fe","#c4f2ff","#e3ffb8","#fff8c4","#f1b8fb",
                     "#baeeff","#eddcca")
names(lightbinfcolors) <- c("lightblueberry","lightcherry","lightapple","lightbarbie","lightfanta",
                            "lightgrape","lightocean","lightmtndew","lightgold","lightorchid",
                            "lightaceblue","lightpoop")


usethis::use_data(x = binfcolors, lightbinfcolors,
                  diaAnnotationColums, diaContamColums, diaMetaColums,
                  norm.methods,
                  internal = T, overwrite = T)
