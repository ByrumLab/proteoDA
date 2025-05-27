### proteoDA parameters
## April 25, 2025

#working_dir ="/home/sbyrum/tests/proteoda2_test4"
working_dir = getwd()

#### Variables for METHODS and PowerPoint REPORT
author = "Vishwajeeth R. Pagala"
author2 = "Stephanie Byrum"
project_name = "Shelat_010725"
instrument = "Bruker timsTOF HT"
search_engine = "DIA-NN (version 2.1.0)"
database = "UniprotKB (57352 proteins)"
db_version = "version 2024/11/27"
organism = "Mus musculus"

### input files: subset quan and protein annotation
input_quan = "quan_results_20250508/publication/uni_prot_quan_rmNA_norm.csv"
## metadata requires column names "sample" and "group"
metadata = "data/Shelatgrp_010725_metadata2.csv"
sample_start = 11  # first column where the samples start, currently 10 for diann, different for spectronaut
anno_start = 1
anno_end = 10

### filter proteins with missing values options
### rule of thumb, require 2/3 of reps to have values 
filt_min_reps = 2           # number of biological replicates that must have a value > 0
filt_min_groups = 1         # number of sample groups that must have the biological reps with values > 0 
group = "group"             # defines the sample metadata column containing group information
require_both_groups = FALSE  # if filtering_per_contrasts, are the biological reps with values > 0 required in both groups

## report directories. For now don't change these or you must update the Project_Summary.qmd
QC_dir = "figures"               # can change this in the GenerateReports.R so doesn't mess up Project_Summary?
DA_dir = "interactive_results"

### limma design and group comparisons
design = ~0 + group # limma no intercept model, comparisons based on contrasts
#design = ~ group   # limma intercept model, everything will be compared to a reference
#design = ~0 + group + (1 |animal)  # Mixed model (group is fix and animal pair factor treated as random factor)
contrasts = "data/contrasts.csv"   # file with sample group comparisons, order of groups determines fold change

# significance thresholds for plotting and tables
bin_size = 500  # depends on number of proteins in data, default is 1000
p.val = 0.05    # significance threshold for plots and tables, default is 0.05
logFC = 1       # significance threshold for plots and tables, default is 1 (2-fold)
stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSDs", "logFC_z_scores", "sig.PVal", "sig.FDR") 

## column names for interactive tables and plots - must match Protein ANNOTATION!!!!!!!
DA_table_cols <- c("uniprot_id","Accession.Number","Protein.Description")
DA_title_col <- "uniprot_id" 
tmp_subdir <- "tmp"

#### QC plot options - Must match sample Metadata columns!!!!!!!!!
plot_labels = "sample"
barplot_grouping_columns <- c("group")
pca_grouping_columns <- c("group", "sample")
den_grouping_columns <- c("group", "sample")
