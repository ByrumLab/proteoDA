### proteoDA parameters
## April 25, 2025

#working_dir ="/home/sbyrum/tests/proteoda2_test4"
working_dir = getwd()

#### Variables for METHODS and PowerPoint REPORT
author = "Vishwajeeth Pagala"
author2 = "Stephanie Byrum and Zuo-Fei Yuan"
project_name = "Durbigrp_071825_MMDIA_SNE"
instrument = "Bruker timsTOF HT"
search_engine = "Spectronaut (version 19)"
database = "uniprot_human_reviewedProteome_20407_20230306 (20407 proteins)"
db_version = "March 2023"
organism = "Homo sapiens"

### input files: subset quan and protein annotation
input_quan = "data/20250910_Durbigrp_071825_Protein_proteoDA_input_test2.csv"
## metadata requires column names "sample" and "group"
metadata = "data/Durbigrp_071825_MMDIA_SampleMetadata.csv"
sample_start = 14  # first column where the samples start, currently 10 for diann, different for spectronaut
anno_start = 1
anno_end = 13

## contrasts_cell adds the interaction effect
contrasts = "data/contrasts_groupmeans.csv"   # file with sample group comparisons, order of groups determines fold change


### filter proteins with missing values options
### rule of thumb, require 2/3 of reps to have values 
filt_min_reps = 2           # number of biological replicates that must have a value > 0
filt_min_groups = 1         # number of sample groups that must have the biological reps with values > 0 
group = "group"             # defines the sample metadata column containing group information
require_both_groups = FALSE  # if filtering_per_contrasts, are the biological reps with values > 0 required in both groups

## report directories. For now don't change these or you must update the Project_Summary.qmd
QC_dir = "figures"               # can change this in the GenerateReports.R so doesn't mess up Project_Summary?
DA_dir = "interactive_results"

######################
# Modeling the four groups (~ 0 + cell_line:treatment) lets you:
#   Remove (control for) each cell line’s DMSO baseline explicitly.
#   Test treatment effects within a cell line (bio vs DMSO).
#   Test whether treatment effects differ between cell lines (difference-of-differences), which is usually the biological question.
#   Adding (1 | batch) (or similar) improves variance estimation when you have repeated measures/blocks, via duplicateCorrelation downstream.
##########################

### limma design and group comparisons
# limma no intercept model, comparisons based on contrasts correcting for cell line batch
design = ~0 + group 

# interaction model to answer "what is treatment effect in A after removing its DMSO baseline,
# compared to treatment effect in B after removing its DMSO baseline? 
#design = ~0 + cell:treatment   

### other examples ---------------
# Random batch, while still modeling the 4 groups explicitly
# design = ~ 0 + cell_line:treatment + (1 | batch)
# design = ~ group   # limma intercept model, everything will be compared to a reference
#design = ~0 + group + (1 |animal)  # Mixed model (group is fix and animal pair factor treated as random factor)

# significance thresholds for plotting and tables
bin_size = 500  # depends on number of proteins in data, default is 1000
p.val = 0.05    # significance threshold for plots and tables, default is 0.05
logFC = 1       # significance threshold for plots and tables, default is 1 (2-fold)
stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSDs", "logFC_z_scores", "sig.PVal", "sig.FDR") 

## column names for interactive tables and plots - must match Protein ANNOTATION!!!!!!!
#DA_table_cols <- c("uniprot_id","Accession.Number","Protein.Description")  # DIANN
DA_table_cols <- c("uniprot_id","PG.ProteinLabel","PG.ProteinDescriptions")  # Spectronaut
DA_title_col <- "uniprot_id" 
tmp_subdir <- "tmp"
ctrl_proteins <-  c( "P46100")        # c( "P58004", "Q12766", "MDM2") highlighted in volcano static plots/powerpoint

#### QC plot options - Must match sample Metadata columns!!!!!!!!!
plot_labels = "sample"
barplot_grouping_columns <- c("group")
pca_grouping_columns <- c("group", "sample","cell")
den_grouping_columns <- c("group", "sample","cell")
