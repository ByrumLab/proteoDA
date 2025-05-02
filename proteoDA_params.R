### proteoDA parameters
## February 25, 2025

working_dir ="/home/sbyrum/projects/Demongrp_032525/proteoda"

#### Variables for METHODS REPORT
author = "Vishwajeeth R. Pagala"
author2 = "Stephanie Byrum"
project_name = "Demontisgrp_032525"
instrument = "Bruker timsTOF HT"
search_engine = "DIA-NN (version 2.1.0)"
database = "UniprotKB (sp_tr), 57352 proteins"
db_version = "version 2024/11/27"
organism = "Mus musculus"

### input files and subset quan and protein annotation
input_quan = "data/uni_prot_quan_rmNA_norm_proteoda.csv"
## metadata requires column names "sample" and "group"
metadata = "data/Demontisgrp_032525_Sample_Metadata.csv"
sample_start = 10
anno_start = 1
anno_end = 9

### filter proteins with missing values options
filt_min_reps = 3
filt_min_groups = 2
require_both_groups = TRUE
group = "group"

## report directories
QC_dir = "figures"
DA_dir = "interactive_results"

### limma design and group comparisons
# design = ~0 + group # no intercept model
design = ~0 + group # animal is fixed
#design = ~0 + group + (1 |animal)  # Mixed model (group is fix and animal random factor)
contrasts = "data/contrasts.csv"

# significance thresholds for plotting and tables
p.val = 0.05
logFC = 1
## column names for interactive tables and plots - must match Protein ANNOTATION!!!!!!!
DA_table_cols <- c("uniprot_id","Accession.Number","Protein.Description")
DA_title_col <- "uniprot_id"
tmp_subdir <- "tmp"

#### QC plot options - Must match sample Metadata columns!!!!!!!!!
plot_labels = "sample"
barplot_grouping_columns <- c("group")
pca_grouping_columns <- c("group", "antibody","age")
den_grouping_columns <- c("group", "antibody", "age")
