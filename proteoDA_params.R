### proteoDA parameters
## March 4, 2025

working_dir ="/research/rgs01/home/clusterHome/sbyrum/projects/proteoda_test/proteoDA"

#### Variables for METHODS REPORT
author = "Vishwajeeth R. Pagala"
author2 = "Zuo-Fei Yuan"
project_name = "PI_DATE"
instrument = "Bruker timsTOF HT"
SearchEngine = "DIA-NN (version 1.8.2)"
database = "Proteome ID: UP000005640, 81,791 proteins"
db_version = "version 2023/04/28"

### input files and subset quan and protein annotation 
input_quan = "data/uni_prot_quan_rmNA_norm_proteoDA.csv"
## metadata requires column names "sample" and "group" 
metadata = "data/Sample_metadata.csv"
sample_start = 10
anno_start = 1
anno_end = 9

### filter proteins with missing values options
filt_min_reps = 2
filt_min_groups = 1
group = "group"

## report directories
QC_dir = "01_QC"
DA_dir = "02_DA"

### limma design and group comparisons 
design = ~0 + group
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
pca_grouping_columns <- c("group","sample")
den_grouping_columns <- c("group", "sample")
