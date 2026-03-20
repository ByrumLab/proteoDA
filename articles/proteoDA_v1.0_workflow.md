# proteoDA_v1.0_workflow

## Installation

To install proteoDA v1.0 from GitHub (with dependencies and vignettes
built):

``` r
# Install proteoDA version 1.0.0
# install.packages("remotes")
remotes::install_github("ByrumLab/proteoDA@v1.0.0",
                        dependencies = TRUE,
                        build_vignettes = TRUE,
                        force = TRUE)
```

### Load the package and view the available vignettes:

``` r
library(proteoDA)
# browseVignettes(package = "proteoDA")
```

## Input File Requirements and example

`proteoDA` ships with small DIA-style test datasets for quick
evaluation. Three data matrix are required as input:

1.  protein intensity data - proteins are in rows and samples are in
    columns
2.  sample metadata - links sample columns in the protein intensity to
    sample group information
3.  protein annotation - includes `uniprot_id` required column, protein
    description, gene names, and any other protein information.

*NOTE:* 1. The protein intensity dataframe must include a column labeled
`uniprot_id`. This becomes the rownames for the lists in R and must have
unique values for each row. 2. The sample metadata csv file should have
a column labeled `sample` that matches the colnames of the samples in
the protein intensity dataframe.

``` r
input_data <- read.csv(system.file("extdata/DIA_data.csv.gz", package = "proteoDA"))
sample_metadata <- read.csv(system.file("extdata/metafile.csv", package = "proteoDA"))

head(input_data)
```

    ##   uniprot_id
    ## 1 A0A087X1C5
    ## 2 A0A0B4J1Y9
    ## 3 A0A0C4DH29
    ## 4     A0AVT1
    ## 5     A0FGR8
    ## 6     A0MZ66
    ##                                                                                                 protein_info
    ## 1         sp|A0A087X1C5|CP2D7_HUMAN Putative cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=5 SV=1
    ## 2 sp|A0A0B4J1Y9|HV372_HUMAN Immunoglobulin heavy variable 3-72 OS=Homo sapiens OX=9606 GN=IGHV3-72 PE=3 SV=1
    ## 3   sp|A0A0C4DH29|HV103_HUMAN Immunoglobulin heavy variable 1-3 OS=Homo sapiens OX=9606 GN=IGHV1-3 PE=3 SV=1
    ## 4 sp|A0AVT1|UBA6_HUMAN Ubiquitin-like modifier-activating enzyme 6 OS=Homo sapiens OX=9606 GN=UBA6 PE=1 SV=1
    ## 5                   sp|A0FGR8|ESYT2_HUMAN Extended sy0ptotagmin-2 OS=Homo sapiens OX=9606 GN=ESYT2 PE=1 SV=1
    ## 6                                 sp|A0MZ66|SHOT1_HUMAN Shootin-1 OS=Homo sapiens OX=9606 GN=SHTN1 PE=1 SV=4
    ##            accession_number molecular_weight Pool_1.mzML Pool_2.mzML
    ## 1 sp|A0A087X1C5|CP2D7_HUMAN           57 kDa    12128420    10100400
    ## 2 sp|A0A0B4J1Y9|HV372_HUMAN           13 kDa    34858280    37054710
    ## 3 sp|A0A0C4DH29|HV103_HUMAN           13 kDa    83996110    99612510
    ## 4      sp|A0AVT1|UBA6_HUMAN          118 kDa    70885410    77754340
    ## 5     sp|A0FGR8|ESYT2_HUMAN          102 kDa    25028960    21329540
    ## 6     sp|A0MZ66|SHOT1_HUMAN           72 kDa    10811940    11644230
    ##   Pool_3.mzML Sample_01_0325.mzML Sample_02_0147.mzML Sample_03_4085.mzML
    ## 1     8047184             6400703             9803543             2147631
    ## 2    39324880            88847090            40182130            39822860
    ## 3   107389100            42280640           179020600           194450400
    ## 4    76753210            73468270           119136300            45476380
    ## 5    26742190            59603060            18563650             6483842
    ## 6    10496870            20095690            11457240             5871754
    ##   Sample_04_0039.mzML Sample_05_5198.mzML Sample_06_7066.mzML
    ## 1              767043             7872038             2384339
    ## 2            39975930            42115980            24683670
    ## 3            49979340            39133770            38161670
    ## 4            69640920            75263440           103112300
    ## 5            14402970            15958740            35201970
    ## 6            30741410            18112560             9614461
    ##   Sample_07_0909.mzML Sample_08_0348.mzML Sample_09_7028.mzML
    ## 1            14989790            17759460             9190121
    ## 2            23234440             1862563            28740020
    ## 3            24654850            16704030            51880990
    ## 4            71260680            60429110            83455310
    ## 5            19612370            77828310            36625060
    ## 6             6803437            10232220             5312412
    ##   Sample_10_2119.mzML Sample_11_6501.mzML Sample_12_9999.mzML
    ## 1             1747821              734365             1101549
    ## 2             9271566            73529600            40445100
    ## 3            22337610           418867800           110062200
    ## 4            27706830           104025000            55297540
    ## 5             5416264            17541150            11326900
    ## 6            14760770            15365450             9556587
    ##   Sample_13_1014.mzML Sample_14_1059.mzML
    ## 1            12225410            11632510
    ## 2            35323730            14474530
    ## 3            74030230            72548530
    ## 4            67913170            87962450
    ## 5            22483620            18911320
    ## 6            16628400            19235090

``` r
head(sample_metadata)
```

    ##      data_column_name sample batch  group
    ## 1         Pool_1.mzML     P1     1   Pool
    ## 2         Pool_2.mzML     P2     1   Pool
    ## 3         Pool_3.mzML     P3     1   Pool
    ## 4 Sample_01_0325.mzML   S325     1 normal
    ## 5 Sample_02_0147.mzML   S147     1 normal
    ## 6 Sample_03_4085.mzML  S4085     1 normal

### Preparing the input files

Most of the protein data files generated from various database search
engines will include the protein annotation along with the protein
intensities for each sample in the same file. In order to create the
`DAList` in `proteoDA`, we need to split the database search results
into the protein annotation and intensity data. Only the columns needed
for analysis should be included.

``` r
annotation_data <- input_data[, 1:4]            # protein-level annotations
intensity_data  <- input_data[, 5:ncol(input_data)]  # quantitative sample columns

# Ensure metadata rownames match intensity column names
rownames(sample_metadata) <- sample_metadata$data_column_name
```

## Create a `DAList` Object

`DAList` is the central container class in proteoDA for holding
intensity, metadata, and downstream analysis results.

``` r
raw <- DAList(
data      = intensity_data,
annotation= annotation_data,
metadata  = sample_metadata
)
```

## Preprocessing

### Filter Samples

The
[`filter_samples()`](https://byrumlab.github.io/proteoDA/reference/filter_samples.md)
allows you to remove any outlier samples, pooled or QC samples not used
in limma group comparisons.

``` r
# samples with the group named "Pool" will be removed 
filtered_samples <- filter_samples(raw, group != "Pool")
```

### How proteoDA v1.0 deals with missing values

Depending on your experimental design, you may wish to treat intensity
values below the detection limit of your instrument as missing data (NA)
or as a value of 0 intensity (protein is truly knocked out). Though this
can be done manually before importing data into a DAList object, we also
provide some simple utility functions to convert missing values to 0 and
vice-verse.

[`missing_to_zero()`](https://byrumlab.github.io/proteoDA/reference/missing_data.md) -
converts missing values (NA by default, or other user-provided values)
to 0 - is mainly for downstream export/compatibility for tools that
require a fully numeric matrix, not for model fitting, such as
performing PCA analysis

[`zero_to_missing()`](https://byrumlab.github.io/proteoDA/reference/missing_data.md) -
converts 0s to NA - limma treats them as missing rather than true
measurements. - NA values are required for
[`filter_proteins_by_group()`](https://byrumlab.github.io/proteoDA/reference/filter_proteins_by_group.md)
and for the missing value heatmap

#### How limma treats NA values

The way NA vs 0 is handled is exactly how limma does it where NA values
in the expression matrix are treated as missing.

For each protein/feature, lmFit(): 1. drops samples that are NA for that
feature only, and 2. fits the linear model using only the remaining
non-missing samples.

This means: - Different proteins can be fit with different numbers of
samples. - If a whole group (e.g. all “cancer” samples) is NA for a
protein, the corresponding design column / contrast becomes
non-estimable and `limma` returns NA for logFC, t-statistic, p-value for
that protein/contrast. - `eBayes()` uses the available residual degrees
of freedom per protein, so proteins with many missing values borrow
information but still have appropriately smaller df.

#### How limma treats zeros

- A value of 0 is just a number to limma.

- On the log2 scale, 0 means “intensity of 1”; on raw scale it means
  “exactly zero”.

- `lmFit()` and `eBayes()` assume these zeros are real observed values
  and

1.  includes them in the mean and variance calculations,
2.  allows them to contribute to fold-changes and t-statistics.

If zeros actually mean “not detected / missing” rather than a true
measured abundance, treating them as real values can:

1.  Pull group means strongly downward,
2.  Inflate log2 fold-changes when one group has many zeros and the
    other doesn’t,
3.  Artificially increase within-group variance for “mixed” features
    (some real intensities, some zeros),
4.  Distort the moderated variance estimates in eBayes().

*NOTE:* NA = “this measurement is not available; ignore it for this
protein’s model fit”. limma never tries to impute internally.

`proteoDA v1.0` does not include an imputation method; however, version
2.0 now includes imputation using the normal distribution as performed
in the Perseus Software (Cox lab).

*Recommendataion:* Convert explicit zeros to NA (recommended before
filtering):

``` r
filtered_samples <- zero_to_missing(filtered_samples)
```

### Filter Proteins by Group

Filter features based on replicate presence (e.g., ≥4 values in ≥1
group):

This function is used to remove proteins from a DAList, filtering out
proteins based on levels of missing values in the “data” data frame of
the `DAList`. The `grouping_column` must be a column in the metadata of
the `DAList`, which lists the group membership for each sample. The
`min_reps` and `min_groups` arguments that determine the number of
replicates/samples per group (`min_reps`) and number of groups
(`min_groups`) in which a protein must have non-missing intensity values
in order to be retained. This function assumes that all missing values
are encoded as NA.

*Rule of Thumb* require at least 2/3 of replicates to have values in at
least one sample group.

``` r
filtered_proteins <- filter_proteins_by_group(
filtered_samples,
min_reps        = 4,
min_groups      = 1,
grouping_column = "group"
)
```

## Normalization and QC

[`write_norm_report()`](https://byrumlab.github.io/proteoDA/reference/write_norm_report.md)
saves a PDF report containing a variety of qualitative plots which give
information about the performance of different normalization methods.
The report is useful for choosing a normalization method to use for
downstream analysis.

### Evaluating Normalization Methods with write_norm_report()

Proteomics intensity data almost always contain systematic,
non-biological biases—arising from sample handling, instrument drift,
batch effects, or other uncontrolled sources. As emphasized in the
`proteiNorm` publication (Graw et al), selecting an appropriate
normalization strategy is a critical preprocessing step as poor
normalization can obscure true biological changes, inflate variance, or
introduce artificial fold-changes.

To support this evaluation in `proteoDA`, the function
[`write_norm_report()`](https://byrumlab.github.io/proteoDA/reference/write_norm_report.md)
provides an automated, reproducible way to compare the performance of
eight widely used normalization methods. Internally, this wrapper runs
the same QC, normalization, and diagnostic steps implemented in the
`proteiNorm` framework, but produces a static PDF or HTML report
suitable for documentation and publication.

#### What normalization methods are evaluated?

[`write_norm_report()`](https://byrumlab.github.io/proteoDA/reference/write_norm_report.md)
applies the following methods to the protein abundance matrix

1.  Log2 transformation (log2)
2.  Median normalization (median)
3.  Mean normalization (mean)
4.  Variance Stabilizing Normalization (vsn)
5.  Quantile normalization (quantile)
6.  Cyclic loess normalization (cycloess)
7.  Robust linear regression (rlr)
8.  Global intensity scaling (gi)

These methods represent a combination of microarray-derived
normalization strategies and proteomics-specific regressions/scaling
approaches, many of which consistently rank among the top performers in
systematic evaluations of label-free proteomics data.

#### What diagnostics does the report generate?

Following the `proteiNorm` workflow, the `proteoDA` report evaluates
each method with a consistent QC panel:

1.  *Total intensity boxplots*: Checks whether global sample intensities
    become comparable after normalization.
2.  *PCA plots*: Visualizes grouping structure, batch effects, and
    separation of biological conditions.
3.  *Missing value heatmap*: Displays patterns of missingness (MAR vs
    MNAR) as clustered heatmap.
4.  *Intragroup variation metrics*

For each normalization method, the report quantifies: - PCV (pooled
coefficient of variation) - PEV (pooled estimate of variance) - PMAD
(pooled median absolute deviation)

These metrics are key evaluation criteria in `proteiNorm` because lower
intragroup variation indicates a more stable, better-performing
normalization method.

5.  *Correlation metrics*:

- Boxplots of within-group correlations
- Full Pearson correlation heatmaps

Higher intragroup correlation suggests that biological replicates align
after normalization.

6.  *Log2 ratio distributions*: Distribution of all pairwise log2
    fold-changes between sample groups, expected to center around zero
    unless strong directed regulation is known, such as a 2-fold spike
    in study.

#### How to interpret the results

The goal is to choose the normalization method that: - Minimizes
within-group variation (PCV/PEV/PMAD) - Maximizes correlation among
replicates - Produces centered log2-ratio distributions - Yields PCA
plots with expected grouping and minimal batch effects - Shows stable,
comparable total intensities

As demonstrated in the `proteiNorm` paper’s case studies (mouse, breast
cancer, spike-in datasets) , different experiments may favor different
methods—VSN and cyclic loess often rank highly, but no single method is
optimal for all designs. write_norm_report() makes these comparisons
reproducible for any proteoDA workflow.

*Reference:* Graw, S., Tang, J., Zafar, M. K., Byrd, A. K., Bolden, C.,
Peterson, E. C., & Byrum, S. D. (2020). proteiNorm – A User-Friendly
Tool for Normalization and Analysis of TMT and Label-Free Protein
Quantification. ACS Omega, 5(40), 25625-25633.
<https://doi.org/10.1021/acsomega.0c02564>

``` r
write_norm_report(
filtered_proteins,
grouping_column = "group",
output_dir = "01_QC_report_v1.0",
filename = "normalization.pdf",
overwrite = TRUE
)

# Available methods: "log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi"

normalized <- normalize_data(filtered_proteins, norm_method = "cycloess")
```

#### Create a QC report including PCA and clustering:

``` r
write_qc_report(
  normalized,
  color_column = "group",
  output_dir = "01_QC_report_v1.0",
  filename = "QC_report.pdf",
  overwrite = TRUE,
  top_proteins = 500,
  standardize  = TRUE,
  pca_axes     = c(1, 2),
  dist_metric  = "euclidean",
  clust_method = "complete"
)
```

### 6. Differential Analysis (No-Intercept Model)

### 6.1 Add limma design and contrasts.

For standard group comparisons, use a no-intercept model: `~ 0 + group`

Contrasts can be added using either the `contrasts_vector` or the
`contrasts_file` option. When using the file option, save the contrasts
in column A as a .csv file.

\[Provide example table here\]

``` r
no_intercept <- add_design(normalized, design_formula = ~0 + group)

# Define contrasts of interest
no_intercept <- add_contrasts(
no_intercept,
contrasts_vector = "cancer_vs_normal = cancer - normal"
)

no_intercept$design$design_matrix
```

    ##                     cancer normal
    ## Sample_01_0325.mzML      0      1
    ## Sample_02_0147.mzML      0      1
    ## Sample_03_4085.mzML      0      1
    ## Sample_04_0039.mzML      0      1
    ## Sample_05_5198.mzML      0      1
    ## Sample_06_7066.mzML      0      1
    ## Sample_07_0909.mzML      0      1
    ## Sample_08_0348.mzML      1      0
    ## Sample_09_7028.mzML      1      0
    ## Sample_10_2119.mzML      1      0
    ## Sample_11_6501.mzML      1      0
    ## Sample_12_9999.mzML      1      0
    ## Sample_13_1014.mzML      1      0
    ## Sample_14_1059.mzML      1      0
    ## attr(,"assign")
    ## [1] 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$group
    ## [1] "contr.treatment"

### 6.2 Fit the linear model

When you call fit_limma_model(no_intercept) in proteoDA, under the hood
it’s doing: fit \<- limma::lmFit(exprs, design) fit \<-
limma::eBayes(fit)

``` r
fit <- fit_limma_model(no_intercept)
```

### 6.3 Extract the differential analysis results:

The `pval_thresh` and `lfc_thresh` parameters will define how
significant features are defined and highlighted in plots and tables.

``` r
results <- extract_DA_results(
fit,
pval_thresh      = 0.05,
lfc_thresh       = 1,
adj_method       = "BH",
extract_intercept = FALSE
)
```

## 7. Reporting Results

### Write statistical tables (CSV, Excel, and summary):

``` r
# in version 2.0, DA_table_cols must be defined in order to write the contrasts sheets. 
DA_table_cols <- c("uniprot_id","Accession.Number","Protein.Description")

write_limma_tables(
results,
output_dir        = "02_DA_results_v1.0",
overwrite         = TRUE,
add_filter        = TRUE,
contrasts_subdir  = NULL,
combined_file_csv = NULL,
spreadsheet_xlsx  = NULL
)
```

### Generate interactive and static volcano plots for each contrast:

``` r
write_limma_plots(
results,
grouping_column = "group",
table_columns   = c("uniprot_id", "protein_info"),
title_column    = "uniprot_id",
output_dir      = "02_DA_results_v1.0",
tmp_subdir      = "tmp",
height          = 1000,
width           = 1000,
overwrite       = TRUE
)
```

## 8. Summary

### This workflow demonstrates the standard proteoDA v1.0 pipeline:

1.  Create a `DAList` from raw data and metadata
2.  Filter samples and proteins
3.  Normalize and QC
4.  Define the limma model with
    [`add_design()`](https://byrumlab.github.io/proteoDA/reference/add_design.md)
    and
    [`add_contrasts()`](https://byrumlab.github.io/proteoDA/reference/add_contrasts.md)
5.  Fit the model with
    [`fit_limma_model()`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md)
6.  Extract results with
    [`extract_DA_results()`](https://byrumlab.github.io/proteoDA/reference/extract_DA_results.md)
7.  Generate reports with
    [`write_limma_tables()`](https://byrumlab.github.io/proteoDA/reference/write_limma_tables.md)
    and
    [`write_limma_plots()`](https://byrumlab.github.io/proteoDA/reference/write_limma_plots.md)

This basic structure extends easily to factorial, paired, or blocked
designs using the same core functions (add_design, add_contrasts, and
fit_limma_model).

## 9. Advanced: 2×2 Factorial Design (Intercept Model with Interaction)

This example demonstrates a 2×2 factorial design with interaction using
the **core v1.0 API** only.  
We assume your metadata has two factors:

- `cell` with levels `WT`, `HET` (genotype)  
- `treat` with levels `Sham`, `Gy` (treatment)

> **Note:** `proteoDA` simplifies design-matrix column names (e.g.,
> `Intercept`, `HET`, `Gy`, `HET.treatGy`) so contrasts and plot labels
> are cleaner.

### 9.1 set sample group reference levels

``` r
## ---- nine_setup_factors, echo=TRUE ------------------------------------------
# Ensure the factors exist in the metadata and set reference levels
stopifnot(all(c("cell", "treat") %in% names(sample_metadata)))

normalized$metadata$cell  <- factor(normalized$metadata$cell,  levels = c("WT","HET"))
normalized$metadata$treat <- factor(normalized$metadata$treat, levels = c("Sham","Gy"))
```

### 9.2 setup contrasts example

``` r
## ---- nine_show_contrasts_table, echo=FALSE, message=FALSE -------------------
# Render a small table of the contrasts we will use (for the vignette)
contrasts_df <- data.frame(
  Label = c(
    "WT_Gy_vs_HET_Gy",
    "WT_Sham_vs_HET_Sham",
    "TreatEffect_WT",
    "TreatEffect_HET",
    "Diff_in_effect_WT_HET",
    "Avg_Treat_Gy_vs_Sham_unscaled"
  ),
  Definition = c(
    "(Intercept + Gy) - (Intercept + HET + Gy + HET.treatGy)",
    "Intercept - (Intercept + HET)",
    "Gy",
    "Gy + HET.treatGy",
    "Gy - (Gy + HET.treatGy)",
    "(Intercept + Gy) + (Intercept + HET + Gy + HET.treatGy) - (Intercept) - (Intercept + HET)"
  ),
  stringsAsFactors = FALSE
)

knitr::kable(contrasts_df, caption = "Example contrasts for 2×2 factorial (intercept model)")
```

| Label                         | Definition                                                                                |
|:------------------------------|:------------------------------------------------------------------------------------------|
| WT_Gy_vs_HET_Gy               | (Intercept + Gy) - (Intercept + HET + Gy + HET.treatGy)                                   |
| WT_Sham_vs_HET_Sham           | Intercept - (Intercept + HET)                                                             |
| TreatEffect_WT                | Gy                                                                                        |
| TreatEffect_HET               | Gy + HET.treatGy                                                                          |
| Diff_in_effect_WT_HET         | Gy - (Gy + HET.treatGy)                                                                   |
| Avg_Treat_Gy_vs_Sham_unscaled | (Intercept + Gy) + (Intercept + HET + Gy + HET.treatGy) - (Intercept) - (Intercept + HET) |

Example contrasts for 2×2 factorial (intercept model)

#### Save the exact same contrasts to a CSV that add_contrasts() can read

#### (a simple one-contrast-per-line plain text file also works; here we use writeLines for clarity)

``` r
## ---- nine_write_contrasts_csv, echo=TRUE ------------------------------------
# Save the exact same contrasts to a CSV that add_contrasts() can read
# (a simple one-contrast-per-line plain text file also works; here we use writeLines for clarity)
fac_con_lines <- c(
  "WT_Gy_vs_HET_Gy = (Intercept + Gy) - (Intercept + HET + Gy + HET.treatGy)",
  "WT_Sham_vs_HET_Sham = Intercept - (Intercept + HET)",
  "TreatEffect_WT = Gy",
  "TreatEffect_HET = Gy + HET.treatGy",
  "Diff_in_effect_WT_HET = Gy - (Gy + HET.treatGy)",
  "Avg_Treat_Gy_vs_Sham_unscaled = (Intercept + Gy) + (Intercept + HET + Gy + HET.treatGy) - (Intercept) - (Intercept + HET)"
)
writeLines(fac_con_lines, "contrasts_factorial_intercept.csv")
```

### 9.3 Add design and contrasts for an interaction model (two cell types and two treatments).

#### The factorial approach formally tests whether treatment effects differ between cell lines.

#### Two-factor encoding: ~ treat1 \* treat2 includes a main effect and an interaction term; the : column quantifies interaction (synergy or repression).

``` r
# Add factorial design with interaction (intercept model)
fac <- add_design(
  normalized,
  design_formula = ~ cell * treat
)

# Add contrasts from the CSV we just wrote
fac <- add_contrasts(
  fac,
  contrasts_file = "contrasts_factorial_intercept.csv"
)

# (Optional sanity check)
colnames(fac$design$design_matrix)
# Expected: "Intercept", "HET", "Gy", "HET.treatGy"
```

### 9.4 Fit and extract results

``` r
## ---- nine_fit_and_extract, echo=TRUE ----------------------------------------
# Fit and extract (v1.0 core functions)
fit_fac <- fit_limma_model(fac)

results_fac <- extract_DA_results(
  fit_fac,
  pval_thresh       = 0.05,
  lfc_thresh        = 1,
  adj_method        = "BH",
  extract_intercept = FALSE   # keep FALSE: we are interested in contrasts, not the Intercept coefficient
)

# Write tables and plots for the factorial contrasts
write_limma_tables(
  results_fac,
  output_dir = "03_DA_factorial",
  overwrite  = TRUE,
  add_filter = TRUE
)

write_limma_plots(
  results_fac,
  grouping_column = "group",                 # or another metadata column for color
  table_columns   = c("uniprot_id","protein_info"),
  title_column    = "uniprot_id",
  output_dir      = "03_DA_factorial",
  tmp_subdir      = "tmp",
  height          = 1000,
  width           = 1000
)
```

### 9.5 Interpretation of the Interaction model `~ cell * treat`

**Interpretation tips**

- `TreatEffect_WT` is the Gy vs Sham effect **within WT** (reference
  genotype).  
- `TreatEffect_HET` is the Gy vs Sham effect **within HET** (WT
  treatment effect + interaction).  
- `Diff_in_effect_WT_HET` compares treatment effects between genotypes
  (the interaction contrast).  
- `Avg_Treat_Gy_vs_Sham_unscaled` is the **sum** of the two
  genotype-specific treatment effects.  
  If you need the literal average, divide the logFC by 2 after the fact
  in a downstream step (see proteoDA version 2.0 tutorial)

## 10. Paired / Blocking Designs (v1.0 Functions)

This section illustrates two approaches to handle non-independent
samples in proteoDA v1.0:

- **Paired fixed-effect model**: explicitly include a `pair_id` term in
  the design (~ 0 + group + pair_id).  
- **Blocking with `duplicateCorrelation`**: set
  `DAList$design$random_factor` to the blocking column name before
  calling
  [`fit_limma_model()`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md).  
  The function then estimates the intra-block correlation internally and
  applies it to the limma fit.

> ✅ Use the fixed-effect approach when the number of pairs is modest
> and you want strict within-pair comparisons.  
> ✅ Use `duplicateCorrelation` when you have many features and repeated
> measures (efficient and robust for omics).

### 10.1 Paired (fixed-effect) model

Design: `~ 0 + group + pair_id` Each pair gets its own coefficient.
contrasts compare groups within pairs.

``` r
## ---- ten_paired_fixed_effect_setup, echo=TRUE --------------------------------
# Requirements in metadata:
#  - a grouping variable used in your comparison, e.g. "group" with levels {Pre, Post}
#  - a pairing variable, e.g. "pair_id" (or "subject"), repeated for matched samples

# Metadata must contain both "group" and "pair_id"
stopifnot(all(c("group","pair_id") %in% names(normalized$metadata)))

# Ensure factors (and set reference levels if needed)
normalized$metadata$group   <- factor(normalized$metadata$group)         # e.g., c("Pre","Post")
normalized$metadata$pair_id <- factor(normalized$metadata$pair_id)

# Add paired fixed-effect design
paired_fe <- add_design(
  normalized,
  design_formula = ~ 0 + group + pair_id
)

# Contrasts: Post vs Pre (labels show up in plots)
paired_fe <- add_contrasts(
  paired_fe,
  contrasts_vector = "Post_vs_Pre = Post - Pre"
)

# Fit and extract
fit_paired_fe <- fit_limma_model(paired_fe)

res_paired_fe <- extract_DA_results(
  fit_paired_fe,
  pval_thresh       = 0.05,
  lfc_thresh        = 0.5,
  adj_method        = "BH",
  extract_intercept = FALSE
)

# Write outputs
write_limma_tables(
  res_paired_fe,
  output_dir = "10A_Paired_FixedEffect",
  overwrite  = TRUE,
  add_filter = TRUE
)

write_limma_plots(
  res_paired_fe,
  grouping_column = "group",
  table_columns   = c("uniprot_id","protein_info"),
  title_column    = "uniprot_id",
  output_dir      = "10A_Paired_FixedEffect",
  tmp_subdir      = "tmp",
  height          = 1000,
  width           = 1000
)
```

#### Notes

1.  This absorbs baseline differences per pair via a fixed effect
    (`pair_id` defined in sample_metadata.csv).
2.  Use when every pair has both conditions (e.g., Pre/Post) and the
    number of pairs is not huge.

### 10.2 Blocking with internal duplicateCorrelation

Design: `~ 0 + group` Random factor: assign the blocking column name to
`DAList$design$random_factor`.

[`fit_limma_model()`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md)
will:

1.  Estimate the consensus correlation via `duplicateCorrelation()`
2.  Inform the user of the estimated value
3.  Apply it automatically during the model fit.

limma estimates an intra-block correlation and uses it in the model
fit—computationally efficient and close in spirit to a random intercept.

``` r
## ---- ten_blocking_dupCor_setup, echo=TRUE -----------------------------------
# Requirements in metadata:
#  - a grouping variable (e.g. "group" with levels {Pre, Post})
#  - a blocking variable (e.g. "pair_id" or "subject"), repeated for correlated samples

stopifnot(all(c("group","pair_id") %in% names(normalized$metadata)))

# Factor coding (reference choice doesn't matter for no-intercept designs)
normalized$metadata$group   <- factor(normalized$metadata$group)
normalized$metadata$pair_id <- factor(normalized$metadata$pair_id)

# Add simple group-means design (no intercept)
blocked <- add_design(
  normalized,
  design_formula = ~ 0 + group
)

# Contrast: Post vs Pre
blocked <- add_contrasts(
  blocked,
  contrasts_vector = "Post_vs_Pre = Post - Pre"
)

# Specify the blocking variable for duplicateCorrelation
blocked$design$random_factor <- "pair_id"

# Fit the model; duplicateCorrelation will run automatically
fit_blocked <- fit_limma_model(blocked)

res_blocked <- extract_DA_results(
  fit_blocked,
  pval_thresh       = 0.05,
  lfc_thresh        = 0.5,
  adj_method        = "BH",
  extract_intercept = FALSE
)

# Write outputs
write_limma_tables(
  res_blocked,
  output_dir = "10B_Blocked_dupCor",
  overwrite  = TRUE,
  add_filter = TRUE
)

write_limma_plots(
  res_blocked,
  grouping_column = "group",
  table_columns   = c("uniprot_id","protein_info"),
  title_column    = "uniprot_id",
  output_dir      = "10B_Blocked_dupCor",
  tmp_subdir      = "tmp",
  height          = 1000,
  width           = 1000
)
```

#### Interpretation

1.  `duplicateCorrelation()` estimates an average
    within-pair/within-subject correlation and passes it to lmFit() so
    the residual variance accounts for dependence.

2.  Great for repeated measures and many features (omics), where a full
    mixed model would be heavy.

3.  [`fit_limma_model()`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md)
    prints an estimated intra-block correlation (e.g., 0.12). If it’s
    low (\< 0.05), the function issues a yellow message suggesting that
    a simple model may suffice.

4.  This approach is efficient for large omics datasets and conceptually
    similar to a random-intercept mixed model.

5.  It provides most of the benefit of modeling correlation within
    subjects without the heavy computation of lme4 or dream.

#### Choosing between the two

1.  Use **paired fixed-effect** when each pair is complete and you want
    strict within-pair comparisons (each pair gets a coefficient).

2.  Use **blocking** (duplicateCorrelation) for scalability and when you
    expect a common within-block correlation across features;
    conceptually similar to a **random intercept**.

## 11. Diagnostics for models

### Checks full rank, condition number, non-estimable contrasts, and (if using blocking) prints the duplicateCorrelation estimate. It works for either path:

1.  paired_fe (paired fixed-effect model)
2.  blocked (blocking with random_factor)

``` r
## ---- diagnostics_model_checks, echo=TRUE, message=FALSE, warning=FALSE -------
# Helper to print a clean header
.hdr <- function(x) cat("\n", paste0("## ", x), "\n", sep="")

# Choose which object to diagnose (set one that exists)
obj <- if (exists("blocked")) blocked else if (exists("paired_fe")) paired_fe else NULL
stopifnot(!is.null(obj))

X <- obj$design$design_matrix
C <- obj$design$contrast_matrix

.hdr("Design matrix summary")
cat("Dimensions (n x p):", paste(dim(X), collapse=" x "), "\n")
cat("Column names:\n"); print(colnames(X))

# Full-rank / rank deficiency check
qrX   <- qr(X)
rankX <- qrX$rank
p     <- ncol(X)
cat("QR rank:", rankX, "  |  Columns:", p, "\n")
if (rankX < p) {
  cat("WARNING: design is NOT full rank. Check factor coding / redundant columns.\n")
} else {
  cat("OK: design is full rank.\n")
}

# Condition number as a rough collinearity diagnostic
k <- tryCatch(kappa(X), error = function(e) NA_real_)
cat("Condition number (kappa):", if (is.na(k)) "NA" else sprintf("%.2f", k), "\n")
if (!is.na(k) && k > 1e6) cat("NOTE: Very high kappa suggests near-collinearity.\n")

# Contrast estimability
.hdr("Contrast estimability (limma::nonEstimable)")
if (!is.null(C)) {
  ne <- limma::nonEstimable(C, X)
  if (is.null(ne)) {
    cat("All contrasts appear estimable given the design.\n")
  } else {
    cat("Non-estimable (or suspect) contrasts detected:\n")
    print(ne)
  }
} else {
  cat("No contrast matrix present.\n")
}

# If using blocking via random_factor, show duplicateCorrelation estimate
.hdr("Blocking / duplicateCorrelation diagnostic")
if (!is.null(obj$design$random_factor)) {
  blk_col <- obj$design$random_factor
  stopifnot(blk_col %in% colnames(obj$metadata))
  block_vec <- obj$metadata[, blk_col, drop = TRUE]
  dc <- limma::duplicateCorrelation(object = obj$data, design = X, block = block_vec)
  cat("Estimated intra-block correlation (duplicateCorrelation):",
      sprintf("%.3f", dc$consensus.correlation), "\n")
} else {
  cat("No random_factor set; duplicateCorrelation not applied.\n")
}

# Quick peek at coefficients and residual df if already fitted
.hdr("Fit overview (if available)")
if (!is.null(obj$eBayes_fit)) {
  efit <- obj$eBayes_fit
  cat("Coefficients matrix dims:", paste(dim(efit$coefficients), collapse=" x "), "\n")
  if (!is.null(efit$df.residual)) cat("Residual df (median across features):",
                                      stats::median(efit$df.residual, na.rm = TRUE), "\n")
  # Show first few coefficient names
  cat("Coefficient names:\n"); print(colnames(efit$coefficients))
} else {
  cat("Model not fitted yet (no eBayes_fit found in object).\n")
}
```

### What it does:

1.  Full rank: uses [`qr()`](https://rdrr.io/r/base/qr.html) to compare
    rank vs. number of columns; warns if deficient.
2.  Collinearity: reports `kappa(X)` (large values ~ near-collinearity).
3.  Estimability: uses
    [`limma::nonEstimable()`](https://rdrr.io/pkg/limma/man/isfullrank.html)
    to flag contrasts that aren’t estimable under the design.
4.  Blocking: if `obj$design$random_factor` is set, re-computes and
    prints the duplicateCorrelation consensus correlation (purely
    diagnostic).
5.  Fit summary: if you’ve already run
    [`fit_limma_model()`](https://byrumlab.github.io/proteoDA/reference/fit_limma_model.md),
    prints coefficient dimensions, residual df, and coefficient names.
