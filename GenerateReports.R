## Generate Reports
source("proteoDA_params.R")
# Generate Methods.html report
# test changing render to ".Rmd" within working directory so output is in that directory
rmarkdown::render("DIA_methods_report.Rmd",
                  knit_root_dir = getwd(),
                  #intermediates_dir = tmp_subdir,
                  output_file = paste0(project_name, "_methods.html"),
                  quiet = T)

## Generate QC_report

QC_file <- "QC_report.html"
print(paste("Generating report:", QC_file))

rmarkdown::render("QC_report.Rmd",
                  knit_root_dir = getwd(),
                  #intermediates_dir = tmp_subdir,
                  output_file = paste0(project_name, "_QC_report.html"),
                  quiet = T)

# Render the Quarto document with a dynamic output file name bash
# quarto render Project_Summary_PowerPoint_20250228_v2.qmd --to pptx # installed quarto using conda

PP_file <- "PowerPoint_report.html"
print(paste("Generating report:", PP_file))

library(yaml)
library(quarto)

# Load YAML variables
vars <- yaml::read_yaml("_variables.yml")

# Define output file path (adjust the path as needed)
output_path <- paste0(vars$project, ".pptx")

# Render the Quarto document with a dynamic output filename
quarto_render(
  input = "Project_Summary_20250425.qmd", # Change to your actual .qmd file
  output_file = output_path
)

# Print message for confirmation
message("Report rendered successfully: ", output_path)

######
## CREATE results directory
######

print("Generating project deliverables directory:")

# Define your results_dir (target root directory)
results_dir <- project_name 

# Define source paths
figures_src <- "figures"
static_plots_src <- "interactive_results/static_plots"
interactive_src <- "interactive_results"

# Define destination paths
qc_figures_dest <- file.path(results_dir, "figures", "QC")
volcano_md_dest <- file.path(results_dir, "figures", "volcano_MD_plots")
interactive_dest <- file.path(results_dir, "interactive_results")

# Create destination directories
dir.create(qc_figures_dest, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_md_dest, recursive = TRUE, showWarnings = FALSE)
dir.create(interactive_dest, recursive = TRUE, showWarnings = FALSE)

# Copy all files from figures/ to figures/QC
figures_files <- list.files(figures_src, full.names = TRUE)
file.copy(from = figures_files, to = qc_figures_dest, overwrite = TRUE)

# Copy all PDF files from interactive_results/static_plots to volcano_MD_plots
volcano_pdfs <- list.files(static_plots_src, pattern = "\\.pdf$", full.names = TRUE)
file.copy(from = volcano_pdfs, to = volcano_md_dest, overwrite = TRUE)

# Copy .html, .xlsx, and DA_summary.csv from interactive_results to interactive_results in results_dir
html_files <- list.files(interactive_src, pattern = "\\.html$", full.names = TRUE)
xlsx_files <- list.files(interactive_src, pattern = "\\.xlsx$", full.names = TRUE)
da_summary <- file.path(interactive_src, "DA_summary.csv")

files_to_copy <- c(html_files, xlsx_files, da_summary)
file.copy(from = files_to_copy, to = interactive_dest, overwrite = TRUE)

cat("Files copied to:", results_dir, "\n")

# ---- ZIP the results_dir ----

# Define zip file name
zipfile_name <- paste0(basename(results_dir), ".zip")

# Zip the entire directory recursively (9 maximum compression, X no extra file metadata)
zip(zipfile_name, files = results_dir, flags = "-r9X", zip = "/usr/bin/zip")

cat("📦 Zipped directory to:", zipfile_name, "\n")