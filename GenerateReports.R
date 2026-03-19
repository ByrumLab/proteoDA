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

