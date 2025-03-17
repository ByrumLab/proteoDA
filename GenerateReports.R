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

rmarkdown::render("QC_report.Rmd",
                  knit_root_dir = getwd(),
                  #intermediates_dir = tmp_subdir,
                  output_file = paste0(project_name, "_QC_report.html"),
                  quiet = T)

## Render the Quarto PowerPoint with a dynamic output file name based on project_name
# quarto render Project_Summary_PowerPoint_20250228_v2.qmd --to pptx # bash command line

library(yaml)
library(quarto)

# load yml variables
vars <- yaml::read_yaml("_variables.yml")

# Define output file name (NO PATHS allowed)
output_name <- paste0(vars$project, ".pptx")

# Render the Quarto PowerPoint with dynamic filename
quarto_render(
  input = "Project_Summary_PowerPoint_20250228_v2.qmd", 
  output_file = output_name
)

message("Report rendered successfully:", output_name)
