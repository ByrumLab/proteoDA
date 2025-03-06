## Generate Reports

source("proteoDA_params.R")

# Generate Methods.html report
# test changing render to ".Rmd" within working directory so output is in that directory
rmarkdown::render("DIA_methods_report.Rmd",
                  knit_root_dir = getwd(),
                  #intermediates_dir = tmp_subdir,
                  output_file = paste0(methods_filename, "_methods.html"),
                  quiet = T)

## Generate QC_report

rmarkdown::render("QC_report.Rmd",
                  knit_root_dir = getwd(),
                  #intermediates_dir = tmp_subdir,
                  output_file = paste0(methods_filename, "_QC_report.html"),
                  quiet = T)
