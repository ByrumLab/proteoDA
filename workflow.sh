module load conda3/202402
conda activate /hpcf/authorized_apps/proteomics_apps/conda/proteoda

Rscript proteoDA_QC_localrun_HPC_params2.R
mv *.png 02_DA/static_plots
Rscript GenerateReports.R
quarto render Project_Summary_PowerPoint_20250228_v2.qmd --to pptx

conda deactivate