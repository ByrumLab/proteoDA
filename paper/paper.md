---
title: 'proteoDA: a package for quantitative proteomics'
tags:
  - R
  - mass spectrometry
  - intensity data
  - normalization
  - linear models
authors:
  - name: Charity L. Washam
    orcid: 0000-0001-5761-9304
    equal-contrib: true
    affiliation: 1 
  - name: Timothy Thurman
    orcid: 0000-0002-9602-6226
    equal-contrib: true 
    affiliation: 1
  - name: Duah Alkam
    orcid: 0000-0002-5965-7694
    affiliation: 1
  - name: Jordan T. Bird
    orcid: 0000-0001-5753-6058
    affiliation: 1
  - name: Allen Gies
    orcid: 0000-0003-2492-0429
    affiliation: 1
  - name: Kalyani Dhusia
    orcid: 0000-0002-8803-1295
    affiliation: 1
  - name: Michael S. Robeson, II
    orcid: 0000-0001-7119-6301
    affiliation: 2
  - name: Stephanie D. Byrum
    orcid: 0000-0002-1783-3610
    corresponding: true 
    affiliation: "1, 3"
affiliations:
 - name: Department of Biochemistry and Molecular Biology, University of Arkansas for Medical Sciences, Little Rock, AR, USA
   index: 1
 - name: Department of Biomedical Imaging, University of Arkansas for Medical Sciences, Little Rock, AR, USA
   index: 2
 - name: Arkansas Children's Research Institute, Little Rock, AR, USA
   index: 3
date: 23 January 2023
bibliography: references.bib
---

# Summary
`proteoDA` is an R package designed to analyze intensity based high resolution mass spectrometry data. The workflow assesses raw protein intensities using proteiNorm [@Graw2021] to evaluate eight normalization methods, provides a quality control report, and provides differential abundance analysis using limma models [@Ritchie2015]. The final results are included in an interactive .html file, which can be open in a web browser to explore the data [@ glimma references].

The R package requires a data frame or matrix containing protein intensities for each sample, an annotation data frame describing the proteins including uniprotID, description, gene symbol, and a metadata data frame containing information about the samples such as the sample name and group. 

The data is imported as a DAList, which is an S3 object containing seven slots in order to hold the data and results. The slots include data, annotation, metadata, design, eBayes_fit, results, and tags. The first three slots are required. The protein annotation must include a column called "uniprot_id". This should be the unique protein accession number from the database search results. 

Once the data is imported into the DAList object, several functions are included to process the data including functions to remove low quality samples, remove proteins with missing values in the majority of samples, normalize the data, and perform differential abundance analysis. An example workflow is provided in the vignette. The functions for proteiNorm are included to evaluate eight normalization methods and provide the final log2 normalized intensities for each sample [@ proteinorm references]. 

After the data is appropriately normalized, the quality control report is generated. The report includes a violin plot of the samples, PCA plot, clustered dendrogram, and a heatmap of proteins with missing values in the samples. 

The next step in the workflow includes defining the limma model and sample group comparisons for the differential abundance analysis [@Ritchie and limma refs]. Several model options are included and are described in the ?add_design function. 

`proteoDA` provides all the functions necessary to evaluate the quality, remove unwanted samples, filter proteins with missing values, normalize the data, perform statistical analysis, and export the results into an interactive html report. 


# Statement of need

`proteoDA` was designed to be used by researchers and students who would like to analyze quantitative proteomics data but have limited experience with mass spectrometry data and R programming. The package is used in the classroom as well as the IDeA National Resource for Quantitative Proteomics workshops for core directors, faculty, and students. 

# Acknowledgements

The development of this R package was supported by the National Institutes of Health National Institute of General Medical Sciences (NIH/NIGMS) grants P20GM121293, R24GM137786, the National Science Foundation Award No. OIA-1946391, and the UAMS Winthrop P. Rockefeller Cancer Institute. 

# References
