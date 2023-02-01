---
title: 'proteoDA: a package for quantitative proteomics'
tags:
  - R
  - mass spectrometry
  - intensity data
  - normalization
  - linear models
authors:
  - name: Timothy J. Thurman
    orcid: 0000-0002-9602-6226
    equal-contrib: true 
    affiliation: 1
  - name: Charity L. Washam
    orcid: 0000-0001-5761-9304
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
`proteoDA` is an R package designed to analyze high resolution, intensity-based mass spectrometry data. `proteoDA` was designed to be streamlined and user-friendly: it uses simple syntax, functions can be easily assembled into pipelines, model specification is simple yet powerful and flexible, and `proteoDA` provides functions to generate rich plots, results, and interactive reports. 

`proteoDA` is built around a `DAList`, a custom R class, which is used to hold the data, statistical design, and results for a differential abundance analysis. Users import protein abundance data, protein annotation data, and sample metadata into a `DAList`, and all further functions operate on that list. Once the data are in `DAList`, `proteoDA` provides functions for further steps of the analysis \autoref{fig:workflow}. 

![A flowchart of the proteoDA workflow.\label{fig:workflow}](proteoDA_flowchart.png) 
`proteoDA` includes functions for 1) handling missing data and filtering samples and proteins based on metadata, annotation information, or missing data thresholds, 2) evaluation of multiple normalization methods using a graphical report based on the `proteiNorm` normalization tool [@Graw2021; @Chawade2014], 3) generation of graphical quality control reports to assess data quality and sample clustering, 4) flexible specification and fitting of differential abundance models (including mixed models), using the R package `limma` to perform model fitting [@Ritchie2015; @Law2020], and 5) generation of tabular results files, as well as interactive and portable HTML result files, using the R package `Glimma` [@Su2017]. `proteoDA` includes a detailed tutorial vignette that explains the analytic pipeline and key functions in further detail.  


# Statement of need

`proteoDA` was designed to help researchers with minimal knowledge of R extract insights from proteomic data. `proteoDA` was originally developed to support the needs of a high-volume proteomics core facility which analyzes hundreds of proteomics projects per year and produces research outputs and reports for customers with a wide range of expertise in proteomics. Not coincidentally, the features that make proteoDA useful for a high-volume core (simple yet flexible syntax, ease of pipeline creation, minimal need for R code for customization, and rich output reports) also make it ideal for users with minimal R experience.

`proteoDA` allows users to quickly assess the quality of a mass spectrometry experiment, normalize the data, control for batch effects, and define flexible linear model designs for a wide variety of proteomic experiments. In addition to providing quantitative analysis for experiments, `proteoDA` is used in the classroom as well as the IDeA National Resource for Quantitative Proteomics workshops for faculty and students.

# Acknowledgements

The development of this R package was supported by the National Institutes of Health National Institute of General Medical Sciences (NIH/NIGMS) grants P20GM121293, R24GM137786, the National Science Foundation Award No. OIA-1946391, and the UAMS Winthrop P. Rockefeller Cancer Institute. 

# References
