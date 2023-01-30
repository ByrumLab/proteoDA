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
---

# Summary

`proteoDA` is an R package for the analysis of intensity-based, high resolution protein mass spectrometry data. `proteoDA` was designed to be streamlined and user-friendly, especially for users who need to analyze proteomic data but have minimal knowledge of R: `proteoDA` uses simple syntax, all functions are designed to be easily assembled into pipelines, model specification is simple yet powerful and flexible, and `proteoDA` includes functions to generate rich plots, results, and interactive reports.

`proteoDA` is built around a `DAList`, a custom R class which is used to hold the data, statistical design, and results for a differential abundance analysis. Users import protein abundance data, protein annotation data, and sample metadata into a `DAList`, and all further functions operate on that list. 


![A flowchart of the proteoDA workflow.\label{fig:workflow}](proteoDA_flowchart.png)

Once the data are in `DAList`, `proteoDA` provides functions for further steps of the analysis \autoref{fig:workflow}. `proteoDA` includes functions for:

* Filtering of samples and proteins based on sample data, annotation data, or missing data thresholds.

* Evaluation of multiple normalization methods using a graphical report based on the `proteiNorm` normalization tool [@Graw2021; @Chawade2014].

* Generation of graphical quality control reports to assess data quality and sample clustering.

* Flexible specification and fitting of differential abundance models (including mixed models), using the R package `limma` to perform model fitting [@Ritchie2015; @Law2020].

* Generation of tabular results files, as well as interactive and portable HTML results files, using the R package `Glimma` [@Su2017].


# Statement of need

XXXPARAGRAPH IDEA: Here it might be nice to have a paragraph that opens with a couple sentences on the importance of proteomics as a field (with some general review citations), maybe some references to some recent papers (maybe even ones from our core?). Would help emphasize the need for the package more. Might be especially nice if we highlight some examples coming from clinicians/clinician-researchers, since they might not be expected to know R. Though I don't know enough about the field to write this section.XXX

`proteoDA` was designed to help researchers with minimal knowledge of R extract insights from their proteomic data. `proteoDA` was originally developed to support the needs of a high-volume proteomics core facility which analyzes hundreds of proteomics projects per year and produces research outputs and reports for customers with a wide range of expertise in proteomics. Not coincidentally, the features that make `proteoDA` useful for a high-volume core (simple yet flexible syntax, ease of pipeline creation, minimal need for R code for customization, and rich output reports) also make it ideal for users with minimal R experience. XXXINCLUDE REFERENCES TO SIMILAR PACKAGES, OR NOT?XXX This package is also used as a teaching tool in the classroom, as well as in IDeA National Resource for Quantitative Proteomics workshops for core directors, faculty, and students. 

# Acknowledgements

The development of this R package was supported by the National Institutes of Health National Institute of General Medical Sciences (NIH/NIGMS) grants P20GM121293, R24GM137786, the National Science Foundation Award No. OIA-1946391, and the UAMS Winthrop P. Rockefeller Cancer Institute. 

# References
