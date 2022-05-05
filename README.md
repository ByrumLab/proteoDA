
<!-- README.md is generated from README.Rmd. Please edit that file -->

# proteomicsDIA

<!-- badges: start -->
<!-- badges: end -->

The `proteomicsDIA` package will (eventually) be our internal R package
that supplies all the functions we need for doing proteomics analysis of
DIA data.

## Installation

The Github repository where the package is stored is private, so
installation is a little more complicated that usual. You’ll need the
`devtools` package installed (which you’ll want to have installed anyway
for developing the package, see below). You’ll also need the personal
access token (PAT) that gives you access to the repository. Email Tim or
Stephanie to get it. Once you have it, you can install the development
version of `proteomicsDIA` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ByrumLab/proteomicsDIA",
                         auth_token = "COPY_PAT_HERE")
```

## Example

The package is getting developed in pipeline order, adding each of the
functions and subfunctions needed for the various steps in the example
workflows (see the progress tracker below). The goal is to get all the
current functionality into the R package before we start adding new
features or making changes.

``` r
# Load the proteomicsDIA package
library(proteomicsDIA)

# Extract Maxquant data
extracted_data <- extract_data("path/to/Samples Report of from Maxquant.csv",
                               pipe = "DIA",
                               enrich = "protein")

# make targets
targets <- make_targets(file = "path/to/metdata.csv",
                        sampleIDs = colnames(extracted_data$data),
                        pipe = "DIA",
                        enrich = "protein")

# Subset targets
sub <- subset_targets(targets = targets, 
                      filter_column = "group",
                      rm.vals = "Pool")

# Process data
norm <- process_data(data = extracted_data$data,
                     targets = sub_higgs$targets,
                     min.reps = 5,
                     min.grps = 3)
# Make the proteinorm report
make_proteinorm_report(normList = norm$normList,
                       groups = norm$targets$group, 
                       enrich = "protein")

# Make the QC report
make_qc_report(normList = norm$normList, 
               norm.method = "vsn",
               groups = norm$targets$group,
               batch = norm$targets$group,
               enrich = "protein")

# Make the design matrix
design <- make_design(targets=norm$targets,
                      group_column = "group",
                      factor_columns = NULL,
                      paired_column = NULL)
# Make the contrasts matrix
contrasts <- make_contrasts(file = "path/to/contrasts/file.csv",
                            design = design$design)

# Functions are documented, check them out:
?extract_data
?make_targets
?subset_targets
?process_data
?make_proteinorm_report
?make_qc_report
?make_design
?make_contrasts

# If you see any typos or things that aren't clear, let me know!
```

## Contributing to package development

I would highly recommend checking out the [R
Packages](https://r-pkgs.org/) book from Hadley Wickham and Jenny Bryan.
It is a pretty comprehensive book on creating R packages, and provides a
nice framework that we can all follow (and the framework which the
package is currently following).

You’ll want to install a few packages for development:

``` r
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "usethis"))
```

You may want to set up your `.Rprofile` to automatically load some of
all of the packages when you start R (see
[here](https://r-pkgs.org/setup.html#personal-startup-configuration) for
more info), though I don’t usually bother and just load them as I need
them.

### Collaboration with Git

We need to figure this out! Deciding on the best path forward with this
depends a little on everyone’s familiarity and level of comfort with
Git. I (Tim) have used Git for the last couple years for all my
projects, though I have rarely used it to collaborate with others and
thus have rarely needed to worry about the optimal collaboration
workflow. We can work together to figure that out more, but in the
meantime here are some resources that might be helpful:

-   [Happy Git with R](https://happygitwithr.com/index.html)- A guide to
    using Git and Github with RStudio.
-   [Chapter 18 of R packages](https://r-pkgs.org/git.html) is on Git.
-   Git GUIs- Personally, I mostly avoid using Git in the command line
    and instead use GUIs. RStudio has a minimal one, explained in the
    first link above. Initially, I mostly used
    [GitKraken](https://www.gitkraken.com/git-client), which is free
    with a student license and has some nice touches (better diff
    viewing, better branch visualization). Lately, I’ve just been using
    the Git functionality in [Visual Studio
    Code](https://code.visualstudio.com/), which also has some nice
    extensions for using Git.

Once we’re all set up and comfortable with Git, we’ll need to figure out
how to manage workflows. That might be a topic for another day.

### Misc ideas

See the GitHub issue tracker for other ideas

-   Add a changelog.md file, to textually track/explain changes (with
    version #s).

-   Check out code coverage and CI tools, to automate testing. Some
    options in usethis. Might want to wait to do much testing until
    after we decide on whether we want to do much re-writing and
    reformatting of the code.

### Progress tracker

I’m working my way through the “top-level” functions in the example DIA
workflow sent by Charity. I’m not doing much on the TMT, LF, or
phosphoTMT side of things, though I have done a little code restyling
for those sometimes.

On the code side, I’ve been doing relatively little. I’ve been doing
some restyling of the code for clarity, and reduced redundancy. I’ve
dropped TODOs through the code to mark ideas for re-doing of code. The
main thing I’ve done is to rework print statements and error messages
(using the `cli` R package) to make these better formatted and more
informative.

On the documentation side, I’m just trying to get some initial
documentation for each main function and all its subfunctions, linking
them together. The goal is to explain all the inputs and outputs to get
a better idea of how they fit together.

On the test side, I started writing some unit tests for the data
extraction functions, but I think I might hold off on that for now. If
we do any changing of the code structure, we may have to re-write the
tests as well. Maybe better to just do it once.

| Function                 | Code               | Document           | Test               |
|--------------------------|--------------------|--------------------|--------------------|
| `extract_data`           | :heavy_check_mark: | :heavy_check_mark: | :heavy_check_mark: |
| `make_targets`           | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `subset_targets`         | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `process_data`           | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_proteinorm_report` | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_qc_report`         | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_design`            | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `make_contrasts`         | :heavy_check_mark: | :heavy_check_mark: | :x:                |
| `run_limma_analysis`     | :x:                | :x:                | :x:                |
| `add_limma_results`      | :x:                | :x:                | :x:                |
