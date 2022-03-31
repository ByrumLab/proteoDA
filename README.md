
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

At the start, `proteomicsDIA` only has one function:
`uselsss_function()`:

``` r
# Load the proteomicsDIA package
library(proteomicsDIA)


# call our only function
useless_function("It doesn't matter what argument you supply")

# Get help on how to use useless_function:
?useless_function
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

-   Add a changelog.md file, to textually track/explain changes (with
    version #s).

-   Functions to make project directories with standard formats? See
    rstantools::rstan_create_package() and usethis::create_project() as
    inspiration. Could have a function that initializes a new RStudio
    project (one for each analysis project), with a pre-defined folder
    structure, preloaded template workflow scripts, etc. Might ease
    things for running projects and generating outputs for both end
    users and for us.

-   Check out code coverage and CI tools, to automate testing. Some
    options in usethis.
