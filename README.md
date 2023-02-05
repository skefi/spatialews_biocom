<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialews\_biocom

<!-- badges: start -->

[![License: GPL (&gt;=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
[![Dependencies](https://img.shields.io/badge/dependencies-2/94-green?style=flat)](#)
<!-- badges: end -->

Research Compendium of the project **{{ PLEASE ADD A FEW WORDS }}**

### How to cite

Please cite this compendium as:

> **{{ PLEASE ADD A CITATION }}**

### Content

This repository is structured as follow:

-   [`data/`](https://github.com/skefi/spatialews_biocom/tree/master/data):
    contains all raw data required to perform analyses

-   [`analyses/`](https://github.com/skefi/spatialews_biocom/tree/master/analyses/):
    contains R scripts to run each step of the workflow

-   [`outputs/`](https://github.com/skefi/spatialews_biocom/tree/master/outputs):
    contains all the results created during the workflow

-   [`figures/`](https://github.com/skefi/spatialews_biocom/tree/master/figures):
    contains all the figures created during the workflow

-   [`R/`](https://github.com/skefi/spatialews_biocom/tree/master/R):
    contains R functions developed especially for this project

-   [`man/`](https://github.com/skefi/spatialews_biocom/tree/master/man):
    contains help files of R functions

-   [`DESCRIPTION`](https://github.com/skefi/spatialews_biocom/tree/master/DESCRIPTION):
    contains project metadata (author, date, dependencies, etc.)

-   [`make.R`](https://github.com/skefi/spatialews_biocom/tree/master/make.R):
    main R script to run the entire project by calling each R script
    stored in the `analyses/` folder

### Usage

Clone the repository, open R/RStudio and run:

    source("make.R")

### Notes

-   All required packages, listed in the `DESCRIPTION` file, will be
    installed (if necessary)
-   All required packages and R functions will be loaded
-   Some analyses listed in the `make.R` might take time
