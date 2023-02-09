
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialews_biocom

<!-- badges: start -->

[![License: GPL (\>=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
[![Dependencies](https://img.shields.io/badge/dependencies-2/94-green?style=flat)](#)
<!-- badges: end -->

Research Compendium of the project **Spatial metrics calculation in the
Biocom data set**

### How to cite

Please cite this compendium as:

> **Sonia Kéfi and Alexandre Génin. Code associated with the paper
> ‘Dryland resilience to aridity is associated with stonger vegetation
> patchiness’**

### Content

This repository is structured as follow:

- [`data/`](https://github.com/skefi/spatialews_biocom/tree/master/raw-data):
  contains raw data required to perform analyses (model and real data)

- [`R/`](https://github.com/skefi/spatialews_biocom/tree/master/R):
  contains R functions allowing to prepare the data and calculate the
  spatial metrics in the model and in the data

- [`outputs/`](https://github.com/skefi/spatialews_biocom/tree/master/outputs):
  contains the results created during the workflow

- [`analyses/`](https://github.com/skefi/spatialews_biocom/tree/master/analyses/):
  contains R scripts to calculate indicator trends and plot figures

- [`figures/`](https://github.com/skefi/spatialews_biocom/tree/master/figures):
  contains the figures created during the workflow

- [`man/`](https://github.com/skefi/spatialews_biocom/tree/master/man):
  contains help files of R functions

- [`DESCRIPTION`](https://github.com/skefi/spatialews_biocom/tree/master/DESCRIPTION):
  contains project metadata (author, date, dependencies, etc.)

- [`make.R`](https://github.com/skefi/spatialews_biocom/tree/master/make.R):
  main R script to run the entire project by calling each R script
  stored in the `analyses/` folder

### Usage

Clone the repository, open spatialews_biocom.Rproj with R/Rstudio.

To calculate the spatial metrics from the the real images (that are in
data/data_CA) and from data generated by the model (that are in
data/data_images_biocom), run:

``` r
source("make.R")
```

This generates result files that are saved in the folder ‘outputs’. In
particular, biocom-grps.rda contains information about the field sites,
indics-data.rda contains a set of spatial metrics measured on each of
the real images, and indics-model.rda contains the same set of spatial
metrics measured on the images generated by the model.

Using these files, code available in the folder ‘analyses’ allows
performing additional analyses and plotting the figures of the paper.

In the folder ‘analyses’, calculate_indicator_trends.R and
calculate_model_indicator_trends.R allow calculating the trends of each
of the indicators along a stress gradient. They generate the files
trends_model.rda, trends_one_group.rda, trends_two_groups.rda and
trends_three_groups.rda which are in the ‘outputs’ folder and are
required to generate the figures of the paper.

The plot_fig.R files allow generating the figures of the paper. The
files clustering_analyses.R and compare_groups.R perform additional
analyses.

### Notes

- All required packages, listed in the `DESCRIPTION` file, will be
  installed (if necessary)
- All required packages and R functions will be loaded
- Some analyses listed in the `make.R` might take time
