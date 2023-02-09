#' spatialews_biocom
#' 
#' @description 
#' Spatial metrics calculation in the Biocom data set
#' 
#' @author Sonia Kefi \email{sonia.kefi@umontpellier.fr}
#' 
#' @date 2023/02/07



## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")


## Run Project ----

targets::tar_make()


# In the R console, type : 
#source("make.R")

# This will run the information in _targets.R file (i.e. prepare the data and calculates all spatial metrics)

# This will generate the fies 'data_biocom.rda', 'biocom-grps.rda', 'data_model.rda', 'indics-data-scaled.rda', 'indics-model.rda' in the folder 'outputs'

# Once this is finished, you can visualize the steps of the analyses by typing in the console:
#library(targets)
#tar_visnetwork()

# To calculate indicator trends in the real or model data, run the codes 'calculate_indicator_trends.R' or 'calculate_model_indicator_trends.R' in the folder 'analyses'. 

# This will generate the files 'trends_one_group.rda', 'trends_two_group.rda', 'trends_three_group.rda', 'trends_model.rda' in the folder 'outputs'. 

# To plot the figures of the paper, run the files 'plot_fig.R' in the folder 'analyses'




