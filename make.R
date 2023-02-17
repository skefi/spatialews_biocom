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


## Load Project Addins (R Functions and Packages) ----

devtools::load_all()


## Global Variables ----

path_data <- here::here("data")
path_model_data <- here::here("data/data_CA")
path_output <- here::here("outputs")
NPERM <- 3 #the value used for the analyses in the paper is 199
N_SHUFFLE <- 3
BOOTN <- 3
ALPHA <- 0.5

# We recommend 199 for NPERM because in the analysis we only use 
# the means of the distribution created. 199 is deemed sufficient to estimate a mean 
# with sufficient accuracy for our purposes. 


## Run Project ----

function_prepare_data(path_data, path_output)
function_compute_indicators(NPERM,path_output)
function_compute_model_indicators(NPERM,path_output,path_model_data)
function_make_branches(NPERM,path_output)
function_indicator_trends(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA) 
function_model_indicator_trends(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA)



####### NOTES #######
# In the R console, type : 
#source("make.R")

# This will run the functions above (i.e. prepare the data and calculates all spatial metrics)

# This will generate the files 'data_biocom.rda', 'biocom-grps.rda', 'data_model.rda', 'indics-data-scaled.rda', 'indics-model.rda', 'trends_one_group.rda', 'trends_two_group.rda', 'trends_three_group.rda', 'trends_model.rda' in the folder 'outputs'. 

# To plot the figures of the paper, run the files 'plot_fig.R' in the folder 'analyses'. 




