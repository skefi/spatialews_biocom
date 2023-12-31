#' spatialews_biocom
#' 
#' @description 
#' Spatial metrics calculation in the Biocom data set
#' 
#' @author Sonia Kefi \email{sonia.kefi@umontpellier.fr}
#' 
#' @date 2023/30/11



## Install Dependencies (listed in DESCRIPTION) ----
#devtools::install_deps(upgrade = "never")

## Load Project Addins (R Functions and Packages) ----
devtools::load_all()

## Global Variables ----
path_model_data <- here::here("data")
path_output <- here::here("outputs")

NPERM <- 199 #the value used for the analyses in the paper is 199
N_SHUFFLE <- NPERM
BOOTN <- 2999
ALPHA <- 0.5

# We recommend 199 for NPERM because in the analysis we only use 
# the means of the distribution created. 199 is deemed sufficient to estimate a mean 
# with sufficient accuracy for our purposes. 

## Run Project ----
function_compute_indicators_data50(NPERM, N_SHUFFLE, path_output)
function_compute_indicators_data200(NPERM,N_SHUFFLE, path_output)
function_compute_indicators_kefimodel_all(NPERM,N_SHUFFLE,path_output,path_model_data)
function_compute_indicators_scanlonmodel(NPERM,N_SHUFFLE,path_output,path_model_data)

function_make_branches_data50(NPERM,path_output)
function_make_branches_data200(NPERM,path_output)

function_indicator_trends_data50(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA)
function_indicator_trends_data200(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA)
function_indicator_trends_kefimodel_all(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA)
function_indicator_trends_scanlonmodel(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA)


####### NOTES #######
# In the R console, type : 
#source("make.R")

# This will run the functions above (i.e. prepare the data and calculates all spatial metrics)

# This will use the files 'data_biocom.rda' and 'data_biocom_200x200.rda' to generate:
# 'biocom-grps.rda' (including information about the groups of sites)
# 'indics-data50_Nperm_199_rev.rda', 'indics-data200_Nperm_199_rev.rda' (which contain all the spatial metrics calculated on the images)
# 'indics-kefimodel_all_Nperm_199_rev.rda', 'indics-scanlonmodel_Nperm_199_rev.rda' (which contain all the spatial metrics calculated on the model simulation outcomes)
# 'trends_one_group_Nperm_199_Bootn_2999_rev.rda', 'trends_two_groups_Nperm_199_Bootn_2999_rev.rda', 'trends_three_groups_Nperm_199_Bootn_2999_rev.rda' 
# (which contain the trends in spatial metrics along the aridity gradient in the data (all or by group), and in the models).

# All files are in the folder 'outputs'. 

# To plot the figures of the paper, run the files 'plot_fig.R' in the folder 'analyses'. 




