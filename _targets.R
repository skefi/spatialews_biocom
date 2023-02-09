#' Targets plan
#' 
#' @author Sonia Kefi \email{sonia.kefi@umontpellier.fr}
#' 
#' @date 2023/02/07


## Attach required packages ----
library(targets)
library(tarchetypes)

## Load Project R Functions ----
source(here::here("R", "functions_helper.R"))
source(here::here("R", "function_prepare_data.R"))
source(here::here("R", "function_compute_indicators.R"))
source(here::here("R", "function_compute_model_indicators.R"))
source(here::here("R", "function_make_branches.R"))

## Analyses pipeline ----
list(
   tar_target(path_data, here::here("data")),
   tar_target(path_model_data, here::here("data/data_CA")),
   tar_target(path_output, here::here("outputs")),
   tar_target(NPERM, 3), #the value used for the analyses in the paper is 199
   tar_target(create_ourdata, function_prepare_data(path_data, path_output)),
   tar_target(compute_indic, function_compute_indicators(NPERM,path_output)),
   tar_target(compute_model_indic, function_compute_model_indicators(NPERM,path_output,path_model_data)),
   tar_target(compute_branches, function_make_branches(NPERM,path_output))
 )

# We recommend 199 for NPERM because in the analysis we only use 
# the means of the distribution created. 199 is deemed sufficient to estimate a mean 
# with sufficient accuracy for our purposes. 

