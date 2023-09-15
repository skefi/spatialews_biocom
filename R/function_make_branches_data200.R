#' Function which includes the groups of dryland sites into the 200 x 200 data set
#'
#' @param NPERM 
#' @param path_output 
#'
#' @return indics
#' @export 
#'
#' @examples function_make_branches_data200(NPERM,path_output)


function_make_branches_data200 <- function(NPERM,path_output){
  
  load(file.path(path_output,"biocom-grps.rda"))
  filename <- paste0("indics-data200_Nperm_", NPERM,"_rev.rda")
  load(file.path(path_output,filename))
  
  # Joining by plotid 
  extr <- biocom[ ,c("plotid", "pretty_grps2", "pretty_grps2c", "pretty_grps3", 
                                "grps2c","grps2","grps3")]
  indics <- plyr::join(indics, extr, by = "plotid")
  
  
  # Save the data into a output_data folder for distribution with the results
  filename <- paste0("indics-data200-grps_Nperm_", NPERM,"_rev.rda")
  save(indics, file = file.path(path_output,filename))
  
} # end function
  