#' Function which classifies sites into branches (or groups of dryland sites)
#'
#' @param NPERM 
#' @param path_output 
#'
#' @return indics
#' @export 
#'
#' @examples function_make_branches(NPERM,path_output)


function_make_branches <- function(NPERM,path_output){
  
  load(file.path(path_output,"data_biocom.rda"))
  filename <- paste0("indics-data_Nperm_", NPERM,".rda")
  load(file.path(path_output,filename))
  
  biocom  <- ourdata$biocom
  
  # Compute and make clusters 
  clust_dat <- ourdata[["biocom"]][ ,c("MF", "imgcover")]
  clust_dat <- apply(clust_dat, 2, function(X) ( X - mean(X) ) / sd(X) )
  
  hcl <- stats::hclust(stats::dist(clust_dat), method = "ward.D")
  grps3 <- stats::cutree(hcl, k = 3) # 3 groups mf+cover
  grps2 <- stats::cutree(hcl, k = 2) # 2 groups mf+cover
  
  # Get pretty names for each branch 
  cov_means <- tapply(clust_dat[ ,"imgcover"], grps3, mean)
  cover_grps <- ifelse(cov_means == max(cov_means), 
                       "high cover", "low cover")
  mf_means <- tapply(clust_dat[ ,"MF"], grps3, mean)
  mf_grps <- ifelse(mf_means == min(mf_means), 
                    "low MF", "high MF")
  grps_names <- paste(cover_grps, mf_grps, sep = " - ")
  pretty_grps3 <- grps_names[grps3]
  
  cov_means <- tapply(clust_dat[ ,"imgcover"], grps2, mean)
  cover_grps <- ifelse(cov_means == max(cov_means), 
                       "high cover", "low cover")
  grps_names <- paste(cover_grps, mf_grps, sep = " - ")
  pretty_grps2 <- grps_names[grps2]
  
  hcl <- stats::hclust(dist(clust_dat[ ,"imgcover", drop = FALSE]), method = "ward.D")
  grps2c <- stats::cutree(hcl, k = 2)
  cov_means <- tapply(clust_dat[ ,"imgcover"], grps2c, mean)
  cover_grps <- ifelse(cov_means == max(cov_means), 
                       "high cover", "low cover")
  pretty_grps2c <- cover_grps[grps2c]
  
  # Add group columns to original dataset + indicators 
  ourdata[["biocom"]][ ,"grps3"] <- grps3
  ourdata[["biocom"]][ ,"grps2"] <- grps2
  ourdata[["biocom"]][ ,"grps2c"] <- grps2c
  ourdata[["biocom"]][ ,"pretty_grps2"] <- pretty_grps2
  ourdata[["biocom"]][ ,"pretty_grps3"] <- pretty_grps3
  ourdata[["biocom"]][ ,"pretty_grps2c"] <- pretty_grps2c
  
  biocom <- cbind(biocom,grps3)
  biocom <- cbind(biocom,grps2)
  biocom <- cbind(biocom,grps2c)
  biocom <- cbind(biocom,pretty_grps3)
  biocom <- cbind(biocom,pretty_grps2)
  biocom <- cbind(biocom,pretty_grps2c)
  
  indics$plotn <- as.factor(indics$plotn)
  
  # Joining by plotid 
  extr <- ourdata[["biocom"]][ ,c("plotid",
                                  "pretty_grps2", "pretty_grps2c", "pretty_grps3", 
                                  "grps2c","grps2","grps3")]
  indics <- plyr::join(indics, extr, by = "plotid")
  
  
  # Save the data into a output_data folder for distribution with the results
  filename <- paste0("indics-data-grps_Nperm_", NPERM,".rda")
  save(indics, file = file.path(path_output,filename))
  
  save(biocom, file = file.path(path_output,"biocom-grps.rda"))
  
} # end function
  