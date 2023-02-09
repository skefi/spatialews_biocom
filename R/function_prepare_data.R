#' This file reads the Biocom data from files and merges them into suitable tables for this analysis. Matrices are read and put in the same order as rows in the biocom database. 
#'
#' @return ourdata.rda
#' @export 
#'
#' @examples function_prepare_data(path_data, path_output)


function_prepare_data <- function(path_data, path_output){
  
  file_biocom <- "biocom_table.txt"
  file_fl <- "biocom_fl.xlsx"
  
  biocom <- read.table(file.path(path_data,file_biocom))
  
  # Read biocom_fl file, which contains
  # information about the image resolution, that we are interested in, 
  # and which is not in biocom_table
  biocom_fl <- readxl::read_excel(file.path(path_data,file_fl))
  
  biocom_fl <- biocom_fl[ ,c("plotn", "Resol")]
  biocom_fl <- as.data.frame(biocom_fl)
  names(biocom_fl) <- tolower(names(biocom_fl))
  
  # Import the resolution data. Note that we set it to NA for all the subplots
  # except the first one as this is the one to which resolution information refers 
  biocom <- plyr::join(biocom, biocom_fl, type = "left", match = "first")
  biocom <- plyr::mutate(biocom, resol = ifelse(duplicated(plotn), NA, resol))
  
  ## Read matrices
  matdir <- file.path(path_data,"data_images_biocom")
  # Scan for matrices in the directory, then load them and convert them to actual 
  # matrix R objects (spw needs matrices to work).
  datalist <- dir(matdir, full.names = TRUE, pattern = "*.txt") 
  matrices <- lapply(datalist, read.table)
  matrices <- lapply(matrices, as.matrix) 
  # Converts matrices to binary matrices (so that the package knows that there are only 2
  # values possible in the matrices); true means 1, false means 0
  matrices <- lapply(matrices, function(x) { x == 1 } )
  # Strip matrix names for easier reading on the terminal
  matrices <- lapply(matrices, function(x) { dimnames(x) <- NULL; x })
  
  # We add names to the list of matrices, so that we know explicitely what plotn/subplot
  # the refer to.
  files <- basename(datalist)
  files2 <- gsub(".txt", "", files)
  names(matrices) <- files2
 
  stopifnot( all(gsub("-.+$", "", names(matrices)) == biocom[ ,"plotn"]) )
  biocom[ ,"file"] <- files2 
  
  # Extract plot number and sort it the same way as the list of matrices, so that 
  # one line in biocom corresponds to one matrix in the list of matrices
  files3 <- as.numeric(gsub("[A-z \\.\\(\\)\\-]", "", files))
  biocom <- lapply(unique(files3), function(x) { subset(biocom, plotn == x) })
  biocom <- do.call(rbind, biocom)
  
  # Compute the mean cover from the matrices and add it to the biocom dataset
  biocom[ ,"imgcover"] <- plyr::laply(matrices, function(o) mean(o))
  
  # Compute number of pixels per image and record it
  biocom[ ,"nbpixels"] <- plyr::laply(matrices, function(o) dim(o)[1]*dim(o)[2])
  
  # Convert the plot number to an unordered factor (instead of a numeric)
  biocom[ ,'plotn'] <- factor(as.character(biocom[ ,"plotn"]), 
                              ordered = FALSE)
  
  # Add a unique id for each of the 345 images
  biocom[ ,"plotid"] <- as.factor(1:nrow(biocom))
  
  ourdata <- list(biocom = biocom,
                  matrices = matrices)
  
  save(ourdata, file = file.path(path_output,"data_biocom.rda"))
  
  return(file.path(path_output,"data_biocom.rda"))
  
} # end function prepare data 
