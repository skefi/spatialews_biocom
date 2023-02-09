#' Function which computes all the indicators
#'
#' @param NPERM 
#' @param path_output 
#' @param path_model_data
#'
#' @return model_indics_df
#' @export 
#'
#' @examples function_compute_model_indicators(NPERM,path_output,path_model_data)


function_compute_model_indicators <- function(NPERM,path_output,path_model_data){
 
  if ( packageVersion("spatialwarnings") < "2.99.99" ) { 
    stop("This code requires spatialwarnings >= 2.99.99. Please update the package.")
  }
  
  # Read the data 
  datalist <- dir(path_model_data, full.names = TRUE, pattern = "*.txt") 
  
  matrices <- lapply(datalist, read.table) # step 2: read the files - result is a data frame
  matrices <- lapply(matrices, as.matrix) # step 3: convert the data frames to matrices
  # step 4: convert to logicals (so that the package knows that there are only 2 values
  # possible in the matrices); true means 1, false means 0
  matrices <- lapply(matrices, function(x){ x==1 }) 
  # Discard column/row names
  matrices <- lapply(matrices, function(o) { 
    colnames(o) <- rownames(o) <- NULL; 
    o
  })
  names(matrices) <- gsub(".txt$", "", gsub("^CA_PotAna_", "", basename(datalist)))
  
  b <- as.numeric(names(matrices))
  aridity <- 1-b
  model_params <- data.frame(b = b, aridity = aridity)
  matrices <- matrices[names(matrices)] 
  
  save(model_params, matrices, file = file.path(path_output,"data_model.rda"))
  
  model_indics <- spatialwarnings::compute_indicator(matrices, 
                                    fun = all_indicf_with_fl, 
                                    taskname = "All model indicators")
  
  # Compute null values 
  model_indics <- spatialwarnings::indictest(model_indics, nulln = NPERM, null_method = "perm")
  
  # Add gradient information to each data.frame
  model_indics_df <- plyr::ldply(seq.int(length(matrices)), function(i) { 
    data.frame(model_params[i, ], model_indics[[i]], check.names = FALSE)
  })
  
  # We do not rely on the z-scores scores computed for each matrix, because the 
  # SD of the null distribution for many indicators depend on the size of the
  # matrix, and we have matrices of different sizes. 
  # 
  # Instead, we compute the mean and sd of the null distribution, and use that 
  # to derive a z-score for the indicator. This overwrites the z-score value 
  # reported by spatialwarnings. 
  # 
  model_indics_df <- plyr::ddply(model_indics_df, ~ indic, function(df) { 
    mu    <- mean(df[ ,"null_mean"], na.rm = TRUE)
    sigma <- sd(df[ ,"null_mean"], na.rm = TRUE)
    df[ ,"z_score"] <- ( df[ ,"value"] - mu ) / sigma
    return(df)
  })
  
  # Standardize indicators value and null. We standardize them separately, which should 
  # not affect the slopes, but the overall mean of a given indicator will be equal to 
  # its null mean, because it will be zero, by definition of the standardization.
  model_indics_df <- plyr::ddply(model_indics_df, ~ indic, plyr::mutate, 
                           value_scaled = stdz(value, na.rm = TRUE), 
                           null_scaled = stdz(null_mean, na.rm = TRUE), 
                           diff = value - null_mean, 
                           diff_scaled = stdz(value - null_mean))
  
  # Add the standardization taking the mean/sd of { obs + null }, to allow their slopes 
  # to be comparable 
  model_indics_df <- plyr::ddply(model_indics_df, ~ indic, plyr::mutate, 
                           value_scaled_orig = stdz(value, 
                                                    mean(c(value, null_mean), na.rm = TRUE), 
                                                    sd(c(value, null_mean), na.rm = TRUE)), 
                           null_mean_scaled_orig = stdz(null_mean, 
                                                    mean(c(value, null_mean), na.rm = TRUE), 
                                                    sd(c(value, null_mean), na.rm = TRUE)))
  
  # Save the results for future reuse
  filename <- paste0("indics-model_Nperm_", NPERM,".rda")
  save(model_indics_df, file = file.path(path_output,filename))
  
  return(file.path(path_output,filename))

} # end function
