#' Function which computes all the spatial metrics in the 50 x 50 m images (three per site in 115 sites of the Biocom data set)
#'
#' @param NPERM
#' @param N_SHUFFLE
#' @param path_output
#'
#' @return indics
#' @export
#'
#' @examples function_compute_indicators_data50(NPERM,N_SHUFFLE,path_output)


function_compute_indicators_data50 <- function(NPERM,N_SHUFFLE,path_output){

  if ( packageVersion("spatialwarnings") < "2.99.99" ) {
    stop("This code requires spatialwarnings >= 2.99.99. Please update the package.")
  }

  load(file.path(path_output,"data_biocom.rda"))
  
  ourdata[["biocom"]][["matrixn"]] <- seq_len(nrow(ourdata[["biocom"]]))

  reset_workers(workers=parallel::detectCores()-1)
  indic_results <- spatialwarnings::compute_indicator(
    ourdata[["matrices"]],
    fun = all_indicf,
    N_SHUFFLE = N_SHUFFLE,
    taskname = "All indicators"
  )

  # Compute null values
  # NOTE: using many workers can use a lot of memory
  # browser()
  # Compute observed values
  # ~1.3m with reltol = 1e-5 and maxit = 1e8 above for the full dataset
  # ~3m with reltol = 1e-6 and maxit = 1e8 above for the full dataset
  # ~15m with reltol = 1e-7 and maxit = 1e8 above for the full dataset
  # (hardcoded parameters in CRAN spatialwarnings)
  reset_workers(workers=parallel::detectCores()-1)
  indic_results <- spatialwarnings::indictest(indic_results,
                                              nulln = NPERM,
                                              null_method = "perm")

  # Compute the flowlength: the flowlength is a bit special because it is
  # only computed for one plot, as it requires information about the slope and
  # resolution of the image, so we do it a bit differently here
  has_info <- which( ! is.na(ourdata[["biocom"]][ , "resol"]) )

  fl_results <- plyr::ldply(has_info, function(i) {
    ic_obs <- spatialwarnings::raw_flowlength_uniform(ourdata[["matrices"]][[i]],
                                                      slope = ourdata[["biocom"]][i,"SLO"],
                                                      cell_size = ourdata[["biocom"]][i,"resol"])

    # We need to compute the nulls by hand because of the issue that spatialwarnings
    # cannot pass vectorized arguments to subfunctions.
    nulls <- future.apply::future_lapply(seq.int(NPERM), function(nrep) {
      m <- ourdata[["matrices"]][[i]]
      m[] <- sample(m)
      spatialwarnings::raw_flowlength_uniform(m,
                                              slope = ourdata[["biocom"]][i,"SLO"],
                                              cell_size = ourdata[["biocom"]][i,"resol"])
    }, future.seed = TRUE)
    nulls <- do.call(c, nulls)
    data.frame(indic = "flowlength",
               ourdata[["biocom"]][i,c("plotn", "plotid", "matrixn")],
               value = ic_obs, null_mean = mean(nulls), null_sd = sd(nulls),
               # We don't use those, so set them to NA
               z_score = NA, pval = NA, null_qsup = NA, null_qinf = NA,
               row.names = NULL)
  })

  # Convert indicator results to data frame and merge FL dataset with the rest of the
  # indicators
  indics_df <- as.data.frame(indic_results)
  indics_df[ ,"plotn"] <- ourdata[["biocom"]][indics_df[ ,"matrixn"], "plotn"]
  indics_df  <- dplyr::left_join( #add plotid image by image:
    indics_df,
    dplyr::select(ourdata[["biocom"]], # for all images
                  matrixn, plotid), # select matrixn and plotid
    by = "matrixn" # join by matrixn
  )
  # adds NA for columns that are missing in fl_results
  stopifnot( length(setdiff(names(fl_results), names(indics_df))) == 0 )
  stopifnot( all(sort(names(fl_results)) == sort(names(indics_df))) )
  indics_df <- rbind(indics_df, fl_results[ ,names(indics_df)])

  indics_df[ ,"plotn"] <- as.character(indics_df[ ,"plotn"])
  indics_df <- indics_df[with(indics_df, order(plotn, indic)), ]

  indics_df <- indics_df[ , ! names(indics_df) %in% c("matrixn")]

  indics_df <- plyr::join(indics_df, ourdata[["biocom"]], by = c("plotid", "plotn"))

  # We don't need matrixn anymore
  indics_df  <- dplyr::select(indics_df, ! matrixn)

  # We do not rely on the z-scores scores computed for each matrix, because the
  # SD of the null distribution for many indicators depends on the size of the
  # matrix, and we have matrices of different sizes.
  #
  # Instead, we compute the mean and sd of the null distribution, and use that
  # to derive a z-score for the indicator. This overwrites the z-score value
  # reported by spatialwarnings.
  #
  indics_df <- plyr::ddply(indics_df, ~ indic, function(df) {
    mu    <- mean(df[ ,"null_mean"], na.rm = TRUE)
    sigma <- sd(df[ ,"null_mean"], na.rm = TRUE)
    df[ ,"z_score"] <- ( df[ ,"value"] - mu ) / sigma
    return(df)
  })

  # Standardize indicators value and null. We standardize them separately, which should
  # not affect the slopes, but the overall mean of a given indicator will be equal to
  # its null mean, because it will be zero, by definition of the standardization.
  indics_df <- plyr::ddply(indics_df, ~ indic, plyr::mutate,
                           value_scaled = stdz(value, na.rm = TRUE),
                           null_scaled = stdz(null_mean, na.rm = TRUE),
                           diff = value - null_mean,
                           diff_scaled = stdz(value - null_mean))

  # Add the standardization taking the mean/sd of { obs + null }, to allow the slopes to
  # be comparable within an indicator between its observed values and the null
  indics_df <- plyr::ddply(indics_df, ~ indic, plyr::mutate,
                           value_scaled_orig = stdz(value,
                                                    mean(c(value, null_mean), na.rm = TRUE),
                                                    sd(c(value, null_mean), na.rm = TRUE)),
                           null_mean_scaled_orig = stdz(null_mean,
                                                        mean(c(value, null_mean), na.rm = TRUE),
                                                        sd(c(value, null_mean), na.rm = TRUE)))

  # Save the results into output_data folder
  indics <- indics_df # rename
  filename <- paste0("indics-data50_Nperm_", NPERM,"_rev.rda")
  save(indics, file = file.path(path_output,filename))

  return(file.path(path_output,filename))

} # end function compute indicators
