#' Function which computes all the spatial metrics on the outcomes of simulations from the model of Kéfi et al. 2007, Theoretical Population Biology
#'
#' @param NPERM
#' @param N_SHUFFLE
#' @param path_output
#' @param path_model_data
#'
#' @return model_indics_df
#' @export
#'
#' @examples function_compute_indicators_kefimodel_all(NPERM,N_SHUFFLE,path_output,path_model_data)


function_compute_indicators_kefimodel_all <- function(NPERM,N_SHUFFLE,path_output,path_model_data){

  if ( packageVersion("spatialwarnings") < "2.99.99" ) {
    stop("This code requires spatialwarnings >= 2.99.99. Please update the package.")
  }

  # Read the data
  kefi_runs <- readRDS(file.path(path_model_data,"aridvege_ca_runs.rds"))
  kefi_runs[[1]]$indics
  kefi_runs[[1]]$mat

  # Extrait les matrices
  matrices <- lapply(kefi_runs, function(o) o[["mat"]])
  model_params <- plyr::ldply(kefi_runs, function(o) o[["covers"]])

  # These constants are here to determine how TPL functions are being fitted.
  # Fitting a TPL requires computing sum_{k=1}^{inf} k^(-a) exp(-b k). This
  # sum is infinite so we need some stopping criterion, here either
  # a relative change in the sum below a threshold (reltol below), or a
  # maximum number of iterations (1e8). These are fixed values in
  # spatialwarnings <= 3.0.3, and can be set via global options in versions
  # above (still in dev on 2023-07-27) as we do below.
  # Decrease reltol and increase maxit below to increase precision of
  # computations, at the cost of higher computation time.
  options(spatialwarnings.constants.reltol = 1e-6)
  options(spatialwarnings.constants.maxit = 1e8)

  reset_workers(workers=parallel::detectCores()-1)
  model_indics <- spatialwarnings::compute_indicator(matrices,
                                                     fun = all_indicf_with_fl,
                                                     N_SHUFFLE = N_SHUFFLE,
                                                     taskname = "All model indicators")

  # Compute null values
  reset_workers(workers=parallel::detectCores()-1)
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
  filename <- paste0("indics-kefimodel_all_Nperm_", NPERM,"_rev.rda")
  save(model_indics_df, file = file.path(path_output,filename))

  return(file.path(path_output,filename))

} # end function
