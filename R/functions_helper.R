#' Standardize metrics
#'
#' @param X 
#' @param mu 
#' @param sigma 
#' @param na.rm 
#'
#' @return ( X - mu ) / sigma
#' @export 
#'
#' @examples stdz()
stdz <- function(X, 
                 mu = mean(X, na.rm = na.rm), 
                 sigma = sd(X, na.rm = na.rm), 
                 na.rm = FALSE) { 
  return(( X - mu ) / sigma)
}



#' Takes the flowlength indicator. flowlength does not only depend on the matrix, but also on external variables (slope, cell_size). There is no way at the moment to pass those matrix-by-matrix, so we read them from attributes attached to the matrix. 
#'
#' @param m 
#'
#' @return flowlength
#' @export
#'
#' @examples fl_indic(m)
fl_indic <- function(m) { 
  slope <- attr(m, "slope") 
  cell_size <- attr(m, "cell_size") 
  return(spatialwarnings::raw_flowlength_uniform(m, slope, cell_size))
}



#' Function that computes indicators
#'
#' @param mat 
#'
#' @return all indic values
#' @export
#'
#' @examples all_indicf(mat)
all_indicf <- function(mat) { 
  #   sink() # close sink so all messages appear directly on the console
  #   cat(".")
  genindics <- c(spatialwarnings::raw_cg_variance(mat, subsize = 4), 
                 spatialwarnings::raw_cg_moran(mat, subsize = 1), 
                 spatialwarnings::raw_sdr(mat, 
                         # These two numbers deserve to be adjusted. E.g. by setting 
                         # the high range to c(0, 1) instead of the default c(0.8, 1)
                         sdr_low_range = c(0, 0.2), 
                         sdr_high_range = c(0, 1)), 
                 cover = mean(mat),
                 cv = sqrt(spatialwarnings::raw_cg_variance(mat, subsize = 4)) / mean(mat)
  )
  psd <- spatialwarnings::patchsizes(mat)
  sizemat <- dim(mat)[1]*dim(mat)[2]
  psd_based_indics <- c(logfmaxpatch = log10(max(psd)/sizemat),
                        # /!\ slope and cutoff are always positive. But their negative 
                        # values are used in the fit, e.g. y ~ x^{-slope}, or 
                        # y ~ x^{-slope]e^{-cutoff}
                        slope    = spatialwarnings::tpl_fit(psd)[['plexpo']],
                        cutoff   = spatialwarnings::tpl_fit(psd)[['cutoff']],
                        # raw_plrange returns its own name for the vector
                        spatialwarnings::raw_plrange(mat))
  return(c(genindics, psd_based_indics))
}



#' Function used for computing model indicators. We define a function that computes all indicators + the flowlength. We had to do that separately before because we had to extract the information on the slope from the biocom database. Here on model data it is different because we have no slope information: we just assume the slope is 8.707794 (which is the average slope in the biocom dataset, and cell size is 1 (i.e. flowlength is defined in terms of cell size). 
#'
#' @param m 
#'
#' @return
#' @export
#'
#' @examples
all_indicf_with_fl <- function(m) { 
  assumed_slope <- 8.707794
  fl <- spatialwarnings::raw_flowlength_uniform(m, slope = assumed_slope, cell_size = 1)
  
  names(fl) <- "flowlength"
  c( all_indicf(m), fl )
}


#' Cleans up data set
#'
#' @param x 
#'
#' @return x
#' @export
#'
#' @examples tidydf(x)
tidydf <- function(x) { 
  as.data.frame( broom::tidy(x) ) 
}



#' Return the sign of something if significant. Used to set colors in graphs.
#'
#' @param x 
#' @param alpha 
#'
#' @return a
#' @export
#'
#' @examples signf(x, alpha)
signf <- function(x, alpha) { 
  a <- ifelse(x > (1-alpha), "pos", ifelse(x < alpha, "neg", "ns"))
  a <- factor(a, levels = c("pos", "neg", "ns"))
}



#' Associates number of stars to a p-value
#'
#' @param X 
#' @param NAval 
#'
#' @return a
#' @export
#'
#' @examples pvalstar(X,NAval)
pvalstar <- function(X, NAval = NA_character_) { 
  if ( length(X) > 1 ) { 
    return( sapply(X, pvalstar) ) 
  }
  
  a <- ""
  if ( is.na(X) ) return( NAval )
  if ( X < 0.10 )  { a <- "." } 
  if ( X < 0.05 )  { a <- "*" } 
  if ( X < 0.01 )  { a <- "**" } 
  if ( X < 0.001 ) { a <- "***" } 
  return(a)
}



#' Get the slope of a set of samples contained in dat, identified by samples. The
# structure of this function (arguments) is dictated by R(s bootstrap boot() function.
#'
#' @param dat 
#' @param samples 
#'
#' @return a
#' @export
#'
#' @examples get_lm(dat, samples)
get_lm <- function(dat, samples) { 
  stopifnot(dat[1, "indic"] %in% names(indic_groups))
  
  # We return NA if we could not make the fit 
  out <- NA_real_
  
  # For FL, regular lm()
  if ( dat[1, "indic"] == "flowlength" ) { 
    # Standard model: 
    a <- coef(lm(indic_value ~ Aridity, data = dat[samples, ]))
    out <- a[2] # Return slope only
  } else { 
    # For other indicators, mixed model. 
    a <- try({ 
      # This was wrong, it was 'tbl' instead of dat[samples, ], which produced wrong 
      # results (same value each time because each time we were using the same 
      # original table). 
      lme4::lmer(indic_value ~ Aridity + (1 | plotn), data = dat[samples, ])
    })
    if ( ! inherits(a, "try-error") ) { 
      out <- a@beta[2] # Return slope only
    } 
  }
  return(out)
}



#' Function that given a dataset, will return BOOTN estimates of the slope along 
# aridity
#'
#' @param tbl 
#' @param ... 
#'
#' @return 
#' @export
#'
#' @examples get_slope_along_aridity(tbl)
get_slope_along_aridity <- function(tbl, ...) {
  # NOTE: boot() with multicore uses a lot of memory, so we use only 4 cpus
  bootstrap_ests <- boot::boot(tbl, get_lm, R = BOOTN, 
                               ncpus = 4, parallel = "multicore")
  data.frame(estimate = bootstrap_ests[["t"]])
}



#' Function that estimates a slope on model data. There is no mixed effect because no need on model data, and we use aridity instead of Aridity because that is what is in the model data.We need to have the two arguments, dat and samples, because those are imposed by the boot() function given in base R. 
#'
#' @param dat 
#' @param samples 
#'
#' @return a[2]
#' @export
#'
#' @examples get_lm_model(dat, samples)
get_lm_model <- function(dat, samples) { 
  a <- coef(lm(indic_value ~ aridity, data = dat[samples, ]))
  a[2] # return the slope
}



#' Function that given a dataset, will return BOOTN estimates of the slope along 
# aridity for the model data
#'
#' @param tbl 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples get_slope_along_aridity_model(tbl)
get_slope_along_aridity_model <- function(tbl, ...) {
  bootstrap_ests <- boot::boot(tbl, get_lm_model, R = BOOTN, 
                               ncpus = min(16, parallel::detectCores()), 
                               parallel = "multicore")
  data.frame(estimate = bootstrap_ests[["t"]])
}

