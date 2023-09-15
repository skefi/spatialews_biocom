#' Function which the trend (slope) of each of the spatial metrics along the aridity gradient in the 50 x 50 images, for the whole data set or per group of sites
#'
#' @param NPERM
#' @param N_SHUFFLE
#' @param path_output
#' @param BOOTN
#' @param ALPHA
#'
#' @return trends
#' @export
#'
#' @examples function_indicator_trends_data50(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA)


function_indicator_trends_data50 <- function(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA){

  #---------------------------------------------------------------------------
  # LOAD data
  #---------------------------------------------------------------------------
  filename <- paste0("indics-data50-grps_Nperm_", NPERM, "_rev.rda")
  load(file.path(path_output,filename))

   # What type of indicator value to use to compute and slopes.
  # Here we use the standardized obs and the standardized null, but using the same mean
  # and s.d., so that we can compare the slopes between obs and null.
  SLOPE_OBS_TYPE <- "value_scaled_orig"
  SLOPE_NULL_TYPE <- "null_mean_scaled_orig"
  # Switch to below to not used the standardize values:
  #SLOPE_OBS_TYPE <- "value"
  #SLOPE_NULL_TYPE <- "null_mean"

  # Groups to which each indicator belongs
  indic_groups <- c(moran    = "Generic",
                    sdr      = "Generic",
                    cv.variance ="Generic",
                    fmaxpatch = "Patch-based",
                    logfmaxpatch = "Patch-based",
                    cutoff = "Patch-based",
                    slope    = "Patch-based",
                    flowlength  = "Hyd.")

  # Order of indicators in plots
  indic_order <- c('logfmaxpatch', 'fmaxpatch', 'slope', 'cutoff', 'flowlength','moran', 'sdr','cv.variance') 
  
   # Format indicator results
  indics_fmt <- tidyr::gather(indics, "indic_value_type", "indic_value",
                       all_of(c("value", "null_mean",
                                "value_scaled_orig",
                                "null_mean_scaled_orig"))) # weird rlang syntax

  # Remove indicators that are not needed
  indics_fmt <- subset(indics_fmt, indic != "cover")

 
  
  #---------------------------------------------------------------------------
  # TRENDS IN THE WHOLE DATA SET
  #---------------------------------------------------------------------------


  indics_fmt[ ,"indic"] <- as.factor(indics_fmt[ ,"indic"])
  indics_fmt[ ,"indic_value_type"] <- as.factor(indics_fmt[ ,"indic_value_type"])
  indics_fmt[,"plotn"] <- as.numeric(indics_fmt[ ,"plotn"])


  trends <- plyr::ddply(indics_fmt, ~ indic + indic_value_type,
                  function(df) {
                    df <- subset(df, is.finite(indic_value))
                    if ( nrow(df) < 3 ) { # not enough data points
                      return( data.frame(estimate = NA) )
                    }
                    get_slope_along_aridity(df)
                  }, .progress = "time")


  # Order indicators for plotting
  trends[ ,"indic_order"] <- factor(trends[ ,"indic"],
                                    levels = rev(indic_order),
                                    ordered = TRUE)

  # NOTE: the as.character() here so that indic_groups is indexed by name and not by
  # the conversion of a factor to integer... (wtf R).
  trends[ ,"indic_group"] <- factor(indic_groups[as.character(trends[ ,"indic"])],
                                    levels = c("Patch-based","Hyd.","Generic"))

  # Compute median slope
  # Here trends contain BOOTN estimates of slopes, effectively a distribution. We
  # compute here the median of that distribution, its proportion below zero (0 means
  # all values below zero = significant, 1 means all values above zero = significant too),
  # and the number of non-NA values (n_estimate, which should be close to BOOTN).
  # n_estimate is not always equal to bootn because for some resamplings the fit can
  # fail, but it's okay I suppose as long as n_estimate is large (>=2000).
  trends_med <- plyr::ddply(trends,
                      ~ indic_value_type + indic + indic_group + indic_order,
                      plyr::summarise,
                      med_estimate = median(estimate, na.rm = TRUE),
                      pval_estimate = mean(estimate < 0, na.rm = TRUE),
                      n_estimate = sum(! is.na(estimate)),
                      q05 = quantile(estimate, .05, na.rm = TRUE),
                      q95 = quantile(estimate, .95, na.rm = TRUE))

  if ( any(trends_med[ ,"n_estimate"] < (BOOTN * 3/4)) ) {
    warning("The number of estimate values in the bootstrap distribution is quite low")
  }

  # We add the trend information to the original data.frame
  trends_all <- plyr::join(trends, trends_med, type = "left")

  # We compute the differences between the trend in observed landscapes and null
  # landscapes. We are effectively computing differences between two distributions,
  # of approx. BOOTN values each.
  # We compute the average difference (diff_obs_m_null), and the probability that a
  # randoml picked value of one is below the other (pnull_inf_obs). The latter represents
  # a probability between zero and one, zero when all null slopes are below their
  # observed counterpart, one when they are all above.
  trends_diffs <-
    plyr::ddply(subset(trends_all, indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)),
          ~ indic_group + indic_order + indic, #+ pretty_grps
          function(df) {
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(pnull_inf_obs = mean(x1 < x2, na.rm = TRUE),
                       diff_obs_m_null = mean(x2 - x1, na.rm = TRUE),
                       n = sum(!is.na(x1)))
          })

  # We compute the differences, but this time we do not average. We obtain a distribution
  # of BOOTN differences between obs and null. This allows keeping the whole distribution,
  # which quantiles are then represented in the graph by stat_pointinterval(). We could
  # manually compute quantiles and represent them in the graph, but stat_pointinterval()
  # does it for us so here we go.
  trends_diffs_distribs <-
    plyr::ddply(subset(trends_all, indic_value_type %in% c(SLOPE_OBS_TYPE,SLOPE_NULL_TYPE)),
          ~ indic_group + indic_order + indic,
          function(df) {
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(difference = x2 - x1)
          })

  # We create a data.frame which contains all the distributions of differences, with
  # added columns with summary stats (pnull_inf_obs, diff_obs_m_null), so that we can
  # plot the distributions, and color them by significance.
  trends_diffs_distribs <- plyr::join(trends_diffs_distribs, trends_diffs,
                                type = "left", match = "first")






  #---------------------------------------------------------------------------
  # TRENDS per 2 branches of cover
  #---------------------------------------------------------------------------


  # Compute trends branch-by-branch
  trends2 <-
    plyr::ddply(indics_fmt, ~ indic + pretty_grps2 + indic_value_type,
          function(df) {
            df <- subset(df, is.finite(indic_value))
            if ( nrow(df) < 3 ) { # not enough data points
              return( data.frame(estimate = NA) )
            }
            get_slope_along_aridity(df)
          }, .progress = "time")


  trends2[ ,"indic_order"] <- factor(trends2[ ,"indic"],
                                     levels = rev(indic_order),
                                     ordered = TRUE)

  trends2[ ,"indic_group"] <- factor(indic_groups[as.character(trends2[ ,"indic"])],
                                     levels = c("Patch-based", "Hyd.","Generic"))

  # Display trends along gradient or raw value
  trends2_med <- plyr::ddply(trends2, ~ indic_value_type + indic_group + pretty_grps2 +
                         indic_order + indic,
                       plyr::summarise,
                       med_estimate = median(estimate, na.rm = TRUE),
                       pval_estimate = mean(estimate < 0, na.rm = TRUE),
                       n_estimate = sum(!is.na(estimate)),
                       q05 = quantile(estimate, .05, na.rm = TRUE),
                       q95 = quantile(estimate, .95, na.rm = TRUE))

  if ( any(trends2_med[ ,"n_estimate"] < (BOOTN * 3/4)) ) {
    warning("The number of estimate values in the bootstrap distribution is quite low")
  }


  trends2_all <- plyr::join(trends2, trends2_med, type = "left")

  # We compute the differences between the trend in observed landscapes and null
  # landscapes
  trends2_diffs <-
    plyr::ddply(subset(trends2_all, indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)),
          ~ indic_group + pretty_grps2 + indic_order + indic,
          function(df) {
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(pnull_inf_obs = mean(x1 < x2, na.rm = TRUE),
                       diff_obs_m_null = mean(x2 - x1, na.rm = TRUE),
                       n = sum(!is.na(x1)))
          })

  trends2_diffs_distribs <-
    plyr::ddply(subset(trends2_all, indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)),
          ~ indic_group + pretty_grps2 + indic_order + indic,
          function(df) {
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(difference = x2 - x1)
          })

  # Join the summary stats to the original distributions, so we can color points,
  # etc in stat_pointinterval()
  trends2_diffs_distribs <- plyr::join(trends2_diffs_distribs, trends2_diffs,
                                 type = "left", match = "first")





  #---------------------------------------------------------------------------
  # TRENDS for 3 branches
  #---------------------------------------------------------------------------

  # Compute trends branch-by-branch
  trends3 <- plyr::ddply(indics_fmt, ~ indic + pretty_grps3 + indic_value_type,
                   function(df) {
                     df <- subset(df, is.finite(indic_value))
                     if ( nrow(df) < 3 ) { # not enough data points
                       return( data.frame(estimate = NA) )
                     }
                     get_slope_along_aridity(df)
                   }, .progress = "time")

  trends3[ ,"indic_order"] <- factor(trends3[ ,"indic"],
                                     levels = rev(indic_order),
                                     ordered = TRUE)

  trends3[ ,"indic_group"] <- factor(indic_groups[as.character(trends3[ ,"indic"])],
                                     levels = c("Patch-based", "Hyd.","Generic"))

  # Display trends along gradient or raw value
  trends3_med <- plyr::ddply(trends3,
                       ~ indic_value_type + indic_group + pretty_grps3 +
                         indic_order + indic,
                       plyr::summarise,
                       med_estimate = median(estimate, na.rm = TRUE),
                       pval_estimate = mean(estimate < 0, na.rm = TRUE),
                       n_estimate = sum(!is.na(estimate)),
                       q05 = quantile(estimate, .05, na.rm = TRUE),
                       q95 = quantile(estimate, .95, na.rm = TRUE))

  if ( any(trends3_med[ ,"n_estimate"] < (BOOTN * 3/4)) ) {
    warning("The number of estimate values in the bootstrap distribution is quite low")
  }

  trends3_all <- plyr::join(trends3, trends3_med, type = "left")

  # We compute the differences between the trend in observed landscapes and null
  # landscapes
  trends3_diffs <-
    plyr::ddply(subset(trends3_all, indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)),
          ~ indic_group + pretty_grps3 + indic_order + indic,
          function(df) {
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(pnull_inf_obs = mean(x1 < x2, na.rm = TRUE),
                       diff_obs_m_null = mean(x2 - x1, na.rm = TRUE),
                       n = sum(!is.na(x1)))
          })

  trends3_diffs_distribs <-
    plyr::ddply(subset(trends3_all, indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)),
          ~ indic_group + pretty_grps3 + indic_order + indic,
          function(df) {
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(difference = x2 - x1)
          })

  # Join the summary stats to the original distributions, so we can color points,
  # etc in stat_pointinterval()
  trends3_diffs_distribs <- plyr::join(trends3_diffs_distribs, trends3_diffs,
                                 type = "left", match = "first")




  #---------------------------------------------------------------------------
  # SAVE outputs
  #---------------------------------------------------------------------------

  filename = paste0("trends_data50_one_group_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
  save(trends, trends_med, trends_all, trends_diffs, trends_diffs_distribs, file = file.path(path_output,filename))

  filename2 = paste0("trends_data50_two_groups_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
  save(trends2, trends2_med, trends2_all, trends2_diffs, trends2_diffs_distribs,
       file = file.path(path_output,filename2))

  filename3 = paste0("trends_data50_three_groups_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
  save(trends3, trends3_med, trends3_all, trends3_diffs, trends3_diffs_distribs,
         file = file.path(path_output,filename3))


  return(file.path(path_output,filename))

} # end function
