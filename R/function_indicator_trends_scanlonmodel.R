#' Function which the trend (slope) of each of the spatial metrics along the aridity gradient in the model
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
#' @examples function_indicator_trends_scanlonmodel(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA)


function_indicator_trends_scanlonmodel <- function(NPERM,N_SHUFFLE,path_output,BOOTN,ALPHA){
  
  #---------------------------------------------------------------------------
  # LOAD data
  #---------------------------------------------------------------------------
  filename <- paste0("indics-scanlonmodel_Nperm_", NPERM,"_rev.rda")
  load(file.path(path_output,filename), verbose = TRUE)
  unique(model_indics_df$indic)
  
  # Format indicator results
  model_indics_fmt <- tidyr::gather(model_indics_df, "indic_value_type", "indic_value", 
                             all_of(c("z_score", "value", "null_mean", 
                                      "value_scaled", "null_scaled", 
                                      "value_scaled_orig", 
                                      "null_mean_scaled_orig"))) # weird rlang syntax
  
  # Remove indicators that are not needed 
  model_indics_fmt <- subset(model_indics_fmt, indic != "cover")
  
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
  
  model_indics_fmt_tot <- model_indics_fmt
  
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  # Model original
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  
  model_indics_fmt <- subset(model_indics_fmt_tot,simu_type=="scanlon_original")
  
  #---------------------------------------------------------------------------
  # TRENDS
  #---------------------------------------------------------------------------
  
  model_indics_fmt[ ,"indic"] <- as.factor(model_indics_fmt[ ,"indic"])
  model_indics_fmt[ ,"indic_value_type"] <- as.factor(model_indics_fmt[ ,"indic_value_type"])
  
  model_indics_fmt$aridity <- 1-model_indics_fmt$ft_star
  
  model_trends <- plyr::ddply(model_indics_fmt, ~ indic + indic_value_type, 
                        function(df) { 
                          df <- subset(df, is.finite(indic_value)) 
                          if ( nrow(df) < 3 ) { # not enough data points
                            return( data.frame(estimate = NA) )
                          }
                          get_slope_along_aridity_model(df)
                        }, .progress = "time")
  
  
  # Order indicators for plotting
  model_trends[ ,"indic_order"] <- factor(model_trends[ ,"indic"], 
                                          levels = rev(indic_order), 
                                          ordered = TRUE)
  
  # NOTE: the as.character() here so that indic_groups is indexed by name and not by 
  # the conversion of a factor to integer... (wtf R).
  model_trends[ ,"indic_group"] <- 
    factor(indic_groups[as.character(model_trends[ ,"indic"])], 
           levels = c("Patch-based","Hyd.","Generic")) #, "Composite"
  
  # Compute median slope
  # Here trends contain BOOTN estimates of slopes, effectively a distribution. We 
  # compute here the median of that distribution, its proportion below zero (0 means 
  # all values below zero = significant, 1 means all values above zero = significant too), 
  # and the number of non-NA values (n_estimate, which should be close to BOOTN). 
  # n_estimate is always equal to bootn here because the fit always succeeds on model 
  # data. 
  model_trends_med <- plyr::ddply(model_trends, 
                            ~ indic_value_type + indic + indic_group + indic_order,
                            plyr::summarise, 
                            med_estimate = median(estimate, na.rm = TRUE), 
                            pval_estimate = mean(estimate < 0, na.rm = TRUE), 
                            n_estimate = sum(! is.na(estimate)))
  
  if ( any(model_trends_med[ ,"n_estimate"] < (BOOTN * 3/4)) ) { 
    warning("The number of estimate values in the bootstrap distribution is quite low")
  }
  
  # We add the trend information to the original data.frame
  model_trends_all <- plyr::join(model_trends, model_trends_med, type = "left")
  
  # We compute the differences between the trend in observed landscapes and null
  # landscapes. We are effectively computing differences between two distributions, 
  # of exactly BOOTN values each. 
  # We compute the average difference (diff_obs_m_null), and the probability that a 
  # randomly picked value of one sample is below the other (pnull_inf_obs). The latter
  # represents a probability between zero and one, zero when all null slopes are below 
  # their observed counterpart, one when they are all above. 
  model_trends_diffs <- 
    plyr::ddply(subset(model_trends_all, 
                 indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)), 
          ~ indic_group + indic_order + indic,  
          function(df) { 
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(pnull_inf_obs = mean(x1 < x2), 
                       diff_obs_m_null = mean(x2 - x1))
          })
  
  # We compute the differences, but this time we do not average. We obtain a distribution
  # of BOOTN differences between obs and null. This allows keeping the whole distribution, 
  # which quantiles are then represented in the graph by stat_pointinterval(). We could
  # manually compute quantiles and represent them in the graph, but stat_pointinterval()
  # does it for us so here we go.
  model_trends_diffs_distribs <- 
    plyr::ddply(subset(model_trends_all, 
                 indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)), 
          ~ indic_group + indic_order + indic, 
          function(df) { 
            x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
            x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
            data.frame(difference = x2 - x1)
          })
  
  # We create a data.frame which contains all the distributions of differences, with 
  # added columns with summary stats (pnull_inf_obs, diff_obs_m_null), so that we can 
  # plot the distributions, and color them by significance.
  model_trends_diffs_distribs <- plyr::join(model_trends_diffs_distribs, model_trends_diffs, 
                                      type = "left", match = "first")
  
  
  #---------------------------------------------------------------------------
  # SAVE outputs
  #---------------------------------------------------------------------------
  
  filename = paste0("trends_scanlonmodel_fac_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
  save(model_trends, model_trends_med, model_trends_diffs_distribs, model_trends_diffs, 
       file = file.path(path_output,filename))
  
  
  
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  # Model no fac
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  
  model_indics_fmt <- subset(model_indics_fmt_tot,simu_type=="scanlon_nospace")
  
  #---------------------------------------------------------------------------
  # TRENDS
  #---------------------------------------------------------------------------
  
  model_indics_fmt[ ,"indic"] <- as.factor(model_indics_fmt[ ,"indic"])
  model_indics_fmt[ ,"indic_value_type"] <- as.factor(model_indics_fmt[ ,"indic_value_type"])
  
  model_indics_fmt$aridity <- 1-model_indics_fmt$ft_star
  
  model_trends <- plyr::ddply(model_indics_fmt, ~ indic + indic_value_type, 
                              function(df) { 
                                df <- subset(df, is.finite(indic_value)) 
                                if ( nrow(df) < 3 ) { # not enough data points
                                  return( data.frame(estimate = NA) )
                                }
                                get_slope_along_aridity_model(df)
                              }, .progress = "time")
  
  
  # Order indicators for plotting
  model_trends[ ,"indic_order"] <- factor(model_trends[ ,"indic"], 
                                          levels = rev(indic_order), 
                                          ordered = TRUE)
  
  # NOTE: the as.character() here so that indic_groups is indexed by name and not by 
  # the conversion of a factor to integer... (wtf R).
  model_trends[ ,"indic_group"] <- 
    factor(indic_groups[as.character(model_trends[ ,"indic"])], 
           levels = c("Patch-based","Hyd.","Generic")) #, "Composite"
  
  # Compute median slope
  # Here trends contain BOOTN estimates of slopes, effectively a distribution. We 
  # compute here the median of that distribution, its proportion below zero (0 means 
  # all values below zero = significant, 1 means all values above zero = significant too), 
  # and the number of non-NA values (n_estimate, which should be close to BOOTN). 
  # n_estimate is always equal to bootn here because the fit always succeeds on model 
  # data. 
  model_trends_med <- plyr::ddply(model_trends, 
                                  ~ indic_value_type + indic + indic_group + indic_order,
                                  plyr::summarise, 
                                  med_estimate = median(estimate, na.rm = TRUE), 
                                  pval_estimate = mean(estimate < 0, na.rm = TRUE), 
                                  n_estimate = sum(! is.na(estimate)))
  
  if ( any(model_trends_med[ ,"n_estimate"] < (BOOTN * 3/4)) ) { 
    warning("The number of estimate values in the bootstrap distribution is quite low")
  }
  
  # We add the trend information to the original data.frame
  model_trends_all <- plyr::join(model_trends, model_trends_med, type = "left")
  
  # We compute the differences between the trend in observed landscapes and null
  # landscapes. We are effectively computing differences between two distributions, 
  # of exactly BOOTN values each. 
  # We compute the average difference (diff_obs_m_null), and the probability that a 
  # randomly picked value of one sample is below the other (pnull_inf_obs). The latter
  # represents a probability between zero and one, zero when all null slopes are below 
  # their observed counterpart, one when they are all above. 
  model_trends_diffs <- 
    plyr::ddply(subset(model_trends_all, 
                       indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)), 
                ~ indic_group + indic_order + indic,  
                function(df) { 
                  x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
                  x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
                  data.frame(pnull_inf_obs = mean(x1 < x2), 
                             diff_obs_m_null = mean(x2 - x1))
                })
  
  # We compute the differences, but this time we do not average. We obtain a distribution
  # of BOOTN differences between obs and null. This allows keeping the whole distribution, 
  # which quantiles are then represented in the graph by stat_pointinterval(). We could
  # manually compute quantiles and represent them in the graph, but stat_pointinterval()
  # does it for us so here we go.
  model_trends_diffs_distribs <- 
    plyr::ddply(subset(model_trends_all, 
                       indic_value_type %in% c(SLOPE_OBS_TYPE, SLOPE_NULL_TYPE)), 
                ~ indic_group + indic_order + indic, 
                function(df) { 
                  x1 <- subset(df, indic_value_type == SLOPE_NULL_TYPE)[ ,"estimate"]
                  x2 <- subset(df, indic_value_type == SLOPE_OBS_TYPE)[ ,"estimate"]
                  data.frame(difference = x2 - x1)
                })
  
  # We create a data.frame which contains all the distributions of differences, with 
  # added columns with summary stats (pnull_inf_obs, diff_obs_m_null), so that we can 
  # plot the distributions, and color them by significance.
  model_trends_diffs_distribs <- plyr::join(model_trends_diffs_distribs, model_trends_diffs, 
                                            type = "left", match = "first")
  
  
  #---------------------------------------------------------------------------
  # SAVE outputs
  #---------------------------------------------------------------------------
  
  filename = paste0("trends_scanlonmodel_nofac_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
  save(model_trends, model_trends_med, model_trends_diffs_distribs, model_trends_diffs, 
       file = file.path(path_output,filename))
  
  

  
  return(file.path(path_output,filename))
  
} # end function  
  
