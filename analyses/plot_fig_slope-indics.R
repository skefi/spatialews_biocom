# Plot slopes of the indicators in the model and in the data (for the whole data set or for groups) ; Fig 2, S3, S10, S14

library(ggplot2)
require(tidybayes)

source(here::here("R", "functions_helper.R"))

#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
# From _targets.R
path_output <- here::here("outputs")
NPERM = 199 #value used for final analyses: 199
BOOTN = 2999 #value used for final analyses: 2999

ALPHA <- 0.05 # Significance level for indicator trends

SLOPE_OBS_TYPE <- "value_scaled_orig"
SLOPE_NULL_TYPE <- "null_mean_scaled_orig"

indic_order <- c('plrange', 'logfmaxpatch', 'slope','cutoff', 'flowlength','moran', 'sdr','cv.variance') 

# Color scale for significance in plots
sign_color_scale <- scale_color_manual(values = c(neg = "#ca0020", 
                                                  pos = "#0571b0", 
                                                  ns  = "#7e7e7e"), 
                                       guide = "none")
sign_fill_scale <- scale_fill_manual(values = c(neg = "#ca0020", 
                                                pos = "#0571b0", 
                                                ns  = "#7e7e7e"), 
                                     guide = "none")


#---------------------------------------------------------------------------
# Load data
#---------------------------------------------------------------------------

filename = paste0("trends_one_group_Nperm_",NPERM,"_Bootn_",BOOTN,".rda")
load(file.path(path_output,filename))

filename2 = paste0("trends_two_groups_Nperm_",NPERM,"_Bootn_",BOOTN,".rda")
load(file.path(path_output,filename2))

filename3 = paste0("trends_three_groups_Nperm_",NPERM,"_Bootn_",BOOTN,".rda")
load(file.path(path_output,filename3))

filename_m = paste0("trends_model_Nperm_",NPERM,"_Bootn_",BOOTN,".rda")
load(file.path(path_output,filename_m))


#---------------------------------------------------------------------------
# Figure 2 main text
#---------------------------------------------------------------------------

plot_trends  <- 
  ggplot(NULL, aes(x = estimate, y = indic_order)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  stat_pointinterval(aes(x = estimate, y = indic_order), 
                     color = '#AAAAAA', point_size = 1, 
                     data = subset(trends_all, indic_value_type == SLOPE_NULL_TYPE)) + 
  stat_pointinterval(aes(x = estimate, y = indic_order, 
                         color = signf(pval_estimate, ALPHA)), 
                     point_size = 2, 
                     data = subset(trends_all, indic_value_type == SLOPE_OBS_TYPE)) + 
  sign_color_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  labs(x = "Estimated trend", 
       y = "Sp. metric",  
       title = "Data (all)") +
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))

model_trends_all <- join(model_trends, model_trends_med, type = "left")

plot_model_trends <- 
  ggplot(NULL, aes(x = estimate, y = indic_order)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  stat_pointinterval(aes(x = estimate, y = indic_order), 
                     color = '#AAAAAA', point_size = 1, 
                     data = subset(model_trends_all,
                                   indic_value_type == SLOPE_NULL_TYPE)) + 
  stat_pointinterval(aes(x = estimate, y = indic_order, 
                         color = signf(pval_estimate, ALPHA)), 
                     point_size = 2, 
                     data = subset(model_trends_all, 
                                   indic_value_type == SLOPE_OBS_TYPE)) +
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  sign_color_scale + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "Model")+ 
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))

plot_trends2 <- 
  ggplot(NULL, aes(x = estimate, y = indic_order)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  stat_pointinterval(aes(x = estimate, y = indic_order), 
                     color = '#AAAAAA', point_size = 1, 
                     data = subset(trends2_all, indic_value_type == SLOPE_NULL_TYPE)) + 
  stat_pointinterval(aes(x = estimate, y = indic_order, 
                         color = signf(pval_estimate, ALPHA)), 
                     point_size = 2, 
                     data = subset(trends2_all, indic_value_type == SLOPE_OBS_TYPE)) + 
  sign_color_scale + 
  facet_grid(indic_group ~ pretty_grps2, space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  scale_fill_manual(values = c("#c2c4d4", "#d40000"), 
                    palette = "RdBu", 
                    na.value = "#7F7F7F", 
                    guide = "none") + 
  labs(x = "Estimated trend", 
       y = "Sp. metric",
       title = "Data (per group)")+ 
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


# 8 x 9.5
# from 05 and 06
top_row <- plot_grid(plot_model_trends, plot_trends, labels= c("A","B"),ncol=2, nrow=1)
bottom_row <- plot_grid(plot_trends2, labels = c('C'), ncol=1,nrow=1)
plot_grid(top_row, bottom_row, ncol = 1)


#---------------------------------------------------------------------------
# Figure S3 trends diff - Plot with mean differences in slopes obs vs. null
# model vs all data
#---------------------------------------------------------------------------

plot_trends_diffs <- 
  ggplot(trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  sign_fill_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="Data (all)")+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


# Plot with differences in trends vs. null
plot_model_trends_diffs <- 
  ggplot(model_trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  sign_fill_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="Model") +
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


plot_trends2_diffs <- 
  ggplot(trends2_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  sign_fill_scale + 
  facet_grid(indic_group ~ pretty_grps2, space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="Data (per group)")

### Figure S3
top_row <- plot_grid(plot_model_trends, plot_trends, labels= c("A","B"),ncol=2, nrow=1)
bottom_row <- plot_grid(plot_model_trends_diffs, plot_trends_diffs, labels= c("C","D"),ncol=2, nrow=1)
plot_grid(top_row, bottom_row, ncol = 1)


#---------------------------------------------------------------------------
# Figure S10 trends diff - 2 groups
#---------------------------------------------------------------------------

plot_grid(plot_trends2, plot_trends2_diffs, labels= c("A","B"), ncol = 1)


#---------------------------------------------------------------------------
# Figure S14 trends diff - 3 groups
#---------------------------------------------------------------------------

plot_trends3 <- 
  ggplot(NULL, aes(x = estimate, y = indic_order)) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  stat_pointinterval(aes(x = estimate, y = indic_order), 
                     color = '#AAAAAA', point_size = 1, 
                     data = subset(trends3_all, indic_value_type == SLOPE_NULL_TYPE)) + 
  stat_pointinterval(aes(x = estimate, y = indic_order, 
                         color = signf(pval_estimate, ALPHA)), 
                     point_size = 2, 
                     data = subset(trends3_all, indic_value_type == SLOPE_OBS_TYPE)) + 
  sign_color_scale + 
  facet_grid(indic_group ~ pretty_grps3, space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "Data (per group)")+ 
  theme(panel.grid.minor.x = element_blank()) 

plot_trends3_diffs <- ggplot(trends3_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  sign_fill_scale + 
  facet_grid(indic_group ~ pretty_grps3, space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "", x = "Est. diff",title="Data (per group)")


#9.5 x 9.5
plot_grid(plot_trends3, plot_trends3_diffs, labels= c("A","B"), ncol = 1)
