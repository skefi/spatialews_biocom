# Plot slopes of the indicators in the model and in the data (for the whole data set or for groups) ; 
# Fig 5, S7, S8, S9, S15, S19

library(ggplot2)
require(tidybayes)
library(plyr) # for join
library(cowplot) # for plot_grid

#source(here::here("R", "functions_helper.R"))

#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
path_output <- here::here("outputs")
NPERM = 199 #value used for final analyses: 199 
BOOTN = 2999 #value used for final analyses: 2999

ALPHA <- 0.05 # Significance level for indicator trends

SLOPE_OBS_TYPE <- "value_scaled_orig"
SLOPE_NULL_TYPE <- "null_mean_scaled_orig"

indic_order <- c('logfmaxpatch', 'slope','cutoff', 'flowlength', 'sdr','cv.variance') 


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
# Load data 50*50 images
#---------------------------------------------------------------------------

filename = paste0("trends_data50_one_group_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename))

filename2 = paste0("trends_data50_two_groups_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename2))

filename3 = paste0("trends_data50_three_groups_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename3))

filename_m = paste0("trends_kefimodel_fac_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename_m))




#---------------------------------------------------------------------------
# Figure 5 main text
#---------------------------------------------------------------------------

trends_all <- subset(trends_all, indic != "fmaxpatch")
trends_all <- subset(trends_all, indic != "moran")

unique(trends_all$indic)

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
  xlim(-6,9)+
  sign_color_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  labs(x = "Estimated trend", 
       y = "Sp. metric",  
       title = "Data (all)") +
  theme(panel.grid.minor.x = element_blank())+
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


model_trends_all <- join(model_trends, model_trends_med, type = "left")
model_trends_all <- subset(model_trends_all, indic != "fmaxpatch")
model_trends_all <- subset(model_trends_all, indic != "moran")
unique(model_trends_all$indic)

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
  xlim(-6,9)+
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  sign_color_scale + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "Model 1")+ 
  theme(panel.grid.minor.x = element_blank())+
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))

trends2_all <- subset(trends2_all, indic != "fmaxpatch")
trends2_all <- subset(trends2_all, indic != "moran")

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
  xlim(-6,9)+
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
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


#### Fig. 5
top_row <- plot_grid(plot_model_trends, plot_trends, labels= c("A","B"),ncol=2, nrow=1)
bottom_row <- plot_grid(plot_trends2, labels = c('C'), ncol=1,nrow=1)
plot_grid(top_row, bottom_row, ncol = 1)
ggsave("./figures/fig5_slope_mod1_data_2groups.pdf", width = 9.5, height = 8)


#---------------------------------------------------------------------------
# Figure S7 trends diff - Plot with mean differences in slopes obs vs. null
# model vs all data
#---------------------------------------------------------------------------

#unique(trends_diffs$indic)
trends_diffs <- subset(trends_diffs, indic != "fmaxpatch")
trends_diffs <- subset(trends_diffs, indic != "moran")

plot_trends_diffs <- 
  ggplot(trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  xlim(-3,8)+
  sign_fill_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="Data (all)")+
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


#unique(model_trends_diffs$indic)
model_trends_diffs <- subset(model_trends_diffs, indic != "fmaxpatch")
model_trends_diffs <- subset(model_trends_diffs, indic != "moran")

# Plot with differences in trends vs. null
plot_model_trends_diffs <- 
  ggplot(model_trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  xlim(-3,8)+
  sign_fill_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="Model 1") +
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


trends2_diffs <- subset(trends2_diffs, indic != "fmaxpatch")
trends2_diffs <- subset(trends2_diffs, indic != "moran")

plot_trends2_diffs <- 
  ggplot(trends2_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  xlim(-2,5)+
  sign_fill_scale + 
  facet_grid(indic_group ~ pretty_grps2, space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="Data (per group)")


#### Figure S7
top_row <- plot_grid(plot_model_trends, plot_trends, labels= c("A","B"),ncol=2, nrow=1)
bottom_row <- plot_grid(plot_model_trends_diffs, plot_trends_diffs, labels= c("C","D"),ncol=2, nrow=1)
plot_grid(top_row, bottom_row, ncol = 1)
ggsave("./figures/figS7_slopes_mod_data_points_diff.pdf", width = 9.5, height = 8)


#---------------------------------------------------------------------------
# Figure S15 trends diff - 2 groups
#---------------------------------------------------------------------------

#### Figure S15
plot_grid(plot_trends2, plot_trends2_diffs, labels= c("A","B"), ncol = 1)
ggsave("./figures/figS15_slopes_2groups_points_diff.pdf", width = 9.5, height = 8)


#---------------------------------------------------------------------------
# Figure S19 trends diff - 3 groups
#---------------------------------------------------------------------------

trends3_all <- subset(trends3_all, indic != "fmaxpatch")
trends3_all <- subset(trends3_all, indic != "moran")

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
  xlim(-6,9)+
  sign_color_scale + 
  facet_grid(indic_group ~ pretty_grps3, space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "Data (per group)")+ 
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank()) 


trends3_diffs <- subset(trends3_diffs, indic != "fmaxpatch")
trends3_diffs <- subset(trends3_diffs, indic != "moran")

plot_trends3_diffs <- ggplot(trends3_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  xlim(-2,5)+
  sign_fill_scale + 
  facet_grid(indic_group ~ pretty_grps3, space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "", x = "Est. diff",title="Data (per group)")


#### Fig S19
plot_grid(plot_trends3, plot_trends3_diffs, labels= c("A","B"), ncol = 1)
ggsave("./figures/figS19_slopes_3groups_points_diffs.pdf", width = 9.5, height = 8)



#---------------------------------------------------------------------------
# Figures S8 and S9
# Figure TPB model no fac
#---------------------------------------------------------------------------

## No fac, glob comp

filename_tg = paste0("trends_kefimodel_nofacglob_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename_tg))

#unique(model_trends$indic)
model_trends_all <- join(model_trends, model_trends_med, type = "left")

model_trends_all <- subset(model_trends_all, indic != "moran")
model_trends_all <- subset(model_trends_all, indic != "fmaxpatch")

plot_tpbmodel_glob_trends <- 
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
  xlim(-6,9)+
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  sign_color_scale + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "no fac, glob")+ 
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))



#unique(model_trends_diffs$indic)
model_trends_diffs <- subset(model_trends_diffs, indic != "moran")
model_trends_diffs <- subset(model_trends_diffs, indic != "fmaxpatch")


# Plot with differences in trends vs. null
plot_tpbmodel_glob_trends_diffs <- 
  ggplot(model_trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  xlim(-3,8)+
  sign_fill_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="no fac, glob") +
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


## No fac local comp

filename_tl = paste0("trends_kefimodel_nofacloc_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename_tl))

model_trends_all <- join(model_trends, model_trends_med, type = "left")

model_trends_all <- subset(model_trends_all, indic != "moran")
model_trends_all <- subset(model_trends_all, indic != "fmaxpatch")

plot_tpbmodel_loc_trends <- 
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
  xlim(-6,9)+
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  sign_color_scale + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "no fac, loc")+ 
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))



#unique(model_trends_diffs$indic)
model_trends_diffs <- subset(model_trends_diffs, indic != "moran")
model_trends_diffs <- subset(model_trends_diffs, indic != "fmaxpatch")

# Plot with differences in trends vs. null
plot_tpbmodel_loc_trends_diffs <- 
  ggplot(model_trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  xlim(-3,8)+
  sign_fill_scale + 
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="no fac, loc") +
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))








#---------------------------------------------------------------------------
# Figure Scanlon 
#---------------------------------------------------------------------------

filename_sf = paste0("trends_scanlonmodel_fac_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename_sf))

#unique(model_trends$indic)
model_trends_all <- join(model_trends, model_trends_med, type = "left")

unique(model_trends_all$indic)
model_trends_all <- subset(model_trends_all, indic != "moran")
model_trends_all <- subset(model_trends_all, indic != "fmaxpatch")


plot_scanlonmodel_fac_trends <- 
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
  xlim(-6,9)+
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  sign_color_scale + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "Model 2")+ 
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))

#unique(model_trends_diffs$indic)
model_trends_diffs <- subset(model_trends_diffs, indic != "moran")
model_trends_diffs <- subset(model_trends_diffs, indic != "fmaxpatch")

# Plot with differences in trends vs. null
plot_scanlonmodel_fac_trends_diffs <- 
  ggplot(model_trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  sign_fill_scale + 
  xlim(-3,8)+
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="Model 2") +
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))


## No fac

filename_sf = paste0("trends_scanlonmodel_nofac_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename_sf))

model_trends_all <- join(model_trends, model_trends_med, type = "left")

model_trends_all <- subset(model_trends_all, indic != "moran")
model_trends_all <- subset(model_trends_all, indic != "fmaxpatch")

plot_scanlonmodel_nofac_trends <- 
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
  xlim(-6,9)+
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  sign_color_scale + 
  labs(x = "Estimated trend", 
       y = "Sp. metric", 
       title = "no fac")+ 
  theme(text = element_text(size=12))+
  theme(panel.grid.minor.x = element_blank())+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))

#unique(model_trends_diffs$indic)
model_trends_diffs <- subset(model_trends_diffs, indic != "moran")
model_trends_diffs <- subset(model_trends_diffs, indic != "fmaxpatch")

# Plot with differences in trends vs. null
plot_scanlonmodel_nofac_trends_diffs <- 
  ggplot(model_trends_diffs, aes(x = diff_obs_m_null, y = indic_order)) + 
  # We want a P-value close to zero when the obs is above the null, so we take 
  # 1 - the P-value to plug into the color scale
  geom_col(aes(fill = signf(1-pnull_inf_obs, ALPHA))) + 
  sign_fill_scale + 
  xlim(-3,8)+
  facet_grid(indic_group ~ ., space = "free_y", 
             scale = "free_y", switch = "y") + 
  theme_minimal() + 
  theme(panel.grid.minor.x = element_blank()) + 
  labs(y = "Sp. metric", x = "Est. diff",title="no fac") +
  theme(text = element_text(size=12))+
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 12))



#---------------------------------------------------------------------------
# Figure 
#---------------------------------------------------------------------------

intermediate_row <- plot_grid(plot_model_trends, plot_tpbmodel_glob_trends, plot_tpbmodel_loc_trends, labels= c("A","B","C"),ncol=3, nrow=1)

bottom_row <- plot_grid(plot_scanlonmodel_fac_trends, plot_scanlonmodel_nofac_trends, labels = c("D","E"), ncol=3,nrow=1)


#### Fig. S8
#plot_grid(top_row, intermediate_row, bottom_row, ncol = 1)
plot_grid(intermediate_row, bottom_row, ncol = 1)
ggsave("./figures/figS8_slopes_mod1_mod2_points.pdf", width = 9.5, height = 8)


#top_row <- plot_grid(plot_trends_diffs, labels= c("A"),ncol=3, nrow=1)
intermediate_row <- plot_grid(plot_model_trends_diffs, plot_tpbmodel_glob_trends_diffs, plot_tpbmodel_loc_trends_diffs, labels= c("A","B","C"),ncol=3, nrow=1)
bottom_row <- plot_grid(plot_scanlonmodel_fac_trends_diffs, plot_scanlonmodel_nofac_trends_diffs, labels = c("D","E"), ncol=3,nrow=1)


#### Fig. S9
#plot_grid(top_row, intermediate_row, bottom_row, ncol = 1)
plot_grid(intermediate_row, bottom_row, ncol = 1)
ggsave("./figures/figS9_slopes_data_mod1_mod2_diffs.pdf", width = 9.5, height = 8)







