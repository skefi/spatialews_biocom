# Plot Figs. S20 and S21

library(ggplot2)
require(tidybayes)
library(plyr) # for join
library(cowplot) # for plot_grid


source(here::here("R", "functions_helper.R"))

#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
# From _targets.R
path_output <- here::here("outputs")
NPERM = 199 #value used for original submission: 199, for revision: 200
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
# Load data 200*200 images
#---------------------------------------------------------------------------

filename = paste0("trends_data200_one_group_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename))

filename2 = paste0("trends_data200_two_groups_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename2))

filename3 = paste0("trends_data200_three_groups_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename3))

filename_m = paste0("trends_kefimodel_fac_Nperm_",NPERM,"_Bootn_",BOOTN,"_rev.rda")
load(file.path(path_output,filename_m))


load(file.path(path_output,"data_biocom_200x200.rda"))

ourdata_200x200[["biocom"]][["matrixn"]] <- seq_len(nrow(ourdata_200x200[["biocom"]]))

load(file.path(path_output,"biocom-grps.rda"))

biocom200 <- ourdata_200x200$biocom

extr <- biocom[ ,c("plotid", "pretty_grps2", "pretty_grps2c", "pretty_grps3", "grps2c","grps2","grps3")]

biocom200_grps <- plyr::join(biocom200, extr, by = "plotid") #51 sites 


grp1 <- biocom200_grps[biocom200_grps$grps2==1,] #low cov, low MF, 32 sites
grp2 <- biocom200_grps[biocom200_grps$grps2==2,] #high, high, 19 sites


#---------------------------------------------------------------------------
# Figure S20
#---------------------------------------------------------------------------

trends_all <- subset(trends_all, indic != "fmaxpatch")
trends_all <- subset(trends_all, indic != "moran")

#unique(trends_all$indic)

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



#### fig S20
plot_grid(plot_model_trends, plot_trends, labels= c("A","B"),ncol=2, nrow=1)
ggsave("./figures/figS20_slopes_mod_data200_2groups.pdf.pdf", width = 9.5, height = 5)


#---------------------------------------------------------------------------
# Figure S21 trends diff - Plot with mean differences in slopes obs vs. null
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



### Figure S21
top_row <- plot_grid(plot_model_trends, plot_trends, labels= c("A","B"),ncol=2, nrow=1)
bottom_row <- plot_grid(plot_model_trends_diffs, plot_trends_diffs, labels= c("C","D"),ncol=2, nrow=1)
plot_grid(top_row, bottom_row, ncol = 1)
ggsave("./figures/figS21_slopes_mod_data200_points_diff.pdf", width = 9.5, height = 8)


