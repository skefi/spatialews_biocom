# Plot indicators along gradients
# Figs. S6, S14, S18

library(ggplot2)
library(tidyr) # for gather
library(gridExtra) # for grid.arrange
library(cowplot) # for plot_grid

#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
NPERM = 199 #value use for final analyses: 199
BOOTN = 2999 #value use for final analyses: 2999
path_output <- here::here("outputs")


ALPHA <- 0.05 # Significance level for indicator trends

indic_order <- c('logfmaxpatch', 'slope','cutoff', 'flowlength', 'sdr','cv.variance') 

#---------------------------------------------------------------------------
# Load data
#---------------------------------------------------------------------------

load(here::here("outputs", "biocom-grps.rda"))
load(file.path(path_output,"data_biocom.rda"))
arid <- biocom

filename <- paste0("indics-data50-grps_Nperm_", NPERM,"_rev.rda")
load(file.path(path_output,filename))

filename_m = paste0("indics-kefimodel_all_Nperm_",NPERM,"_rev.rda")
load(file.path(path_output,filename_m))


#---------------------------------------------------------------------------
# Fig S2: Trends in model
#---------------------------------------------------------------------------

# Format indicator results
model_indics_fmt <- tidyr::gather(model_indics_df, "indic_value_type", "indic_value", 
                                  all_of(c("z_score", "value", "null_mean", 
                                           "value_scaled", "null_scaled", 
                                           "value_scaled_orig", 
                                           "null_mean_scaled_orig"))) 


unique(model_indics_fmt$indic_value_type)
unique(model_indics_fmt$indic)

# Remove indicators that are not needed 
model_indics_fmt <- subset(model_indics_fmt, indic != "cover")
model_indics_fmt <- subset(model_indics_fmt, indic != "fmaxpatch")
model_indics_fmt <- subset(model_indics_fmt, indic != "moran")
model_indics_fmt$aridity <- 1-model_indics_fmt$b
model_indics_fmt <- subset(model_indics_fmt, simu_type == "original")

# We only work on aridity values above 0.3 
model_indics_fmt <- subset(model_indics_fmt, aridity > 0.3)

# change order
model_indics_fmt[ ,"indic_order"] <- 
  factor(model_indics_fmt[ ,"indic"], 
         levels = rev(c('flowlength','cutoff','cv.variance','slope','sdr','logfmaxpatch')), 
         ordered = TRUE)


facet_names <- list(
  'logfmaxpatch'="log10(fmaxpatch)",
  'sdr'="sdr",
  'slope'="slope",
  'cv.variance'="cv",
  'cutoff'="cutoff",
  'flowlength'="flowlength"
)

facet_labeller <- function(variable,value){
  return(facet_names[value])
}


#9.5 x 5.8
fig_mod <- ggplot(subset(model_indics_fmt, indic_value_type %in% c("value", "null_mean")), 
                  aes(x = aridity, y = indic_value)) + 
  geom_point(aes(color = indic_value_type), 
             alpha = .8) + #0.8 
  #geom_smooth(aes(x=aridity,y=indic_value,color=indic_value_type),method=lm,formula=y~x)+ #
  #scale_color_discrete(name = "Displayed value") + 
  #facet_wrap( ~ indic_order, labeller = facet_labeller, 
  facet_wrap( ~ indic_order, labeller = facet_labeller, 
              ncol=2, scales = "free_y") + 
  scale_color_manual(values = c("value"="black","null_mean"="gray"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "Aridity", 
       y = "Sp. metric")+
  theme(legend.position="none")




#---------------------------------------------------------------------------
# Fig S2: Trends in data - all 
#---------------------------------------------------------------------------

# Format indicator results
indics_fmt_dat <- tidyr::gather(indics, "indic_value_type", "indic_value", all_of(c("z_score", "value", "null_mean"))) 

unique(indics_fmt_dat$indic_value_type)
unique(indics_fmt_dat$indic)

# Remove mean which is not an indicator 
indics_fmt_dat <- subset(indics_fmt_dat, indic != "cover")
indics_fmt_dat <- subset(indics_fmt_dat, indic != "fmaxpatch")
indics_fmt_dat <- subset(indics_fmt_dat, indic != "moran")

# change order
indics_fmt_dat[ ,"indic_order"] <- 
  factor(indics_fmt_dat[ ,"indic"], 
         levels = rev(c('flowlength','cutoff','cv.variance','slope','sdr','logfmaxpatch')), 
         ordered = TRUE)


#9.5 x 5.8
fig_dat <- ggplot(subset(indics_fmt_dat, indic_value_type %in% c("value", "null_mean")), aes(x = Aridity, y = indic_value)) + 
  geom_point(aes(color = indic_value_type), alpha = .8) + #0.8
  #geom_smooth(aes(x=Aridity,y=indic_value,color=indic_value_type),method=lm,formula=y~x)+
  #facet_wrap( ~ indic_order, labeller = facet_labeller, ncol=2, scales = "free_y") + 
  facet_wrap( ~ indic_order, labeller = facet_labeller, ncol=2, scales = "free_y") + 
  
  scale_color_manual(values = c("value"="black","null_mean"="gray"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "Aridity", 
       y = "Sp. metric")+
  theme(legend.position="none")


## Fig S6
# figS6_rawtrends_model_data.pdf
#8x9.5
plot_grid(fig_mod, fig_dat, labels = c('A', 'B'), ncol=2,nrow=1)


 #---------------------------------------------------------------------------
# Fig. S15 (S8): Trends in data - 2 groups 
#---------------------------------------------------------------------------

branch1 <- subset(indics_fmt_dat, grps2==1)
branch2 <- subset(indics_fmt_dat, grps2==2)

plot1 <- ggplot(subset(branch1, indic_value_type %in% c("value","null_mean")),
                aes(x = Aridity, y = indic_value)) + 
  geom_point(aes(color = indic_value_type), 
             alpha = .8) + 
  geom_smooth(aes(x=Aridity,y=indic_value,color=indic_value_type),method=lm,formula=y~x)+
  facet_wrap( ~ indic_order, labeller = facet_labeller, ncol=2, scales = "free_y") + 
  scale_color_manual(values = c("value"="black","null_mean"="gray"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "Aridity", 
       y = "Sp. metric",
       title = "Low-low")+
  theme(legend.position="none")

plot2 <- ggplot(subset(branch2, indic_value_type %in% c("value", "null_mean")), 
                aes(x = Aridity, y = indic_value)) + 
  geom_smooth(aes(x=Aridity,y=indic_value,color=indic_value_type),method=lm,formula=y~x)+
  geom_point(aes(color = indic_value_type), 
             alpha = .8) + 
  facet_wrap( ~ indic_order, labeller = facet_labeller, ncol=2, scales = "free_y") + 
  scale_color_manual(values = c("value"="black","null_mean"="gray"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "Aridity", 
       y = "Sp. metric",
       title = "High-high")+
  theme(legend.position="none")

# figS15_rawtrends_data_2groups.pdf
#8 x 9.5
lay <- rbind(c(1,2))
grid.arrange(plot1, plot2, layout_matrix = lay)



##---------------------------------------------------------------------------
### Fig. S19 (S12): plot raw data per 3 branch
##---------------------------------------------------------------------------

branch31 <- subset(indics_fmt_dat, grps3==1)
branch32 <- subset(indics_fmt_dat, grps3==2)
branch33 <- subset(indics_fmt_dat, grps3==3)


plot31 <- ggplot(subset(branch31, indic_value_type %in% c("value","null_mean")),
                 aes(x = Aridity, y = indic_value)) + 
  geom_point(aes(color = indic_value_type), 
             alpha = .8) + 
  geom_smooth(aes(x=Aridity,y=indic_value,color=indic_value_type),method=lm,formula=y~x)+
  facet_wrap( ~ indic_order, labeller = facet_labeller, ncol=2, scales = "free_y") + 
  scale_color_manual(values = c("value"="black","null_mean"="gray"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "Aridity", 
       y = "Sp. metric",
       title = "Low-low")+
  theme(legend.position="none")

plot32 <- ggplot(subset(branch32, indic_value_type %in% c("value", "null_mean")), 
                 aes(x = Aridity, y = indic_value)) + 
  geom_smooth(aes(x=Aridity,y=indic_value,color=indic_value_type),method=lm,formula=y~x)+
  geom_point(aes(color = indic_value_type), 
             alpha = .8) + 
  facet_wrap( ~ indic_order, labeller = facet_labeller, ncol=2, scales = "free_y") + 
  scale_color_manual(values = c("value"="black","null_mean"="gray"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "Aridity", 
       y = "Sp. metric",
       title = "Low-high")+
  theme(legend.position="none")

plot33 <- ggplot(subset(branch33, indic_value_type %in% c("value", "null_mean")), 
                 aes(x = Aridity, y = indic_value)) + 
  geom_smooth(aes(x=Aridity,y=indic_value,color=indic_value_type),method=lm,formula=y~x)+
  geom_point(aes(color = indic_value_type), 
             alpha = .8) + 
  facet_wrap( ~ indic_order, labeller = facet_labeller, ncol=2, scales = "free_y") + 
  scale_color_manual(values = c("value"="black","null_mean"="gray"))+
  theme_minimal() +
  theme(text = element_text(size=12))+
  labs(x = "Aridity", 
       y = "Sp. metric",
       title = "High-high")+
  theme(legend.position="none")


# figS19_rawtrends_data_3groups.pdf
#8 x 11
lay <- rbind(c(1,2,3))
grid.arrange(plot31, plot32, plot33, layout_matrix = lay)



