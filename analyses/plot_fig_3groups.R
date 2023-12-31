# plots Fig. S16 with 3 groups

library(ggplot2)
library(tidyr)
library(ade4)
library(randomForest)
library(gmodels)
library(patchwork)
library(plyr)
library(gridExtra)

#devtools::install_github('alexgenin/rollply')
library(rollply)


#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
NPERM = 199 
path_output <- here::here("outputs")

#---------------------------------------------------------------------------
# Load data
#---------------------------------------------------------------------------

load(here::here("outputs", "biocom-grps.rda"))
arid <- biocom

filename <- paste0("indics-data50-grps_Nperm_", NPERM,"_rev.rda")
load(file.path(path_output,filename))


#-------------------------------------------------------------------------------
#### Visualize the split of the data in 3 branches
# 3 branches of cover + MF
#-------------------------------------------------------------------------------

rolling_means2 <- rollply(arid, ~ Aridity | pretty_grps3, 
                          wdw.size = 0.05, grid_npts = 128, 
                          summarise, 
                          mean.cover = mean(imgcover), 
                          mean.mf = mean(MF), 
                          n = length(imgcover)
)

# Keep only areas with a lot of points
rolling_means2 <- subset(rolling_means2, n > 20)

interesting_thresholds <-   
  # End of the two upper branches = 0.75
  # End of the lower branch on low aridity = 0.58
  # End of the lower branch on high aridity (end of the data really) ~= 0.92
  geom_vline(aes(xintercept = x), 
             linetype = "dashed", 
             data = data.frame(x = c(0.7, 0.54, 0.82)))

# panel A
fig.branches3.cover = ggplot(NULL, aes(x = Aridity, y = imgcover, color = pretty_grps3)) + 
  geom_point(data = arid,alpha=.6,size=1) + 
  scale_colour_manual(values=c("#66C2A5", "#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "cover")+
  geom_point(aes(x = Aridity, y = mean.cover, size = n), size = 1.5, alpha=0.8, data = rolling_means2)+
  interesting_thresholds


# panel B
fig.gaussian.fit3.cover = ggplot(arid, aes(x = imgcover)) + 
  geom_density(aes(fill = pretty_grps3), 
               bw = .04,alpha= .7) +
  geom_point(aes(color = pretty_grps3), 
             y = rnorm(nrow(arid), 0, 0.02), alpha = .4) +
  scale_colour_manual(values=c("#66C2A5", "#8DA0CB","#E1BE6A"))+
  scale_fill_manual(values=c("#66C2A5", "#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "cover", 
       y = "density")+
  theme(legend.position="none")

# panel C
fig.boxplot3.cover = ggplot(arid, aes(x = pretty_grps3, y = imgcover, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5", "#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "cover (***)")+
  theme(legend.position="none")


# with random effect
library(lme4)

data <- subset(indics, indics$indic=="cutoff")
data$grps3<-as.factor(data$grps3)
#data$value<-log10(data$value)
mod <- lmer(data = data, value ~ grps3 + ( 1 |plotn)) 
#summary(mod)
null <- lmer(data = data, value ~ ( 1 |plotn)) 
anovtable <- as.data.frame(anova(mod,null))
comp <- anovtable$`Pr(>Chisq)`[2]
#value = fixef(mod)["Aridity"]
ifelse(comp < 0.05, "p < 0.05","not significant")
comp

data <- subset(indics, indics$indic=="flowlength")
data$grps3<-as.factor(data$grps3)
mod <- lm(data = data, value ~ grps3) 
summary(mod)


# panel D
fig.branches3.mf = ggplot(NULL, aes(x = Aridity, y = MF, color = pretty_grps3)) + 
  geom_point(data = arid,alpha=.6,size=1) + 
  scale_color_manual(values=c("#66C2A5", "#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "MF")+
  geom_point(aes(x = Aridity, y = mean.mf, size = n), size = 1.5, alpha=0.8, data = rolling_means2)+
  interesting_thresholds



# panel E
fig.gaussian.fit3.mf = ggplot(arid, aes(x = MF)) + 
  geom_density(aes(fill = pretty_grps3), 
               bw = .14,alpha= .7) + 
  geom_point(aes(color = pretty_grps3), 
             y = rnorm(nrow(arid), 0, 0.02), alpha = .4) + 
  scale_color_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "MF", 
       y = "density")+
  theme(legend.position="none")


# panel F
fig.boxplot3.mf = ggplot(arid, aes(x = pretty_grps3, y = MF, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  #geom_point(aes(color = classif_mean))+
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #annotate(geom="text", x=0.8, y=0.9, label="data colored by type cover", color="black")+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "MF (***)")+
  theme(legend.position="none")


# without random effect
#data <- subset(indics, indics$indic=="MF")
#data$value<-log10(data$value)
arid$grps3<-as.factor(arid$grps3)
mod <- lm(data = arid, MF ~ grps3) 
summary(mod)

#anova(lm(data$value ~ data$grps2))

pairwise.t.test(arid$MF, arid$pretty_grps3, p.adj = "bonferroni")
pairwise.t.test(arid$imgcover, arid$pretty_grps3, p.adj = "bonferroni")


#-------------------------------------------------------------------------------
## With all envi var 


fig.boxplot3.arid = ggplot(arid, aes(x = pretty_grps3, y = Aridity, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  #geom_point(aes(color = classif_mean))+
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #annotate(geom="text", x=0.8, y=0.9, label="data colored by type cover", color="black")+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Aridity (***)")+
  theme(legend.position="none")

fig.boxplot3.prod = ggplot(arid, aes(x = pretty_grps3, y = prod, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  #geom_point(aes(color = classif_mean))+
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #annotate(geom="text", x=0.8, y=0.9, label="data colored by type cover", color="black")+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Productivity (***)")+
  theme(legend.position="none")

fig.boxplot3.sand = ggplot(arid, aes(x = pretty_grps3, y = Sand, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  #geom_point(aes(color = classif_mean))+
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #annotate(geom="text", x=0.8, y=0.9, label="data colored by type cover", color="black")+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Sand (***)")+
  theme(legend.position="none")

fig.boxplot3.sr = ggplot(arid, aes(x = pretty_grps3, y = sr, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  #geom_point(aes(color = classif_mean))+
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #annotate(geom="text", x=0.8, y=0.9, label="data colored by type cover", color="black")+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Sp. richness (**)")+
  theme(legend.position="none")


top_row <- plot_grid(fig.branches3.cover, fig.branches3.mf, fig.gaussian.fit3.cover, fig.gaussian.fit3.mf, labels= c("A","B","C","D"),ncol=2, nrow=2)

bottom_row2bis <- plot_grid(fig.boxplot3.cover, fig.boxplot3.mf, fig.boxplot3.prod,fig.boxplot3.arid, fig.boxplot3.sand, fig.boxplot3.sr, labels = c('E', '', '','', '',''), ncol=6,nrow=1)


#-----------------------------------------------------------------------------
## For revision - without log for maxpatch
#-------------------------------------------------------------------------------

indics_c <- indics
unique(indics_c$indic)
indics_c <- subset(indics_c, indic != "cover")
indics_c <- subset(indics_c, indic != "logfmaxpatch")
indics_c <- subset(indics_c, indic != "moran")

# transform indics_c in long format
indics_c$indic <- as.factor(indics_c$indic)
indics_sub <- indics_c[,c("indic","value","plotid","pretty_grps3")]
indics_wide <- spread(indics_sub, indic, value)
indics_wide$pretty_grps3 <- as.factor(indics_wide$pretty_grps3)

#fmaxpatch
fig.boxplot.fmaxpatch = ggplot(indics_wide, aes(x = pretty_grps3, y = fmaxpatch, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+  
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "fmaxpatch (***; log scale)")+
  scale_y_continuous(trans='log10')+
  theme(legend.position="none")

#slope
fig.boxplot.slope = ggplot(indics_wide, aes(x = pretty_grps3, y = slope, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "slope (***)")+
  theme(legend.position="none")

#cutoff
fig.boxplot.cutoff = ggplot(indics_wide, aes(x = pretty_grps3, y = cutoff, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "cutoff (***)")+
  theme(legend.position="none")

#cv
fig.boxplot.cv = ggplot(indics_wide, aes(x = pretty_grps3, y = cv.variance, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "cv (***)")+
  theme(legend.position="none")

#sdr
fig.boxplot.sdr = ggplot(indics_wide, aes(x = pretty_grps3, y = sdr, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "sdr (NS)")+
  theme(legend.position="none")

#fl
fig.boxplot.fl = ggplot(indics_wide, aes(x = pretty_grps3, y = flowlength, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "flowlength (***)")+
  theme(legend.position="none")

fig.leg = ggplot(indics_wide, aes(x = pretty_grps3, y = flowlength, fill=pretty_grps3)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "flowlength (***)")+
  theme(legend.position="bottom",legend.title=element_blank())

legend <- get_legend(fig.leg)

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

third_row <- plot_grid(fig.boxplot.fmaxpatch, fig.boxplot.slope, fig.boxplot.cutoff, labels = c('', '', ''), ncol=3, nrow=1) 

fourth_row <- plot_grid(fig.boxplot.cv, fig.boxplot.sdr,fig.boxplot.fl, labels = c('', '', ''), ncol=3,nrow=1)  

## Fig S16 no log
# fig3_spmetrics_boxplots_nolog.pdf
# 5 x 9.5
#plot_grid(top_row, bottom_row2, ncol = 1)
#grid.arrange(top_row, bottom_row2, legend, ncol = 1,nrow=3, heights = c(1,1,0.2))

plot_grid(top_row, bottom_row2bis,third_row,fourth_row,ncol = 1)

lay <- rbind(c(1),
             c(1),
             c(2),
             c(3),
             c(4))

grid.arrange(top_row, bottom_row2bis,third_row,fourth_row, layout_matrix = lay)

lay <- rbind(c(1),
             c(1),
             c(2),
             c(3),
             c(4),
             c(5))

grid.arrange(top_row, bottom_row2bis,third_row,fourth_row,legend,layout_matrix = lay)

grid.arrange(top_row, bottom_row2bis,third_row,fourth_row,legend,nrow=5,ncol=1,heights = c(2,1,1,1,0.1))

#figS16_3groups_all.pdf
#9.5x9.5


