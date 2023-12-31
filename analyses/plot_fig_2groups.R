# Version for revision PNAS
# Figures 2, 3, S10, S12, S13

#devtools::install_github('alexgenin/rollply')
library(rollply)
library(ggplot2)
library(cowplot)
library(dplyr)
library(plyr)
library(tidyr)
library(gridExtra)


#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
NPERM = 199 #value use for first sub: 199, for revision: 200
path_output <- here::here("outputs")

#---------------------------------------------------------------------------
# Load data
#---------------------------------------------------------------------------

load(here::here("outputs", "biocom-grps.rda"))
arid <- biocom

filename <- paste0("indics-data50-grps_Nperm_", NPERM,"_rev.rda")
load(file.path(path_output,filename))


#---------------------------------------------------------------------------
# Figure 2 : 2 groups
#---------------------------------------------------------------------------

rolling_means2 <- rollply(arid, ~ Aridity | pretty_grps2, 
                          wdw.size = 0.07, grid_npts = 128, 
                          #wdw.size = 0.1, grid_npts = 128, 
                          summarise, 
                          mean.cover = mean(imgcover), 
                          mean.mf = mean(MF), 
                          n = length(imgcover)
)

# Keep only areas with a lot of points
rolling_means2 <- subset(rolling_means2, n > 30)

fig.branches2.cover = ggplot(NULL, aes(x = Aridity, y = imgcover, color = pretty_grps2)) + 
  geom_point(data = arid,alpha=.6,size=1) + 
  #scale_color_manual(values = c("#5ab4ac","#d8b365"))+ # colors from the 1st submission
  scale_color_manual(values = c("#40B0A6","#E1BE6A"))+ 
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "cover")+
  geom_point(aes(x = Aridity, y = mean.cover, size = n), size = 2, alpha=0.9, data = rolling_means2)


fig.branches2.mf = ggplot(NULL, aes(x = Aridity, y = MF, color = pretty_grps2)) + 
  geom_point(data = arid,alpha=.6,size=1) + 
  scale_color_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "MF")+
  geom_point(aes(x = Aridity, y = mean.mf, size = n), size = 2, alpha=0.9, data = rolling_means2)

fig.gaussian.fit2.cover = ggplot(arid, aes(x = imgcover)) + 
  geom_density(aes(fill = pretty_grps2), 
               bw = .04,alpha= .7) +
  geom_point(aes(color = pretty_grps2), 
             y = rnorm(nrow(arid), 0, 0.02), alpha = .4) +
  scale_colour_manual(values=c("#40B0A6","#E1BE6A"))+
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "cover", 
       y = "density")+
  theme(legend.position="none")

fig.gaussian.fit2.mf = ggplot(arid, aes(x = MF)) + 
  geom_density(aes(fill = pretty_grps2), 
               bw = .14,alpha= .7) + 
  geom_point(aes(color = pretty_grps2), 
             y = rnorm(nrow(arid), 0, 0.02), alpha = .4) + 
  scale_color_manual(values=c("#40B0A6","#E1BE6A"))+
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "MF", 
       y = "density")+
  theme(legend.position="none")

# for the legend
fig.leg = ggplot(arid, aes(x = MF)) + 
  geom_density(aes(fill = pretty_grps2), 
               bw = .14,alpha= .7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(legend.position="bottom",legend.title=element_blank())

legend <- get_legend(fig.leg)


# Figure from 1st submission; now top row of S13
#5 x 9.5
toprowS13 <- plot_grid(fig.branches2.cover, fig.branches2.mf, fig.gaussian.fit2.cover, fig.gaussian.fit2.mf, labels= c("A","B","C","D"),ncol=2, nrow=2)



#---------------------------------------------------------------------------
# Figure 2 revised
# use potential analysis 
#---------------------------------------------------------------------------

load(here::here("outputs", "data_biocom.rda"))

matrices <- ourdata[["matrices"]]
biocom <- ourdata[["biocom"]]

# 
# Do the density analysis/detection of maxima/minima for cover
# 
# 
# horizontal resolution 
hres <- 0.01
# bins (actually points along the gradient because bins overlap if the window size is 
# greater than hres/2 below)
bins <- seq(min(biocom[ ,"Aridity"]), 
            max(biocom[ ,"Aridity"]), 
            by = hres)
# window size
window_width <- hres * 15

# y values at which the density will be evaluated
yvalues <- seq(min(biocom[ ,"imgcover"]), 
               max(biocom[ ,"imgcover"]), 
               l = 1024)

# 
window_analysis <- llply(bins, function(bin) { 
  # Compute density of points in window
  ok <- subset(biocom, Aridity > (bin - window_width) & Aridity < (bin + window_width))
  dens <- with(ok, density(imgcover, bw = bw.nrd(imgcover), 
                           from = min(yvalues), to = max(yvalues), n = 2048))
  dens <- approx(dens[["x"]], dens[["y"]], yvalues)[["y"]]
  
  # Get extrema and classify them into minima/maxima
  ddens <- c(NA, diff(dens))
  dddens <- c(NA, diff(ddens))
  # We use dens > 1e-10 to avoid extrema due to noise in the tails of the densities, 
  # use a higher value if needed to clean up even more here
  extrem <- c(NA, diff(sign(ddens)) != 0 & dens > 1e-10)
  extrem_type <- ifelse(ddens < 0, 1, 2)
  #plot(yvalues, dens, pch = 20)
  #points(yvalues, dens, col = extrem * (extrem_type + 2), cex = 2)
  
  extrems <- yvalues[extrem & ! is.na(extrem)]
  extrems_types <- ifelse(ddens < 0, "maximum", "minimum")[extrem & ! is.na(extrem)]

  list(values = data.frame(x = bin, y = yvalues, z = dens), 
       extrems = data.frame(x = bin, y = extrems, type = extrems_types))
}, .progress = "time")

# Extract results
grid <- ldply(window_analysis, function(o) o[["values"]])
extrems <- ldply(window_analysis, function(o) o[["extrems"]])
extrems$type2 <- extrems$type

extrems$type2[extrems$type =="maximum" & extrems$y > 0.5] <- "max1" 
extrems$type2[extrems$type =="maximum" & extrems$y < 0.5] <- "max2" 

# plot
fig.pot.cover <- 
  ggplot(NULL, aes(x = x, y = y)) + 
  geom_raster(aes(alpha = z), 
              data = subset(grid, z > 0.1))+ 
              #color = "black") + 
  geom_point(aes(color = type2), data = extrems) + 
  scale_color_manual(values = c("#40B0A6","#E1BE6A","white")) + 
  scale_alpha(guide = "none") + 
  theme_minimal() +
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "cover")
#  labs(x = "Aridity", y = "Cover", title = "Density analysis for cover") 


extrems_sub <- subset(extrems, extrems$type=="maximum")
  
  
fig.branches2.cover.bis = 
    ggplot(NULL) + 
    geom_point(aes(x = Aridity, y = imgcover, color = pretty_grps2), data = arid,alpha=.8,size=1) + 
    #scale_color_manual(values = c("#5ab4ac","#d8b365"))+ # colors from the 1st submission
    scale_color_manual(values = c("#40B0A6","#E1BE6A","black"))+ 
    theme_minimal() + 
    theme(legend.position="none",
          text = element_text(size=12))+
    labs(x = "aridity", 
         y = "cover")+
    #geom_point(aes(x = Aridity, y = mean.cover, size = n), size = 2, alpha=0.9, data = rolling_means2)
      #ggplot(NULL) +
      geom_point(aes(x = x, y = y, size = n, color = type), size = 2, alpha=0.9, data = extrems_sub)


# 
# 
# Do the same density analysis/detection of maxima/minima for MF
# 
# 
hres <- 0.01
bins <- seq(min(biocom[ ,"Aridity"]), 
            max(biocom[ ,"Aridity"]), 
            by = hres)
window_width <- hres * 15
yvalues <- seq(min(biocom[ ,"MF"]), 
               max(biocom[ ,"MF"]), 
               l = 1024)

window_analysis <- llply(bins, function(bin) { 
  ok <- subset(biocom, Aridity > (bin - window_width) & Aridity < (bin + window_width))
  dens <- with(ok, density(MF, bw = bw.nrd(MF), 
                           from = min(yvalues), to = max(yvalues), n = 2048))
  dens <- approx(dens[["x"]], dens[["y"]], yvalues)[["y"]]
  ddens <- c(NA, diff(dens))
  dddens <- c(NA, diff(ddens))
  extrem <- c(NA, diff(sign(ddens)) != 0 & dens > 0.5)
  extrem_type <- ifelse(ddens < 0, 1, 2)
  #plot(yvalues, dens, pch = 20)
  #points(yvalues, dens, col = extrem * (extrem_type + 2), cex = 2)
  
  # Extract extrema 
  extrems <- yvalues[extrem & ! is.na(extrem)]
  extrems_types <- ifelse(ddens < 0, "maximum", "minimum")[extrem & ! is.na(extrem)]
  
  list(values = data.frame(nvals = nrow(ok), 
                           x = bin, y = yvalues, z = dens), 
       extrems = data.frame(nvals = nrow(ok), 
                            x = bin, y = extrems, type = extrems_types))
}, .progress = "time")

grid <- ldply(window_analysis, function(o) o[["values"]])
extrems <- ldply(window_analysis, function(o) o[["extrems"]])

extrems$type2 <- extrems$type

extrems$type2[extrems$type =="maximum" & extrems$y > 0] <- "max1" 
extrems$type2[extrems$type =="maximum" & extrems$y < 0] <- "max2" 

fig.pot.mf <- 
  ggplot(NULL, aes(x = x, y = y)) + 
  geom_raster(aes(alpha = z), 
              data = subset(grid, z > 0.1))+ 
              #color = "black") + 
  geom_point(aes(color = type2), data = extrems) + 
  scale_color_manual(values = c("#40B0A6","#E1BE6A","white")) + 
  scale_alpha(guide = "none") + 
  theme_minimal() +
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "MF")

extrems_sub <- subset(extrems, extrems$type=="maximum")

fig.branches2.mf.bis = 
  ggplot(NULL, aes(x = Aridity, y = MF, color = pretty_grps2)) + 
  geom_point(data = arid,alpha=.6,size=1) + 
  scale_color_manual(values=c("#40B0A6","#E1BE6A","black"))+
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "MF")+
    geom_point(aes(x = x, y = y, size = n, color = type), size = 2, alpha=0.9, data = extrems_sub)



#plot_grid(fig.branches2.cover.bis, fig.branches2.mf.bis, fig.gaussian.fit2.cover, fig.gaussian.fit2.mf, labels= c("A","B","C","D"),ncol=2, nrow=2)


#### Fig 2
#fig2_2pot-cover-MF_densities_rev2
#landscape, 6x9.5
fig2 <- plot_grid(fig.pot.cover, fig.pot.mf, fig.gaussian.fit2.cover, fig.gaussian.fit2.mf, labels= c("A","B","C","D"),ncol=2, nrow=2)
grid.arrange(fig2, legend,nrow=2,ncol=1,heights = c(1,0.1))
#ggsave("./figures/fig2_2pot-cover-MF_densities.pdf", width = 9.5, height = 6)






#---------------------------------------------------------------------------
# Figure 3 : Difference in sp metrics
#---------------------------------------------------------------------------

indics_c <- indics
unique(indics_c$indic)
indics_c <- subset(indics_c, indic != "cover")
indics_c <- subset(indics_c, indic != "fmaxpatch")
indics_c <- subset(indics_c, indic != "moran")

indics_patch <- subset(indics_c, indic == "cutoff"  | indic == "logfmaxpatch" | indic == "slope")
indics_patch$indic <- as.factor(indics_patch$indic)
indics_patch$prop_order = factor(indics_patch$indic,levels=c("logfmaxpatch","slope","cutoff"),ordered=TRUE) 

indics_patch$pretty_grps2 <- as.factor(indics_patch$pretty_grps2)

facet_names3 <- list(
  'logfmaxpatch'="log(fmaxpatch) (***)",
  'slope'="slope (***)",
  'cutoff'="cutoff (***)"
)

facet_labeller3 <- function(variable,value){
  return(facet_names3[value])
}

row3 <- ggplot(indics_patch,aes(x=pretty_grps2, y=value,fill=pretty_grps2,alpha=0.7)) +
  facet_wrap(~prop_order,
             labeller = facet_labeller3,
             scale="free_y",nrow=1) +
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"),name=NULL)+ #name removes the name of the legend
  theme_minimal() +
  theme(legend.position="none",
        text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x= "",y="")+
  geom_boxplot()

indics_ews <- subset(indics_c, indic == "cv.variance" | indic == "sdr" | indic == "flowlength") 
indics_ews$indic <- as.factor(indics_ews$indic)
indics_ews$prop_order = factor(indics_ews$indic,levels=c("cv.variance","sdr","flowlength"),ordered=TRUE) #,"skewness"

facet_names2 <- list(
  'cv.variance'="cv (***)",
  'sdr'="sdr (NS)",
  'flowlength'="flowlength (***)"
)

facet_labeller2 <- function(variable,value){
  return(facet_names2[value])
}

row2 <- ggplot(indics_ews,aes(x=pretty_grps2, y=value, fill=pretty_grps2,alpha=0.7)) + 
  facet_wrap(~prop_order,
             labeller = facet_labeller2,
             scale="free_y",nrow=1) +
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"),name=NULL)+
  theme_minimal() +
  theme(legend.position="bottom",
        text = element_text(size=12))+
  scale_alpha(guide = 'none')+ #removes alpha legend
  theme(axis.text.x = element_blank())+
  labs(x= "",y="")+
  geom_boxplot()


#### Fig 3 old
# fig3_spmetrics_boxplots.pdf
# 5x9.5
#plot_grid(row3, row2, ncol=1, nrow=2)


#--------------------------------------------------------------------------# Figure 3 revised : no log for fmaxpatch
#--------------------------------------------------------------------------

indics_c <- indics
unique(indics_c$indic)
indics_c <- subset(indics_c, indic != "cover")
indics_c <- subset(indics_c, indic != "logfmaxpatch")
indics_c <- subset(indics_c, indic != "moran")

# transform indics_c in long format
indics_c$indic <- as.factor(indics_c$indic)
indics_sub <- indics_c[,c("indic","value","plotid","pretty_grps2")]
indics_wide <- spread(indics_sub, indic, value)
indics_wide$pretty_grps2 <- as.factor(indics_wide$pretty_grps2)

#fmaxpatch
fig.boxplot.fmaxpatch = ggplot(indics_wide, aes(x = pretty_grps2, y = fmaxpatch, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "fmaxpatch (***; log scale)")+
  scale_y_continuous(trans='log10')+
  theme(legend.position="none")

#slope
fig.boxplot.slope = ggplot(indics_wide, aes(x = pretty_grps2, y = slope, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "slope (***)")+
  theme(legend.position="none")

#cutoff
fig.boxplot.cutoff = ggplot(indics_wide, aes(x = pretty_grps2, y = cutoff, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "cutoff (***)")+
  theme(legend.position="none")

#cv
fig.boxplot.cv = ggplot(indics_wide, aes(x = pretty_grps2, y = cv.variance, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "cv (***)")+
  theme(legend.position="none")

#moran
#fig.boxplot.moran = ggplot(indics_wide, aes(x = pretty_grps2, y = moran, fill=pretty_grps2)) + 
 # geom_boxplot(alpha=.7) + 
#  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
#  theme_minimal() + 
#  theme(text = element_text(size=12))+
  #theme(axis.text.x = element_text(angle=45))+
#  theme(axis.text.x = element_blank())+
#  labs(x = " ", 
#       y = "moran (NS)")+
#  theme(legend.position="none")

#sdr
fig.boxplot.sdr = ggplot(indics_wide, aes(x = pretty_grps2, y = sdr, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "sdr (NS)")+
  theme(legend.position="none")

#fl
fig.boxplot.fl = ggplot(indics_wide, aes(x = pretty_grps2, y = flowlength, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "flowlength (***)")+
  theme(legend.position="none")


top_row <- plot_grid(fig.boxplot.fmaxpatch, fig.boxplot.slope, fig.boxplot.cutoff, ncol=3, nrow=1)

bottom_row2 <- plot_grid(fig.boxplot.cv, fig.boxplot.sdr,fig.boxplot.fl, ncol=3,nrow=1)  


#### Fig 3 no log
# fig3_spmetrics_boxplots_nolog.pdf
# lanbdscape, 6 x 9.5
fig3 <- plot_grid(top_row, bottom_row2, ncol = 1)
grid.arrange(fig3, legend,nrow=2,ncol=1,heights = c(1,0.1))
#ggsave("./figures/fig2_2pot-cover-MF_densities.pdf", width = 9.5, height = 6)


#---------------------------------------------------------------------------
# Figure S10 : 2 groups of cover
#---------------------------------------------------------------------------

fig.density2d.cover = ggplot(arid)+
  geom_point(aes(x=Aridity,y=imgcover, color = pretty_grps2c))+
  scale_color_manual(values = c("#40B0A6","#E1BE6A"))+ 
  stat_density2d(aes(x=Aridity,y=imgcover),color="black")+
  theme_minimal()+
  labs(x = "aridity", 
       y = "cover")+
  theme(legend.position="none")

rolling_means2c <- rollply(arid, ~ Aridity | pretty_grps2c, 
                          wdw.size = 0.07, grid_npts = 128, 
                          summarise, 
                          mean.cover = mean(imgcover), 
                          n = length(imgcover)
)

# Keep only areas with a lot of points
rolling_means2c <- subset(rolling_means2c, n > 30)

fig.branches2c.cover = ggplot(NULL, aes(x = Aridity, y = imgcover, color = pretty_grps2c)) + 
  geom_point(data = arid,alpha=.6,size=1) + 
  scale_color_manual(values = c("#40B0A6","#E1BE6A"))+ 
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=12))+
  labs(x = "aridity", 
       y = "cover")+
  geom_point(aes(x = Aridity, y = mean.cover, size = n), size = 2, alpha=0.4, data = rolling_means2c,color="black")

fig.gaussian.fit2c.cover = ggplot(arid, aes(x = imgcover)) + 
  geom_density(aes(fill = pretty_grps2c), 
               bw = .04,alpha= .7) +
  geom_point(aes(color = pretty_grps2c), 
             y = rnorm(nrow(arid), 0, 0.02), alpha = .4) +
  scale_colour_manual(values=c("#40B0A6","#E1BE6A"))+
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  labs(x = "cover", 
       y = "density")+
  theme(legend.position="none")

#boxplot
fig.boxplot.cover = ggplot(arid, aes(x = pretty_grps2c, y = imgcover, fill=pretty_grps2c)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "cover (***)")+
  theme(legend.position="none",legend.title=element_blank())


#### Fig S10
# figS10_2groupscover.pdf
# landscape, 6 x 9.5
top_row <- plot_grid(fig.density2d.cover, fig.branches2c.cover, fig.gaussian.fit2c.cover, fig.boxplot.cover, labels= c("A","B","C","D"),ncol=2, nrow=2)

grid.arrange(top_row, legend,nrow=2,ncol=1,heights = c(1,0.1))



#---------------------------------------------------------------------------
# Figure S12 : Vege type per branch
#---------------------------------------------------------------------------

arid.branch1 = arid[arid$grps2==1,]
arid.branch1$veg_type=as.factor(arid.branch1$veg_type)
arid.branch2 = arid[arid$grps2==2,]
arid.branch2$veg_type=as.factor(arid.branch2$veg_type)

fig.branch1.vegtype = ggplot(arid.branch1)+
  geom_point(aes(x=Aridity,y=imgcover,color=veg_type))+
  #scale_color_manual(values = c("olivedrab4", "olivedrab2","yellow3"))+ 
  scale_color_manual(values=c("#407a78","#feb624","#d44206"))+
  theme_minimal()+
  ggtitle("low cover, low MF")+
  labs(x= "aridity",y="cover")+
  theme(legend.position="none")

fig.branch2.vegtype = ggplot(arid.branch2)+
  geom_point(aes(x=Aridity,y=imgcover,color=veg_type))+
  #scale_color_manual(values = c("olivedrab4", "olivedrab2","yellow3"))+ 
  scale_color_manual(values=c("#407a78","#feb624","#d44206"))+
  theme_minimal()+
  ggtitle("high cover, high MF")+
  labs(x= "aridity",y="cover")+
  theme(legend.position="none")


#---------------------
# calculate prop

### we want to know the prop of each vege type on each branch
df <- subset(arid[,c("veg_type","grps2")])
df$veg_type[df$veg_type==1] <- "grassland" 
df$veg_type[df$veg_type==2] <- "shrubland" 
df$veg_type[df$veg_type==4] <- "savanna" 

# count the nb of veg_type in each branch
tab <- data.frame(branch=1,veg_type="grassland", nb=length(df$veg_type[df$veg_type=="grassland"&df$grps2==1]),
                  frac=length(df$veg_type[df$veg_type=="grassland"&df$grps2==1])/length(df$veg_type[df$veg_type=="grassland"]),
                  frac2=length(df$veg_type[df$veg_type=="grassland"&df$grps2==1])/length(df$veg_type[df$grps2==1]))

tab2 <- data.frame(branch=2,veg_type="grassland", nb=length(df$veg_type[df$veg_type=="grassland"&df$grps2==2]),
                   frac=length(df$veg_type[df$veg_type=="grassland"&df$grps2==2])/length(df$veg_type[df$veg_type=="grassland"]),
                   frac2=length(df$veg_type[df$veg_type=="grassland"&df$grps2==2])/length(df$veg_type[df$grps2==2]))


tab<-rbind(tab,tab2)

tab1 <- data.frame(branch=1,veg_type="shrubland", nb=length(df$veg_type[df$veg_type=="shrubland"&df$grps2==1]),
                   frac=length(df$veg_type[df$veg_type=="shrubland"&df$grps2==1])/length(df$veg_type[df$veg_type=="shrubland"]),
                   frac2=length(df$veg_type[df$veg_type=="shrubland"&df$grps2==1])/length(df$veg_type[df$grps2==1]))

tab2 <- data.frame(branch=2,veg_type="shrubland", nb=length(df$veg_type[df$veg_type=="shrubland"&df$grps2==2]),
                   frac=length(df$veg_type[df$veg_type=="shrubland"&df$grps2==2])/length(df$veg_type[df$veg_type=="shrubland"]),
                   frac2=length(df$veg_type[df$veg_type=="shrubland"&df$grps2==2])/length(df$veg_type[df$grps2==2]))


tab<-rbind(tab,tab1,tab2)

tab1 <- data.frame(branch=1,veg_type="savanna", nb=length(df$veg_type[df$veg_type=="savanna"&df$grps2==1]),
                   frac=length(df$veg_type[df$veg_type=="savanna"&df$grps2==1])/length(df$veg_type[df$veg_type=="savanna"]),
                   frac2=length(df$veg_type[df$veg_type=="savanna"&df$grps2==1])/length(df$veg_type[df$cgrps2==1]))

tab2 <- data.frame(branch=2,veg_type="savanna", nb=length(df$veg_type[df$veg_type=="savanna"&df$grps2==2]),
                   frac=length(df$veg_type[df$veg_type=="savanna"&df$grps2==2])/length(df$veg_type[df$veg_type=="savanna"]),
                   frac2=length(df$veg_type[df$veg_type=="savanna"&df$grps2==2])/length(df$veg_type[df$grps2==2]))

tab<-rbind(tab,tab1,tab2)

tab$veg_type=as.factor(tab$veg_type)
tab$branch[tab$branch==1] <- "low" 
tab$branch[tab$branch==2] <- "high" 
tab$branch=as.factor(tab$branch)

fig.frac <- ggplot(tab, aes(fill=veg_type, y=nb, x=branch)) + 
  theme_minimal() + 
  labs(x = "group of sites", 
       y = "number of sites")+
  theme(axis.title=element_text(size=10))+
  scale_fill_manual(values=c("#407a78","#d44206","#feb624"))+
  geom_bar(position="stack", stat="identity")+
  labs(fill = "vegetation type")


#### Fig S12
plot_grid(fig.branch2.vegtype, fig.branch1.vegtype,fig.frac, labels=c("A","B","C"),ncol=2,nrow=2)
ggsave("./figures/figS12_vege_type.pdf.pdf", width = 9.5, height = 6)



#---------------------------------------------------------------------------
# Figure S13 : Envi var 
#---------------------------------------------------------------------------

fig.boxplot2.cover = ggplot(arid, aes(x = pretty_grps2, y = imgcover, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "cover (***)")+
  theme(legend.position="none")


fig.boxplot2.mf = ggplot(arid, aes(x = pretty_grps2, y = MF, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
   theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "MF (***)")+
  theme(legend.position="none")

fig.boxplot2.arid = ggplot(arid, aes(x = pretty_grps2, y = Aridity, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Aridity (***)")+
  theme(legend.position="none")

fig.boxplot2.prod = ggplot(arid, aes(x = pretty_grps2, y = prod, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Productivity (***)")+
  theme(legend.position="none")

fig.boxplot2.sand = ggplot(arid, aes(x = pretty_grps2, y = Sand, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Sand (***)")+
  theme(legend.position="none")

fig.boxplot2.sr = ggplot(arid, aes(x = pretty_grps2, y = sr, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#40B0A6","#E1BE6A"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
   theme(axis.text.x = element_blank())+
  labs(x = " ", 
       y = "Sp. richness (***)")+
  theme(legend.position="none")


top_row <- plot_grid(fig.branches2.cover, fig.branches2.mf, fig.gaussian.fit2.cover, fig.gaussian.fit2.mf, labels= c("A","B","C","D"),ncol=2, nrow=2)

bottom_row2 <- plot_grid(fig.boxplot2.cover, fig.boxplot2.mf, fig.boxplot2.prod,fig.boxplot2.arid, fig.boxplot2.sand, fig.boxplot2.sr, labels = c('E', '', '','', '',''), ncol=6,nrow=1)


#### Fig S13
#figS13_cover_mf_2branches.pdf
# 8 x 9.5
#plot_grid(top_row, bottom_row2, ncol = 1)
grid.arrange(top_row, bottom_row2, legend, nrow=3,ncol=1,heights = c(2,1,0.1))

