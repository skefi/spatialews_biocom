# Allows to plot the figures 2, 3, S7 and S8



#devtools::install_github('alexgenin/rollply')
library(rollply)
library(ggplot2)

#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
# From _targets.R
NPERM = 3 #value use for final analyses: 199
path_output <- here::here("outputs")

#---------------------------------------------------------------------------
# Load data
#---------------------------------------------------------------------------

load(here::here("outputs", "biocom-grps.rda"))
#load(file.path(path_output,"data_biocom.rda"))
arid <- biocom

filename <- paste0("indics-data-grps_Nperm_", NPERM,".rda")
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
  scale_color_manual(values = c("#5ab4ac","#d8b365"))+ 
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=10))+
  labs(x = "aridity", 
       y = "cover")+
  geom_point(aes(x = Aridity, y = mean.cover, size = n), size = 2, alpha=0.9, data = rolling_means2)

fig.branches2.mf = ggplot(NULL, aes(x = Aridity, y = MF, color = pretty_grps2)) + 
  geom_point(data = arid,alpha=.6,size=1) + 
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(legend.position="none",
        text = element_text(size=10))+
  labs(x = "aridity", 
       y = "MF")+
  geom_point(aes(x = Aridity, y = mean.mf, size = n), size = 2, alpha=0.9, data = rolling_means2)

fig.gaussian.fit2.cover = ggplot(arid, aes(x = imgcover)) + 
  geom_density(aes(fill = pretty_grps2), 
               bw = .04,alpha= .7) +
  geom_point(aes(color = pretty_grps2), 
             y = rnorm(nrow(arid), 0, 0.02), alpha = .4) +
  scale_colour_manual(values=c("#5ab4ac","#d8b365"))+
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  labs(x = "cover", 
       y = "density")+
  theme(legend.position="none")

fig.gaussian.fit2.mf = ggplot(arid, aes(x = MF)) + 
  geom_density(aes(fill = pretty_grps2), 
               bw = .14,alpha= .7) + 
  geom_point(aes(color = pretty_grps2), 
             y = rnorm(nrow(arid), 0, 0.02), alpha = .4) + 
  scale_color_manual(values=c("#5ab4ac","#d8b365"))+
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  labs(x = "MF", 
       y = "density")+
  theme(legend.position="none")

#5 x 9.5
plot_grid(fig.branches2.cover, fig.branches2.mf, fig.gaussian.fit2.cover, fig.gaussian.fit2.mf, labels= c("A","B","C","D"),ncol=2, nrow=2)



#---------------------------------------------------------------------------
# Figure 3 : Difference in sp metrics
#---------------------------------------------------------------------------

indics_c <- indics
unique(indics_c$indic)
indics_c <- subset(indics_c, indic != "cover")
indics_c <- subset(indics_c, indic != "skewness")

indics_patch <- subset(indics_c, indic == "cutoff"  | indic == "logfmaxpatch" | indic == "plrange" | indic == "slope")
indics_patch$indic <- as.factor(indics_patch$indic)
indics_patch$prop_order = factor(indics_patch$indic,levels=c("plrange","logfmaxpatch","slope","cutoff"),ordered=TRUE) 

facet_names3 <- list(
  'plrange'="plr (***)",
  'logfmaxpatch'="logfmaxpatch (***)",
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
  scale_fill_manual(values=c("#5ab4ac","#d8b365"),name=NULL)+ #name removes the name of the legend
  theme_minimal() +
  theme(legend.position="none",
        text = element_text(size=12))+
  theme(axis.text.x = element_blank())+
  labs(x= "",y="")+
  geom_boxplot()

indics_ews <- subset(indics_c, indic == "moran"  | indic == "cv.variance" | indic == "sdr" | indic == "flowlength") 
indics_ews$indic <- as.factor(indics_ews$indic)
indics_ews$prop_order = factor(indics_ews$indic,levels=c("cv.variance","moran","sdr","flowlength"),ordered=TRUE) #,"skewness"

facet_names2 <- list(
  'cv.variance'="cv (***)",
  'moran'="moran (NS)",
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
  scale_fill_manual(values=c("#5ab4ac","#d8b365"),name=NULL)+
  theme_minimal() +
  theme(legend.position="bottom",
        text = element_text(size=12))+
  scale_alpha(guide = 'none')+ #removes alpha legend
  theme(axis.text.x = element_blank())+
  labs(x= "",y="")+
  geom_boxplot()

# Fig 2 5x9.5
plot_grid(row3, row2, ncol=1, nrow=2)



#---------------------------------------------------------------------------
# Figure S7 : Vege type per branch
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
  ggtitle("Low branch (low, low)")+
  labs(x= "Aridity",y="cover")+
  theme(legend.position="none")

fig.branch2.vegtype = ggplot(arid.branch2)+
  geom_point(aes(x=Aridity,y=imgcover,color=veg_type))+
  #scale_color_manual(values = c("olivedrab4", "olivedrab2","yellow3"))+ 
  scale_color_manual(values=c("#407a78","#feb624","#d44206"))+
  theme_minimal()+
  ggtitle("High branch (high, high)")+
  labs(x= "Aridity",y="cover")+
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
  geom_bar(position="stack", stat="identity")

plot_grid(fig.branch2.vegtype, fig.branch1.vegtype,fig.frac, labels=c("A","B","C"),ncol=2,nrow=2)


#---------------------------------------------------------------------------
# Figure S8 : Envi var 
#---------------------------------------------------------------------------

fig.boxplot2.cover = ggplot(arid, aes(x = pretty_grps2, y = imgcover, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  #theme(axis.text.x = element_text(angle=45))+
  theme(axis.text.x = element_blank())+
  labs(x = "branch", 
       y = "cover (***)")+
  theme(legend.position="none")


fig.boxplot2.mf = ggplot(arid, aes(x = pretty_grps2, y = MF, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
   theme(axis.text.x = element_blank())+
  labs(x = "branch", 
       y = "MF (***)")+
  theme(legend.position="none")

fig.boxplot2.arid = ggplot(arid, aes(x = pretty_grps2, y = Aridity, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_blank())+
  labs(x = "branch", 
       y = "Aridity (***)")+
  theme(legend.position="none")

fig.boxplot2.prod = ggplot(arid, aes(x = pretty_grps2, y = prod, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_blank())+
  labs(x = "branch", 
       y = "Productivity (***)")+
  theme(legend.position="none")

fig.boxplot2.sand = ggplot(arid, aes(x = pretty_grps2, y = Sand, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
  theme(axis.text.x = element_blank())+
  labs(x = "branch", 
       y = "Sand (***)")+
  theme(legend.position="none")

fig.boxplot2.sr = ggplot(arid, aes(x = pretty_grps2, y = sr, fill=pretty_grps2)) + 
  geom_boxplot(alpha=.7) + 
  scale_fill_manual(values=c("#5ab4ac","#d8b365"))+
  theme_minimal() + 
  theme(text = element_text(size=10))+
   theme(axis.text.x = element_blank())+
  labs(x = "branch", 
       y = "Sp. richness (***)")+
  theme(legend.position="none")


top_row <- plot_grid(fig.branches2.cover, fig.branches2.mf, fig.gaussian.fit2.cover, fig.gaussian.fit2.mf, labels= c("A","B","C","D"),ncol=2, nrow=2)

bottom_row2 <- plot_grid(fig.boxplot2.cover, fig.boxplot2.mf, fig.boxplot2.prod,fig.boxplot2.arid, fig.boxplot2.sand, fig.boxplot2.sr, labels = c('E', '', '','', '',''), ncol=6,nrow=1)

# 8 x 9.5
plot_grid(top_row, bottom_row2, ncol = 1)

