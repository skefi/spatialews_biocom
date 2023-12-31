# Perform and plot PCA analysis of data
# Version for revision PNAS
# Figures S4, S5


library(tidyr)
library(umap)
library(GGally) #for ggpairs
library(ade4)
library(plotly)
library(ggfortify)


#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------

NPERM = 199 #value use for final analyses: 199
BOOTN = 2999 #value use for final analyses: 2999
path_output <- here::here("outputs")
ALPHA <- 0.05 # Significance level for indicator trends

indic_order <- c('logfmaxpatch', 'slope','cutoff', 'flowlength','moran', 'sdr','cv.variance') 

#---------------------------------------------------------------------------
# Load data
#---------------------------------------------------------------------------

load(here::here("outputs", "biocom-grps.rda"))
load(file.path(path_output,"data_biocom.rda"))
arid <- biocom

filename <- paste0("indics-data50-grps_Nperm_", NPERM,"_rev.rda")
load(file.path(path_output,filename))
#unique(indics$indic)


#---------------------------------------------------------------------------
# PCA of data
#---------------------------------------------------------------------------

# transform indics into long format 
indics_sub <- subset(indics, select=c(indic, value, plotid,pretty_grps2))
indics_wide <- spread(indics_sub,indic,value)
indics_wide$pretty_grps2 <- as.factor(indics_wide$pretty_grps2)

indics_wide_sub <- subset(indics_wide, select=-fmaxpatch)
indics_wide_sub <- subset(indics_wide_sub, select=-cover)

arid.pca <- indics_wide_sub[,3:length(indics_wide_sub)]


#### Fig. S5
ggpairs(arid.pca)
ggsave("./figures/figS5_pairs.pdf", width = 9.5, height = 9.5)



#indics_wide_sub <- subset(indics_wide_sub, select=-flowlength)
indics_wide_sub <- subset(indics_wide_sub, select=-moran)
indics_wide_sub <- na.omit(indics_wide_sub)

arid.pca <- indics_wide_sub[,3:length(indics_wide_sub)]

res.pca <- dudi.pca(arid.pca,scannf=FALSE,nf=2)
summary(res.pca)
res.pca$co # coordinates of the 3 variables on the 2 axis
inertia.pca <- inertia.dudi(res.pca,row.inertia=TRUE)
scatter(res.pca,posieig="none",clab.row=0)

## Using plotly

pca_res <- prcomp(arid.pca, scale. = TRUE)
PCAvalues <- data.frame(grps2=indics_wide_sub$pretty_grps2,pca_res$x) # extract PCA axes
PCAloadings <- data.frame(Variables=rownames(pca_res$rotation),pca_res$rotation)

ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = grps2)) +
  geom_point(size = 3) +
  scale_color_manual(values=c("#40B0A6","#E1BE6A"))+
  geom_segment(data = PCAloadings, aes(x = 0, y = 0, xend = (PC1*4),yend = (PC2*4)),
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black") +
  annotate("text", x = (PCAloadings$PC1*5.5), y = (PCAloadings$PC2*4.5),
           label = PCAloadings$Variables)+
  theme_minimal()+
  theme(text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "top")


#### Fig. S4
ggsave("./figures/figS4_pca_data.pdf", width = 7, height = 7)

  
