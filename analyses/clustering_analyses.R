# Calculates the trend (slope) of each of the spatial metrics along the aridity gradient

library(ggplot2)
#https://cran.r-project.org/web/packages/dendextend/vignettes/FAQ.html
library(dendextend)
library(mclust)

source(here::here("R", "functions_helper.R"))


#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------

# From _targets.R
NPERM = 199 #value use for final analyses: 199
path_output <- here::here("outputs")


#---------------------------------------------------------------------------
# LOAD data
#---------------------------------------------------------------------------

filename <- paste0("indics-data-grps_Nperm_", NPERM,".rda")
load(file.path(path_output,filename))

load(file.path(path_output,"data_biocom.rda"))

load(file.path(path_output,"biocom-grps.rda"))
arid <- biocom


#---------------------------------------------------------------------------
# Visualize distributions
#---------------------------------------------------------------------------

clust_dat <- ourdata[["biocom"]][ ,c("MF", "imgcover")]

GGally::ggpairs(clust_dat, 
                aes(alpha = .3))

#fig.gaussian.cover = 
  ggplot(arid, aes(x = imgcover)) + 
  geom_density(fill="grey",bw = .04,alpha= .7) +
  theme_minimal() + 
  theme(text = element_text(size=10))+
  labs(x = "cover", 
       y = "density")+
  theme(legend.position="none")

#fig.gaussian.MF = 
  ggplot(arid, aes(x = MF)) + 
  geom_density(fill="grey",bw = .14,alpha= .7) +
  theme_minimal() + 
  theme(text = element_text(size=10))+
  labs(x = "MF", 
       y = "density")+
  theme(legend.position="none")


#---------------------------------------------------------------------------
# Hierarchical clustering on standardized data 
#---------------------------------------------------------------------------

clust_dat_scaled <- apply(clust_dat, 2, function(X) ( X - mean(X) ) / sd(X) )
hcl <- hclust(dist(clust_dat_scaled), method = "ward.D") # We use the ward distance because it is the one that works well for "globular" clusters

plot(hcl)

# The hierchichal clustering suggests that there are two groups in the data. 
# Let's visualize them
classif_hclust <- cutree(hcl, k = 2)
envdat_clust <- data.frame(clust_dat_scaled, 
                           classif = as.factor(classif_hclust))
GGally::ggpairs(envdat_clust, 
                aes(color = classif, alpha = .3))


#---------------------------------------------------------------------------
# Elbow rule
#---------------------------------------------------------------------------
#Let's try the simple "elbow rule": adding a group to the classification should only provide a reduced improvement in accuracy. Our metric of accuracy will be the percent of variance explained by the classification, i.e. the sum of square of each group to its centroid. 

maxgroups <- 20
elbow_data <- lapply(seq.int(maxgroups), function(n) { 
  classif <- cutree(hcl, k = n)
  # Compute distance to centroids
  within_groups_ss <- unlist(lapply(unique(classif), function(groupn) { 
    dists_to_centroid <- apply(clust_dat_scaled[classif == groupn, ], 2, 
                               function(X) { ( X - mean(X) )^2 })
    return(sum(dists_to_centroid))
  }))
  
  total_ss <- sum(apply(clust_dat_scaled, 2, function(X) { ( X - mean(X) )^2 }))
  
  # Percentage explained 
  data.frame(ngroup = n, 
             perc_var = 1 - sum(within_groups_ss) / total_ss)
})
elbow_data <- do.call(rbind, elbow_data)

elbow_data[ ,"improvement"] <- c(NA, diff(elbow_data[, "perc_var"]))

ggplot(elbow_data, aes(x = ngroup, y = perc_var)) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(breaks = seq(1, maxgroups)) + 
  labs(x = "Number of groups", 
       y = "Percent variance explained (within group sum of square / total SS)")

ggplot(elbow_data, aes(x = ngroup, y = improvement)) + 
  geom_point() + 
  geom_line() + 
  scale_x_continuous(breaks = seq(1, maxgroups)) + 
  labs(x = "Number of groups", 
       y = "Difference in percent variance explained as we add group number")

# After we go from one group to two, we get an improvement of about 45%, then it is under 10%. It can be an argument to say that we have two groups.


#---------------------------------------------------------------------------
# Gaussian mixture
#---------------------------------------------------------------------------
#We can use an information criterion to test how many groups there are in the  data. This requires defining a likelihood, which is often done by assuming that  the points are drawn from a mixture of normal distributions. Mclust provides  what's necessary to fit the distributions. 

clustfits <- mclust::mclustBIC(clust_dat_scaled) #, verbose = (.PROGRESS == "time") )
summary(clustfits)
figclust <- plot(clustfits)

## extract the case for 2 clusters:
mod <- Mclust(clust_dat_scaled, modelName="VEI")
summary(mod, parameters = TRUE)
plot(mod, what = "classification")

length(which(arid$grps2!=mod$classification)) #14
length(which(arid$grps2==mod$classification)) #331

#Here the lines are different types of mixtures. For example, EEE means that scale, variance and covariance of all the (multidimensional) normal distributions in the mixture are equal. 
# Regardless of the type of distribution, their is a huge increase of fit quality when we reach 2 groups, then the increase is not as important. The best models based on BIC are within the 3 and 4 groups, and they are a bit better than models with two groups. Let's look at the classification it gives us. 

# try G=2, G=3
classif_mclust2 <- mclust::Mclust(clust_dat_scaled, G = 2, modelNames = "VVV")$classification
envdat_mclust2 <- data.frame(clust_dat_scaled, 
                             classif = as.factor(classif_mclust2))
GGally::ggpairs(envdat_mclust2, 
                aes(color = classif, alpha = .3))

classif_mclust3 <- mclust::Mclust(clust_dat_scaled, G = 3, modelNames = "VVV")$classification
envdat_mclust3 <- data.frame(clust_dat_scaled, 
                             classif = as.factor(classif_mclust3))
GGally::ggpairs(envdat_mclust3, 
                aes(color = classif, alpha = .3))

# One reason why Mclust always gives more group that intuition is because the form  of the distributions is assumed to be normal. This means that all deviations from a mixture of normal distribution (which is bound to happen a lot in messy datasets) is accomodated by adding distributions. 




#---------------------------------------------------------------------------
# Color the leafs of the tree by the final group chosen (grps2 and grps3) to see the agreement between approaches
#---------------------------------------------------------------------------

dend <- as.dendrogram(hclust(dist(clust_dat_scaled), method = "ward.D"))
od <- order.dendrogram(dend)

# Hierachical clustering
dend1 <- color_branches(dend, k = 2,col=c("#d8b365","#5ab4ac"),groupLabels=TRUE)
dend1 <- set(dend1, "labels_cex", 0.01)
h2groups <- plot(dend1, main = "2 groups")

dend2 <- color_branches(dend, k = 3,col=c("#FC8D62","#8DA0CB","#66C2A5"),groupLabels=TRUE)
dend2 <- set(dend2, "labels_cex", 0.01)
h3groups <- plot(dend2, main = "3 groups")

# By the group eventually chosen, i.e. grps2 and grps3
cluster2 <- c()
for(i in 1: length(order.dendrogram(dend))){
  cluster2[i] <- arid$grps2[od[i]]
}
dend3 <- color_branches(dend, clusters=cluster2)
plot(dend3, main = "2 groups")

#par(mfrow = c(1, 2))

## By mclust
#classif_mclust2 <- mclust::Mclust(clust_dat_scaled, G = 2, modelNames = "VVV")$classification
#od <- order.dendrogram(dend)
cluster3 <- c()
for(i in 1: length(order.dendrogram(dend))){
  cluster3[i] <- classif_mclust2[od[i]]
}
dend4 <- color_branches(dend, clusters=cluster3,col=c("#d8b365","#5ab4ac"),groupLabels=TRUE)
dend4 <- set(dend4, "labels_cex", 0.01)
plot(dend4, main = "2 groups (Gaussian)")
#dend3 <- color_branches(dend, clusters=classif_mclust2)
#plot(dend3, main = "2 groups (gaussian)")

