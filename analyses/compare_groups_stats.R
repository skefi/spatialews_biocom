# Statistically compare the 2 and 3 groups of sites

library(lme4)
library(plyr)


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
# Effect of aridity on spatial metrics 
#---------------------------------------------------------------------------

## Role of aridity
data <- subset(indics, indics$indic=="cover")

mod <- lmer(data = data, value ~ Aridity + ( 1 |plotn)) 
#summary(mod)
null <- lmer(data = data, value ~ ( 1 |plotn)) 
anovtable <- as.data.frame(anova(mod,null))
anovtable
comp <- anovtable$`Pr(>Chisq)`[2]
#value = fixef(mod)["Aridity"]
ifelse(comp < 0.05, "p < 0.05","not significant")
comp


# for fl
data <- subset(indics, indics$indic=="flowlength")
mod <- lm(data = data, value ~ Aridity) 
summary(mod)


#---------------------------------------------------------------------------
# 2 groups: trends of cover and MF with aridity in each of the two groups
#---------------------------------------------------------------------------

dat.1 = arid[arid$grps2==1,]
dat.2 = arid[arid$grps2==2,]

lm_mf1 <- lm(formula = MF ~ Aridity, data = dat.1)
summary(lm_mf1)

lm_c1 <- lm(formula = imgcover ~ Aridity, data = dat.1)
summary(lm_c1)

lm_mf2 <- lm(formula = MF ~ Aridity, data = dat.2)
summary(lm_mf2)

lm_c2 <- lm(formula = imgcover ~ Aridity, data = dat.2)
summary(lm_c2)


#---------------------------------------------------------------------------
# 2 groups t-tests: compare envi variables between groups
#---------------------------------------------------------------------------

aridity.1 = arid[arid$grps2==1,c("Aridity")]
aridity.2 = arid[arid$grps2==2,c("Aridity")]
t.test(aridity.1,aridity.2) # t = 9.2969, df = 254.92, p-value < 2.2e-16

mf.1 = arid[arid$grps2==1,c("MF")]
mf.2 = arid[arid$grps2==2,c("MF")]
t.test(mf.1,mf.2) # p-value < 2.2e-16

sr.1 = arid[arid$grps2==1,c("sr")]
sr.2 = arid[arid$grps2==2,c("sr")]
t.test(sr.1,sr.2) # p-value = 0.0009092

prod.1 = arid[arid$grps2==1,c("prod")]
prod.2 = arid[arid$grps2==2,c("prod")]
t.test(prod.1,prod.2) # p-value < 2.2e-16

sand.1 = arid[arid$grps2==1,c("Sand")]
sand.2 = arid[arid$grps2==2,c("Sand")]
t.test(sand.1,sand.2) # p-value < 2.2e-16

cover.1 = arid[arid$grps2==1,c("imgcover")]
cover.2 = arid[arid$grps2==2,c("imgcover")]
t.test(cover.1,cover.2) # p-value < 2.2e-16


## ANOVA + bonferroni test
signi.comp2 = function(sub.table){
  temp = anova(lm(sub.table$value ~ sub.table$pretty_grps2))
  p = temp$"Pr(>F)"[1]
  #https://www.r-bloggers.com/r-tutorial-series-anova-pairwise-comparison-methods/
  #d = diff(temp$estimate)/abs(temp$estimate[1])
  temp2 = pairwise.t.test(sub.table$value, sub.table$pretty_grps2, p.adj = "bonferroni")
  #The Bonferroni adjustment simply divides the Type I error rate (.05) by the number of tests (in this case, three). Hence, this method is often considered overly conservative. The Bonferroni adjustment can be made using p.adj = “bonferroni” in the pairwise.t.test() function.
  #The Holm adjustment sequentially compares the lowest p-value with a Type I error rate that is reduced for each consecutive test. In our case, this means that our first p-value is tested at the .05/3 level (.017), second at the .05/2 level (.025), and third at the .05/1 level (.05). This method is generally considered superior to the Bonferroni adjustment and can be employed using p.adj = “holm” in the pairwise.t.test() function.
  d = temp2$"p.value"
  d = t(as.vector(d)) # 1 :1vs2, 2: 1vs3, 3: NA, 4: 2vs3  
  d = as.data.frame(d)
  res = cbind(p,d)
  return(res)
}

test <- ddply(indics,~indic,signi.comp2) 
test 


## lm
data <- subset(indics, indics$indic=="flowlength")
mod <- lm(data = data, value ~ Aridity) 
#mod <- lm(data = data, value ~ Aridity) 
summary(mod)
anova(lm(data$value ~ data$grps2))


#---------------------------------------------------------------------------
# 2 groups : mixed models, compare spatial spatial metrics between the 2 groups
#---------------------------------------------------------------------------

data <- subset(indics, indics$indic=="logfmaxpatch")
data$grps2<-as.factor(data$grps2)
mod <- lmer(data = data, value ~ grps2 + ( 1 |plotn)) 
#summary(mod)
null <- lmer(data = data, value ~ ( 1 |plotn)) 
anovtable <- as.data.frame(anova(mod,null))
anovtable
comp <- anovtable$`Pr(>Chisq)`[2]
#value = fixef(mod)["Aridity"]
ifelse(comp < 0.05, "p < 0.05","not significant")
comp

# for fl
data <- subset(indics, indics$indic=="flowlength")
#data$value<-log10(data$value)
data$grps2<-as.factor(data$grps2)
mod <- lm(data = data, value ~ grps2) 
summary(mod)

anova(lm(data$value ~ data$grps2))



#---------------------------------------------------------------------------
# 3 groups ANOVA
#---------------------------------------------------------------------------

signi.comp3 = function(sub.table){
  temp = anova(lm(sub.table$value ~ sub.table$pretty_grps3))
  p = temp$"Pr(>F)"[1]
  #https://www.r-bloggers.com/r-tutorial-series-anova-pairwise-comparison-methods/
  #d = diff(temp$estimate)/abs(temp$estimate[1])
  temp2 = pairwise.t.test(sub.table$value, sub.table$pretty_grps3, p.adj = "bonferroni")
  #The Bonferroni adjustment simply divides the Type I error rate (.05) by the number of tests (in this case, three). Hence, this method is often considered overly conservative. The Bonferroni adjustment can be made using p.adj = “bonferroni” in the pairwise.t.test() function.
  #The Holm adjustment sequentially compares the lowest p-value with a Type I error rate that is reduced for each consecutive test. In our case, this means that our first p-value is tested at the .05/3 level (.017), second at the .05/2 level (.025), and third at the .05/1 level (.05). This method is generally considered superior to the Bonferroni adjustment and can be employed using p.adj = “holm” in the pairwise.t.test() function.
  d = temp2$"p.value"
  d = t(as.vector(d)) # 1 :1vs2, 2: 1vs3, 3: NA, 4: 2vs3  
  d = as.data.frame(d)
  res = cbind(p,d)
  return(res)
}

test <- ddply(indics,~indic,signi.comp3) 
test 


#---------------------------------------------------------------------------
# 3 groups mixed models
#---------------------------------------------------------------------------

data <- subset(indics, indics$indic=="logfmaxpatch")
data$grps3<-as.factor(data$grps3)
mod <- lmer(data = data, value ~ grps3 + ( 1 |plotn)) 
#summary(mod)
null <- lmer(data = data, value ~ ( 1 |plotn)) 
anovtable <- as.data.frame(anova(mod,null))
anovtable
comp <- anovtable$`Pr(>Chisq)`[2]
ifelse(comp < 0.05, "p < 0.05","not significant")
comp

# for fl
data <- subset(indics, indics$indic=="flowlength")
#data$value<-log10(data$value)
data$grps3<-as.factor(data$grps3)
mod <- lm(data = data, value ~ grps3) 
summary(mod)

anova(lm(data$value ~ data$grps3))
