# Plot PSD along gradients
# Fig 4  

library(tidyr) # for gather
library(spatialwarnings)
library(dplyr)

#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
NPERM = 199 #value use for final analyses: 199
BOOTN = 2999 #value use for final analyses: 2999
path_output <- here::here("outputs")
ALPHA <- 0.05 # Significance level for indicator trends

# These constants are here to determine how TPL functions are being fitted.
# Fitting a TPL requires computing sum_{k=1}^{inf} k^(-a) exp(-b k). This
# sum is infinite so we need some stopping criterion, here either
# a relative change in the sum below a threshold (reltol below), or a
# maximum number of iterations (1e8). These are fixed values in
# spatialwarnings <= 3.0.3, and can be set via global options in versions
# above (still in dev on 2023-07-27) as we do below.
# Decrease reltol and increase maxit below to increase precision of
# computations, at the cost of higher computation time.
options(spatialwarnings.constants.reltol = 1e-6)
options(spatialwarnings.constants.maxit = 1e8)


#---------------------------------------------------------------------------
# Function needed 
#---------------------------------------------------------------------------

n_random_psds  <- function(mat, nsim = 5) {
  
  #psd_real <- patchsizes(to_logical(mat))
  psd_real <- patchsizes(mat)
  
  mats_randoms <- replicate(
    nsim,
    expr = {
      #to_logical(spatialwarnings:::shuffle_matrix(mat))
      matrix(as.logical(spatialwarnings:::shuffle_matrix(mat)),nrow=dim(mat)[1])
      #(spatialwarnings:::shuffle_matrix(mat))
    },
    simplify = FALSE)
  
  psd_randoms <- lapply(
    mats_randoms,
    patchsizes
  )
  
  cumpsd_randoms <- lapply(
    psd_randoms,
    function(psd) {
      return(spatialwarnings:::cumpsd(psd))
    }
  )
  
  cumpsd_randoms <- lapply(
    seq_along(cumpsd_randoms),
    function(i, dflist) {
      dflist[[i]] %>%
        mutate(matn = i)
    },
    dflist = cumpsd_randoms
  )
  
  out  <- do.call(bind_rows, cumpsd_randoms) %>%
    mutate(psdtype = "rand")
  
  return(out)
  
}


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Load data
#--------------------------------------------------------------------------#--------------------------------------------------------------------------

load(here::here("outputs", "biocom-grps.rda"))
load(file.path(path_output,"data_biocom.rda"))
arid <- biocom

filename <- paste0("indics-data50-grps_Nperm_", NPERM,"_rev.rda")
load(file.path(path_output,filename))


#---------------------------------------------------------------------------
# Short code version (after having chosen the sites)
#---------------------------------------------------------------------------

## Sites kept: healthier: 148-b for grasslands 
## degraded: 192-c for grasslands 

indics$veg_type <- as.factor(indics$veg_type)

### Upper branch 
#--------------

# We want sites on the upper branch that are grassland (veg_type==1) vs shrublands (veg_type==2), i.e. grps2 = 2

## Grassland
indics_high_g <- subset(indics,grps2==2 & veg_type==1)
indics_high_g <- subset(indics_high_g,indic=="slope" | indic=="cutoff")
indics_high_g <- indics_high_g[,c("Aridity","indic","value","file")]
indics_high_g <- spread(indics_high_g,indic,value)

# subset of sites
chosen_h_g <- subset(indics_high_g,file=="148-b") 
pics_h_g <- ourdata$matrices[[chosen_h_g$file]]
#display_matrix(pics_h_g,palette="Paired")

fig.img1 <- display_matrix(pics_h_g,palette="Greys")+
  theme(legend.position = "none")+
  scale_fill_grey(start = 1, end = 0)+
  theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) 
  
mat <- pics_h_g
psd_h_g <- patchdistr_sews(mat, best_by = "AIC", fit_lnorm = FALSE)
plot_psd_h_g <- plot_distr(psd_h_g,best_only = TRUE,plrange=FALSE) +
                          theme(legend.position = "none",
                                text = element_text(size=12)) + 
                          scale_x_continuous(trans = "log10", limits = c(1, 10000), breaks = 10^seq(log10(1), log10(10000), by = 1))

CORES <- min(7, parallel::detectCores()-1)
plan(multisession, workers = CORES) # turn on parallel computation
random <- n_random_psds(mat, nsim = 10)
plan(sequential) # turn off parallel computation

plot_psdn_h_g   <- plot_psd_h_g +
  ggtitle("healthier")+
  geom_line(data = random, aes(x = patchsize, y = y, group = matn), color = "lightgrey",show.legend = FALSE)


### Lower branch 
#--------------

## Grassland
indics_low_g <- subset(indics,grps2==1 & veg_type==1)
indics_low_g <- subset(indics_low_g,indic=="slope" | indic=="cutoff")
indics_low_g <- indics_low_g[,c("Aridity","indic","value","file")]
indics_low_g <- spread(indics_low_g,indic,value)

# subset of sites
chosen_l_g <- subset(indics_low_g,file=="192-c") 
pics_l_g <- ourdata$matrices[[chosen_l_g$file]]
fig.img2 <-   display_matrix(pics_l_g,palette="Greys")+
  theme(legend.position = "none")+
  scale_fill_grey(start = 1, end = 0)+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 

mat2 <- pics_l_g
psd_l_g <- patchdistr_sews(mat2, best_by = "AIC", fit_lnorm = FALSE)
plot_psd_l_g <- plot_distr(psd_l_g,best_only = TRUE,plrange=FALSE) +
  theme(legend.position = "none",
        text = element_text(size=12)) + 
  scale_x_continuous(trans = "log10", limits = c(1, 10000), breaks = 10^seq(log10(1), log10(10000), by = 1))

CORES <- min(7, parallel::detectCores()-1)
plan(multisession, workers = CORES) # turn on parallel computation
random2 <- n_random_psds(mat2, nsim = 10)
plan(sequential) # turn off parallel computation

plot_psdn_l_g   <- plot_psd_l_g +
  ggtitle("degraded")+
  geom_line(data = random2, aes(x = patchsize, y = y, group = matn), color = "lightgrey",show.legend = FALSE)

#plot_grid(plot_psdn_h_g, plot_psdn_l_g, ncol=2, nrow=1) 
#pics_final <- ourdata$matrices[c("148-b","192-c")]
#display_matrix(pics_final,palette="Paired")


blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

first_row <- plot_grid(blankPlot, fig.img1, blankPlot, fig.img2, labels = c('A', '', 'B',""), ncol=4, nrow=1) 

second_row <- plot_grid(plot_psdn_h_g, plot_psdn_l_g, ncol=2, nrow=1) 


#### Fig 4
grid.arrange(first_row, second_row, ncol = 1,nrow=2, heights = c(0.6,1))
#ggsave("./figures/fig4_PSD_grass-high-low.pdf", width = 9.5, height = 5)

chosen_final_h <- subset(indics,file=="148-b")
#Aridity 0.645
#cover  6.813029e-01  
#cutoff  1.490116e-08  
#fmaxpatch  5.498201e-01  
#perco  1.000000e+00  
#slope  1.468552e+00  

chosen_final_l <- subset(indics,file=="192-c") 
# for 192-c
#Aridity 0.904
#cover  0.233972711  
#cutoff  0.007929005  
#fmaxpatch  0.007512017  
#perco  0.000000000  
#slope  0.433799393  



