# Plot PSD along gradients
# Fig S2

library(tidyr) # for gather

#source(here::here("R", "functions_helper.R"))

#---------------------------------------------------------------------------
# Parameters needed 
#---------------------------------------------------------------------------
# From _targets.R
NPERM = 200 #value use for final analyses: 199
BOOTN = 2999 #value use for final analyses: 2999
path_output <- here::here("outputs")
ALPHA <- 0.05 # Significance level for indicator trends

#indic_order <- c('plrange', 'logfmaxpatch', 'slope','cutoff', 'flowlength','moran', 'sdr','cv.variance') 


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# Load model data
#--------------------------------------------------------------------------#--------------------------------------------------------------------------


load(here::here("outputs", "biocom-grps.rda"))
load(file.path(path_output,"data_biocom.rda"))
arid <- biocom

#filename <- paste0("indics-data-grps_Nperm_", NPERM,".rda")
# For revision:
filename <- paste0("indics-data50-grps_Nperm_", NPERM,"_rev.rda")
load(file.path(path_output,filename))

#filename_m <- paste0("indics-model_Nperm_", NPERM,".rda")
# For revision:
filename_m = paste0("indics-model_Nperm_",NPERM,"_rev.rda")
load(file.path(path_output,filename_m))




#-------------------------------------------------------------------------------
### Build the data set
#-------------------------------------------------------------------------------
location_data <- '/Users/soniakefi/Dropbox/MyProjects/2022_Scaling/TPBmodel/Matlab_2011/data/facilitation_f09'

# Go to the folder where the data are stored
# Produce a character vector of all the names of the data files in the selected
# directory (214 entries corresponding to all the landscapes of the data set)

### data set 1
#datalist <- dir('/Users/soniakefi/Dropbox/MyProjects/2019_EWS-Biocom/tutorial/data_CA_full_2011', full.names = TRUE, pattern = "*.txt") 
## Read the files - results in a data frame of 214 elements
#matrices <- lapply(datalist, read.table) 
#matrices <- lapply(matrices, as.matrix) 
#matrices <- lapply(matrices, function(x) { x == 1 }) 
#files <- dir("/Users/soniakefi/Dropbox/MyProjects/2019_EWS-Biocom/tutorial/data_CA_full_2011") 
#files2 <- gsub(".txt", "", files)
#files3 <- gsub("CA_PotAna_", "", files2)
#names(matrices) <- files3

### data set 2
datalist <- dir(location_data, full.names=TRUE, pattern="*.txt")

# Read the files - results in a data frame of 214 elements
matrices <- lapply(datalist, read.table) 

# Convert the data frame into matrices
# The resulting matrices contain values: 1 (vegetation), 2 (empty but
# recolonizable) and 3 (degraded)
matrices <- lapply(matrices, as.matrix) 
#str(matrices)

# Convert to logicals (so that the package knows that there are only 2 values
# possible in the matrices); true means 1 (vegetation), false means 0 (degraded
# or recolonizable)
matrices <- lapply(matrices, function(x) { x == 1 }) 

# Produce a character vector of the names of files in the current directory
files <- dir(location_data) 
files2 <- gsub(".txt", "", files)
files3 <- gsub("CA_b", "", files2)
names(matrices) <- files3

# Extract parameter b from the file name, which is the stress level
b <- as.numeric(files3)
aridity <- 1 - b

# Start creating the result file
new.parameters <- data.frame(b = b, aridity = aridity)
matrices <- matrices[order(files3)] 

ourdata <- list(new.parameters = new.parameters, matrices = matrices)




#display_matrix(mat,along = b,palette="Paired")



#-------------------------------------------------------------------------------
### PSD
#-------------------------------------------------------------------------------

#' Returns 1 if pl better, 0 if tpl better, NA if impossible fit
#'
#' @param mat 
#'
#' @return best fit
#' @export 
#'
#' @examples best_pl_or_tpl(mat)

best_pl_or_tpl <- function(mat) {
  res <- spatialwarnings::patchdistr_sews(mat, best_by = "AIC", fit_lnorm = TRUE)
  res <- as.data.frame(res)
  if (nrow(res) == 1) {
    return(as.numeric(NA))
  }
  aic_pl <- subset(res, type == "pl")[ ,"AIC"]
  aic_tpl <- subset(res, type == "tpl")[ ,"AIC"]
  as.numeric(aic_pl < aic_tpl)
}


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

########### PL

mat = matrices[[16]]
new.parameters$b[16]
display_matrix(mat,palette="Paired")

psd_indic <- patchdistr_sews(mat)
plot_psd <- plot_distr(psd_indic,best_only = TRUE,plrange=FALSE)

#****
CORES <- min(7, parallel::detectCores()-1)
plan(multisession, workers = CORES) # turn on parallel computation
#random <- n_random_psds(to_logical(mat), nsim = 10)
random <- n_random_psds(mat, nsim = 10)
plan(sequential) # turn off parallel computation
#****

plot_psdn   <- plot_psd +
  #scale_x_log10(limits = c(1,1e4))+
  geom_line(data = random, aes(x = patchsize, y = y, group = matn), color = "lightgrey",show.legend = FALSE)
#theme_minimal()+
#xlim(log10(1),log10(1000))

plot_psdn

cover <- mean(mat)
psd <- patchsizes(mat)
sizemat <- dim(mat)[1]*dim(mat)[2]
logmaxpatch = log10(max(psd)) 
logfmaxpatch = log10(max(psd)/sizemat)
maxpatch = (max(psd))
fmaxpatch = (max(psd)/sizemat)
# /!\ slope and cutoff are always positive. But their negative 
# values are used in the fit, e.g. y ~ x^{-slope}, or 
# y ~ x^{-slope]e^{-cutoff}
slope    = psd_indic$psd_type$plexpo[2]
cutoff   = psd_indic$psd_type$cutoff[2]
slope_pl    = psd_indic$psd_type$plexpo[1]
#ks_dist  = get_ks_distance(mat), 
plr = raw_plrange(mat)
#?spatialwarnings::percolation
perco <- percolation(mat)
#pl_or_tpl = best_pl_or_tpl(mat)  # Returns 1 if pl


########### TPL

mat = matrices[[10]]
new.parameters$b[10]

psd_indic2 <- patchdistr_sews(mat)
plot_psd2 <- plot_distr(psd_indic2,best_only = TRUE,plrange=FALSE)

#****
CORES <- min(7, parallel::detectCores()-1)
plan(multisession, workers = CORES) # turn on parallel computation
#random <- n_random_psds(to_logical(mat), nsim = 10)
random2 <- n_random_psds(mat, nsim = 10)
plan(sequential) # turn off parallel computation
#****

plot_psd2n   <- plot_psd2 +
  #scale_x_log10(limits = c(1,1e4))+
  geom_line(data = random2, aes(x = patchsize, y = y, group = matn), color = "lightgrey")
#theme_minimal()+
#xlim(log10(1),log10(1000))

plot_psd2n

cover <- mean(mat)
slope    = psd_indic2$psd_type$plexpo[2]
cutoff   = psd_indic2$psd_type$cutoff[2]
slope_pl    = psd_indic2$psd_type$plexpo[1]

########### TPL

mat = matrices[[4]]
new.parameters$b[4]

psd_indic3 <- patchdistr_sews(mat)
plot_psd3 <- plot_distr(psd_indic3,best_only = TRUE,plrange=FALSE)

plot_psd4 <- plot_distr(psd_indic3,best_only = FALSE,plrange=FALSE)


#****
CORES <- min(7, parallel::detectCores()-1)
plan(multisession, workers = CORES) # turn on parallel computation
#random <- n_random_psds(to_logical(mat), nsim = 10)
random3 <- n_random_psds(mat, nsim = 10)
plan(sequential) # turn off parallel computation
#****

plot_psd3n   <- plot_psd3 +
  #scale_x_log10(limits = c(1,1e4))+
  geom_line(data = random3, aes(x = patchsize, y = y, group = matn), color = "lightgrey")
#theme_minimal()+
#xlim(log10(1),log10(1000))

plot_psd4n   <- plot_psd4 +
  #scale_x_log10(limits = c(1,1e4))+
  geom_line(data = random3, aes(x = patchsize, y = y, group = matn), color = "lightgrey")+
  geom_point(aes(x=39.4,y=0))

plot_psd3n

cover <- mean(mat)
sizemat <- dim(mat)[1]*dim(mat)[2]
psd <- patchsizes(mat)
sizemat <- dim(mat)[1]*dim(mat)[2]
logmaxpatch = log10(max(psd)) 
logfmaxpatch = log10(max(psd)/sizemat)
maxpatch = (max(psd))
fmaxpatch = (max(psd)/sizemat)
slope    = psd_indic3$psd_type$plexpo[2]
cutoff   = psd_indic3$psd_type$cutoff[2]
slope_pl    = psd_indic3$psd_type$plexpo[1]

plot_grid(plot_psdn, plot_psd2n, plot_psd3n, labels= c("A","B","C"),ncol=3, nrow=1)




#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------#--------------------------------------------------------------------------









# extract needed results
indics_df.sub <- indics_df[indics_df$indic=="slope_pl",]
res.temp$slope_pl <- indics_df.sub$value
res.temp$slope_pl_null <- indics_df.sub$null_mean
res.temp$slope_pl_diff <- indics_df.sub$diff

indics_df.sub <- indics_df[indics_df$indic=="slope_tpl",]
res.temp$slope_tpl <- indics_df.sub$value
res.temp$slope_tpl_null <- indics_df.sub$null_mean
res.temp$slope_tpl_diff <- indics_df.sub$diff

indics_df.sub <- indics_df[indics_df$indic=="cutoff",]
res.temp$cutoff <- indics_df.sub$value
res.temp$cutoff_null <- indics_df.sub$null_mean
res.temp$cutoff_diff <- indics_df.sub$diff

indics_df.sub <- indics_df[indics_df$indic=="pl_or_tpl",]
res.temp$best <- indics_df.sub$value
res.temp$best_null <- indics_df.sub$null_mean
res.temp$best_diff <- indics_df.sub$diff

indics_df.sub <- indics_df[indics_df$indic=="perco",]
res.temp$perco <- indics_df.sub$value
res.temp$perco_null <- indics_df.sub$null_mean
res.temp$perco_diff <- indics_df.sub$diff

# merge with rest of results
res <- rbind(res,res.temp)
rm(res.temp) # clear res.temp


