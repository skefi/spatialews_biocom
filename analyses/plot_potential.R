# 
# Make small potential graphs to show on the pictures
# 

library(ggplot2)

# The values of aridity for which we want to draw the potential
X_VALUES <- c(0.55, 0.7, 0.85)


### For the potential in cover 

# Computing the potential requires a number of points. We take the points that fall withing this width around the X value considered. Must be of length 1
X_WIDTH <- 0.1

# Computing the potential requires computing the density, then transforming it. Here we set the bandwidth to compute the density. Play with different values to see what you like 
Y_BW <- 0.06

# The variable to compute potential on 
Y_VAR <- "imgcover"

# Load biocom data 
load(file.path(path_output,"data_biocom.rda"))

# For each x, values compute the potential 
pot_estimates <- plyr::ldply(X_VALUES, function(xref) { 
  
# Select values that fall within the range we are interested in. Here the width is divided by two because we want values on both side xref
potential_data <- subset(ourdata[["biocom"]], 
                           abs(Aridity - xref) < (X_WIDTH/2))
  
if ( nrow(potential_data) == 0 ) { 
    stop(sprintf("No data left to compute potential for %s and width %s", 
                 xref, X_WIDTH))
}
  
  # Estimate density 
  dens <- density(potential_data[ ,Y_VAR], bw = Y_BW, n = 2048)
  
  pot_result <- data.frame(dens[["x"]], dens[["y"]])
  names(pot_result) <- c(Y_VAR, "dens")
  
  # Convert to potential with the formula 
  # (Livina, Kwasniok, and Lenton 2010 Climate of the past; page 2)
  # sigma here is the "noise level". We just set it to unity as we are only interested
  # in the relative values of the potential. It does not matter as long as the noise 
  # level can be considered the same along the x axis (aridity), which is a reasonable
  # assumption. I don't think there is any way to estimate sigma from spatial data.
  sigma <- 1
  pot_result[ ,"potential"] <- - (sigma^2/2) * log(pot_result[["dens"]])
  
  data.frame(aridity = xref, pot_result)
})

# plot the potentials 
ggplot(pot_estimates) + 
  geom_line(aes_string(x = Y_VAR, y = "potential"),linewidth=1) + # to adjust for Y_VAR
  facet_wrap( ~ aridity, labeller = label_both) + 
  #theme_minimal() + # consider theme_void to remove everything
  theme_void()+
  # Consider adjusting y scale so that the potential is not as compressed
  theme(strip.text.x = element_blank())+
  coord_cartesian(ylim = c(min(pot_estimates[ ,"potential"]), 1)) +
  labs(x = Y_VAR, y = "Potential (U)")


### For the potential in MF

# Computing the potential requires a number of points. We take the points that fall 
# withing this width around the X value considered. Must be of length 1
X_WIDTH <- 0.1

# Computing the potential requires computing the density, then transforming it. Here 
# we set the bandwidth to compute the density. Play with different values to see 
# what you like 
Y_BW <- 0.18

# The variable to compute potential on 
Y_VAR <- "MF"

# For each x, values compute the potential 
pot_estimates <- plyr::ldply(X_VALUES, function(xref) { 
  
  # Select values that fall within the range we are interested in. Here the width 
  # is divided by two because we want values on both side xref
  potential_data <- subset(ourdata[["biocom"]], 
                           abs(Aridity - xref) < (X_WIDTH/2))
  
  if ( nrow(potential_data) == 0 ) { 
    stop(sprintf("No data left to compute potential for %s and width %s", 
                 xref, X_WIDTH))
  }
  
  # Estimate density 
  dens <- density(potential_data[ ,Y_VAR], bw = Y_BW, n = 2048)
  
  pot_result <- data.frame(dens[["x"]], dens[["y"]])
  names(pot_result) <- c(Y_VAR, "dens")
  
  # Convert to potential with the formula 
  # (Livina, Kwasniok, and Lenton 2010 Climate of the past; page 2)
  # sigma here is the "noise level". We just set it to unity as we are only interested
  # in the relative values of the potential. It does not matter as long as the noise 
  # level can be considered the same along the x axis (aridity), which is a reasonable
  # assumption. I don't think there is any way to estimate sigma from spatial data.
  sigma <- 1
  pot_result[ ,"potential"] <- - (sigma^2/2) * log(pot_result[["dens"]])
  
  data.frame(aridity = xref, pot_result)
})

# Plot the potentials
ggplot(pot_estimates) + 
  geom_line(aes_string(x = Y_VAR, y = "potential")) + # to adjust for Y_VAR
  facet_wrap( ~ aridity, labeller = label_both) + 
  theme_minimal() + # consider theme_void to remove everything
  # Consider adjusting y scale so that the potential is not as compressed
  coord_cartesian(ylim = c(min(pot_estimates[ ,"potential"]), 1)) + 
  labs(x = Y_VAR, y = "Potential (U)")



