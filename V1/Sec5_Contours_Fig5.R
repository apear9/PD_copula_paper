# This script is a part of the analysis of the Danish fire insurance data
# in Section 5, and it produces Fig. 5.

# Danish insurance data
library(fitdistrplus)
library(copula)
library(tidyverse)

# set seed
set.seed(1072025) # date of the first ver. of script: 1 July 2025

# Load data
data("danishmulti")

# Get subset with > 0 profit losses
dsub <- filter(danishmulti, Profits > 0) # n=616

# Create matrices with only the required variables
X1 <- with(dsub, Building + Contents)
X2 <- dsub$Profits
X <- cbind(X1, X2)

# Convert data to 'copula data'
U <- pobs(X)

# Compute empirical copula
ugrd <- as.matrix(expand.grid(u=seq(0,1,0.002),v=seq(0,1,0.002)))
ecop <- copula::empCopula(U)
eval <- copula::pCopula(ugrd, ecop)
edf <- as.data.frame(cbind(ugrd, eval))
names(edf)[3] <- "Empirical"

# Kendall's tau
kest <- cor(U[,1], U[,2], method = "kendall")

# Initialise standard copula families
clcop <- rotCopula(claytonCopula()) # survival of Clayton
jocop <- joeCopula() # Joe
gacop <- normalCopula() # Gaussian
gmcop <- gumbelCopula() # Gumbel
frcop <- frankCopula() # Frank
hrcop <- huslerReissCopula() # Husler-Reiss
twcop <- tawnCopula() # Tawn
glcop <- galambosCopula() # Galambos

# Fit them
clfit <- fitCopula(clcop, U, method = "itau")
gmfit <- fitCopula(gmcop, U, method = "itau")
frfit <- fitCopula(frcop, U, method = "itau")
jofit <- fitCopula(jocop, U, method = "itau")
gafit <- fitCopula(gacop, U, method = "itau")
hrfit <- fitCopula(hrcop, U, method = "itau")
twfit <- fitCopula(twcop, U, method = "itau")
glfit <- fitCopula(glcop, U, method = "itau")

# Define functions to fit PD copula
phi <- function(u, lambda){
  1/(lambda * (lambda + 1)) * (u^(lambda + 1) - u + lambda * (1 - u))
}
rootfunc <- function(x, tt, lambda){
  x^(lambda + 1) - (lambda + 1) * x + lambda - lambda * (lambda + 1) * tt
}
findroot <- function(tt, lambda){
  
  vals <- rep(0, length(tt))
  for(i in 1:length(tt)){
    if(tt[i] == 0){
      vals[i] <- 1
    } else {
      vals[i] <- uniroot(rootfunc, c(0, 1), lambda = lambda, tt = tt[i])$root
    }

  }
  return(vals)
  
}
PDC <- function(u, v, lambda, density = F, maxtol = 1e6){
  
  if(lambda > -1){
    phi0 <- 1/(lambda + 1)
  } else {
    phi0 <- Inf
  }
  tt <- phi(u, lambda) + phi(v, lambda)
  ind1 <- tt <= phi0
  ind2 <- tt < maxtol
  vals <- rep(0, length(u))
  vals[ind1 & ind2] <- findroot(tt[ind1 & ind2], lambda)
  vals[!ind2] <- 0
  out <- list(C = vals)
  if(density == T){
    dens <- -lambda * vals ^ (lambda - 1) * (u^lambda - 1) * (v^lambda - 1) / (vals ^ lambda - 1)^3
    dens[!ind1] <- NA
    out$dens <- dens
  }
  return(out)
  
}
Kt <- function(lambda){
  
  if(lambda == 0){
    return(3 - 4 * log(2))
  } 
  if(lambda == -1){
    return(7 - 2 * pi ^2 / 3)
  }
  integral_part <- integrate(function(x) (1 - x)/(x^lambda - 1), 0, 1)$value
  1 + 2/(lambda + 1) + 4 * lambda / (lambda + 1) * integral_part
  
}

# Estimate lambda
lambda_fitted <- uniroot(function(x) Kt(x) - kest, interval = c(-1, 0))$root

# Compute PD copula output
pdval <- PDC(ugrd[,1], ugrd[,2], lambda_fitted, T)

# Examine contours
glval <- pCopula(ugrd, glfit@copula)
twval <- pCopula(ugrd, twfit@copula)

# Contours
edf_contour <- edf %>%
  mutate(
    `Power-divergence` = pdval$C,
    Galambos = glval, 
    Tawn = twval
  )
edf_contour_long <- edf_contour %>%
  gather(Copula, Value, -c(u, v)) %>%
  mutate(Copula = factor(Copula, levels = c("Empirical", "Power-divergence", "Galambos", "Tawn"), ordered = T))

# Plot
contours_plot <- ggplot() +
  geom_contour(
    data = filter(edf_contour_long, Copula == "Empirical"), aes(x = u, y = v, z = Value),
    col = 'black',
    breaks = seq(.05, 1-.05, .1),
    linewidth = 1.2
  ) +
  geom_contour(
    data = filter(edf_contour_long, Copula != "Empirical"), aes(x = u, y = v, z = Value, col = Copula),
    breaks = seq(.05, 1-.05, .1),
    linewidth = .9
  ) +
  scale_colour_manual(values = c("Power-divergence" = "red", "Galambos" = "darkgreen", "Tawn" = "blue")) +
  coord_equal() +
  labs(x = expression(u[1]), y = expression(u[2]), title = "(a) Contours", colour = "") +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = "bottom")

# Similar thing for the copula density of the CR copula
dens_df <- edf %>%
  mutate(Density = pdval$dens)

dens_plot <- ggplot() +
  geom_raster(data = dens_df, aes(x = u, y = v, fill = Density^0.25)) +
  scale_fill_distiller(palette = "Spectral") +
  geom_point(data = data.frame(U), aes(x = X1, y = X2), alpha = 0.25) +
  theme_bw() +
  coord_equal() +
  labs(
    title = "(b) Density",
    x = expression(u[1]), 
    y = expression(u[2]), 
    fill = expression(c[-0.647](u[1],u[2])^(1/4))
  ) +
  theme(legend.position = "bottom", text = element_text(size = 16))

contours_plot
dens_plot
