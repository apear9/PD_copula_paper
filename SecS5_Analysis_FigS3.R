# This script is part of the analysis of the Danish fire insurance data
# presented in Sec 5 of the main article and Sec S5 of the Supplement.

# Danish insurance data
library(fitdistrplus)
library(copula)
library(tidyverse)

# set seed
set.seed(30062025)

# Load data
data("danishmulti")

# Get subset with > 0 profit losses
dsub <- filter(danishmulti, Profits > 0) # n=616

# Create matrix with only the data we need
X1 <- with(dsub, Building + Contents)
X2 <- dsub$Profits
X <- cbind(X1, X2)
# Convert to 'copula data'
U <- pobs(X)

# Number of observations
nobs <- nrow(X)

# Kendall's tau
kest <- cor(U[,1], U[,2], method = "kendall")

# Initialise standard copula families
twcop <- tawnCopula()
gbcop <- galambosCopula()

# Fit them
twfit <- fitCopula(twcop, U, method = "itau")
gbfit <- fitCopula(gbcop, U, method = "itau")

# Define functions for calculating (pseudo)inverse of phi
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
phi_inv <- function(tt, lambda){
  if(lambda > -1){
    phi0 <- 1/(lambda + 1)
  }
  if(lambda <= -1){
    ind <- !is.infinite(tt)
  } else {
    ind <- tt < phi0 & tt != 0
  }
  vals <- rep(0, length(tt))
  vals[ind] <- findroot(tt[ind], lambda)
  vals[tt == 0] <- 1
  return(vals)
}

# Define function for simulating from PD copula
simcop <- function(n, lambda, tol = 1e-4){
  if(lambda > -1){
    phi0 <- 1/(lambda + 1)
  } else {
    phi0 <- Inf
  }
  # Independent uniform variates
  uu <- runif(n)
  tt <- runif(n)
  # Vector to store information
  vv <- numeric(n)
  # May need these values later
  if(lambda > -1) alt <- phi_inv(phi0 - phi(uu, lambda), lambda)
  # Get the inner part
  inner_inner_part <- (1 + (uu^lambda - 1)/tt)^(1/lambda)
  inner_part <- phi(inner_inner_part, lambda) - phi(uu, lambda)
  # Check for nan-s 
  inner_inner_part[is.nan(inner_inner_part)] <- -9999
  # Check for negatives
  ind <- inner_inner_part > 0
  # Calculate v for the selected values of inner_part
  vv[ind] <- phi_inv(inner_part[ind], lambda)
  # For the other parts, use alt
  if(lambda > -1) vv[!ind] <- alt[!ind]
  # Return result
  return(cbind(uu, vv))
  
}

# Define kendall's tau for the PD copula
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

# Estimate lambda for the PD copula
lambda_fitted <- uniroot(function(x) Kt(x) - kest, interval = c(-1, 0))$root

# Simulate  contours
twsim <- rCopula(nobs, twfit@copula)
gbsim <- rCopula(nobs, gbfit@copula)

# From PD copula
pdsim <- simcop(nobs, lambda_fitted)

# Contours
sims <- rbind(
  U,
  pdsim,
  twsim,
  gbsim
) %>% 
  as.data.frame %>%
  mutate(
    Copula = rep(c("Observed data", "Power-divergence", "Galambos", "Tawn"), each = nobs)
  ) %>%
  mutate(Copula = factor(Copula, ordered = T, levels = c("Observed data", "Power-divergence", "Galambos", "Tawn")))

# This produces Fig. S3
ggplot(data = sims, aes(x = X1, y = X2)) +
  geom_point(alpha = 0.5) +
  coord_equal() +
  facet_wrap(~Copula, nrow = 2) + 
  labs(x = expression(u[1]), y = expression(u[2])) +
  theme_bw() +
  theme(
    text = element_text(size = 16), 
    legend.position = "bottom", 
    axis.text.x = element_text(angle = 45, hjust= 1),
    axis.text.y = element_text(angle = 45, vjust= -1)
  )