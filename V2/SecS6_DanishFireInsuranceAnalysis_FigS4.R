# This script is required for the analysis of Danish insurance data
# presented in Section S6 of the Supplement.

# Packages
library(fitdistrplus)
library(copula)
library(tidyverse)

# set seed (date that I wrote the first version of this script)
set.seed(19062025)

# Load data
data("danishmulti")

# Get subset with > 0 profit losses
dsub <- danishmulti %>% filter(Profits > 0)

# Create matrix containing only the relevant data
X1 <- with(dsub, Building + Contents)
X2 <- dsub$Profits
X <- cbind(X1, X2)

# Convert the data to 'copula data'
U <- pobs(X)

# Prepare data for plotting
Udf <- data.frame(U) %>%
  mutate(Units = "Copula data")
Xdf <- data.frame(X) %>% 
  mutate(Units = "Original data (MDK)")
alldf <- bind_rows(Udf, Xdf) %>%
  mutate(Units = factor(Units, levels = c("Original data (MDK)", "Copula data"), ordered = T)) %>%
  rename(TotalMaterial = X1, Profits = X2)

# This code makes Fig. S4
ggplot(data = alldf, aes(x = TotalMaterial, y = Profits)) +
  geom_point(alpha = 0.5, size = 2) +
  theme_bw() +
  labs(x = "Building + Contents", y = "Profits") +
  facet_wrap(~Units, nrow = 1, scales = "free") +
  theme(text = element_text(size = 16)) +
  NULL

# Initialise standard copula families
clcop <- rotCopula(claytonCopula()) # Survival copula of Clayton copula
jocop <- joeCopula() # Joe copula
gacop <- normalCopula() # Gaussian copula
gmcop <- gumbelCopula() # Gumbel copula
frcop <- frankCopula() # Frank copula
hrcop <- huslerReissCopula() # Husler-Reiss copula
twcop <- tawnCopula() # Tawn copula
glcop <- galambosCopula() # Galambos copula

# Fit all copulas to the data
clfit <- fitCopula(clcop, U, method = "mpl") # 145.6, param = 1.202
gmfit <- fitCopula(gmcop, U, method = "mpl") # 127.9, param = 1.603
frfit <- fitCopula(frcop, U, method = "mpl") # 94.9, param = 3.683
jofit <- fitCopula(jocop, U, method = "mpl") # 147.6, param = 2.047
gafit <- fitCopula(gacop, U, method = "mpl") # 84.24, param = 0.4946
hrfit <- fitCopula(hrcop, U, method = "mpl") # 120.8, param = 1.293
twfit <- fitCopula(twcop, U, method = "mpl") # 128.1, param = 0.8969
glfit <- fitCopula(glcop, U, method = "mpl") # 125.4, param = 0.883

# Fit PD copula by ml
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
PDC <- function(u, v, lambda, density = F, maxtol = 1e6){
  
  if(lambda > -1){
    phi0 <- 1/(lambda + 1)
  } else {
    phi0 <- Inf
  }
  tt <- phi(u, lambda) + phi(v, lambda)
  ind1 <- tt < phi0
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
# Function to simulate from PD copula
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
# ML functions
llik_PD <- function(x, X, ...){
  
  lambda <- x
  uu <- X[,1]
  vv <- X[,2]
  out <- PDC(uu, vv, lambda, T, ...)
  llik <- sum(log(out$dens))
  return(llik)
  
}
fit_lambda <- function(X, init, ...){
  ML_out <- optim(init, llik_PD, X = X, method = "Brent", ..., control = list(fnscale = -1))
  lambda_fitted <- list(
    con = ML_out$con, 
    par = ML_out$par,
    val = ML_out$value
  )
  return(lambda_fitted)
}

fit_lambda(U, -0.8, lower = -0.99, upper = -0.6) # 155.6789, param = -0.789