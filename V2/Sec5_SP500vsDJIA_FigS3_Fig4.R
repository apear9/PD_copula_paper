# Functions needed to evaluate PD copula, simulate from it, etc.
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

# Packages
library(tidyverse)
library(copula)
library(forecast)
library(xtable)

# Load data
SP500 <- read_csv("SP500.csv") %>% filter(observation_date >= "2022-10-01" & observation_date <= "2026-09-30")
DJIA <- read_csv("DJIA.csv") %>% filter(observation_date >= "2022-10-01" & observation_date <= "2026-09-30")
IDX <- left_join(SP500, DJIA, by = "observation_date")

# Check NAs
IDX <- na.omit(IDX)

# Fit ARIMA models and get residuals
SP500_model <- auto.arima(IDX$SP500, stepwise = F)
DJIA_model <- auto.arima(IDX$DJIA, stepwise = F)

# ACF/PACF plots for residuals (Fig S3 of Supplementary Material)
par(mfrow=c(2, 2), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
acf(SP500_model$residuals, main = "ACF: S&P 500 residuals")
pacf(SP500_model$residuals, main = "PACF: S&P 500 residuals")
acf(DJIA_model$residuals, main = "ACF: DJIA residuals")
pacf(DJIA_model$residuals, main = "PACF: DJIA residuals")
par(mfrow=c(1,1))

# Now check the copula data and fit copula modls
XY <- cbind(SP500_model$residuals, DJIA_model$residuals)
UV <- pobs(XY) 

# Plot time series and copula data based on ARIMA residuals (Fig 4 of main text)
par(mfrow=c(1, 3), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
plot(ts(IDX$DJIA), xlab = "Days", ylab = "DJIA", main = "(a) DJIA")
plot(ts(IDX$SP500), xlab = "Days", ylab = "S&P500", main = "(b) S&P500")
plot(UV, asp = 1, xlab = "S&P500", ylab = "DJIA", main = "(c) Residuals: DJIA vs S&P500"); rect(0, 0, 1, 1, lty = 2)
par(mfrow=c(1, 1))

# Fit models
clcop <- claytonCopula() # Clayton
jocop <- joeCopula() # Joe
gacop <- normalCopula() # Gaussian
gmcop <- gumbelCopula() # Gumbel
frcop <- frankCopula() # Frank
hrcop <- huslerReissCopula() # Husler-Reiss
glcop <- galambosCopula() # Galambos

# Fit copulas
clfit <- fitCopula(clcop, UV, method = "mpl")
scfit <- fitCopula(rotCopula(clcop), UV, method = "mpl")
jofit <- fitCopula(jocop, UV, method = "mpl")
sjfit <- fitCopula(rotCopula(jocop), UV, method = "mpl")
gmfit <- fitCopula(gmcop, UV, method = "mpl")
sgfit <- fitCopula(rotCopula(gmcop), UV, method = "mpl")
frfit <- fitCopula(frcop, UV, method = "mpl")
hrfit <- fitCopula(hrcop, UV, method = "mpl")
glfit <- fitCopula(glcop, UV, method = "mpl")
gafit <- fitCopula(gacop, UV, method = "mpl")

# Fitting PD model
PDfit <- fit_lambda(UV, -1.5, maxtol = 1e16, lower = -10, upper = -1)

# Compile parameter estimates and log-likelihood values
pms <- c(clfit@estimate, scfit@estimate, frfit@estimate,
         gmfit@estimate, sgfit@estimate, jofit@estimate, 
         sjfit@estimate, PDfit$par,
         glfit@estimate, hrfit@estimate, 
         gafit@estimate)
lls <- c(clfit@loglik, scfit@loglik, frfit@loglik,
         gmfit@loglik, sgfit@loglik, jofit@loglik, 
         sjfit@loglik, PDfit$val,
         glfit@loglik, hrfit@loglik, 
         gafit@loglik)
to_print <- data.frame(
  Model = c("Clayton", "S-Clayton", "Frank", "Gumbel", "S-Gumbel", "Joe", "S-Joe", "Power-divergence", "Galambos", "Husler-Reiss", "Gaussian"), 
  Parms = pms,
  LL = lls
)

# Print a table (rough version of Table 1)
print(xtable(to_print)) 