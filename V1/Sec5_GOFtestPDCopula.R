# This script performs the parametric-bootstrap goodness-of-fit test
# for the PD copula using the algorithm given in Appendix A of 
# Genest et al. (2009). Goodness-of-fit tests for copulas: A review and power 
# study. Insurance: Mathematics and Economics, 44, 199-213.
# I have annotated the code to show how each code block
# corresponds to the steps outlined in Appendix A of Genest et al. (2009).
# This script takes about 8 minutes on my computer, which is just a laptop
# with 2Ghz processor, running R 4.4.2. 

# Packages
library(copula)
library(fitdistrplus)
library(tidyverse)

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
# Kendall's tau
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

# Set seed 
set.seed(29062025)  # date of first ver. of script: 29 Jun 2025

# Load and process data
data("danishmulti")
dsub <- filter(danishmulti, Profits > 0) # n=616

# Create matrix with only the relevant data
X1 <- with(dsub, Building + Contents)
X2 <- dsub$Profits
X <- cbind(X1, X2)

# Compute copula data
U <- pobs(X)

# Compute number of samples
n <- nrow(X)

## Step (1) in Appendix A of Genest et al. (2009). Insur: Math & Econ, 44, 199-213.

# Estimate empirical copula
ecop <- empCopula(U)

# Estimate parameter of PD copula (use Kendall's tau)
kest <- cor(U[,1], U[,2], method = "kendall")
lambda_fitted <- uniroot(function(x) Kt(x) - kest, interval = c(-1, 0))$root
lambda_fitted

## Step (2): Get test statistic
Sn <- sum((pCopula(U, ecop)-PDC(U[,1], U[,2], lambda = lambda_fitted, F)$C)^2)

## Step (3): Loop over k = 1, ..., N
N <- 5e3
Sstars <- numeric(N)
t1 <- Sys.time()
for(k in 1:N){
  
  if(k %% 500 == 0){
    t2 <- Sys.time()
    print(t2-t1)
    t1 <- Sys.time()
  }
  
  # (a) Get random sample from copula and compute rank vectors
  Ystarnk <- simcop(n, lambda_fitted)
  Rstarnk <- pobs(Ystarnk)
  
  # (b) Adjust the ranks, get empirical copula, and estimate lambda
  Ustarnk <- Rstarnk # the pobs() function computes rank/(n+1) in Point (3)(b)
  Cstarnk <- empCopula(Ustarnk)
  Ktaunk <- cor(Ustarnk[,1], Ustarnk[,2], method = "kendall")
  lambdank <- uniroot(function(x) Kt(x) - Ktaunk, interval = c(-1.5, 0.5))$root
  
  # (c) Estimate Sstarnk
  Cthetastarnk <- PDC(Ustarnk[,1], Ustarnk[,2], lambda = lambdank, F)$C
  Sstars[k] <-sum((pCopula(Ustarnk, Cstarnk)-Cthetastarnk)^2)
  
}

# Compute approx p value
pval <- mean(Sstars > Sn)
Sn # statistic
pval # p-value

# # Uncomment to save output
# save.image("./PD-GOF.Rdata")