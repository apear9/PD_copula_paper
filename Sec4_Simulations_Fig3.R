# This script simulates random variates from the power-divergence copula
# using the conditional distribution method. It produces Fig. 3 in Section 4
# of the main article.

# Packages needed
library(pracma)
library(tidyverse)

# Compute copula for general lambda
phi <- function(u, lambda){
  if(lambda == 0){
    return(1 - u +  u * log(u))
  }
  if(lambda == -1){
    return(u - 1 - log(u))
  }
  return(1/(lambda * (lambda + 1)) * (u^(lambda + 1) - u + lambda * (1 - u)))
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
  # Two special cases
  if(lambda == 0){
    
    return(vals)
  }
  if(lambda == -1){
    
    return(vals)
  }
  # The more general case
  if(lambda > -1){
    phi0 <- 1/(lambda + 1)
  }
  if(lambda < -1){
    ind <- !is.infinite(tt)
  } else {
    ind <- tt < phi0 & tt != 0
  }
  vals <- rep(0, length(tt))
  vals[ind] <- findroot(tt[ind], lambda)
  vals[tt == 0] <- 1
  return(vals)
}

# Function to simulate
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

# Plot in ggplot2
set.seed(1)
nsim <- 1000
df10 <- data.frame(simcop(nsim, 10))
df3 <- data.frame(simcop(nsim, 3))
df2 <- data.frame(simcop(nsim, 2))
df1 <- data.frame(simcop(nsim, 1))
dfn05 <- data.frame(simcop(nsim, -0.5))
dfn2 <- data.frame(simcop(nsim, -2))
dfn3 <- data.frame(simcop(nsim, -3))
dfn10 <- data.frame(simcop(nsim, -10))
df_all <- bind_rows(df10, df3, df2, df1, dfn05, dfn2, dfn3, dfn10)
df_all <- df_all %>% 
  mutate(
    lambda = rep(c(10, 3, 2, 1, -0.5, -2, -3, -10), each = nsim)
  ) %>%
  mutate(
    lambda = factor(
      paste0("lambda==",lambda), 
      levels = paste0("lambda==",sort(unique(lambda))), 
      ordered = T
    )
  )

# This produces Fig. 3 in the main article
sims_plot <- ggplot(data = df_all, aes(x = uu, y = vv)) +
  geom_point(alpha = 0.2) +
  facet_wrap(~lambda, nrow = 2, labeller = label_parsed) + 
  coord_equal() +
  theme_bw() +
  labs(x = expression(u[1]), y = expression(u[2])) +
  theme(text = element_text(size = 20))
sims_plot