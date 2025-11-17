# This script shows how to numerically calculate the power-divergence
# copula and copula density (if it exists). The script also 
# produces Fig. 1 and Fig. 2 in Section 4 of the main article.

# Packages
library(tidyverse)

# Compute copula for general lambda
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

# Prepare dataframes for plotting
df <- expand.grid(u = seq(0, 1, 0.005), v = seq(0, 1, 0.005))
lambda_ns2 <- with(df, PDC(u, v, -sqrt(2), T))
lambda_s2 <- with(df, PDC(u, v, sqrt(2), T))
df$`lambda==-sqrt(2)` <- lambda_ns2$C
df$`lambda==sqrt(2)` <- lambda_s2$C
df_long <- df %>% gather(Lambda, Copula, -c(u, v))
dfd <- expand.grid(u = seq(0, 1, 0.005), v = seq(0, 1, 0.005))
dfd$`lambda==-sqrt(2)` <- lambda_ns2$dens
dfd$`lambda==sqrt(2)` <- lambda_s2$dens
dfd_long <- dfd %>% gather(Lambda, Density, -c(u, v))

# Plot copula and 'densities' (lambda = 1 case is not strictly speaking a density)
copula_plot <- ggplot() +
  geom_raster(data = df_long, aes(x = u, y = v, fill = Copula)) +
  geom_contour(data = df_long, aes(x = u, y = v, z = Copula), col = 'black') +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(~ factor(Lambda, levels = c("lambda==-sqrt(2)", "lambda==sqrt(2)"), ordered = T), labeller = label_parsed, nrow = 1) +
  coord_equal() +
  labs(fill = expression(C[lambda](u[1],u[2])), x = expression(u[1]), y = expression(u[2])) +
  theme(text = element_text(size = 15))
density_plot <- ggplot() +
  geom_raster(data = dfd_long, aes(x = u, y = v, fill = Density^0.25)) +
  geom_contour(data = dfd_long, aes(x = u, y = v, z = Density^0.25), breaks = seq(0, 3, 0.25), col = 'black') +
  scale_fill_distiller(palette = "Spectral") +
  facet_wrap(~ factor(Lambda, levels = c("lambda==-sqrt(2)", "lambda==sqrt(2)"), ordered = T), labeller = label_parsed, nrow = 1) +
  coord_equal() +
  labs(fill = expression(bgroup("(", over(partialdiff^2*C[lambda], partialdiff*u[1]*partialdiff*u[2]), ")")^over(1,4)), x = expression(u[1]), y = expression(u[2])) +
  theme(text = element_text(size = 15))

copula_plot
density_plot