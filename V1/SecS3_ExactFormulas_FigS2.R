# This script generates Fig. S2 (more or less) of the Supplement
# (The plots were arranged together after the fact in an image-editing software)

# Packages
library(tidyverse)
library(gridExtra)

# Functions to evaluate the copula and copula density when lambda = -2
cop <- function(u, v){
  
  tt <- u^-1 + u + v^-1 + v - 4
  val <- 1 + 0.5 * (tt) - sqrt(0.25 * (tt)^2 + u^-1 + u + v^-1 + v - 4)
  return(val)
  
}
copd <- function(u, v){
  C <- cop(u, v)
  2 * C^-3 * (u^-2 - 1) * (v^-2 - 1) / (C^-2 - 1)^3
}

# Construct a data frame with the outputs
df <- expand.grid(u = seq(0, 1, 0.002), v = seq(0, 1, 0.002))
df$C <- with(df, cop(u,v))
df$c <- with(df, copd(u,v))

# Make the plots for lambda = -2
cop_plot_n2 <- ggplot() +
  geom_raster(data = df, aes(x = u, y = v, fill = C)) +
  geom_contour(data = df, aes(x = u, y = v, z = C), col = "black") +
  scale_fill_distiller(palette = "Spectral") + 
  coord_equal() +
  labs(fill = expression(C[-2](u[1],u[2])), x = expression(u[1]), y = expression(u[2])) +
  theme( text= element_text(size = 15)) 
dens_plot_n2 <- ggplot() +
  geom_raster(data = df, aes(x = u, y = v, fill = (c)^(1/4))) +
  geom_contour(data = df, aes(x = u, y = v, z = (c)^(1/4)), col = "black") +
  scale_fill_distiller(palette = "Spectral", breaks = seq(0, 4.5, 0.25)) +
  coord_equal() +
  labs(fill = expression(bgroup("(", over(partialdiff^2*C[-2], partialdiff*u[1]*partialdiff*u[2]), ")")^over(1,4)), x = expression(u[1]), y = expression(u[2])) +
  theme( text = element_text(size = 15))

# Functions to evaluate copula and copula density when lambda = -0.5
cop <- function(u, v){
  
  dom <- 2 - 2  * (sqrt(u) + sqrt(v)) + u + v
  
  ind <- 1 * ((dom >= 0) & (dom < 1))
  
  val <- ind * (3 - 2 * sqrt(dom) - 2 * (sqrt(u) + sqrt(v)) + u + v)
  
  return(val)
  
}
copd <- function(u, v){
  
  dom <- 2 - 2  * (sqrt(u) + sqrt(v)) + u + v
  
  ind <- 1 * ((dom >= 0) & (dom < 1))
  
  val <- ind * (1 - 1/sqrt(u)) * (1 - 1/sqrt(v))/(2 - 2 * (sqrt(u) + sqrt(v)) + u + v)^(3/2)
  val[!ind] <- NA
    
  return(val)
  
}

# Create a data frame with the outputs
df <- expand.grid(u = seq(0, 1, 0.002), v = seq(0, 1, 0.002))
df$C <- with(df, cop(u,v))
df$c <- with(df, copd(u,v))

# Plot the copula and copula density
cop_plot_n0.5 <- ggplot() +
  geom_raster(data = df, aes(x = u, y = v, fill = C)) +
  geom_contour(data = df, aes(x = u, y = v, z = C), col = "black") +
  scale_fill_distiller(palette = "Spectral") + 
  coord_equal() +
  labs(fill = expression(C[-0.5](u[1],u[2])), x = expression(u[1]), y = expression(u[2])) +
  theme( text= element_text(size = 15))
dens_plot_n0.5 <- ggplot() +
  geom_raster(data = df, aes(x = u, y = v, fill = (c)^(1/4))) +
  geom_contour(data = df, aes(x = u, y = v, z = (c)^(1/4)), col = "black", breaks = seq(0, 4.5, 0.25)) +
  scale_fill_distiller(palette = "Spectral") + 
  coord_equal() +
  labs(fill = expression(bgroup("(", over(partialdiff^2*C[-0.5], partialdiff*u[1]*partialdiff*u[2]), ")")^over(1,4)), x = expression(u[1]), y = expression(u[2])) +
  theme( text = element_text(size = 15))

# Functions to evaluate copula and "copula density" when lambda = 1
cop <- function(u, v){
  
  val <- pmax(0, 1 - sqrt(u^2 - 2 * u + v^2 - 2 * v + 2))
  
  return(val)
  
}
copd <- function(u, v){
  
  dom <- u^2 - 2 * u + v^2 - 2 * v + 2
  
  ind <- 1 * ((dom >= 0) & (dom < 1))
  
  val <- ind * (u-1) * (v-1) / (u^2 - 2 * u + v^2 - 2 * v + 2)^(3/2)
  val[!ind] <- NA
  
  return(val)
  
}

# Create a data frame with the outputs
df1 <- expand.grid(u = seq(0, 1, 0.002), v = seq(0, 1, 0.002))
df1$C <- with(df1, cop(u,v))
df1$c <- with(df1, copd(u,v))

# Plot copula and "copula density" 
cop_plot_1 <- ggplot() +
  geom_raster(data = df1, aes(x = u, y = v, fill = C)) +
  geom_contour(data = df1, aes(x = u, y = v, z = C), col = "black") +
  scale_fill_distiller(palette = "Spectral") + 
  coord_equal() +
  labs(fill = expression(C[1](u[1],u[2])*"   "), x = expression(u[1]), y = expression(u[2])) +
  theme( text = element_text(size = 15))
dens_plot_1 <- ggplot() +
  geom_raster(data = df1, aes(x = u, y = v, fill = (c)^(1/4))) +
  geom_contour(data = df1, aes(x = u, y = v, z = (c)^(1/4)), col = "black", breaks = seq(0, 4.5, 0.25)) +
  scale_fill_distiller(palette = "Spectral") + 
  coord_equal() +
  labs(fill = expression(bgroup("(", over(partialdiff^2*C[1], partialdiff*u[1]*partialdiff*u[2]), ")")^over(1,4)), x = expression(u[1]), y = expression(u[2])) +
  theme( text = element_text(size = 15))

# Print the plots
cop_plot_n2
cop_plot_n0.5
cop_plot_1
dens_plot_n2
dens_plot_n0.5
dens_plot_1
