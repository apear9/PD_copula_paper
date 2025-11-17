# This script produces the example of the piecewise phi function and its
# inverse as seen in Section S2 of the Supplement. This script also produces
# Fig. S1.

# Packages
library(tidyverse)
library(gridExtra)

# Functions
phipw <- function(x){
  
  (0.25 - exp(2) + exp(1/x)) * (x < 0.5) + (x-1)^2 * (x >= 0.5)
  
}
phipwinv <- function(t){
  1/(log(t - 0.25 + exp(2))) * (0.25 <= t) + (1 - sqrt(t)) * (t < 0.25) + 0
}

# Create data frame for plotting
phidf <- data.frame(x = seq(0.4, 1, 0.01)) %>%
  mutate(phi = phipw(x))
phinvdf <- data.frame(t = seq(0, 4, 0.01)) %>%
  mutate(phinv = phipwinv(t))
generator <- ggplot(data = phidf, aes(x = x, y = phi)) +
  geom_line() +
  theme_bw() +
  labs(x = "x", y = expression(phi[pw](x)), title = "(a) Piecewise generator") +
  theme(text = element_text(size = 16))
inverse <- ggplot(data = phinvdf, aes(x = t, y = phinv)) +
  geom_line() +
  theme_bw() +
  labs(x = "t", y = expression(phi[pw]^-1*(t)), title = "(b) Piecewise inverse") +
  theme(text = element_text(size = 16))

# Print plots
generator
inverse