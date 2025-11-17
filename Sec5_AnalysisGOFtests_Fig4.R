# This script is required for the analysis of Danish insurance data
# presented in Section 5 of the main article.
# It performs goodness-of-fit testing for all copulas except the 
# power-divergence copula.

# Packages
library(fitdistrplus)
library(copula)
library(tidyverse)

# set seed (date that I wrote the first version of this script)
set.seed(19062025)

# Load data
data("danishmulti")

# Get subset with > 0 profit losses
dsub <- filter(danishmulti, Profits > 0) # n=616

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

# This code makes Fig. 4 (main paper)
ggplot(data = alldf, aes(x = TotalMaterial, y = Profits)) +
  geom_point(alpha = 0.5, size = 2) +
  theme_bw() +
  labs(x = "Building + Contents", y = "Profits") +
  facet_wrap(~Units, nrow = 1, scales = "free") +
  theme(text = element_text(size = 16)) +
  NULL

# Calculate the empirical copula for the Danish fire insurance data
ugrd <- as.matrix(expand.grid(u=seq(0,1,0.002),v=seq(0,1,0.002)))
ecop <- copula::empCopula(U)
eval <- copula::pCopula(ugrd, ecop)
edf <- as.data.frame(cbind(ugrd, eval))
names(edf)[3] <- "Empirical"

# Get the sample version of Kendall's tau
kest <- cor(U[,1], U[,2], method = "kendall")

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
clfit <- fitCopula(clcop, U, method = "itau")
gmfit <- fitCopula(gmcop, U, method = "itau")
frfit <- fitCopula(frcop, U, method = "itau")
jofit <- fitCopula(jocop, U, method = "itau")
gafit <- fitCopula(gacop, U, method = "itau")
hrfit <- fitCopula(hrcop, U, method = "itau")
twfit <- fitCopula(twcop, U, method = "itau")
glfit <- fitCopula(glcop, U, method = "itau")

# Goodness of fit testing
set.seed(123)
clgof <- gofCopula(clcop, U, estim.method = "itau")
clgof
gmgof <- gofCopula(gmcop, U, estim.method = "itau")
gmgof
frgof <- gofCopula(frcop, U, estim.method = "itau")
frgof
jogof <- gofCopula(jocop, U, estim.method = "itau")
jogof
gagof <- gofCopula(gacop, U, estim.method = "itau")
gagof
hrgof <- gofEVCopula(hrcop, U, method = "itau")
hrgof
twgof <- gofEVCopula(twcop, U, method = "itau")
twgof
glgof <- gofEVCopula(glcop, U, method = "itau")
glgof
