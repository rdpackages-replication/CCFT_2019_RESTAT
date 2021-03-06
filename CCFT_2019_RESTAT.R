################################################################################
## Regression Discontinuity Designs Using Covariates
## Author: Calonico, Cattaneo, Farrell and Titiunik
## Last update: 07-JUL-2021
################################################################################
## WEBSITE: https://rdpackages.github.io/
## RDROBUST: install.packages(rdrobust)
################################################################################
## NOTE: if you are using rdrobust version 2020 or newer, the option 
## masspoints="off" and stdvars="on" may be needed in order to replicate the 
## results in the paper.
## For example, line 35:
##    tmp = rdrobust(y, x, c=59.1968)
## should be replaced by:
##    tmp = rdrobust(y, x, c=59.1968, masspoints="off", stdvars="on")
################################################################################

rm(list=ls(all=TRUE))
library(rdrobust)

################################################################################
  ## Head Start Data
################################################################################
data <- read.csv("headstart.csv")
attach(data)

y <- mort_age59_related_postHS
x <- povrate60
z <- cbind(census1960_pop, census1960_pctsch1417, census1960_pctsch534,
  census1960_pctsch25plus, census1960_pop1417, census1960_pop534,
  census1960_pop25plus, census1960_pcturban, census1960_pctblack)
out <- matrix(NA,9,3)

## rho unrestricted; MSE-optimal bandwidths w/o covs; RD w/o covs
rd <- rdrobust(y, x, c=59.1968)
h  <- rd$bws[1,1]
b  <- rd$bws[2,1]
IL <- rd$ci[3,2] - rd$ci[3,1]

out[1,1] <- round(rd$coef[1],3)
out[2,1] <- paste("[", round(rd$ci[3,1],3),",", round(rd$ci[3,2],3), "]", sep="")
out[3,1] <- ""
out[4,1] <- round(rd$pv[3],3)

## rho unrestricted; MSE-optimal bandwidths w/o covs; RD w/ covs
rd <- rdrobust(y, x, c=59.1968, covs=z, h=h, b=b)
ILch <- ((rd$ci[3,2] - rd$ci[3,1])/IL - 1)* 100
out[1,2] <- round(rd$coef[1],3)
out[2,2] <- paste("[", round(rd$ci[3,1],3),",", round(rd$ci[3,2],3), "]", sep="")
out[3,2] <- round(ILch,3)
out[4,2] <- round(rd$pv[3],3)

## rho unrestricted; MSE-optimal bandwidths w/ covs; RD w/ covs
rd <- rdrobust(y, x, c=59.1968, covs=z)
ILch <- ((rd$ci[3,2] - rd$ci[3,1])/IL - 1)* 100
out[1,3] <- round(rd$coef[1],3)
out[2,3] <- paste("[", round(rd$ci[3,1],3),",", round(rd$ci[3,2],3), "]", sep="")
out[3,3] <- round(ILch,3)
out[4,3] <- round(rd$pv[3],3)

## rho=1; MSE-optimal bandwidths w/o covs; RD w/o covs
rd <- rdrobust (y, x, c=59.1968, rho=1)
h  <- rd$bws[1,1]
b  <- rd$bws[2,1]
IL <- rd$ci[3,2] - rd$ci[3,1]
out[5,1] <- paste("[", round(rd$ci[3,1],3),",", round(rd$ci[3,2],3), "]", sep="")
out[6,1] <- ""
out[7,1] <- round(rd$pv[3],3)
out[8,1] <- paste(round(h,3),"|",round(b,3))
out[9,1] <- paste(round(rd$N_h[1],3),"|",round(rd$N_h[2],3))

## rho=1; MSE-optimal bandwidths w/o covs; RD w/ covs
rd <- rdrobust (y, x, c=59.1968, covs=z, h=h, b=b)
ILch <- ((rd$ci[3,2] - rd$ci[3,1])/IL - 1)* 100
out[5,2] <- paste("[", round(rd$ci[3,1],3),",", round(rd$ci[3,2],3), "]", sep="")
out[6,2] <- round(ILch,3)
out[7,2] <- round(rd$pv[3],3)
out[8,2] <- paste(round(h,3),"|",round(b,3))
out[9,2] <- paste(round(rd$N_h[1],3),"|",round(rd$N_h[2],3))

## rho=1; MSE-optimal bandwidths w/ covs; RD w/ covs
rd <- rdrobust (y, x, c=59.1968, covs=z, rho=1)
ILch <- ((rd$ci[3,2] - rd$ci[3,1])/IL - 1)* 100
out[5,3] <- paste("[", round(rd$ci[3,1],3),",", round(rd$ci[3,2],3), "]", sep="")
out[6,3] <- round(ILch,3)
out[7,3] <- round(rd$pv[3],3)
out[8,3] <- paste(round(rd$bws[1,1],3),"|",round(rd$bws[2,1],3))
out[9,3] <- paste(round(rd$N_h[1],3),"|",round(rd$N_h[2],3))


rownames(out) <- c("RD treatment effect", "Robust 95% CI", "CI length change (%)", "Robust p-value", 
                   "Robust 95% CI", "CI length change (%)", "Robust p-value", "h|b", "n-|n+")

colnames(out) <- c("Standard", "Cov-adjusted", "Cov-adjusted")

