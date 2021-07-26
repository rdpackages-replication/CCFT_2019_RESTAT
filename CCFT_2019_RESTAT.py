################################################################################
## Regression Discontinuity Designs Using Covariates
## Author: Calonico, Cattaneo, Farrell and Titiunik
## Python code created by Ricardo Masini
## Last update: 26-JUL-2021
################################################################################
## WEBSITE: https://rdpackages.github.io/
## RDROBUST: pip install rdrobust
################################################################################
## NOTE: if you are using rdrobust version 2020 or newer, the option 
## masspoints="off" and stdvars=True may be needed in order to replicate the 
## results in the paper.
## For example, line 52:
##    tmp = rdrobust(y, x, c=59.1968)
## should be replaced by:
##    tmp = rdrobust(y, x, c=59.1968, masspoints="off", stdvars=True)
################################################################################

from rdrobust import rdrobust
import pandas  as pd

# Load Data
df = pd.read_csv('headstart.csv')

# Variable selection
y = df['mort_age59_related_postHS']
x = df['povrate60']
z = df[['census1960_pop', 
       'census1960_pctsch1417',
       'census1960_pctsch534',
       'census1960_pctsch25plus',
       'census1960_pop1417',
       'census1960_pop534',
       'census1960_pop25plus', 
       'census1960_pcturban', 
       'census1960_pctblack']]

out = pd.DataFrame(columns = ["Standard", "Cov-adjusted", "Cov-adjusted"],
                   index = pd.Index(["RD treatment effect",
                                      "Robust 95% CI",
                                      "CI length change (%)",
                                      "Robust p-value",
                                      "Robust 95% CI",
                                      "CI length change (%)",
                                      "Robust p-value",
                                      "h|b",
                                      "n-|n+"]))

## rho unrestricted; MSE-optimal bandwidths w/o covs; RD w/o covs
rd = rdrobust(y, x, c=59.1968)
h  = rd.bws['left'][0]
b  = rd.bws['left'][1]
IL = rd.ci.loc['Robust'][1] - rd.ci.loc['Robust'][0]
out.iat[0,0] = round(rd.coef.loc['Conventional'].values[0],3)
out.iat[1,0] = "[" + str(round(rd.ci.loc['Robust'][0],3)) + "," + str(round(rd.ci.loc['Robust'][1],3)) + "]"
out.iat[2,0] = ""
out.iat[3,0] = round(rd.pv.iloc[2].values[0],3)

covs_drop = True

## rho unrestricted; MSE-optimal bandwidths w/o covs; RD w/ covs
rd = rdrobust(y, x, c=59.1968, h=h, b=b, covs = z, covs_drop = covs_drop)
ILch = ((rd.ci.loc['Robust'][1] - rd.ci.loc['Robust'][0])/IL - 1)* 100
out.iat[0,1] = round(rd.coef.loc['Conventional'].values[0],3)
out.iat[1,1] = "[" + str(round(rd.ci.loc['Robust'][0],3)) + "," + str(round(rd.ci.loc['Robust'][1],3)) + "]"
out.iat[2,1] = round(ILch,3)
out.iat[3,1] = round(rd.pv.iloc[2].values[0],3)

## rho unrestricted; MSE-optimal bandwidths w/ covs; RD w/ covs
rd = rdrobust(y, x, c=59.1968, covs=z, covs_drop = covs_drop)
ILch = ((rd.ci.loc['Robust'][1] - rd.ci.loc['Robust'][0])/IL - 1)* 100
out.iat[0,2] = round(rd.coef.loc['Conventional'].values[0],3)
out.iat[1,2] = "[" + str(round(rd.ci.loc['Robust'][0],3)) + "," + str(round(rd.ci.loc['Robust'][1],3)) + "]"
out.iat[2,2] = round(ILch,3)
out.iat[3,2] = round(rd.pv.iloc[2].values[0],3)

## rho=1; MSE-optimal bandwidths w/o covs; RD w/o covs
rd = rdrobust (y, x, c=59.1968, rho=1)
h  = rd.bws['left'][0]
b  = rd.bws['left'][1]
IL = rd.ci.loc['Robust'][1] - rd.ci.loc['Robust'][0]
out.iat[4,0] =  "[" + str(round(rd.ci.loc['Robust'][0],3)) + "," + str(round(rd.ci.loc['Robust'][1],3)) + "]"
out.iat[5,0] = ""
out.iat[6,0] = round(rd.pv.iloc[2].values[0],3)
out.iat[7,0] = str(round(h,3)) + "|" + str(round(b,3))
out.iat[8,0] = str(round(rd.N_h[0],3)) + "|" + str(round(rd.N_h[1],3))

## rho=1; MSE-optimal bandwidths w/o covs; RD w/ covs
rd = rdrobust (y, x, c=59.1968, covs=z, h=h, b=b,covs_drop = covs_drop)
ILch = ((rd.ci.loc['Robust'][1] - rd.ci.loc['Robust'][0])/IL - 1)* 100
out.iat[4,1] = "[" + str(round(rd.ci.loc['Robust'][0],3)) + "," + str(round(rd.ci.loc['Robust'][1],3)) + "]"
out.iat[5,1] = round(ILch,3)
out.iat[6,1] = round(rd.pv.iloc[2].values[0],3)
out.iat[7,1] = str(round(h,3)) + "|" + str(round(b,3))
out.iat[8,1] = str(round(rd.N_h[0],3)) + "|" + str(round(rd.N_h[1],3))

## rho=1; MSE-optimal bandwidths w/ covs; RD w/ covs
rd = rdrobust (y, x, c=59.1968, covs=z, rho=1, covs_drop = covs_drop)
h  = rd.bws['left'][0]
b  = rd.bws['left'][1]
ILch = ((rd.ci.loc['Robust'][1] - rd.ci.loc['Robust'][0])/IL - 1)* 100
out.iat[4,2] =  "[" + str(round(rd.ci.loc['Robust'][0],3)) + "," + str(round(rd.ci.loc['Robust'][1],3)) + "]"
out.iat[5,2] = round(ILch,3)
out.iat[6,2] = round(rd.pv.iloc[2].values[0],3)
out.iat[7,2] = str(round(h,3)) + "|" + str(round(b,3))
out.iat[8,2] = str(round(rd.N_h[0],3)) + "|" + str(round(rd.N_h[1],3))

print(out)


