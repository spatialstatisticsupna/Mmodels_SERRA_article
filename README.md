# Bayesian inference in multivariate spatio-temporal areal models using INLA: analysis of gender-based violence in small areas
This repository contains the R code to fit the models described in the paper entitled _"Bayesian inference in multivariate spatio-temporal areal models using INLA: analysis of gender-based violence in small areas"_ [(Vicente et al., 2020)](https://doi.org/10.1007/s00477-020-01808-x).

## Table of contents
- [R code](#R-code)

- [References](#References)

# R code
An updated version of the R code to fit most of the multivariate spatio-temporal models described in the paper (compatible with R-INLA 21.11.22 or newer versions) can be found [here](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/). In addition, a new reparameterization of the covariance matrix of the M-models based on its Barlett decomposition is defined for models fitted with INLA. 

The original code to reproduce the results of the paper using R-4.0.4 and the R-INLA version 21.02.23 is also available at: [Version v1](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/releases/tag/v1).

# Acknowledgements
This work has been supported by the Spanish Ministry of Economy, Industry, and Competitiveness (project MTM2017-82553-R, AEI/FEDER, UE), and partially funded by la Caixa Foundation (ID 1000010434), Caja Navarra Foundation and UNED Pamplona, under agreement LCF/PR/PR15/51100007.

# References
Vicente, G., Goicoa, T., and Ugarte, M.D. (2020). Bayesian inference in multivariate spatio-temporal areal models using INLA: analysis of gender-based violence in small areas. _Stochastic Environmental Research and Risk Assessment_, vol. 34, 1421-1440, 2020, [DOI:10.1007/s00477-020-01808-x](https://doi.org/10.1007/s00477-020-01808-x).
