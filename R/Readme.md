# R code

This folder contains the necessary R functions to fit to fit most of the multivariate spatio-temporal models described in Vicente et al. (2020) and reproduce results.

The original code [(Version v1)](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/releases/tag/v1) has been updated to make it compatible with R-INLA 21.11.22 or newer versions.

In addition, the Barlett decomposition of the between-crime spatial and temporal covariance matrices is used for models fitted with INLA.

The ```dataMmodel.RData``` file contains the following R objects:

- ```data_UP```: contains the data set used. It is a list with two `data.frame` objects (corresponding to **rape** and **dowry death** data) with the following variables,
	- **dist**: Districts
	- **state**: State (Uttar Pradesh)
	- **year**: Years (from 2001 to 2014)
	- **pop**: Female population between 15 and 49 years (linear interpolation)
	- **obs**: Observed number of cases (rapes or dowry deaths)
	- **exp**: Expected number of cases (rapes or dowry deaths)
	- **smr**: Standardized incidence/mortality ratio

- ```carto_india```: `sf` object with the cartography of the 33 states (without islands) of India

- ```carto_UP```: `sf` object with the cartography of the 70 districts (year 2001) of Uttar Pradesh


The file [reproduce_paper_Mmodels.R](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/reproduce_paper_Mmodels.R) permits to reproduce similar results to the ones given in the paper.

The [functions](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/functions) folder contains the necessary functions to fit M-models using INLA and WinBUGS.

The file [run_Mmodels_inla.R](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/run_Mmodels_inla.R) allows to fit fixed-effect M-models using INLA.

The rest of the files allows to fit the fixed/random-effect M-models using WinBUGS.
