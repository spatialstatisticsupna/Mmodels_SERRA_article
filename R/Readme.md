# R code

This folder contains the necessary R functions to fit the multivariate spatio-temporal models described in Vicente et al. (2020) and reproduce results.

The ```dataMmodel.RData``` file contains the following R objects:

- ```data```: contains the data set used. It is a dataframe with the following variables,
	- **dist**: Districts
	- **state**: Satte (Uttar Pradesh)
	- **year**: Year (2001:2014)
	- **rape**: Observed number of rapes
	- **dowry**: Observed number of dowry deaths
	- **pop**: Female population between 15 and 49 years (linear interpolation)
	- **e_rape**: Expected number of rapes
	- **e_dowry**: Expected number of dowry deaths
	- **smr_rape**: Standardized incidence ratio (SMR) of rapes
	- **smr_dowry**: Standardized mortality ratio (SMR) of dowry deaths
	- **ID_area**: Area Identifiers (Districts)
	- **ID_year**: Time Identifiers (Year)
	- **ID_area_year**: Area-time Identifiers (Districts-Year)


- ```carto_india```: SpatialPolygonDataFrame object with the cartography of the 33 states (without islands) of India

- ```carto_up```: SpatialPolygonDataFrame object with the cartography of the 70 districts (year 2001) of Uttar Pradesh

- ```Mzz```: is the upper triangular matrix of the Cholesky decomposition of the inverse of the temporal neighborhood matrix


The file [reproduce_paper_Mmodels.R](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/reproduce_paper_Mmodels.R) permit to reproduce the results given in the paper.

The [functions](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/functions) folder contains the necessary functions to fit M-models using INLA and WinBUGS.

The [run_Mmodels_inla.R](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/run_Mmodels_inla.R) files allow you to adjust the fixed/random-effect M-models using INLA.

The [run_winbugs_bym_fe.R](https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/run_winbugs_bym_fe.R) files allow you to adjust the fixed-effect M-models with a BYM model for the spatial random effect, using MCMC.
Similarly for the other models considered in the paper.

 
