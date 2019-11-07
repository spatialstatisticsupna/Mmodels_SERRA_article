# Multivariate spatio-temporal areal models: a comparison of INLA and WinBus for the analysis of gender-based violence in small areas
This repository contains the R code to fit the models described in the paper entitled _"Multivariate spatio-temporal areal models: a comparison of INLA and WinBus for the analysis of gender-based violence in small areas"_ (Vicente et al., 2019).

## Table of contents
- [Figures-Tables](#Figures-Tables)

- [R code](#R-code)

- [References](#References)

# Figures-Tables

## Figures

- [Figure 1.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_1.pdf): Map of the administrative division of Uttar Pradesh into districts and its location in India (top right corner)

- [Figure 2.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_2.pdf): Evolution of the crude rates (per 100000 women) of rapes and dowry deaths in Uttar Pradesh in the period 2001-2014

- [Figure 3.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_3.pdf): Dispersion plots of the final relative risks for rapes and dowry deaths obtained with the Type II interaction RE M-model with in INLA (y-axis) vs. WinBUGS (x-axis), using the iCAR (first row), pCAR (second row), LCAR (third row) and the BYM (last row) spatial priors

- [Figure 4.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_4.pdf): Posterior mean of the district-specific spatial risk, $\exp(\theta_{ij})$ (left column), and the exceedence probabilities, i.e., P(exp(theta_ij)>1|O) (right column), for rape (top) and dowry deaths (bottom)

- [Figure 5.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_5.pdf): Temporal pattern of incidence risks (posterior means of exp(gamma_tj) ) for rape  and dowry deaths in Uttar Pradesh

- [Figure 6.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_6.pdf): Map of estimated incidence risks for rape (top) and  posterior probabilities that the relative risk is greater than one (P(R_itj >1|O) ) (bottom) in Uttar Pradesh

- [Figure 7.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_7.pdf): Map of estimated incidence risks for dowry deaths (top) and  posterior probabilities that the relative risk is greater than one (P(R_itj >1|O) ) (bottom) in Uttar Pradesh 

- [Figure 8.](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/figures/figure_8.pdf): Temporal evolution of final risk estimates for rape and dowry deaths in some districts in Uttar Pradesh: Ghazlabad, Kheri, Mainpuri, Sant Kabir Nagar, and Varanasi

## Tables

- **Table 1.** Descriptive statistics. Minimum (min), first quartile (q1), mean, third quartile (q3), maximum (max), standard desviation (sd), and coefficient of variation (cv) of the number of rapes and dowry deaths in the districts of Uttar Pradesh per year

|	   |     |    | Rapes |   |     |    |    |     |    |Dowry deaths | | | | |
| Year | min | q1 | mean | q3 | max | sd | cv | min | q1 | mean | q3 | max | sd | cv |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2001 |  1 | 13.0 | 27.9 | 41.0 | 93 | 21.5 | 0.8 |  4 | 18.0 | 31.6 | 43.8 | 88 | 19.0 | 0.6 |
| 2002 |  0 | 9.0 | 20.2 | 30.8 | 73 | 14.7 | 0.7 |  3 | 14.2 | 27.0 | 34.8 | 83 | 18.1 | 0.7 |
| 2003 |  0 | 5.0 | 13.0 | 19.5 | 47 | 11.1 | 0.9 |  3 | 10.2 | 18.9 | 24.0 | 55 | 11.5 | 0.6 |
| 2004 |  3 | 9.2 | 19.9 | 25.8 | 72 | 15.0 | 0.8 |  3 | 14.2 | 24.4 | 29.0 | 71 | 15.3 | 0.6 |
| 2005 |  1 | 7.0 | 17.3 | 24.0 | 61 | 14.2 | 0.8 |  1 | 12.2 | 22.3 | 26.8 | 70 | 13.9 | 0.6 |
| 2006 |  2 | 9.0 | 18.8 | 26.0 | 51 | 12.1 | 0.6 |  7 | 14.2 | 25.7 | 34.8 | 67 | 14.4 | 0.6 |
| 2007 |  1 | 10.0 | 23.5 | 32.5 | 82 | 16.6 | 0.7 |  4 | 16.0 | 29.6 | 36.8 | 78 | 17.3 | 0.6 |
| 2008 |  2 | 12.0 | 26.7 | 35.8 | 82 | 19.0 | 0.7 |  5 | 17.2 | 32.0 | 38.8 | 88 | 18.7 | 0.6 |
| 2009 |  3 | 13.0 | 25.1 | 35.2 | 77 | 17.5 | 0.7 |  8 | 19.2 | 31.9 | 40.8 | 83 | 18.0 | 0.6 |
| 2010 |  1 | 10.2 | 21.9 | 26.0 | 75 | 17.4 | 0.8 |  5 | 18.2 | 31.4 | 40.0 | 95 | 19.7 | 0.6 |
| 2011 |  2 | 14.2 | 29.1 | 39.0 | 89 | 20.6 | 0.7 |  6 | 17.0 | 33.2 | 41.8 | 95 | 18.7 | 0.6 |
| 2012 |  4 | 15.0 | 28.0 | 35.8 | 86 | 17.4 | 0.6 |  5 | 19.0 | 32.0 | 40.8 | 97 | 17.9 | 0.6 |
| 2013 |  5 | 23.2 | 43.5 | 53.8 | 119 | 28.5 | 0.7 |  5 | 19.0 | 33.3 | 41.2 | 98 | 19.5 | 0.6 |
| 2014 |  5 | 23.0 | 49.5 | 69.0 | 164 | 31.7 | 0.6 |  6 | 23.2 | 35.3 | 46.8 | 98 | 18.4 | 0.5 |


- **Table 2.** Correlations between spatial (by year) and temporal patterns (by district) of rape and dowry deaths

| Correlation | min | q.25 | median | mean | q.75 | max | sd | cv |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| spatial patterns | 0.319 | 0.371 | 0.449 | 0.449 | 0.538 | 0.621 | 0.099 | 0.220 |
| temporal trends | -0.369 | 0.142 | 0.396 | 0.378 | 0.630 | 0.865 | 0.300 | 0.793 |


- **Table 3.** Model selection criteria, DIC, WAIC and LS, for different models

| Theta |				| Type 		|  DIC | WAIC | LS 	 |
| :--- 	| :--- 			| :--- 		| ---: | ---: | ---: |
| iCAR 	| FE M-models	| Additive	| 14160.929 | 14413.494 | 7210.957 |
|   	|  				| Type I 	| 12608.167 | 12522.138 | 6608.222 |
| 		|  				| Type II 	| 12355.856 | 12379.212 | 6338.481 |
| 		|  				| Type III 	| 12663.282 | 12707.279 | 6624.992 |
| 		|  				| Type IV 	| 12405.457 | 12479.757 | 6370.050 |
|		| RE M-models 	| Additive 	| 14161.084 | 14413.314 | 7210.853 |
| 		|  				| Type I 	| 12607.161 | 12521.936 | 6607.393 |
| 		|  				| Type II 	| 12356.652 | 12387.969 | 6338.562 |
| 		|  				| Type III 	| 12661.840 | 12710.519 | 6623.163 |
| 		|  				| Type IV 	| 12403.473 | 12472.729 | 6369.541 |
| pCAR 	| FE M-models 	| Additive 	| 14161.376 | 14415.178 | 7211.860 |
| 		|  				| Type I 	| 12606.321 | 12507.739 | 6607.746 |
| 		|  				| Type II 	| 12356.132 | 12373.483 | 6338.431 |
| 		|  				| Type III 	| 12660.066 | 12693.335 | 6622.705 |
| 		|  				| Type IV 	| 12403.443 | 12476.556 | 6369.834 |
| 		| RE M-models 	| Additive 	| 14161.125 | 14414.223 | 7211.355 |
| 		|  				| Type I 	| 12607.587 | 12522.777 | 6608.033 |
| 		|  				| Type II 	| 12362.337 | 12399.167 | 6342.436 |
| 		|  				| Type III 	| 12660.802 | 12699.595 | 6622.373 |
| 		|  				| Type IV 	| 12393.019 | 12441.440 | 6365.348 |
| LCAR 	| FE M-models 	| Additive 	| 14160.912 | 14413.994 | 7211.237 |
| 		|  				| Type I 	| 12608.498 | 12529.463 | 6609.148 |
| 		|  				| Type II 	| 12358.354 | 12392.021 | 6339.574 |
| 		|  				| Type III 	| 12663.113 | 12715.302 | 6623.690 |
| 		|  				| Type IV 	| 12396.433 | 12455.568 | 6366.125 |
| RE M-models			| Additive 	| 14159.901 | 14412.886 | 7210.686 |
| 		|  				| Type I 	| 12609.721 | 12522.181 | 6609.431 |
| 		|  				| Type II 	| 12354.983 | 12374.041 | 6337.927 |
| 		|  				| Type III 	| 12657.593 | 12696.315 | 6621.781 |
| 		|  				| Type IV 	| 12404.071 | 12479.757 | 6369.211 |
| BYM 	| FE M-models 	| Additive 	| 14160.500 | 14413.932 | 7211.210 |
| 		|  				| Type I 	| 12608.295 | 12541.938 | 6609.303 |
| 		|  				| Type II 	| 12353.168 | 12375.668 | 6337.332 |
| 		|  				| Type III 	| 12664.510 | 12722.878 | 6623.983 |
| 		|  				| Type IV 	| 12400.066 | 12463.970 | 6368.017 |
| 		| RE M-models 	| Additive 	| 14161.078 | 14413.632 | 7211.020 |
| 		|  				| Type I 	| 12607.490 | 12522.248 | 6607.727 |
| 		|  				| Type II 	| 12354.795 | 12380.401 | 6337.246 |
| 		|  				| Type III 	| 12663.906 | 12707.357 | 6625.187 |
| 		|  				| Type IV 	| 12402.443 | 12473.332 | 6368.702 |


## Appendix

- **Table A.1.** District identifiers (ID) of Uttar Pradesh

| ID 	| Dist 			| ID 	| Dist 					| ID 	| Dist  |
| ---: 	| :--- 			| ---: 	| :--- 					| ---: 	| :---  |
|1		| Agra 			| 25 	| Fatehpur 				| 49 	| Mainpuri |
|2 		| Aligarh 		| 26 	| Firozabad 			| 50 	| Mathura |
|3 		| Allahabad 	| 27 	| Gautam Buddha Nagar	| 51 	| Mau |
|4 		| Ambedkar Nagar| 28 	| Ghaziabad 			| 52 	| Meerut |
|5 		| Auraiya 		| 29 	| Ghazipur 				| 53 	| Mirzapur |
|6 		| Azamgarh 		| 30 	| Gonda 				| 54 	| Moradabad |
|7 		| Baghpat 		| 31 	| Gorakhpur 			| 55 	| Muzaffarnagar |
|8 		| Bahraich 		| 32 	| Hamirpur 				| 56 	| Pilibhit |
|9 		| Ballia 		| 33 	| Hardoi 				| 57 	| Pratapgarh |
|10 	| Balrampur 	| 34 	| Hathras				| 58 	| Rae Bareli |
|11 	| Banda 		| 35 	| Jalaun 				| 59 	| Rampur |
|12 	| Barabanki 	| 36 	| Jaunpur 				| 60 	| Saharanpur |
|13 	| Bareilly 		| 37 	| Jhansi 				| 61 	| Sant Kabir Nagar |
|14 	| Basti 		| 38 	| Jyotiba Phule Nagar 	| 62 	| Sant Ravidas Nagar Bhadohi |
|15 	| Bijnor 		| 39 	| Kannauj 				| 63 	| Shahjahanpur |
|16 	| Budaun 		| 40 	| Kanpur Dehat 			| 64 	| Shrawasti |
|17 	| Bulandshahr 	| 41 	| Kanpur Nagar 			| 65 	| Siddharthnagar |
|18 	| Chandauli 	| 42 	| Kaushambi 			| 66 	| Sitapur |
|19 	| Chitrakoot 	| 43 	| Kheri 				| 67 	| Sonbhadra |
|20 	| Deoria 		| 44 	| Kushinagar 			| 68 	| Sultanpur |
|21 	| Etah 			| 45 	| Lalitpur 				| 69 	| Unnao |
|22 	| Etawah 		| 46 	| Lucknow 				| 70 	| Varanasi |
|23 	| Faizabad		| 47 	| Mahoba 				|  		|  |
|24 	| Farrukhabad 	| 48 	| Mahrajganj 			|  		|  |


- **Table A.2.** Posterior means, standard deviations, and 95\% credible intervals for the crime-specific intercepts (alpha_j, j=1,2) of the models with a spatio-temporal Type II interaction term

|		|				|		|	 	 |    | FE 	  |  	  | 	 | 	  | RE    |  |
|		|				|		|	mean | sd | q.025 | q.975 | mean | sd | q.025 | q.975 |
| :---	| :---			| :---	| ---: 	|---: |---: |---: |---: |---: |---: |---: |
| iCAR 	| Rape 			| MCMC  | -0.187 | 0.008 | -0.203 | -0.170 | -0.186 | 0.008 | -0.201 | -0.169 |
|      	|      			| INLA  | -0.186 | 0.008 | -0.203 | -0.170 | -0.185 | 0.008 | -0.202 | -0.169  |
|      	| Dowry deaths 	| MCMC 	| -0.061 | 0.007 | -0.075 | -0.047 | -0.061 | 0.007 | -0.076 | -0.046 |
|      	|              	| INLA 	| -0.061 | 0.007 | -0.075 | -0.047 | -0.061 | 0.007 | -0.074 | -0.047 |
| pCAR 	| Rape 			| MCMC  | -0.186 | 0.008 | -0.203 | -0.171 | -0.186 | 0.009 | -0.204 | -0.169 |
|      	|      			| INLA 	| -0.187 | 0.008 | -0.203 | -0.170 | -0.185 | 0.008 | -0.202 | -0.169 |
|      	| Dowry deaths 	| MCMC 	| -0.061 | 0.007 | -0.074 | -0.048 | -0.061 | 0.007 | -0.075 | -0.048 |
|      	| 				| INLA 	| -0.061 | 0.007 | -0.075 | -0.048 | -0.062 | 0.007 | -0.076 | -0.049 |
| LCAR 	| Rape 			| MCMC 	| -0.187 | 0.008 | -0.203 | -0.170 | -0.186 | 0.009 | -0.203 | -0.169 |
|      	|      			| INLA 	| -0.185 | 0.008 | -0.201 | -0.169 | -0.186 | 0.008 | -0.203 | -0.170 |
|      	| Dowry deaths 	| MCMC 	| -0.061 | 0.007 | -0.074 | -0.048 | -0.061 | 0.007 | -0.075 | -0.047 |
|      	| 				| INLA 	| -0.060 | 0.007 | -0.074 | -0.047 | -0.061 | 0.007 | -0.075 | -0.048 |
| BYM  	| Rape 			| MCMC 	| -0.187 | 0.008 | -0.203 | -0.172 | -0.186 | 0.008 | -0.203 | -0.170 |
|      	| 				| INLA 	| -0.187 | 0.008 | -0.204 | -0.171 | -0.185 | 0.008 | -0.201 | -0.169 |
|      	| Dowry deaths	| MCMC 	| -0.061 | 0.007 | -0.074 | -0.048 | -0.061 | 0.007 | -0.075 | -0.048 |
|      	| 				| INLA 	| -0.062 | 0.007 | -0.075 | -0.048 | -0.061 | 0.007 | -0.075 | -0.048 |


- **Table A.3.** Posterior means, standard deviations, and 95\% credible intervals for the hyperparameters of the models with a spatio-temporal Type II interaction term
 
 |   |       |           |      | INLA  |      |      | MCMC  |       |
|   | Model | Parameter | mean | q.025 | q.975| mean | q.025 | q.975 |
| :---	| :---		  | :---	      | ---:  |---:   |---:   |  ---: |---:   |---:   |
| iCAR  | FE M-models | sigma_delta_1 | 0.212 | 0.187 | 0.233 | 0.210 | 0.190 | 0.232 |
|       | 	          | sigma_delta_2 | 0.093 | 0.080 | 0.107 | 0.093 | 0.080 | 0.108 |
|       | RE M-models | sigma_theta   | 0.605 | 0.268 | 1.343 | 0.669 | 0.266 | 1.748 |
|       | 	          | sigma_gamma   | 0.222 | 0.095 | 0.488 | 0.242 | 0.094 | 0.626 |
|       | 	          | sigma_delta_1 | 0.210 | 0.190 | 0.237 | 0.210 | 0.190 | 0.233 |
|       | 	          | sigma_delta_2 | 0.092 | 0.080 | 0.109 | 0.094 | 0.080 | 0.108 |
| pCAR  | FE M-models | sigma_delta_1 | 0.210 | 0.191 | 0.234 | 0.209 | 0.189 | 0.231 |
|       | 	          | sigma_delta_2 | 0.092 | 0.081 | 0.104 | 0.093 | 0.080 | 0.108 |
|       | 	          | rho_1         | 0.928 | 0.773 | 0.992 | 0.965 | 0.856 | 0.999 |
|       | 	          | rho_2         | 0.985 | 0.948 | 0.999 | 0.968 | 0.833 | 0.999 |
|       | RE M-models | sigma_theta   | 0.756 | 0.346 | 1.570 | 0.674 | 0.280 | 1.822 |
|       |             | sigma_gamma   | 0.292 | 0.198 | 0.397 | 0.256 | 0.089 | 0.713 |
|       | 	          | sigma_delta_1 | 0.214 | 0.189 | 0.252 | 0.210 | 0.188 | 0.232 |
|       | 	          | sigma_delta_2 | 0.094 | 0.079 | 0.117 | 0.093 | 0.080 | 0.107 |
|       |            	| rho_1         | 0.911 | 0.731 | 0.990 | 0.961 | 0.827 | 0.999 |
|       | 	          | rho_2         | 0.921 | 0.734 | 0.984 | 0.967 | 0.841 | 0.999 |
| LCAR	| FE M-models | sigma_delta_1 | 0.206 | 0.189 | 0.230 | 0.209 | 0.188 | 0.232 |
|       | 	          | sigma_delta_2 | 0.092 | 0.080 | 0.107 | 0.093 | 0.080 | 0.106 |
|       | 	          | lambda_1      | 0.830 | 0.558 | 0.980 | 0.860 | 0.568 | 0.996 |
|       | 	          | lambda_2      | 0.920 | 0.746 | 0.992 | 0.868 | 0.588 | 0.996 |
|       | RE M-models | sigma_theta   | 0.528 | 0.231 | 1.108 | 0.647 | 0.258 | 1.782 |
|       | 	          | sigma_gamma   | 0.184 | 0.142 | 0.261 | 0.258 | 0.094 | 0.741 |
|       | 	          | sigma_delta_1 | 0.212 | 0.199 | 0.226 | 0.209 | 0.187 | 0.230 |
|       | 	          | sigma_delta_2 | 0.093 | 0.082 | 0.104 | 0.093 | 0.080 | 0.108 |
|       | 	          | lambda_1      | 0.730 | 0.409 | 0.942 | 0.844 | 0.560 | 0.995 |
|       | 	          | lambda_2      | 0.945 | 0.808 | 0.996 | 0.853 | 0.548 | 0.996 |
| BYM2  | FE M-models | sigma_delta_1 | 0.210 | 0.200 | 0.219 | 0.210 | 0.189 | 0.232 |
|       | 	          | sigma_delta_2 | 0.094 | 0.087 | 0.101 | 0.093 | 0.081 | 0.107 |
|       | RE M-models | sigma_theta_s | 0.591 | 0.261 | 1.285 | 0.681 | 0.265 | 1.831 |
|       | 	          | sigma_theta_h | 0.042 | 0.017 | 0.079 | 0.053 | 0.002 | 0.199 |
|       |             | sigma_gamma   | 0.209 | 0.113 | 0.372 | 0.246 | 0.089 | 0.629 |
|       |        	    | sigma_delta_1 | 0.211 | 0.190 | 0.233 | 0.210 | 0.189 | 0.231 |
|       | 	          | sigma_delta_2 | 0.093 | 0.080 | 0.109 | 0.093 | 0.080 | 0.107 |



# R code
R code to fit the multivariate spatio-temporal models described in the paper, and to reproduce the results, has been included [here](https://github.com/spatialstatisticsupna/Mmodels_comparison_article/blob/master/R/).

# References
Vicente, G., Goicoa, T., and Ugarte, M.D. (2019). Multivariate spatio-temporal areal models: a comparison of INLA and WinBus for the analysis of gender-based violence in small areas. 
```diff
- PONER REFERENCIA BIEN
```