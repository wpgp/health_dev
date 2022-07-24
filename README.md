health_dev: Subnational reproductive, maternal, newborn, child and
adolescent health and development atlas for India, version 1.0
================
2022-07-24

------------------------------------------------------------------------

**Table 1.** Files and their descriptions within the health_dev GitHub
repository for the paper Subnational reproductive, maternal, newborn,
child, and adolescent health and development atlas for India.

| Name       | Type   | Description                                                                                                                                                                                                                                                                                                                                                                                                  |
|------------|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| out        | Folder | Folder contain the prediction and uncertainty gridded datasets (raster files) produced from the prediction R script and the out- of-sample cross validation summary statistics (csv files) from the validation R script.                                                                                                                                                                                     |
| rda        | Folder | Folder to contain INLA objects and the model summary statistics (saved as rda files) produced from the modelling R script. Files within this folder will be required to run the prediction and validation R scripts.                                                                                                                                                                                         |
| shp        | Folder | This folder contains the shapefiles required to run all R scripts in this repository. These should be the administrative boundaries of the study area as polygons and the location of the clusters in the study area as points (lat/lon) These shapefiles can be obtained from the DHS program at www.dhsprogram.com.                                                                                        |
| tif        | Folder | This folder contains the raster files for all geospatial covariates. Files within this folder are required to run the prediction R script. Examples of geospatial covariate datasets can be found at www.hub.worldpop.org/project/categories?id=14.                                                                                                                                                          |
| covariates | csv    | This file contains a demo of the format of the data extracted from geospatial covariates considered when modelling the health and development indicators. This file is required to run all R scripts in this repository. Examples of geospatial covariate datasets can be found from <https://hub.worldpop.org/project/categories?id=14>.                                                                    |
| indicators | csv    | This file contains a demo of health and development indicators to model. This file is required to run all R scripts in this repository. The indicators were extracted from the India NFHS-4 (National Family Health Survey 4) 2015-16 DHS (Demographic Health Survey) (1-3) database, which are publicly available after registration onto the Measure DHS website (www.dhsprogram.com).                     |
| modelling  | R      | R script for modelling the health and development indicators. The files required to run this script are the covariates and indicator csv files and the files in the shp folder. This script outputs an INLA object and the model summary statistics (both saved as rda files). Further description of the methodology is given in the sections below.                                                        |
| prediction | R      | R script for predicting the health and development indicators. The files required to run this script are the covariates and indicators csv files, the files in the shp folder, the files in the tif folder, and the files in the rda folder. This script outputs a prediction gridded dataset (tif file) and an uncertainty gridded dataset (tif file) for target indicator and are saved to the out folder. |
| validation | R      | R script for out-of-sample (k- fold) validation for the models of the health and development indicators. The files required to run this script are the covariates and indicators csv files, the files in the shp folder, and the files in the rda folder. This script outputs k- fold summary statistics as csv files. Further description of the methodology is given in the sections below                 |

------------------------------------------------------------------------

# Script for modelling the health and development indicators - modelling.R

The geospatial covariate selection is two-staged. In the first stage, we
check for multicollinearity amongst the geospatial covariates. In the
second stage, we employ the back-ward stepwise model selection method.

To check for multicollinearity, a Pearson correlation matrix for the
geospatial covariates is created and any pairs with a Pearson
correlation coefficient
![r\>0.8](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r%3E0.8 "r>0.8")
are flagged. The flagged covariates are then individually fitted in
non-Bayesian binomial generalised linear models (GLMs). The Bayesian
information criteria (BIC) of the models are then calculated. The
covariate in the model with a lower BIC is retained while the covariate
in the model with the greater BIC is omitted for the target indicator.
To further ensure that multicollinearity is not a problem between the
remaining geospatial covariates, variance inflation factors (VIFs) are
calculated. If any covariate returns a VIF \> 4, it is omitted.

After checking for multicollinearity, a backward model selection
algorithm is used to select the best (sub)set of geospatial covariates
for the target indicator. The algorithm is as follows. The remaining
geospatial covariates are fitted in a non-Bayesian binomial GLM and the
BIC is calculated. A covariate is removed from the model and the BIC is
recalculated. If the recalculated BIC is less than the previously
calculated BIC, this subset of covariates is preferred. These steps are
performed iteratively until the recalculated BIC is not less than the
BIC calculated from the previous iteration. At this point, the best
(sub)set of geospatial covariates have been attained and they will be
used when constructing the Bayesian point-referenced spatial binomial
GLM in INLA.

The constructed Bayesian point-referenced spatial binomial GLM is given
as follows.

![
\\begin{align\*}
Y(\\mathbf{s}\_i)\|m(\\mathbf{s}\_i) &\\sim    \\textrm{Binomial}(m(\\mathbf{s}\_i), p(\\mathbf{s}\_i)), \\nonumber \\\\
\\textrm{logit}(p(\\mathbf{s}\_i)) &= \\mathbf{x}(\\mathbf{s}\_i)\\boldsymbol{\\beta} + \\omega(\\mathbf{s}\_i) + \\epsilon(\\mathbf{s}\_i).
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0AY%28%5Cmathbf%7Bs%7D_i%29%7Cm%28%5Cmathbf%7Bs%7D_i%29%20%26%5Csim%20%20%20%20%5Ctextrm%7BBinomial%7D%28m%28%5Cmathbf%7Bs%7D_i%29%2C%20p%28%5Cmathbf%7Bs%7D_i%29%29%2C%20%5Cnonumber%20%5C%5C%0A%5Ctextrm%7Blogit%7D%28p%28%5Cmathbf%7Bs%7D_i%29%29%20%26%3D%20%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_i%29%5Cboldsymbol%7B%5Cbeta%7D%20%2B%20%5Comega%28%5Cmathbf%7Bs%7D_i%29%20%2B%20%5Cepsilon%28%5Cmathbf%7Bs%7D_i%29.%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
Y(\mathbf{s}_i)|m(\mathbf{s}_i) &\sim    \textrm{Binomial}(m(\mathbf{s}_i), p(\mathbf{s}_i)), \nonumber \\
\textrm{logit}(p(\mathbf{s}_i)) &= \mathbf{x}(\mathbf{s}_i)\boldsymbol{\beta} + \omega(\mathbf{s}_i) + \epsilon(\mathbf{s}_i).
\end{align*}
")

![
\\begin{align\*}
&\\omega(\\mathbf{s}\_i) \\sim N_n(\\boldsymbol{0}, \\Sigma\_\\omega), \\\\
&\\Sigma\_\\omega = \\sigma^2\_{\\omega}\\exp(-\\phi D).
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0A%26%5Comega%28%5Cmathbf%7Bs%7D_i%29%20%5Csim%20N_n%28%5Cboldsymbol%7B0%7D%2C%20%5CSigma_%5Comega%29%2C%20%5C%5C%0A%26%5CSigma_%5Comega%20%3D%20%5Csigma%5E2_%7B%5Comega%7D%5Cexp%28-%5Cphi%20D%29.%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
&\omega(\mathbf{s}_i) \sim N_n(\boldsymbol{0}, \Sigma_\omega), \\
&\Sigma_\omega = \sigma^2_{\omega}\exp(-\phi D).
\end{align*}
")

![
\\epsilon(\\mathbf{s}\_i) \\sim N(0, \\sigma^2\_\\epsilon)
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cepsilon%28%5Cmathbf%7Bs%7D_i%29%20%5Csim%20N%280%2C%20%5Csigma%5E2_%5Cepsilon%29%0A "
\epsilon(\mathbf{s}_i) \sim N(0, \sigma^2_\epsilon)
")

The number of occurrence of events of the target indicator
![Y(\\mathbf{s}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%28%5Cmathbf%7Bs%7D_i%29 "Y(\mathbf{s}_i)")
within cluster locations
![\\mathbf{s}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bs%7D_i "\mathbf{s}_i")
for
![i = 1, \\dots, n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%20%3D%201%2C%20%5Cdots%2C%20n "i = 1, \dots, n")
follows a Binomial distribution with the total number of surveys
conducted within the cluster locations
![m(\\mathbf{s}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%28%5Cmathbf%7Bs%7D_i%29 "m(\mathbf{s}_i)")
and the proportion of events happening in the cluster
![p(\\mathbf{s}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%28%5Cmathbf%7Bs%7D_i%29 "p(\mathbf{s}_i)").
With a logit link,
![p(\\mathbf{s}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%28%5Cmathbf%7Bs%7D_i%29 "p(\mathbf{s}_i)")
is calculated with a linear combination of the fixed effects
![\\mathbf{x}(\\mathbf{s}\_i)\\boldsymbol{\\beta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_i%29%5Cboldsymbol%7B%5Cbeta%7D "\mathbf{x}(\mathbf{s}_i)\boldsymbol{\beta}"),
spatial random effects
![\\omega(\\mathbf{s}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Comega%28%5Cmathbf%7Bs%7D_i%29 "\omega(\mathbf{s}_i)")
and independent identical (iid) random effects
![\\epsilon(\\mathbf{s}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%28%5Cmathbf%7Bs%7D_i%29 "\epsilon(\mathbf{s}_i)").

The fixed effects are given by the geospatial covariates
![\\mathbf{x}(\\mathbf{s}\_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bx%7D%28%5Cmathbf%7Bs%7D_i%29 "\mathbf{x}(\mathbf{s}_i)")
selected from the backward model selection algorithm mentioned above and
![\\boldsymbol{\\beta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cbeta%7D "\boldsymbol{\beta}")
is a vector of regression coefficients to be estimated. The spatial
random effects follow a multivariate normal distribution with zero-mean
and some covariance matrix
![\\Sigma\_\\omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CSigma_%5Comega "\Sigma_\omega").
In this study, elements of the covariance matrix are calculated with the
exponential covariance function. The exponential covariance function is
calculated with the spatial variance
![\\sigma^2\_\\omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%5E2_%5Comega "\sigma^2_\omega"),
the spatial decay parameter
![\\phi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cphi "\phi")
and the
![n\\times n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5Ctimes%20n "n\times n")
Euclidean distance matrix
![D](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D "D")
between the cluster locations. The parameters
![\\sigma^2\_\\omega](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%5E2_%5Comega "\sigma^2_\omega")
and
![\\phi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cphi "\phi")
are unknown and are to be estimated in INLA. The iid random effects
follow a normal distribution with a mean of zero and an unknown variance
![\\sigma^2\_\\epsilon](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma%5E2_%5Cepsilon "\sigma^2_\epsilon")
which will be estimated along with the other parameters mentioned above.

Additional components must be constructed before fitting the model in
INLA. First a mesh of the study domain is constructed with the shape
file and coordinates within the target indicator file. Using this mesh
object, a stochastic partial differential equation (SPDE) object is
defined with functions in INLA where the priors of the spatial decay
parameter and spatial variance parameter is defined. With the mesh
object, INLA stack “A” matrices are created and stacked with the INLA
stack functions. Finally, these components, along with the model are
fitted into the INLA function.

# Script for predicting the health and development indicators – prediction.R

The prediction R script loads the generates posterior samples from the
INLA object (saved from the modelling R script). Then it reads in the
raster files corresponding to the geospatial covariates of the model for
the target indicators and compiles it as a prediction data frame.
Finally, the predicted values are computed from the prediction data
frame, INLA mesh objects and INLA posterior sample objects, and are
slotted to the cells in the raster file – producing the high-resolution
(5x5km) prediction and uncertainty gridded datasets / surfaces as tif
files.

# Script for validating the models constructed for health and development indicators – validation.R

The validation R script accesses the performance of the model
constructed for the target indicator from the modelling R script with
k-fold cross validations and compute evaluation metrics. The k-fold
cross validation functions by first partitioning the dataset into k
parts, then training the model with k-1 parts of the dataset and testing
the trained model with the kth part of the dataset. The model is the
Bayesian point-referenced spatial generalized linear model constructed
in the modelling R script (i.e., with the same (sub)set of geospatial
covariates) for the target indicator. For each fold, the following
evaluation metrics are calculated:

![
\\begin{align\*}
\\rho(\\hat{\\mathbf{p}}\_i, \\mathbf{p}), \\\\
\\\\
\\sqrt{\\frac{1}{n\_{\\rm test}}\\sum^{n\_{\\rm test}}\_{i=1}(\\hat{p}\_i - p_i)^2}, \\\\
\\\\
\\frac{1}{n\_{\\rm test}}\\sum^{n\_{\\rm test}}\_{i=1}\|\\hat{p}\_i - p_i\|, \\\\
\\\\
\\bigg( \\frac{\\sum\_{i=1}^{n\_{\\rm test}}(\\hat{p}\_i - p_i)}{\\sum\_{j=1}^{n\_{\\rm test}}p_j}  \\bigg) \\times 100, \\\\
\\\\
\\bigg(\\frac{1}{n\_{\\rm test}}\\sum^{n\_{\\rm test}}\_{i=1}C_i\\bigg)\\times 100,
\\end{align\*}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Balign%2A%7D%0A%5Crho%28%5Chat%7B%5Cmathbf%7Bp%7D%7D_i%2C%20%5Cmathbf%7Bp%7D%29%2C%20%5C%5C%0A%5C%5C%0A%5Csqrt%7B%5Cfrac%7B1%7D%7Bn_%7B%5Crm%20test%7D%7D%5Csum%5E%7Bn_%7B%5Crm%20test%7D%7D_%7Bi%3D1%7D%28%5Chat%7Bp%7D_i%20-%20p_i%29%5E2%7D%2C%20%5C%5C%0A%5C%5C%0A%5Cfrac%7B1%7D%7Bn_%7B%5Crm%20test%7D%7D%5Csum%5E%7Bn_%7B%5Crm%20test%7D%7D_%7Bi%3D1%7D%7C%5Chat%7Bp%7D_i%20-%20p_i%7C%2C%20%5C%5C%0A%5C%5C%0A%5Cbigg%28%20%5Cfrac%7B%5Csum_%7Bi%3D1%7D%5E%7Bn_%7B%5Crm%20test%7D%7D%28%5Chat%7Bp%7D_i%20-%20p_i%29%7D%7B%5Csum_%7Bj%3D1%7D%5E%7Bn_%7B%5Crm%20test%7D%7Dp_j%7D%20%20%5Cbigg%29%20%5Ctimes%20100%2C%20%5C%5C%0A%5C%5C%0A%5Cbigg%28%5Cfrac%7B1%7D%7Bn_%7B%5Crm%20test%7D%7D%5Csum%5E%7Bn_%7B%5Crm%20test%7D%7D_%7Bi%3D1%7DC_i%5Cbigg%29%5Ctimes%20100%2C%0A%5Cend%7Balign%2A%7D%0A "
\begin{align*}
\rho(\hat{\mathbf{p}}_i, \mathbf{p}), \\
\\
\sqrt{\frac{1}{n_{\rm test}}\sum^{n_{\rm test}}_{i=1}(\hat{p}_i - p_i)^2}, \\
\\
\frac{1}{n_{\rm test}}\sum^{n_{\rm test}}_{i=1}|\hat{p}_i - p_i|, \\
\\
\bigg( \frac{\sum_{i=1}^{n_{\rm test}}(\hat{p}_i - p_i)}{\sum_{j=1}^{n_{\rm test}}p_j}  \bigg) \times 100, \\
\\
\bigg(\frac{1}{n_{\rm test}}\sum^{n_{\rm test}}_{i=1}C_i\bigg)\times 100,
\end{align*}
")

the Pearson’s correlation coefficient, the root mean squared error, the
mean absolute error, percentage bias, and the coverage rate. In the
evaluate metrics above,
![p_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_i "p_i")
is used to denote the observed values – i.e., the proportions of the
target indicators partitioned for testing – and
![\\hat{p}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bp%7D_i "\hat{p}_i")
is used to denote the predicted mean values from the Bayesian
point-referenced spatial binomial generalized linear model.

The notation
![\\rho(\\cdot)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Crho%28%5Ccdot%29 "\rho(\cdot)")
is used to the denote the Pearson’s correlation coefficient where
explicitly it is calculated with the covariance of the observed and
predicted values, and the standard deviation of the observed and
predicted values

![
\\rho(\\hat{\\mathbf{p}}\_i, \\mathbf{p}) = \\frac{\\textrm{cov}(\\hat{\\mathbf{p}}\_i, \\mathbf{p})}{\\sigma\_{\\bf \\hat{p}}\\sigma\_{\\bf p}}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Crho%28%5Chat%7B%5Cmathbf%7Bp%7D%7D_i%2C%20%5Cmathbf%7Bp%7D%29%20%3D%20%5Cfrac%7B%5Ctextrm%7Bcov%7D%28%5Chat%7B%5Cmathbf%7Bp%7D%7D_i%2C%20%5Cmathbf%7Bp%7D%29%7D%7B%5Csigma_%7B%5Cbf%20%5Chat%7Bp%7D%7D%5Csigma_%7B%5Cbf%20p%7D%7D%0A "
\rho(\hat{\mathbf{p}}_i, \mathbf{p}) = \frac{\textrm{cov}(\hat{\mathbf{p}}_i, \mathbf{p})}{\sigma_{\bf \hat{p}}\sigma_{\bf p}}
")

Here, note that the vectors
![\\hat{\\mathbf{p}} = (\\hat{p}\_1, \\dots, \\hat{p}\_n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Cmathbf%7Bp%7D%7D%20%3D%20%28%5Chat%7Bp%7D_1%2C%20%5Cdots%2C%20%5Chat%7Bp%7D_n%29 "\hat{\mathbf{p}} = (\hat{p}_1, \dots, \hat{p}_n)")
and
![\\mathbf{p} = (p_1, \\dots, p_n)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bp%7D%20%3D%20%28p_1%2C%20%5Cdots%2C%20p_n%29 "\mathbf{p} = (p_1, \dots, p_n)")
where
![n\_{\\rm test}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_%7B%5Crm%20test%7D "n_{\rm test}")
is the number of observations partitioned for testing. Better predictive
performance is reflected from a greater Pearson’s correlation
coefficient. The root mean squared error (RMSE), mean absolute error
(MAE) and percentage bias have straightforward calculations that does
not require additional explanation. Better predictive performance is
reflected from smaller RMSE, MAE and percentage bias values. The
coverage rate which ranges from 0 to 100. First,
![C_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_i "C_i")
in the equation is defined as follow

![
C_i = 
\\begin{cases}
1 \\hspace{0.5cm} \\textrm{if} \\hspace{0.3cm} \\hat{p}\_{i\\ 0.025q} \< p_i \< \\hat{p}\_{i\\ 0.975q}, \\\\
0 \\hspace{0.5cm} \\textrm{otherwise},
\\end{cases}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0AC_i%20%3D%20%0A%5Cbegin%7Bcases%7D%0A1%20%5Chspace%7B0.5cm%7D%20%5Ctextrm%7Bif%7D%20%5Chspace%7B0.3cm%7D%20%5Chat%7Bp%7D_%7Bi%5C%200.025q%7D%20%3C%20p_i%20%3C%20%5Chat%7Bp%7D_%7Bi%5C%200.975q%7D%2C%20%5C%5C%0A0%20%5Chspace%7B0.5cm%7D%20%5Ctextrm%7Botherwise%7D%2C%0A%5Cend%7Bcases%7D%0A "
C_i = 
\begin{cases}
1 \hspace{0.5cm} \textrm{if} \hspace{0.3cm} \hat{p}_{i\ 0.025q} < p_i < \hat{p}_{i\ 0.975q}, \\
0 \hspace{0.5cm} \textrm{otherwise},
\end{cases}
")

where
![\\hat{p}\_{i\\ 0.025q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bp%7D_%7Bi%5C%200.025q%7D "\hat{p}_{i\ 0.025q}")
and
![\\hat{p}\_{i\\ 0.975q}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bp%7D_%7Bi%5C%200.975q%7D "\hat{p}_{i\ 0.975q}")
represents the ith 0.025 quantile and 0.0975 quantile predicted value.
To put it simply,
![C_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_i "C_i")
is either 1 or 0, for
![i = 1,\\dots,n\_{\\rm test}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%20%3D%201%2C%5Cdots%2Cn_%7B%5Crm%20test%7D "i = 1,\dots,n_{\rm test}"),
depending on some condition. This condition is if the observed value is
within the 0.025 quantile and 0.0975 quantile of the predicted value,
![C_i = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_i%20%3D%201 "C_i = 1"),
otherwise
![C_i = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C_i%20%3D%200 "C_i = 0").
Better predictive performance is reflected from a higher coverage rate.

The validation R script returns csv files with the evaluation metrics
calculated for each fold for the model of the target indicator being
validated.

------------------------------------------------------------------------

# Acknowledement

The work is funded by the Children’s Investment Foundation Fund (CIFF).
The authors acknowledge the support of the PMO Team at WorldPop and
would like to thank EME and India Programme Team at CIFF for their
inputs and continuous support, and all staff at CIFF who provided
feedback at each stage of this work. Moreover, the authors would like to
thank the DHS Program staff for their input on the construction of some
of the indicators.

# Suggested citation

Chan, H.M.T, Dreoni, I., Tejedor-Garavito, N., Kerr D., Bonnie, A.,
Tatem A.J. and Pezzulo, C. 2022. health_dev: Subnational reproductive,
maternal, newborn, child and adolescent health and development atlas for
India, version 1.0. WorldPop, University of Southampton. doi: XXX.

# Reference

1.  International Institute for Population Sciences - IIPS/India and
    ICF. \[Producers\]. 2017. National Family Health Survey NFHS-4,
    \[Datasets IABR74DT.dta; IACR74DT.dta; IAHR74DT.dta; IAIR74DT.dta;
    IAKR74DT.dta; IAMR74DT.dta; IAPR74DT.dta; IAGE71FL.shp\], 2015-16:
    India. Mumbai: IIPS. ICF \[Distributor\], 2017. 6 International
    Institute for Population Sciences - IIPS/India and ICF. 2017.
    National Family Health Survey NFHS-4, 2015-16: India. Mumbai: IIPS.
    (www.dhsprogram.com)
2.  International Institute for Population Sciences (IIPS), I. and ICF.,
    India National Family Health Survey NFHS-4 2015-16. Mumbai, India:
    IIPS and ICF. Available at
    <http://dhsprogram.com/pubs/pdf/FR339/FR339.pdf>. 2017
3.  The DHS Program Code Share Project, Code Library, DHS Program. DHS
    Program Github site. <https://github.com/DHSProgram>., in DHS
    Program Github site. 2022.

------------------------------------------------------------------------
