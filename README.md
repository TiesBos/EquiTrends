
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EquiTrends

<!-- badges: start -->

[![R-CMD-check](https://github.com/TiesBos/EquiTrends/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TiesBos/EquiTrends/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Codecov test
coverage](https://codecov.io/gh/TiesBos/EquiTrends/graph/badge.svg)](https://app.codecov.io/gh/TiesBos/EquiTrends)
<!-- badges: end -->

`EquiTrends` is an R package for equivalence testing in the context of
Difference-in-Differences estimation. It allows users to test if
pre-treatment trends in the treated group are “equivalent” to those in
the control group. Here, “equivalence” means that rejection of the null
hypothesis implies that a function of the pre-treatment placebo effects
(maximum absolute, average or root mean squared value) does not exceed a
pre-specified threshold below which trend differences are considered
negligible. The package is based on the theory developed in Dette &
Schumann ([2024](https://doi.org/10.1080/07350015.2024.2308121)).

The package contains the functions `maxEquivTest` to perform the testing
procedure surrounding the maximum placebo coefficient (see equation
(3.1) of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121))),
`meanEquivTest` to perform the testing procedure surrounding the mean
placebo coefficient (see equation (3.2) of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121))) and
`rmsEquivTest` to perform the testing procedure surrounding the root
mean squared placebo coefficient (see equation (3.3) and (3.4) of Dette
& Schumann ([2024](https://doi.org/10.1080/07350015.2024.2308121))).
Furthermore, the package contains the function `sim_paneldata` to
simulate a paneldataset for such testing purposes.

## Installation

You can install the development version of `EquiTrends` from
[GitHub](https://github.com/TiesBos/EquiTrends) with:

``` r
# install.packages("devtools")
devtools::install_github("TiesBos/EquiTrends")
```

## Data Simulation

The `EquiTrends` package contains a function to simulate panel data,
tailored to the Difference-in-Differences framework. The function
`sim_paneldata` simulates a panel dataset with a given number of
individuals $N$ (`N`), number of periods $T+1$ (in the setting of this
package, indicating the number of pre-treatment periods. In
`sim_paneldata` $T+1$ is referred to as `tt`), number of covariates $p$
(`p`), and treatment effects. Typically, period $T+1$ is referred to as
the “base period”. The function also allows for the simulation of
heterogeneity in treatment effects (specified through `eta`) and time
fixed effects (through `lambda`). Furthermore, the function allows for
heteroscedasticty (specified through the binary variable `het`), serial
correlation (through the AR(1) coefficient `phi`:
$u_{i,t} = \phi u_{i,t-1} + v_{i,t}$ where $v_{i,t}$ follows an i.i.d.
$N(0,\sigma^2)$ distribution and $\sigma$ is specified through `sd`),
and clustering in the model errors $u_{i,t}$. The function returns a
data frame with the following columns: `ID` (the cross-sectional
individual identifier), `period` (the time identifier), `Y` (the
dependent variable), `G` (a binary vector indicating if an individual
receives treatment, indicated by 1, or not, indicated by 0), and `X_1`,
`X_2`, …, `X_p` (additional control variables). The construction of the
dependent variable follows the two-way fixed effect model, similar to
the model in equation (2.5) of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)):

$$Y_{i,t} =  \eta_i + \lambda_t + \sum_{l=1}^{T}{\beta_l}G_iD_l(t) + X_{1, i, t}\gamma_1+ \dots + X_{p,i,t}\gamma_p +u_{i,t} \quad \text{with} \  \ i=1,...,N, \ \ t=1,...,T+1$$

where $D_l(t)$ is a dummy variable that equals 1 if $t=l$ and 0
otherwise. The error-terms $u_{i,t}$ are generated through a normal
distribution with mean 0 and a variance-covariance structure depending
on the user-specified parameters. In the following, the $\beta_l$
coefficients are referred to as placebo coefficients, since they
represent the difference in pre-trends between the treatment and control
group before treatment has been assigned.

An example of the `sim_paneldata` function is provided below:

``` r
library(EquiTrends)

# Simulate a panel dataset with 500 individuals, 5 periods, 2 additional 
# regressors, and a binary treatment variable without heteroscedasticity, 
# serial correlation, and clustering. Furthermore, there are no fixed effects or 
# pre-trends in the model (since all values in beta are 0).
sim_data <- sim_paneldata(N = 500, tt = 5, p = 2, beta = rep(0, 5), 
                          gamma = rep(1, 2), het = 0, phi = 0, sd = 1, 
                          burnins = 50)
head(sim_data)
#>   ID period          Y G        X_1          X_2
#> 1  1      1 -0.8123777 0 -0.4407095 -0.655157012
#> 2  1      2 -1.8888861 0 -0.2212108 -0.349846262
#> 3  1      3 -2.2912561 0 -0.9741446 -0.000781637
#> 4  1      4 -0.5314161 0  0.2259398 -1.557426790
#> 5  1      5 -1.5528134 0 -0.1413597 -1.590501621
#> 6  2      1  1.5202663 1 -0.1386675  1.074761245
```

## Testing for Equivalence of Pre-Trends

The `EquiTrends` package contains functions to test for equivalence of
pre-trends in Difference-in-Differences estimation. The functions
`rmsEquivTest`, `meanEquivTest`, and `maxEquivTest` are used to test for
equivalence of pre-trends in Difference-in-Differences estimation using
the placebo coefficients $\beta_{l} \ (l=1,...,T)$ estimates. The
functions are based on the work of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)).

### The `rmsEquivTest` function

`rmsEquivTest` implements the equivalence testing procedure surrounding
the root mean squared placebo coefficient as described in section 4.2.3
of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)). The function
tests the null hypothesis that the root mean squared placebo coefficient
is larger than or equal to a user-specified equivalence threshold
$\delta$. That is, if

$$\beta_{RMS} = \sqrt{\frac{1}{T}\sum_{l=1}^{T} \beta_l^2},$$

the tested hypotheses can be represented as

$$H_0: \beta_{RMS} \geq \delta \quad \text{vs.} \quad H_1: \beta_{RMS} < \delta.$$

The null and alternative hypothesis can therefore be seen as
non-negligible and negligible differences in pre-trends, respectively.
The function returns an object of class `rmsEquivTest` containing

- `placebo_coefficients`: A numeric vector of the estimated placebo
  coefficients,
- `rms_placebo_coefs`: The root mean squared value of the placebo
  coefficients,
- `significance_level`: The significance level of the test,
- `base_period`: The base period used in the testing procedure,
- `num_individuals`: The number of cross-sectional individuals in the
  panel used for testing,
- `num_periods`: The number of pre-treatment periods in the panel used
  for testing (if the panel is unbalanced, `num_periods` represents the
  range in the number of time periods covered by different individuals),
- `num_observations`: The total number of observations in the panel used
  for testing,
- `is_panel_balanced`: A logical value indicating whether the used panel
  is balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.
- If `equiv_threshold_specified = TRUE`, then additionally:
  - `rms_critical_value`: The critical value at the chosen significance
    level,
  - `reject_null_hypothesis`: A logical value indicating whether to
    reject the null hypothesis,
  - `equiv_threshold`: The equivalence threshold specified.
- If `equiv_threshold_specified = FALSE`, then additionally:
  - `minimum_equiv_threshold`: The minimum equivalence threshold for
    which the null hypothesis of non-negligible trend-differences can be
    rejected.

One should note that rows containing `NA` values are removed from the
panel before the testing procedure is performed.

Please be aware that the equivalence test based on the root mean squared
placebo coefficient applies a randomization technique (as described by
Dette & Schumann (2024)), leading to a stochastic critical value and
minimum equivalence threshold. Therefore, the results may vary between
different runs of the function.

``` r
# Perform the equivalence test using an equivalence threshold of 1 with periods 
# 1-4 as pre-treatment periods based on the RMS testing procedure:
#  - option 1: using column names in the panel
# One can use the names of the columns in the panel to specify the variables:
rmsEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", X = c("X_1", "X_2"),
             data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
             base_period = 4)
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Root Mean Squared Placebo Effect 
#> Significance level: 0.05 
#> Alternative hypothesis: the root mean squared placebo effect does not exceed the equivalence threshold of 1 .
#> ---
#> RMS Placebo Effect   Simulated Crit. Val.    Reject H0 
#> 0.1835               0.9558                  TRUE      
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel: 
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

``` r
#  - option 2: using column numbers in the panel 
# Alternatively, one can use the column numbers in the panel to specify the variables:
rmsEquivTest(Y = 3, ID = 1, G = 4, period = 2, X = c(5, 6),
             data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
             base_period = 4)
             
#  - option 3: using separate variables 
# One can also use the variables directly without specifying the data variable:
data_Y <- sim_data$Y
data_ID <- sim_data$ID
data_G <- sim_data$G
data_period <- sim_data$period
data_X <- cbind(sim_data$X_1, sim_data$X_2)

rmsEquivTest(Y = data_Y, ID = data_ID, G = data_G, period = data_period, X = data_X,
             equiv_threshold = 1, pretreatment_period = 1:4,
             base_period = 4)
```

The testing procedures can also be performed without specifying the
equivalence threshold. Then, the minimum equivalence threshold is
returned for which the null hypothesis of non-negligible
trend-differences can be rejected. Again, the three possible ways of
entering the data as above can be used.

``` r
rmsEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", X = c("X_1", "X_2"),
             data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:4,
             base_period = 4)
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Root Mean Squared Placebo Effect 
#> Significance level: 0.05 
#> Alternative hypothesis: the root mean squared placebo effect does not exceed the equivalence threshold.
#> ---
#> RMS Placebo Effect   Min. Equiv. Threshold 
#> 0.1835               0.2558                
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel: 
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

Finally, one should note that the test procedure also works for
unbalanced panels.

``` r
# To illustrate this, we generate an unbalanced panel dataset by randomly selecting
# 70% of the observations from the balanced panel dataset:

random_indices <- sample(nrow(sim_data), 0.7*nrow(sim_data))
unbalanced_sim_data <- sim_data[random_indices, ]
#  With Equivalence Threshold:
rmsEquivTest(Y = 3, ID = 1, G = 4, period = 2, X = c(5, 6),
             data = unbalanced_sim_data, equiv_threshold = 1, 
             pretreatment_period = 1:4, base_period = 4)

#  Without Equivalence Threshold:
rmsEquivTest(Y = 3, ID = 1, G = 4, period = 2, X = c(5, 6),
             data = unbalanced_sim_data, equiv_threshold = NULL, 
             pretreatment_period = 1:4, base_period = 4)
```

### The `maxEquivTest` function

The `maxEquivTest` function tests the null hypothesis that the maximum
placebo coefficient is larger than or equal to a user-specified
equivalence threshold $\delta$. That is, if

$$\lVert\beta\rVert_\infty = \max_{l=1,...T} |\beta_l|,$$

the tested hypotheses can be represented as

$$H_0: \lVert\beta\rVert_\infty \geq \delta \quad \text{vs.} \quad H_1: \lVert\beta\rVert_\infty < \delta.$$

The null and alternative hypothesis can therefore be seen as
non-negligible and negligible differences in pre-trends, respectively.

The function `maxEquivTest` contains three testing procedures for this
test, as described in Section 4.2.1. of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)). The function
allows for the testing of the equivalence of pre-trends using a
bootstrap for spherical errors (`type = "Boot"`), a wild bootstrap for
clustered standard errors (`type = "Wild"`), and an Intersection Union
approach (`type = "IU"`) that rejects the null if all estimates for
$\beta_1,...,\beta_{T}$ are smaller than their individual critical
values. The function returns an object of class `maxEquivTestBoot` if
`type = "Boot"` or `type = "Wild"` or `maxEquivTestIU` if `type = "IU"`.
If no type is specified, `maxEquivTest` applies the Intersection Union
procedure for efficiency reasons.

#### Implemention of the `maxEquivTest` function with `type = "IU"`

Examples of implementing the Intersection unit test with different
possible variance-covariance matrices (required to perform the test) are
provided below (for more information on the possible variance-covariance
matrices, see the documentation of the `maxEquivTest` function). If an
equivalence threshold is supplied, the function will test the previous
hypothesis. If no equivalence threshold is supplied, the function finds
the minimum equivalence threshold for which the null of non-negligible
trend-differences can be reject using the Intersection Union test. The
function returns an object of class `maxEquivTestIU` containing the
following information:

- `placebo_coefficients`: A numeric vector of the estimated placebo
  coefficients,
- `abs_placebo_coefficients`: A numeric vector with the absolute values
  of estimated placebo coefficients,
- `placebo_coefficients_se`: A numeric vector with the standard errors
  of the placebo coefficients,
- `significance_level`: The chosen significance level of the test,
- `base_period`: The base period used in the testing procedure,
- `placebo_names`: The names corresponding to the placebo coefficients,
- `num_individuals`: The number of cross-sectional individuals in the
  panel used for testing,
- `num_periods`: The number of pre-treatment periods in the panel used
  for testing (if the panel is unbalanced, `num_periods` represents the
  range in the number of time periods covered by different individuals),
- `num_observations`: The number of observations in the panel used for
  testing,
- `is_panel_balanced`: A logical value indicating whether the panel data
  is balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.
- If `equiv_threshold_specified = TRUE`, then additionally:
  - `IU_critical_values`: A numeric vector with the individual critical
    values for each of the placebo coefficients,
  - `reject_null_hypothesis`: A logical value indicating whether the
    null hypothesis of negligible pre-trend differences can be rejected
    at the specified significance level,
  - `equiv_threshold`: The equivalence threshold employed.
- If `equiv_threshold_specified = FALSE`, then additionally:
  - `minimum_equiv_thresholds`: A numeric vector including for each
    placebo coefficient the minimum equivalence threshold for which the
    null hypothesis of negligible pre-trend differences can be rejected
    for the corresponding placebo coefficient individually,
  - `minimum_equiv_threshold`: A numeric scalar minimum equivalence
    threshold for which the null hypothesis of negligible pre-trend
    differences can be rejected for all placebo coefficients.

One should note that rows containing `NA` values are removed from the
panel before the testing procedure is performed.

``` r
# Perform the test with equivalent threshold specified as 1 based on 
# pre-treatment periods 1-4 and homoscedastic error-terms:
  # To select variables, one can use the column names / numbers in the panel data
maxEquivTest(Y = "Y", ID = "ID", G = "G", period = 2, X= c(5,6),
              data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
              base_period = 4, type = "IU")
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Intersection Union 
#> Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of 1 .
#> Reject null hypothesis: TRUE 
#> ( Critical values are printed for the significance level: 0.05 )
#> ---
#> Abs. Estimate    Std. Error  Critical Value 
#> 0.09848          0.1221          0.7992        
#> 0.27253          0.1221          0.7992        
#> 0.13041          0.1220          0.7993        
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel: 
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000

  # Alternatively, one can enter the variables separately:
data_Y <- sim_data$Y
data_ID <- sim_data$ID
data_G <- sim_data$G
data_period <- sim_data$period
data_X <- sim_data[, c(5, 6)]
maxEquivTest(Y = data_Y, ID = data_ID, G = data_G, period = data_period, X = data_X,
             equiv_threshold = 1, pretreatment_period = 1:4,
             base_period = 4, type = "IU")
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Intersection Union 
#> Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of 1 .
#> Reject null hypothesis: TRUE 
#> ( Critical values are printed for the significance level: 0.05 )
#> ---
#> Abs. Estimate    Std. Error  Critical Value 
#> 0.09848          0.1221          0.7992        
#> 0.27253          0.1221          0.7992        
#> 0.13041          0.1220          0.7993        
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel: 
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

``` r
# Perform the test without specifying the equivalence threshold with heteroscedastic 
# and autocorrelation robust variance-covariance matrix estimator:
maxEquivTest(Y = 3, ID = 1, G = 4, period = 2, 
             data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:4,
             base_period = 4, type = "IU", vcov = "HAC")
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Intersection Union 
#> Significance level: 0.05 
#> Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold.
#> Minimum equivalence threshold to accept the alternative: 0.4974 
#> ---
#>  Estimate    Std. Error   Minimum Equivalence Threshold 
#> 0.1028       0.2172      0.4489    
#> 0.1558       0.2088      0.4974    
#> 0.1257       0.2131      0.4711    
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel: 
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

``` r
# Perform the test without specifying the equivalence threshold with a custom
# variance-covariance matrix estimator:
vcov_func <- function(x) {plm::vcovHC(x, method = "white1", type = "HC2")}
maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
             data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
             base_period = 4, type = "IU", vcov = vcov_func)
 
# Perform the test using clustered standard errors based on a vector indicating 
# the cluster. For instance, two clusters with the following rule: all
# individuals with an ID below 250 are in the same cluster.
cluster_ind <- ifelse(sim_data$ID < 250, 1, 2)
maxEquivTest(Y = data_Y, ID = data_ID, G = data_G, period = data_period, X = data_X,
               equiv_threshold = 1, pretreatment_period = 1:4,
               base_period = 4, type = "IU", vcov = "CL", cluster = cluster_ind)
#> Registered S3 method overwritten by 'clubSandwich':
#>   method    from    
#>   bread.mlm sandwich
```

Note that the testing procedure can also handle unbalanced panels.

``` r
# To illustrate this, we generate an unbalanced panel dataset by randomly selecting
# 70% of the observations from the balanced panel dataset:
random_indices <- sample(nrow(sim_data), 0.7*nrow(sim_data))
unbalanced_sim_data <- sim_data[random_indices, ]
maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", X = c(5, 6),
              data = unbalanced_sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
              base_period = 4, type = "IU", vcov = "HAC")
```

##### Implementation of the bootstrap approaches

Examples of implementing the bootstrap based test are provided below.
For both `type = "Boot"` and `type = "Wild"`, an equivalence threshold
is required to perform the test. Furthermore, both testing procedures
return an object of class “maxEquivTestBoot” containing

- `placebo_coefficients`: A numeric vector of the estimated placebo
  coefficients,
- `abs_placebo_coefficients`: A numeric vector with the absolute values
  of estimated placebo coefficients,
- `max_abs_coefficient`: The maximum absolute estimated placebo
  coefficient,
- `B`: The number of bootstrap samples used to find the critical value,
- `significance_level`: The chosen significance level of the test,
- `base_period`: The base period used in the testing procedure,
- `placebo_names`: The names corresponding to the placebo coefficients,
- `num_individuals`: The number of cross-sectional individuals in the
  panel used for testing,
- `num_periods`: The number of pre-treatment periods in the panel used
  for testing (if the panel is unbalanced, `num_periods` represents the
  range in the number of time periods covered by different individuals),
- `num_observations`: The total number of observations in the panel used
  for testing,
- `is_panel_balanced`: A logical value indicating whether the panel data
  is balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.
- If `equiv_threshold_specified = TRUE`, then additionally:
  - `bootstrap_critical_value`: The by bootstrap found critical value
    for the equivalence test based on the maximum absolute placebo
    coefficient,
  - `reject_null_hypothesis`: A logical value indicating whether the
    null hypothesis of negligible pre-trend differences can be rejected
    at the specified significance level,
- If `equiv_threshold_specified = FALSE`, then additionally:
  - `minimum_equiv_threshold`: The minimum equivalence threshold for
    which the null hypothesis of negligible pre-trend differences can be
    rejected for the bootstrap procedure.

One should note that rows containing `NA` values are removed from the
panel before the testing procedure is performed.

On top of that, please be aware that the bootstrap procedures for the
equivalence test based on the maximum absolute placebo coefficient apply
a bootstrap procedure (as described by Dette & Schumann (2024)), leading
to a stochastic critical value and minimum equivalence threshold.
Therefore, the results may vary slightly between different runs of the
function.

The bootstrap for spherical errors with 1000 bootstrap iterations:

``` r
# Perform the test with equivalence threshold specified as 1 based on 
# pre-treatment periods 1:4 (with base period 4) with the general bootstrap procedure:
maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
             data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
             base_period = 4, type = "Boot", B = 1000)
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Bootstrap for Spherical Errors  (Based on 1000 bootstrap samples)
#> Significance level: 0.05 
#> Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of 1 .
#> ---
#> Max. Abs. Coefficient    Bootstrap Critical Value    Reject H0 
#> 0.1558                   0.6586                      TRUE      
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel:
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

The Wild boostrap with 1000 bootstrap iterations:

``` r
# Perform the test with the equivalence threshold specified as 1 based on 
# pre-treatment periods 1:4 (with base period 4) with the wild bootstrap procedure:
maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
             data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
             base_period = 4, type = "Wild")
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Cluster Wild Bootstrap (Based on 1000 bootstrap samples)
#> Significance level: 0.05 
#> Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of 1 .
#> ---
#> Max. Abs. Coefficient    Bootstrap Critical Value    Reject H0 
#> 0.1558                   0.6642                      TRUE      
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel:
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

The bootstrap procedures can handle unspecified equivalence thresholds:

``` r
maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
             data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:4,
             base_period = 4, type = "Boot", B = 1000)
maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
             data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:4,
             base_period = 4, type = "Wild", B = 1000)
```

The bootstrap procedures can handle unbalanced panels:

``` r
 maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
             data = unbalanced_sim_data, equiv_threshold = 1, 
             pretreatment_period = 1:4,
             base_period = 4, type = "Boot", B = 1000)
 maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
             data = unbalanced_sim_data, equiv_threshold = 1, 
             pretreatment_period = 1:4,
             base_period = 4, type = "Wild", B = 1000) 
```

### The `meanEquivTest` function

The `meanEquivTest` implements the equivalence testing procedure
surrounding the mean placebo coefficient, as described in Section 4.2.2.
of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)). The function
tests the null hypothesis that the absolute mean placebo coefficient is
larger than or equal to a user-specified equivalence threshold,
$\delta$. That is, if

$$\bar{\beta} = \frac{1}{T}\sum_{l=1}^{T} \beta_l,$$

the tested hypotheses can be represented as

$$H_0: |\bar{\beta}| \geq \delta \quad \text{vs.} \quad H_1: |\bar{\beta}| < \delta.$$

The null and alternative hypothesis can therefore be seen as
non-negligible and negligible differences in pre-trends, respectively.
Implementation of the test is similar to the `maxEquivTest` function in
terms of the possible variance-covariance matrices (for more information
on the possible variance-covariance matrices, see the documentation of
the `meanEquivTest` function). The function returns an object of class
`meanEquivTest` containing

- `placebo_coefficients`: A numeric vector of the estimated placebo
  coefficients,
- `abs_mean_placebo_coefs`: The absolute value of the mean of the
  placebo coefficients,
- `var_mean_placebo_coef`: The estimated variance of the mean placebo
  coefficient,
- `significance_level`: The significance level of the test,
- `base_period`: The base period used in the testing procedure,
- `num_individuals`: The number of cross-sectional individuals in the
  panel used for testing,
- `num_periods`: The number of pre-treatment periods in the panel used
  for testing (if the panel is unbalanced, `num_periods` represents the
  range in the number of time periods covered by different individuals)
- `num_observations`: The total number of observations in the panel used
  for testing,
- `is_panel_balanced`: A logical value indicating whether the panel is
  balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.
- If `equiv_threshold_specified = TRUE`, then additionally:
  - `mean_critical_value`: The critical value at the chosen significance
    level,
  - `p_value`: The p-value of the test,
  - `reject_null_hypothesis`: A logical value indicating whether to
    reject the null hypothesis,
  - `equiv_threshold`: The equivalence threshold specified.
- If `equiv_threshold_specified = FALSE`, then additionally:
  - `minimum_equiv_threshold`: The minimum equivalence threshold for
    which the null hypothesis of non-negligible (based on the
    equivalence threshold) trend-differences can be rejected.

One should note that rows containing `NA` values are removed from the
panel before the testing procedure is performed.

``` r
# Perform the test with equivalent threshold specified as 1 based on 
# pre-treatment periods 1-4 and assuming homoscedastic error-terms:
  # To select variables, one can use the column names / column numbers in the panel data:
  meanEquivTest(Y = "Y", ID = "ID", G = "G", period = 2, X = c(5, 6),
                data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
                base_period = 4)
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Mean Placebo Effect 
#> Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold of 1 .
#> ---
#> Abs. Mean Placebo Effect Std. Error  p-value Reject H0 
#> 0.1671                   0.09965     <2e-16  TRUE      
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel: 
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

``` r
  # Alternatively, one can use separate variables:
  data_Y <- sim_data$Y
  data_ID <- sim_data$ID
  data_G <- sim_data$G
  data_period <- sim_data$period
  data_X <- sim_data[, c(5, 6)]
  meanEquivTest(Y = data_Y, ID = data_ID, G = data_G, period = data_period, X = data_X,
                equiv_threshold = 1, pretreatment_period = 1:4,
                base_period = 4)
```

``` r
# Perform the test with a heteroscedastic and autocorrelation robust 
# variance-covariance matrix estimator, and without specifying the equivalence threshold:
meanEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", X = c(5, 6),
              data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:4,
              base_period = 4, vcov = "HAC")
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Mean Placebo Effect 
#> Significance level: 0.05 
#> Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold.
#> ---
#> Abs. Mean Placebo Effect Std. Error  Min. Equiv. Threshold 
#> 0.1671                   0.09691     0.3265                
#> ---
#> No. placebo coefficients estimated: 3 
#> Base period: 4 
#>  
#> Balanced Panel: 
#>  + No. pre-treatment periods: 4 
#>  + No. individuals: 500 
#>  + Total no. observations: 2000
```

``` r
# Perform the test with an equivalence threshold of 1 and a custom
# variance-covariance matrix estimator:
vcov_func <- function(x) {plm::vcovHC(x, method = "white1", type = "HC2")}
meanEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
              data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
              base_period = 4, vcov = vcov_func)
               
# Perform the test using clustered standard errors based on a vector indicating 
# the cluster. For instance, two clusters with the following rule: all
# individuals with an ID below 250 are in the same cluster:
cluster_ind <- ifelse(sim_data$ID < 250, 1, 2)
meanEquivTest(Y = data_Y, ID = data_ID, G = data_G, period = data_period, X = data_X,
               equiv_threshold = 1, pretreatment_period = 1:4,
               base_period = 4, vcov = "CL", cluster = cluster_ind)
```

Note that the testing procedure can also handle unbalanced panels:

``` r
# Finally, one should note that the test procedure also works for unbalanced panels.
# To illustrate this, we generate an unbalanced panel dataset by randomly selecting
# 70% of the observations from the balanced panel dataset:
random_indices <- sample(nrow(sim_data), 0.7*nrow(sim_data))
unbalanced_sim_data <- sim_data[random_indices, ]
meanEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", X = c(5, 6),
              data = unbalanced_sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
              base_period = 4, vcov = "HAC")
```

## References

Dette H., & Schumann M. (2024). “Testing for Equivalence of Pre-Trends
in Difference-in-Differences Estimation.” *Journal of Business &
Economic Statistics*, 1–13. DOI:
[10.1080/07350015.2024.2308121](https://doi.org/10.1080/07350015.2024.2308121)
