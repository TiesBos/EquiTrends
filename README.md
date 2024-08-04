
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EquiTrends

<!-- badges: start -->

[![R-CMD-check](https://github.com/TiesBos/EquiTrends/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TiesBos/EquiTrends/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Testing for parallel trends is crucial in the Difference-in-Difference
framework. The goal of EquiTrends is to provide a set of functions to
test for equivalence of pre-trends in difference-in-differences
estimation, based on the placebo coefficients (used to compare the
trends of the treated and control group in the pre-treatment period to
some base period, generally the final period in the pre-treatment
period). The procedures follow the work of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)).

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

You can install the development version of EquiTrends from
[GitHub](https://github.com/TiesBos/EquiTrends) with:

``` r
# install.packages("devtools")
devtools::install_github("TiesBos/EquiTrends")
```

## Data Simulation

The `EquiTrends` package contains a function to simulate panel data,
tailored to the Difference-in-Difference framework. The function
`sim_paneldata` simulates a panel dataset with a given number of
individuals $N$ (`N`), number of periods $T+1$ (in the setting of this
package, indicating the number of pre-treatment periods. In code $T+1$
is referred to as `tt`), number of covariates $p$ (`p`), and treatment
effects. Generally, period $T+1$ is referred to as the . The function
returns a data frame with the following columns: `ID` (the
cross-sectional individual identifier), `period` (the time identifier),
`Y` (the dependent variable), `G` (a binary vector indicating if an
individual receives treatment, indicated by 1, or not, indicated by 0),
`X_1`, `X_2`, …, `X_p` (additional control variables). The function also
allows for the simulation of heterogeneity in treatment effects
(specified through `alpha`), time fixed effects (through `lambda`),
heteroscedasticty (specified through the binary variable `het`), serial
correlation (through the AR(1) coefficient `phi`), and clustering of the
standard errors. The construction of the dependent variable follows the
two-way fixed effect model, similar to the model in equation (2.5) of
Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)):
$$Y_{i,t} =  \alpha_i + \lambda_t + \sum_{l=1}^{T}{\beta_l}G_iD_l(t) + X_{1, i, t}\gamma_1+ \dots + X_{p,i,t}\gamma_p +u_{i,t} \quad \text{with} \; i=1,...,N, t=1,...,T+1$$
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
#>   ID period          Y G         X_1         X_2
#> 1  1      1  0.8924617 0  0.75300314 -1.17823063
#> 2  1      2  0.1826494 0  0.03631835 -0.07145882
#> 3  1      3 -1.6622285 0 -1.06100626  0.22117210
#> 4  1      4  3.4030413 0  0.64505208  1.01616032
#> 5  1      5  1.6724796 0  0.36741150  0.49046374
#> 6  2      1 -0.2539209 1 -0.85272295  0.35119328
```

## Testing for Equivalence of Pre-Trends

The `EquiTrends` package contains functions to test for equivalence of
pre-trends in difference-in-differences estimation. The functions
`maxEquivTest`, `meanEquivTest`, and `rmsEquivTest` are used to test for
equivalence of pre-trends in difference-in-differences estimation using
the placebo coefficients $\beta_{l} \; (l=1,...,T)$ estimates. The
functions are based on the work of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)).

### The `maxEquivTest` function

The maxEquivTest implements the equivalence testing procedure
surrounding the maximum absolute placebo coefficient. The function tests
the null hypothesis that the maximum placebo coefficient is larger than
or equal to a user-specified equivalence threshold, $\delta$, indicating
what negligible is to the user. That is, if
$$\lVert\beta\rVert_\infty = \max_{l=1,...T} |\beta_l|,$$ the testing
procedure can be represented as

$$H_0: \lVert\beta\rVert_\infty \geq \delta \quad \text{vs.} \quad H_1: \lVert\beta\rVert_\infty < \delta.$$
The null and alternative hypothesis can therefore be seen as
non-negligible and negligible differences in pre-trends, respectively.

The function `maxEquivTest` contains three testing procedures for this
test, as described in Section 4.2.1. of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)). The function
allows for the testing of the equivalence of pre-trends using a
bootstrap for spherical errors (`type = "Boot"`), a wild bootstrap for
clustered standard errors (`type = "Wild"`), and an Intersection Union
approach (`type = "IU"`) that rejects the null if all
$\beta_1,...,\beta_{T}$ are smaller than their individual critical
values. The function returns an object of class `maxEquivTestBoot` if
`type = "Boot"` or `type = "Wild"` or `maxEquivTestIU` if `type = "IU"`.
If no type is specified, `maxEquivTest` applies the Intersection Union
procedure for efficiency reasons.

#### Implemention of the `maxEquivTest` function with `type = "IU"`

Examples of implementing the Intersection unit test with different
possible variance-covariance matrices (required to perform the test) is
provided below (for more information on the possible variance-covariance
matrices, see the documentation of the `maxEquivTest` function). If an
equivalence threshold is supplied, the function will test the previous
hypothesis. If no equivalence threshold is supplied, the function finds
the critical value for the test at the specified significance level. The
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
  panel,
- `num_periods`: The number of periods in the panel,
- `num_observations`: The number of observations in the panel,
- `is_panel_balanced`: A logical value indicating whether the panel data
  is balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.
- Additionally, if `equiv_threshold_specified = TRUE`:
  - `IU_critical_values`: A numeric vector with the individual critical
    values for each of the placebo coefficients,
  - `reject_null_hypothesis`: A logical value indicating whether the
    null hypothesis of negligible pre-trend differences can be rejected
    at the specified significance level `alpha`,
  - `equiv_threshold`: The equivalence threshold employed.
- Additionally, if `equiv_threshold_specified = FALSE`:
  - `minimum_equiv_thresholds`: A numeric vector including for each
    placebo coefficient the minimum equivalence threshold for which the
    null hypothesis of negligible pre-trend differences can be rejected
    for the corresponding placebo coefficient individually,
  - `minimum_equiv_threshold`: A numeric scalar minimum equivalence
    threshold for which the null hypothesis of negligible pre-trend
    differences can be rejected for all placebo coefficients
    individually.

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
#> 0.02959          0.005772        0.9905        
#> 0.10194          0.005775        0.9905        
#> 0.19062          0.005779        0.9905        
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
#> 0.02959          0.005772        0.9905        
#> 0.10194          0.005775        0.9905        
#> 0.19062          0.005779        0.9905        
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
#> Minimum equivalence threshold to accept the alternative: 0.5312 
#> ---
#>  Estimate    Std. Error   Minimum Equivalence Threshold 
#> 0.207513     0.010036    0.2240    
#> 0.001381     0.010020    0.0126    
#> 0.515274     0.009705    0.5312    
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

# Note that the testing procedure can also handle unbalanced panels. 
# Finally, one should note that the test procedure also works for unbalanced panels.
# To illustrate this, we generate an unbalanced panel dataset by randomly selecting
# 70% of the observations from the balanced panel dataset:
random_indeces <- sample(nrow(sim_data), 0.7*nrow(sim_data))
unbalanced_sim_data <- sim_data[random_indeces, ]
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
- `bootstrap_critical_value`: The by bootstrap found critical value for
  the equivalence test based on the maximum absolute placebo
  coefficient,
- `reject_null_hypothesis`: A logical value indicating whether the null
  hypothesis of negligible pre-trend differences can be rejected at the
  specified significance level `alpha`,
- `B`: The number of bootstrap samples used to find the critical value,
- `significance_level`: The chosen significance level of the test
  `alpha`,
- `base_period`: The base period used in the testing procedure,
- `placebo_names`: The names corresponding to the placebo coefficients,
- `num_individuals`: The number of cross-sectional individuals in the
  panel,
- `num_periods`: The number of periods in the panel,
- `num_observations`: The total number of observations in the panel,
- `is_panel_balanced`: A logical value indicating whether the panel data
  is balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.

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
#> 0.5153                   0.6311                      TRUE      
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
#> 0.5153                   0.6332                      TRUE      
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
 # The bootstrap procedures can handle unbalanced panels:
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

The meanEquivTest implements the equivalence testing procedure
surrounding the mean placebo coefficient, as described in Section 4.2.2.
of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)). The function
tests the null hypothesis that the absolute mean placebo coefficient is
larger than or equal to a user-specified equivalence threshold,
$\delta$, indicating the what negligible is to the user. That is, if

$$\bar{\beta} = \frac{1}{T}\sum_{l=1}^{T} \beta_l,$$

the testing procedure can be represented as

$$H_0: |\bar{\beta}| \geq \delta \quad \text{vs.} \quad H_1: |\bar{\beta}| < \delta.$$
The null and alternative hypothesis can therefore be seen as
non-negligible and negligible differences in pre-trends, respectively.
Implementation of the test is similar to the `maxEquivTest` function, in
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
  panel,
- `num_periods`: The number of periods in the panel,
- `num_observations`: The total number of observations in the panel,
- `is_panel_balanced`: A logical value indicating whether the panel is
  balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.
- If `equiv_threshold_specified = TRUE`, then additionally:
  - `mean_critical_value`: The critical value at the alpha level,
  - `p_value`: The p-value of the test,
  - `reject_null_hypothesis`: A logical value indicating whether to
    reject the null hypothesis,
  - `equiv_threshold`: The equivalence threshold specified.
- If `equiv_threshold_specified = FALSE`, then additionally:
  - `minimum_equiv_threshold`: The minimum equivalence threshold for
    which the null hypothesis of non-negligible (based on the
    equivalence threshold) trend-differences can be rejected.

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
#> 0.1074                   0.01415     <2e-16  TRUE      
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
#> 0.1074                   0.01395     0.1303                
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

# Note that the testing procedure can also handle unbalanced panels. 
# Finally, one should note that the test procedure also works for unbalanced panels.
# To illustrate this, we generate an unbalanced panel dataset by randomly selecting
# 70% of the observations from the balanced panel dataset:
random_indeces <- sample(nrow(sim_data), 0.7*nrow(sim_data))
unbalanced_sim_data <- sim_data[random_indeces, ]
meanEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", X = c(5, 6),
              data = unbalanced_sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
              base_period = 4, vcov = "HAC")
```

### The `rmsEquivTest` function

The rmsEquivTest implements the equivalence testing procedure
surrounding the root mean squared placebo coefficient as described in
section 4.2.3 of Dette & Schumann
([2024](https://doi.org/10.1080/07350015.2024.2308121)). The function
tests the null hypothesis that the root mean squared placebo coefficient
is larger to a user-specified equivalence threshold, $\delta$,
indicating the what negligible is to the user. That is, if

$$\beta_{RMS} = \sqrt{\frac{1}{T}\sum_{l=1}^{T} \beta_l^2},$$

the testing procedure can be represented as

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
  panel,
- `num_periods`: The number of periods in the panel,
- `num_observations`: The total number of observations in the panel,
- `is_panel_balanced`: A logical value indicating whether the panel data
  is balanced,
- `equiv_threshold_specified`: A logical value indicating whether an
  equivalence threshold was specified.
- If `equiv_threshold_specified = TRUE`, then additionally:
  - `rms_critical_value`: The critical value at the alpha level,
  - `reject_null_hypothesis`: A logical value indicating whether to
    reject the null hypothesis,
  - `equiv_threshold`: The equivalence threshold specified.
- If `equiv_threshold_specified = FALSE`, then additionally:
  - `minimum_equiv_threshold`: The minimum equivalence threshold for
    which the null hypothesis of non-negligible (based on the
    equivalence threshold) trend-differences can be rejected.

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
#> Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold of 1 .
#> ---
#> RMS Placebo Effect   Simulated Crit. Val.    Reject H0 
#> 0.126                0.9807                  TRUE      
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

``` r
# The testing procedures can also be performed without specifying the 
# equivalence threshold specified. Then, the minimum equivalence threshold is returned
# for which the null hypothesis of non-negligible trend-differences can be rejected.
# Again, the three possible ways of entering the data as above can be used:
rmsEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", X = c("X_1", "X_2"),
             data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:4,
             base_period = 4)
#> 
#>                ==================================================
#>                Equivalence Tests for Pre-trends in DiD Estimation
#>                ==================================================
#> Type: Root Mean Squared Placebo Effect 
#> Significance level: 0.05 
#> Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold.
#> ---
#> RMS Placebo Effect   Min. Equiv. Threshold 
#> 0.126                0.2448                
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
# Finally, one should note that the test procedure also works for unbalanced panels.
# To illustrate this, we generate an unbalanced panel dataset by randomly selecting
# 70% of the observations from the balanced panel dataset:

random_indeces <- sample(nrow(sim_data), 0.7*nrow(sim_data))
unbalanced_sim_data <- sim_data[random_indeces, ]
#  With Equivalence Threshold:
rmsEquivTest(Y = 3, ID = 1, G = 4, period = 2, X = c(5, 6),
             data = unbalanced_sim_data, equiv_threshold = 1, 
             pretreatment_period = 1:4, base_period = 4)

#  Without Equivalence Threshold:
rmsEquivTest(Y = 3, ID = 1, G = 4, period = 2, X = c(5, 6),
             data = unbalanced_sim_data, equiv_threshold = NULL, 
             pretreatment_period = 1:4, base_period = 4)
```

## References

Dette H., & Schumann M. (2024). “Testing for Equivalence of Pre-Trends
in Difference-in-Differences Estimation.” *Journal of Business &
Economic Statistics*, 1–13. DOI:
[10.1080/07350015.2024.2308121](https://doi.org/10.1080/07350015.2024.2308121)
