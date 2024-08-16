# ----------- The Mean Test Function -------------------------------------------
#' An internal function of the EquiTrends Mean Equivalence Testing procedure
#' 
#' @description This is a supporting function of the \code{meanEquivTest} function. It calculates the placebo coefficients and the absolute value of the mean of the placebo coefficients. It then calculates the critical value and p-values if an equivalence threshold is supplied for the test, according to Dette & Schumann (2024). If equivalence threshold is not supplied, it calculates the minimum equivalence threshold for which the null of non-negligible pre-trend differences can be rejected.
#'
#' @param data The data.frame object containing the data for the test. Should be of the form what is returned by the \link[EquiTrends]{EquiTrends_dataconstr} function.
#' @param equiv_threshold The equivalence threshold for the test. If NULL, the minimum equivalence threshold for which the null hypothesis can be rejected is calculated.
#' @param vcov The variance-covariance matrix estimator. See \link[EquiTrends]{meanEquivTest} for more information.
#' @param cluster The cluster variable for the cluster-robust variance-covariance matrix estimator. See \link[EquiTrends]{meanEquivTest} for more information.
#' @param alpha The significance level for the test. Only required if no equivalence threshold is supplied.
#' @param n The number of cross-sectional individuals in the data.
#' @param no_periods The number of periods in the data.
#' @param base_period The base period for the test. Must be one of the unique periods in the data.
#' @param is_panel_balanced A logical value indicating whether the panel data is balanced.
#'
#' @references 
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \doi{10.1080/07350015.2024.2308121}
#'
#' @return
#' #' An object of class "meanEquivTest" containing:
#' \item{\code{placebo_coefficients}}{A numeric vector of the estimated placebo coefficients,}
#' \item{\code{abs_mean_placebo_coefs}}{the absolute value of the mean of the placebo coefficients,}
#' \item{\code{var_mean_placebo_coef}}{the estimated variance of the mean placebo coefficient,}
#' \item{\code{significance_level}}{the significance level of the test,}
#' \item{\code{num_individuals}}{the number of cross-sectional individuals in the data,}
#' \item{\code{num_periods}}{the number of periods in the data,}
#' \item{\code{base_period}}{the base period in the data,}
#' \item{\code{num_observations}}{the total number of observations in the data,}
#' \item{\code{equiv_threshold_specified}}{a logical value indicating whether an equivalence threshold was specified.}
#' \item{\code{is_panel_balanced}}{a logical value indicating whether the panel data is balanced.}
#'
#' If \code{is.null(equiv_threshold)}, then additionally \code{minimum_equiv_threshold}: the minimum equivalence threshold for which the null hypothesis of non-negligible (based on the equivalence threshold) trend-differnces can be rejected. 
#' 
#' if \code{!(is.null(equiv_threshold))}, then additionally
#' \itemize{
#' \item \code{mean_critical_value}: the critical value at the alpha level,
#' \item \code{p_value}: the p-value of the test,
#' \item \code{reject_null_hypothesis}: A logical value indicating whether to reject the null hypothesis,
#' \item \code{equiv_threshold}: the equivalence threshold specified.
#' }
#' 
meanTest_func <- function(data, equiv_threshold, vcov, cluster, alpha, n, no_periods, base_period, is_panel_balanced){
  # Construct the formula for the plm() function
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm_formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Run the two-way fixed effects model:
  plm_twfe <- plm::plm(plm_formula, data=data, effect="twoways", 
                       model="within", index=c("ID","period"))
  
  # Extract the estimated coefficients:
  betas <- plm_twfe$coefficients
  placebo_names <- base::grep("placebo_",base::names(betas),value=TRUE)
  betas_placebo <- c(betas[placebo_names])
  
  # Calculating the mean of the absolute placebo coefficients:
  no_placebos <- length(placebo_names)
  one_vec <- rep(1, no_placebos)
  demean_vec <- one_vec/no_placebos
  mean_placebo <- abs(as.numeric(t(demean_vec)%*%betas_placebo))
 
  # Calculate the variance-covariance matrix:
  if(is.null(vcov)){
    vcov_mat <- plm_twfe$vcov
  } else if(!is.function(vcov) && vcov == "HC"){
    vcov_mat <- plm::vcovHC(plm_twfe, type="HC1", method = "white1")
  } else if(!is.function(vcov) && vcov == "HAC"){
    vcov_mat <- plm::vcovHC(plm_twfe, type="HC3", method = "arellano")
  } else if(!is.function(vcov) && vcov == "CL"){
    if(is.null(cluster)){
      vcov_mat <- clubSandwich::vcovCR(plm_twfe, cluster="ID", type="CR0")
    } else {
      vcov_mat <- clubSandwich::vcovCR(plm_twfe, cluster=data[,"cluster"], type="CR0")
    }
  } else if(is.function(vcov)) {
    vcov_mat <- vcov(plm_twfe)
    # Verifying it gives a square matrix:
    if(!is.matrix(vcov_mat) || nrow(vcov_mat) != ncol(vcov_mat)){
      stop("The vcov argument is invalid.")
    }
  } else {
    stop("The vcov argument is invalid.")
  }
  
  # The estimated variance of the mean placebo coefficient estimator:
  placebo_vcov <- vcov_mat[placebo_names, placebo_names]
  mean_placebo_var <- as.numeric(t(demean_vec)%*%placebo_vcov%*%demean_vec)
  
  # If delta is not NULL, we calculate the critical value / p-value. Else, the minimum delta 
  # for which H0 can be rejected.
  if(!is.null(equiv_threshold)){
    # The Critical value at the alpha level:
    crit_val <- VGAM::qfoldnorm(alpha, mean=equiv_threshold, sd=sqrt(mean_placebo_var))
    
    # The p-value
    p_value <- VGAM::pfoldnorm(mean_placebo, mean=equiv_threshold, sd=sqrt(mean_placebo_var))
    
    # reject H0:
    reject_H0 <- mean_placebo < crit_val
    
    results_list <- list(placebo_coefficients = betas_placebo, abs_mean_placebo_coefs = mean_placebo,
                         var_mean_placebo_coef = mean_placebo_var, 
                         mean_critical_value = crit_val, p_value = p_value, 
                         reject_null_hypothesis = reject_H0,
                         equiv_threshold = equiv_threshold,
                         significance_level = alpha,
                         base_period = base_period, 
                         equiv_threshold_specified = TRUE,
                         num_individuals = n, num_periods = no_periods, 
                         num_observations = nrow(data), is_panel_balanced = is_panel_balanced)
    
  } else {
    
    minimum_equiv_threshold <- meanTest_optim_func(mean_placebo, sqrt(mean_placebo_var), alpha)
    
    results_list <- list(placebo_coefficients = betas_placebo, abs_mean_placebo_coefs = mean_placebo,
                         var_mean_placebo_coef = mean_placebo_var,
                         minimum_equiv_threshold = minimum_equiv_threshold, significance_level = alpha,
                         base_period = base_period, 
                         equiv_threshold_specified = FALSE, 
                         num_individuals = n, num_periods = no_periods,
                         num_observations = nrow(data), is_panel_balanced = is_panel_balanced)
    
  }
  class(results_list) <- "meanEquivTest"
  return(results_list)
}

# ----------- Finding the Minimum Delta ----------------------------------------
# Function to minimize:
meanTest_obj_func <- function(coef, mean, sd, alpha){
  return(1e100*(VGAM::pfoldnorm(coef, mean, sd) - alpha)^2)
}


#' @title Finding the minimum equivalence threshold for the mean equivalence test
#' @description \code{meanTest_optim_func} solves the optimization problem to find the minimum equivalence threshold for which one can reject the null hypothesis of non-negligible pre-trend differences at a given significance level for the equivalence test based on the mean placebo coefficient. 
#'
#' @param coef The estimated absolute value of the mean placebo coefficients
#' @param sd The estimated standard deviation of the mean of the placebo coefficients
#' @param alpha The significance level
#'
#' @return
#' The minimum equivalence threshold for which the null hypothesis of non-negligible differences can be rejected for the equivalence test based on the mean placebo coefficient.
meanTest_optim_func <- function(coef, sd, alpha){
  obj_wrapper <- function(x) meanTest_obj_func(coef=coef, mean=x, sd=sd, alpha=alpha)
  
  
  result <- nloptr::nloptr(x0 = coef,
                           eval_f = obj_wrapper,
                           eval_grad_f = NULL,
                           lb = max(0, coef - 4*sd),
                           ub = coef + 4*sd,
                           eval_g_ineq = NULL,
                           eval_jac_g_ineq = NULL,
                           eval_g_eq = NULL,
                           eval_jac_g_eq = NULL,
                           opts = list("algorithm" = "NLOPT_LN_COBYLA",
                                       maxeval=20000000, xtol_rel = 1e-25))
  return(result$solution)
  
}