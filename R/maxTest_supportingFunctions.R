############ Supporting Functions of maxTest ###################################
# ----------- Intersection Union Approach --------------------------------------
#' An internal function of the EquiTrends Maximum Equivalence Testing procedure using the Intersection Union approach. 
#'
#' @description This is a supporting function of the \code{maxEquivTest} function. It calculates the placebo coefficients and the absolute value of the placebo coefficients. It then calculates the critical value and p-values if an equivalence threshold is supplied for the test, according to Dette & Schumann (2024). If no equivalence threshold is supplied, it calculates the minimum equivalence threshold for which the null of non-negligible pre-trend differences can be rejected.
#' 
#'
#' @param data The data.frame object containing the data for the test. Should be of the form what is returned by the \link[EquiTrends]{EquiTrends_dataconstr} function.
#' @param equiv_threshold The equivalence threshold for the test. If NULL, the minimum equivalence threshold for which the null hypothesis of non-negligible can be rejected is calculated.
#' @param vcov The variance-covariance matrix estimator. See \link[EquiTrends]{maxEquivTest} for more information.
#' @param cluster The cluster variable for the cluster-robust variance-covariance matrix estimator. See \link[EquiTrends]{maxEquivTest} for more information.
#' @param alpha The significance level for the test. 
#' @param n The number of cross-sectional individuals in the data.
#' @param no_periods The number of periods in the data.
#' @param base_period The base period for the test. Must be one of the unique periods in the data.
#' @param is_panel_balanced A logical value indicating whether the panel data is balanced.
#'
#' @references 
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \doi{10.1080/07350015.2024.2308121}
#'
#' @return
#' An object of class "maxEquivTestIU" containing:
#' \item{\code{placebo_coefficients}}{A numeric vector of the estimated placebo coefficients,}
#' \item{\code{abs_placebo_coefficients}}{a numeric vector with the absolute values of estimated placebo coefficients,}
#' \item{\code{placebo_coefficient_se}}{a numeric vector with the standard errors of the placebo coefficients,}
#' \item{\code{significance_level}}{the chosen significance level of the test,}
#' \item{\code{num_individuals}}{the number of cross-sectional individuals (n),}
#' \item{\code{num_periods}}{the number of periods (T),}
#' \item{\code{num_observations}}{the total number of observations (N),}
#' \item{\code{base_period}}{the base period in the data,}
#' \item{\code{placebo_names}}{the names corresponding to the placebo coefficients,}
#' \item{\code{equiv_threshold_specified}}{a logical value indicating whether an equivalence threshold was specified.}
#' \item{\code{is_panel_balanced}}{a logical value indicating whether the panel data is balanced.}
#'  Additionally, if \code{!(is.null(equiv_threshold))}\itemize{
#'  \item{\code{IU_critical_values}: a numeric vector with the individual critical values for each of the placebo coefficients,}
#'  \item{\code{reject_null_hypothesis}: a logical value indicating whether the null hypothesis of negligible pre-trend differences can be rejected at the specified significance level \code{alpha},}
#'  \item{\code{equiv_threshold}: the equivalence threshold employed,}
#' }
#' if \code{is.null(equiv_threshold)}\itemize{
#' \item{\code{minimum_equiv_thresholds}: a numeric vector including for each placebo coefficient the minimum equivalence threshold for which the null hypothesis of negligible pre-trend differences can be rejected for the corresponding placebo coefficient individually,}
#' \item{\code{minimum_equiv_threshold}: a numeric scalar minimum equivalence threshold for which the null hypothesis of negligible pre-trend differences can be rejected for all placebo coefficients individually.}
#' }
#'
maxTestIU_func <- function(data, equiv_threshold, vcov, cluster, alpha, n, no_periods, base_period, is_panel_balanced){
  # Construct the formula for the plm() function
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm_formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Run the two-way fixed effects model:
  IU_twfe <- plm::plm(plm_formula, data=data, effect="twoways", 
                      model="within", index=c("ID","period"))
  
  # Extract the estimated coefficients:
  betas <- IU_twfe$coefficients
  placebo_names <- base::grep("placebo_",base::names(betas),value=TRUE)
  betas_placebo <- c(betas[placebo_names])
  abs_betas_placebo <- abs(betas_placebo)
  
  # Extract the Variance-Covariance Matrix based on the user input:
  if(is.null(vcov)){
    vcov_mat <- IU_twfe$vcov
  } else if(!is.function(vcov) && vcov == "HC"){
    vcov_mat <- plm::vcovHC(IU_twfe, type="HC1", method = "white1")
  } else if(!is.function(vcov) && vcov == "HAC"){
    vcov_mat <- plm::vcovHC(IU_twfe, type="HC3", method = "arellano")
  } else if(!is.function(vcov) && vcov == "CL"){
    if(is.null(cluster)){
      vcov_mat <- clubSandwich::vcovCR(IU_twfe, cluster="ID", type="CR0")
    } else {
      vcov_mat <- clubSandwich::vcovCR(IU_twfe, cluster=data[,"cluster"], type="CR0")
    }
  } else {
    if(!is.matrix(vcov(IU_twfe)) || nrow(vcov(IU_twfe)) != ncol(vcov(IU_twfe))){
      stop("The vcov argument is invalid")
    }
    vcov_mat <- vcov(IU_twfe)
  }

  # Extracting the variances of the slope coefficients
  subvcov_mat <- as.matrix(vcov_mat[placebo_names, placebo_names, drop = FALSE])
  beta_var <- diag(subvcov_mat)
  
  # Calculating the standard errors
  beta_se <- sqrt(beta_var)
  
  if(!is.null(equiv_threshold)){
    
    # Critical Values method
    crit_values <- c(base::sapply(beta_se, 
                                  FUN = function(x){return(VGAM::qfoldnorm(alpha, mean=equiv_threshold, sd=x))}))
    reject_vec <- abs(abs_betas_placebo) < crit_values

    reject_H0 <- all(reject_vec==TRUE)
    
    results_list <- structure(list(placebo_coefficients = betas_placebo, 
                                   abs_placebo_coefficients = abs_betas_placebo, 
                                   placebo_coefficients_se = beta_se, 
                                   IU_critical_values = crit_values,
                                   reject_null_hypothesis = reject_H0, equiv_threshold = equiv_threshold,
                                   significance_level = alpha,
                                   base_period = base_period, placebo_coef_names = placebo_names,
                                   equiv_threshold_specified = TRUE,
                                   num_individuals = n,
                                   num_periods = no_periods, num_observations = nrow(data),
                                   is_panel_balanced = is_panel_balanced), class = "maxEquivTestIU")
    
    return(results_list)
  } else{
    
    # Calculating the minimum deltas
    minimum_equiv_thresholds <- base::mapply(FUN = function(x,y){maxTestIU_optim_func(coef=x, sd=y, alpha=alpha)}, abs_betas_placebo, beta_se)
    
    # Then, the minimum delta is the maximum value over all these values:
    minimum_equiv_threshold <- max(minimum_equiv_thresholds)
    
    results_list <- structure(list(placebo_coefficients = betas_placebo,
                                   abs_placebo_coefficients = abs_betas_placebo, 
                                   placebo_coefficients_se = beta_se,
                                   minimum_equiv_threshold = minimum_equiv_threshold,
                                   minimum_equiv_thresholds = minimum_equiv_thresholds,
                                   significance_level = alpha,
                                   base_period = base_period,
                                   placebo_coef_names = placebo_names,
                                   equiv_threshold_specified = FALSE,
                                   num_individuals = n, num_periods = no_periods,
                                   num_observations = nrow(data),
                                   is_panel_balanced = is_panel_balanced), class = "maxEquivTestIU")
    
    return(results_list)
  }
}


# Functions to help find the minimum delta such that the p-value is alpha
maxTestIU_obj_func <- function(coef, mean, sd, alpha){
  return(1e20*(VGAM::pfoldnorm(coef, mean, sd) - alpha)^2)
}

#' @title Finding the minimum equivalence threshold for the equivalence test based on the IU procedure for the maximum placebo coefficient.
#' @description \code{maxTestIU_optim_func} solves the optimization problem to find the minimum equivalence threshold for which one can reject the null hypothesis of non-negligible pre-trend differences at a given significance level for the equivalence test based on the maximum placebo coefficient, especially for the Intersection Union type. 
#'
#' @param coef The estimated absolute value of the mean placebo coefficients
#' @param sd The estimated standard deviation of the mean of the placebo coefficients
#' @param alpha The significance level
#'
#' @return
#' The minimum equivalence threshold for which the null hypothesis of non-negligible differences can be rejected for the equivalence test based on the mean placebo coefficient.

maxTestIU_optim_func <- function(coef, sd, alpha){
  obj_wrapper <- function(x) maxTestIU_obj_func(coef=coef, mean=x, sd=sd, alpha=alpha)
  
  result <- stats::optimize(obj_wrapper, c(max(0, coef - 10*sd), 10*sd), tol = 1e-20)
  return(result$minimum)
  
}


# ----------- The Bootstrap Approaches -----------------------------------------
#' An internal function of the EquiTrends Maximum Equivalence Testing procedure using the Bootstrap approaches.
#'
#' @description This is a supporting function of the \code{maxEquivTest} function. It calculates the placebo coefficients and the absolute value of the placebo coefficients. It then calculates the critical value by bootstrap if an equivalence threshold is supplied for the test, according to Dette & Schumann (2024). 
#'
#' @param data The data.frame object containing the data for the test. Should be of the form what is returned by the \link[EquiTrends]{EquiTrends_dataconstr} function.
#' @param equiv_threshold The equivalence threshold for the test. 
#' @param alpha The significance level for the test.
#' @param n The number of cross-sectional individuals in the data.
#' @param B The number of bootstrap replications.
#' @param no_periods The number of periods in the data.
#' @param base_period The base period for the test. Must be one of the unique periods in the data.
#' @param type The type of bootstrap to be used. Must be one of "Boot" or "Wild".
#' @param original_names The original names of the control variables in the data.
#' @param is_panel_balanced A logical value indicating whether the panel data is balanced.
#'
#' @references
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \doi{10.1080/07350015.2024.2308121}
#'
#' @return
#' an object of class "maxEquivTestBoot"  with
#'  \item{\code{placebo_coefficients}}{A numeric vector of the estimated placebo coefficients,}
#'  \item{\code{abs_placebo_coefficients}}{a numeric vector with the absolute values of estimated placebo coefficients,}
#'  \item{\code{max_abs_coefficient}}{the maximum absolute estimated placebo coefficient,}
#'  \item{\code{bootstrap_critica_value}}{the by bootstrap found critical value for the equivalence test based on the maximum absolute placebo coefficient,}
#'  \item{\code{reject_null_hypothesis}}{a logical value indicating whether the null hypothesis of negligible pre-trend differences can be rejected at the specified significance level \code{alpha},}
#'  \item{\code{B}}{the number of bootstrap samples used to find the critical value,}
#'  \item{\code{significance_level}}{the chosen significance level of the test \code{alpha},}
#'  \item{\code{num_individuals}}{the number of cross-sectional individuals (n),}
#'  \item{\code{num_periods}}{the number of periods (T),}
#'  \item{\code{num_observations}}{the total number of observations (N),}
#'  \item{\code{base_period}}{the base period in the data,}
#'  \item{\code{placebo_names}}{the names corresponding to the placebo coefficients,}
#'  \item{\code{equiv_threshold_specified}}{a logical value indicating whether an equivalence threshold was specified.}
#'  \item{\code{is_panel_balanced}}{a logical value indicating whether the panel data is balanced.}
#'
maxTestBoot_func <- function(data, equiv_threshold, alpha, n, B, no_periods, 
                                 base_period, type, original_names, is_panel_balanced){
  D <- as.matrix(stats::model.matrix(~0+factor(period), data=data))
  # Between transformation on D to obtain WD:
  WD <- matrix_between_transformation(D, matrix(data$ID))
  
  # Check for possible multicolinearity:
  if(qr(WD)$rank < ncol(WD)){
    WD <- remove_multicollinearity(WD, asmatrix = TRUE)$df
  }

  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  X_df <- data[, c(placebo_names, X_names)]
  if(length(c(placebo_names, X_names)) == 1){
    X_df <- as.data.frame(X_df)
  }
  colnames(X_df) <- c(placebo_names, original_names)
  
  # Create a model.matrix:
  model_data <- stats::model.matrix(~.+0, data = X_df)
  model_matrix <- as.matrix(model_data)
  
  # Double demean the data:
  X <- double_demean(x = model_matrix, individual = matrix(data$ID), time = matrix(data$period), WD = WD)

  # Check multicolinearity:
  if(qr(X)$rank < ncol(X)){
    new_X <- remove_multicollinearity(X)
    X <- new_X$df
    removed_ind <- new_X$problematic_vars
    if(any(removed_ind %in% 1:length(placebo_names))){
      removed_placebos <- removed_ind[removed_ind %in% 1:length(placebo_names)]
      period_removed_placebo <- placebo_names[removed_placebos]
      period_removed_placebo <- sub("placebo_", "", period_removed_placebo)
      warning(paste("The placebo corresponding to period(s) ", paste(period_removed_placebo, collapse = ", "), " removed due to multicolinearity."))
      removed_ind <- removed_ind[!(removed_ind %in% 1:length(placebo_names))]
      placebo_names <- placebo_names[-removed_placebos]
    }
    removed_names <- colnames(model_matrix)[removed_ind]
    warning(paste("The following control variables were removed due to multicolinearity: ", removed_names))
  }
  
  Y <- as.matrix(double_demean(matrix(data$Y),individual = matrix(data$ID), time = matrix(data$period), WD = WD))

  # The unconstrained coefficient is:
  unconstrained_coefs <- solve(t(X)%*%X, t(X)%*%Y)
  
  placebo_coefs <- unconstrained_coefs[1:length(placebo_names)]
  names(placebo_coefs) <- placebo_names
  
  
  # its maximum absolute entry is:
  max_unconstr_coef <- max(abs(unconstrained_coefs[1:length(placebo_names)]))
  
  if(!is.null(equiv_threshold)){

    if(max_unconstr_coef >= equiv_threshold){
      constrained_coefs <- unconstrained_coefs
    } else {
      constrained_coefs <- boot_optimization_function(x=X, y=Y, no_placebos = length(placebo_names), 
                                                      equiv_threshold = equiv_threshold,
                                                      start_val = unconstrained_coefs)
    }
    
    if(type == "Boot"){

      # Calculate the variance based on the constrained coefficient:
      resid_variance <- sigma_hathat_c(parameter = constrained_coefs,  
                                       x=X, y=Y,
                                       time = data$period,
                                       ID = data$ID)
  
      # Run the Bootstrap:
      bootstrap_maxcoefs <- maxTestBoot_bootstrap(Xb = X%*%constrained_coefs, X=X, B=B,
                                                  variance = resid_variance, ID = data$ID,
                                                  period = data$period, WD=WD, 
                                                  no_placebos = length(placebo_names))
  
    } else {
      u_ddot <- Y - X%*%constrained_coefs
      
      # Run the Wild Bootstrap:
      bootstrap_maxcoefs <- maxTestBoot_wildbootstrap(Xb = X%*%constrained_coefs, X=X, B=B,
                                                       u_ddot = u_ddot,
                                                       ID = data$ID, period = data$period,
                                                       no_placebos = length(placebo_names),
                                                       WD = WD)
    }
    
    # Find the critical value at the alpha level:
    boot_crit_value <- unname(stats::quantile(bootstrap_maxcoefs, probs = alpha))
    
    # Reject Or Not:
    reject_H0 <- (max(abs(unconstrained_coefs[1:length(placebo_names)])) < boot_crit_value)
  
    results_list <- structure(list(placebo_coefficients = placebo_coefs,
                                   abs_placebo_coefficients = abs(placebo_coefs),
                                   max_abs_coefficient =max(abs(unconstrained_coefs[1:length(placebo_names)])),
                                   bootstrap_critical_value = boot_crit_value,
                                   reject_null_hypothesis = reject_H0, 
                                   equiv_threshold = equiv_threshold,
                                   B = B, 
                                   significance_level = alpha, wild = (type=="Wild"),
                                   base_period = base_period,
                                   equiv_threshold_specified = !is.null(equiv_threshold),
                                   num_individuals = n, num_periods = no_periods, num_observations = nrow(data),
                                   is_panel_balanced = is_panel_balanced), class = "maxEquivTestBoot")
   } else {
    # Find the maximum placebo variance
     placebo_variance <- sigma_hathat_c(parameter = unconstrained_coefs, x=X, y=Y, ID = data$ID, time = data$period)
     vcov_mat <- as.matrix(solve(t(X)%*%X)*placebo_variance)
     vcov_mat_placebos <- as.matrix(vcov_mat[1:length(placebo_names), 1:length(placebo_names)])
     max_sd <- as.numeric(max(sqrt(diag(vcov_mat_placebos))))
     
    # Find the minimum delta for which the null hypothesis can be rejected:
    min_equiv_thr <- min_delta(data = data, equiv_threshold = equiv_threshold, alpha = alpha, n = n, B = B,
                               no_periods = no_periods, base_period = base_period, type = type,
                               original_names = original_names, is_panel_balanced = is_panel_balanced,
                               max_abs_coef = max_unconstr_coef, max_sd = max_sd)


    results_list <- structure(list(placebo_coefficients = placebo_coefs,
                                   abs_placebo_coefficients = abs(placebo_coefs),
                                   max_abs_coefficient = max(abs(unconstrained_coefs[1:length(placebo_names)])),
                                   minimum_equiv_threshold = min_equiv_thr,
                                   significance_level = alpha,
                                   B = B, wild = (type=="Wild"),
                                   base_period = base_period,
                                   equiv_threshold_specified = !is.null(equiv_threshold),
                                   num_individuals = n, num_periods = no_periods, num_observations = nrow(data),
                                   is_panel_balanced = is_panel_balanced), class = "maxEquivTestBoot")
  }
  return(results_list)
}


#   ---- Supporting functions for the Bootstrap Approach ----

# Removing multicolinear columns:
remove_multicollinearity <- function(df, asmatrix = FALSE) {
  # Create the design matrix
  mat <- as.matrix(df)
  
  # Perform QR decomposition
  qr_mat <- qr(mat)
  
  # Identify the problematic variables
  problematic_vars <- qr_mat$pivot[seq(from = qr_mat$rank + 1, to = ncol(mat))]
  
  df <- df[,-problematic_vars]
  if(asmatrix){
    df <- as.matrix(df)
  }
  return(list(df = df, problematic_vars = problematic_vars))
}



# The objective function of the constrained estimator:
boot_objective_function <- function(parameter, x, y, no_placebos, equiv_threshold){
  Xb <- as.numeric(x%*%parameter)
  MSE <- mean((y-Xb)^2)
  return(MSE)
}


# The constraint function for the optimization. Note that we only penalize the 
# coefficients corresponding to the placebos!
boot_constraint_function <- function(parameter, x, y, no_placebos, equiv_threshold){
  value <- -max(abs(parameter[1:no_placebos])) + equiv_threshold
  return(value)
}


#' @title Finding the restricted placebo coefficients for the maximum equivalence test based on the bootstrap approaches
#' @description \code{boot_optimization_function} solves the optimization problem to find the restricted placebo coefficients, according to Dette & Schumann (2024). 
#'
#' @param x The double demeaned independent variables.
#' @param y The double demeaned dependent variable.
#' @param no_placebos The number of placebo coefficients.
#' @param equiv_threshold The equivalence threshold for the test.
#' @param start_val The starting values for the optimization.
#'
#' @references 
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \doi{10.1080/07350015.2024.2308121}
#'
#' @return
#' A numeric vector containing the restricted placebo coefficients
boot_optimization_function <- function(x, y, no_placebos, equiv_threshold, start_val){
  constrained_optimum <- nloptr::nloptr(x0 = start_val,
                                        eval_f = boot_objective_function,
                                        eval_grad_f = NULL,
                                        lb = rep(-Inf, length(start_val)),
                                        ub = rep(Inf, length(start_val)),
                                        eval_g_ineq = boot_constraint_function,
                                        eval_jac_g_ineq = NULL,
                                        eval_g_eq = NULL,
                                        eval_jac_g_eq = NULL,
                                        opts = list("algorithm" = "NLOPT_LN_COBYLA", 
                                                    maxeval=2000000, xtol_rel = 1e-6),
                                        x = x, y=y, no_placebos = no_placebos, equiv_threshold=equiv_threshold)

  return(constrained_optimum$solution)
}

# Constrained Variance:
#' Calculating the constrained variance of the residuals for the Boostrap approaches in the EquiTrends Maximum Equivalence Testing procedure, according to Dette & Schumann (2024).
#'
#' @param parameter The constrained coefficients.
#' @param x The double demeaned independent variables.
#' @param y The double demeaned dependent variable.
#' @param ID The ID variable.
#' @param time The time variable.
#'
#' @references
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \doi{10.1080/07350015.2024.2308121}
#'
#' @return
#' The estimated constrained variance of the residuals.
sigma_hathat_c <- function(parameter, x, y, ID, time){
  Xb <- as.numeric(x%*%parameter)
  
  residuals_demeaned <- y - Xb
  
  N <- length(ID)
  n <- length(unique(ID))
  no_periods <- length(unique(time))
  p <- ncol(x)
  
  df <- N - p - n - no_periods + 1
  
  c_sigma_hathat <- as.numeric(t(residuals_demeaned)%*%residuals_demeaned)/df
  
  return(c_sigma_hathat)
}

# Minimum equivalence threshold for the Bootstrap approaches:

min_delta <- function(data, equiv_threshold, alpha, n, B, no_periods, 
                      base_period, type, original_names, is_panel_balanced, 
                      max_abs_coef, max_sd){
  
  # Using a wrapper function that 
  wrapper_func <- function(x){
    bootstrapTest_result <- maxTestBoot_func(data = data, equiv_threshold = x, alpha = alpha, n = n, B = B, 
                                             no_periods = no_periods, base_period = base_period, type = type, 
                                             original_names = original_names, is_panel_balanced = is_panel_balanced)
    value <- ifelse(bootstrapTest_result$reject_null_hypothesis, -exp(-x), exp(-x))
    return(value)
  }
  
  min_equiv_threshold <- stats::optimize(f = wrapper_func, interval = c(max_abs_coef, max_abs_coef + 15*max_sd))$minimum
  
  return(min_equiv_threshold)
  
}


# --- maxTest Error Function ---
#' @title Additional input checks for the maxEquivTest function
#' 
#' @description This function checks additonal inputs specific to the maxEquivTest function. 
#'
#' @param type the type of test for the maximum absolute placebo coefficient to be conducted; must be one of "IU", "Boot" or "Wild".
#' @param equiv_threshold the equivalence threshold for the test. Must be a numeric scalar or NULL.
#' @param vcov the variance-covariance matrix estimator. See \code{\link[EquiTrends]{maxEquivTest}} for more information.
#' @param B the number of bootstrap iterations. Must be a numeric integer scalar.
#'
#' @return
#' A list with two elements: \code{error} a logical value indicating whether an error was found, and \code{message} a character string with the error message. If no error was found, \code{error} is \code{FALSE} and \code{message} is empty.
maxTest_error <- function(type, equiv_threshold, vcov, B){
  
  if(length(type)!=1 && !identical(type, c("IU", "Boot", "Wild"))){
    return(list(error=TRUE, message = "type is not valid"))
  }
  
  if(length(type)==1 && !(type %in% c("IU", "Boot", "Wild"))){
    return(list(error=TRUE, message = "type is not valid"))
  }
  
  if(type != "IU"){
    if(!is.numeric(B) || B<= 0 || B != round(B) || length(B) != 1){
      return(list(error=TRUE, message = "B must be a strictly positive integer scalar"))
    }
  }
  
  # If type = IU, vcov must be correctly specified:
  if(type == "IU" && !is.null(vcov)){
    if(!is.function(vcov) && !(vcov %in% c("HC", "HAC", "CL"))){
      return(list(error=TRUE, message = "vcov is not valid"))
    }
  }
  
  
  return(list(error=FALSE))
}










