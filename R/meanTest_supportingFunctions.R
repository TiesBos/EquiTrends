# ----------- The Mean Test Function -------------------------------------------
meanTest.func <- function(data, equiv_threshold, vcov, cluster, alpha, n, no_periods, base_period){
  # Construct the formula for the plm() function
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm_formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Run the two-way fixed effects model:
  plm_twfe <- plm::plm(plm_formula, data=data, effect="twoways", 
                       model="within", index=c("ID","period"))
  
  # Extract the estimated coefficients:
  betas <- plm_twfe$coefficients
  betas_placebo <- c(betas[placebo_names])
  
  # Calculating the mean of the absolute placebo coefficients:
  no_placebos <- length(placebo_names)
  one_vec <- rep(1, no_placebos)
  demean_vec <- one_vec/no_placebos
  mean_placebo <- abs(as.numeric(t(demean_vec)%*%betas_placebo))
 
  # Calculate the variance-covariance matrix:
  if(is.null(vcov)){
    vcov_mat <- plm_twfe$vcov
  } else if(vcov == "HC"){
    vcov_mat <- sandwich::vcovHC(plm_twfe, type="HC1", method = "white1")
  } else if(vcov == "HAC"){
    vcov_mat <- sandwich::vcovHC(plm_twfe, type="HC3", method = "arellano")
  } else if(vcov == "CL"){
    if(is.null(cluster)){
      vcov_mat <- clubSandwich::vcovCR(plm_twfe, cluster="ID", type="CR0")
    } else {
      vcov_mat <- clubSandwich::vcovCR(plm_twfe, cluster=data[,"cluster"], type="CR0")
    }
  } else {
    vcov_mat <- vcov(plm_twfe)
    # Verifying it gives a square matrix:
    if(!is.matrix(vcov_mat) || nrow(vcov_mat) != ncol(vcov_mat)){
      stop("The vcov argument is invalid.")
    }
  }
  
  # The estimated variance of the mean placebo coefficient estimator:
  placebo_vcov <- vcov_mat[placebo_names, placebo_names]
  sum_variance <- as.numeric(t(one_vec)%*%placebo_vcov%*%one_vec)
  mean_placebo_var <- sum_variance/n
  
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
                         num_individuals = n,
                         num_periods = no_periods, base_period = base_period, 
                         equiv_threshold_specified = TRUE)
    
  } else {
    
    minimum_equiv_threshold <- meanTest_optim_func(mean_placebo, sqrt(mean_placebo_var), 0.05)
    
    results_list <- list(placebo_coefficients = betas_placebo, abs_mean_placebo_coefs = mean_placebo,
                         var_mean_placebo_coef = mean_placebo_var,
                         minimum_equiv_threshold = minimum_equiv_threshold, significance_level = alpha,
                         num_individuals = n,
                         num_periods = no_periods, base_period = base_period, 
                         equiv_threshold_specified = FALSE)
    
  }
  class(results_list) <- "meanEquivTest"
  return(results_list)
}

# ----------- Finding the Minimum Delta ----------------------------------------
# Function to minimize:
meanTest_obj_func <- function(coef, mean, sd, alpha){
  return(1e100*(VGAM::pfoldnorm(coef, mean, sd) - alpha)^2)
}


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