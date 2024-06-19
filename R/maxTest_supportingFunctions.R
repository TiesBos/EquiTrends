############ Supporting Functions of maxTest ###################################
# ----------- Intersection Union Approach --------------------------------------
maxTestIU_func <- function(data, equiv_threshold, vcov, cluster, alpha, n, no_periods, base_period){
  # Construct the formula for the plm() function
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm_formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Run the two-way fixed effects model:
  IU_twfe <- plm::plm(plm_formula, data=data, effect="twoways", 
                      model="within", index=c("ID","period"))
  
  # Extract the estimated coefficients:
  betas <- IU_twfe$coefficients
  betas_placebo <- c(betas[placebo_names])
  abs_betas_placebo <- abs(betas_placebo)
  
  # Extract the Variance-Covariance Matrix based on the user input:
  if(is.null(vcov)){
    vcov_mat <- IU_twfe$vcov
  } else if(vcov == "HC"){
    vcov_mat <- sandwich::vcovHC(IU_twfe, type="HC1", method = "white1")
  } else if(vcov == "HAC"){
    vcov_mat <- sandwich::vcovHC(IU_twfe, type="HC3", method = "arellano")
  } else if(vcov == "CL"){
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
  subvcov_mat <- vcov_mat[placebo_names, placebo_names]
  beta_var <- diag(subvcov_mat)
  
  # Calculating the standard errors
  beta_se <- as.matrix(sqrt(beta_var/n))
  
  if(!is.null(equiv_threshold)){
    
    # Critical Values method
    crit_values <- c(base::sapply(beta_se, 
                                  FUN = function(x){return(VGAM::qfoldnorm(alpha, mean=equiv_threshold, sd=x))}))
    reject_vec <- abs(abs_betas_placebo) < crit_values
    
    # P-values calculated
    p_values <- c(base::mapply(FUN = function(x, y){return(VGAM::pfoldnorm(x, mean=equiv_threshold, sd=y))},
                               abs_betas_placebo, beta_se))
    
    reject_H0 <- all(reject_vec==TRUE)
    
    results_list <- structure(list(placebo_coefficients = betas_placebo, 
                                   abs_placebo_coefficients = abs_betas_placebo, 
                                   placebo_coefficients_se = beta_se, 
                                   IU_critical_values = crit_values, p_values = p_values,
                                   reject_null_hypothesis = reject_H0, equiv_threshold = equiv_threshold,
                                   significance_level = alpha, num_individuals = n,
                                   num_periods = no_periods, 
                                   base_period = base_period, placebo_coef_names = placebo_names,
                                   equiv_threshold_specified = TRUE), class = "maxEquivTestIU")
    
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
                                   num_individuals = n,
                                   num_periods = no_periods, base_period = base_period,
                                   placebo_coef_names = placebo_names,
                                   equiv_threshold_specified = FALSE), class = "maxEquivTestIU")
    
    return(results_list)
  }
}

# Functions to help find the minimum delta such that the p-value is alpha
maxTestIU_obj_func <- function(coef, mean, sd, alpha){
  return(1e20*(VGAM::pfoldnorm(coef, mean, sd) - alpha)^2)
}

maxTestIU_optim_func <- function(coef, sd, alpha){
  obj_wrapper <- function(x) maxTestIU_obj_func(coef=coef, mean=x, sd=sd, alpha=alpha)
  
  
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

# ----------- The Bootstrap Approaches -----------------------------------------
maxTestBoot_func <- function(data, equiv_threshold, alpha, n, B, no_periods, 
                                 base_period, type){
  # Obtain the double demeaned data:
  dd_data <- double_demean(data)
  
  # find the double-demeaned independent and dependent variable:
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  X <- as.matrix(dd_data[, c(placebo_names, X_names)])
  Y <- as.matrix(dd_data$Y)
  
  # The unconstrained coefficient is:
  unconstrained_coefs <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  # its maximum absolute entry is:
  max_unconstr_coef <- max(abs(unconstrained_coefs[1:length(placebo_names)]))
  
  if(max_unconstr_coef >= equiv_threshold){
    constrained_coefs <- unconstrained_coefs
  } else {
    constrained_coefs <- boot_optimization_function(x=X, y=Y, no_placebos = length(placebo_names), 
                                                    equiv_threshold = equiv_threshold,
                                                    start_val = unconstrained_coefs)$solution
  }
  
  if(type == "Boot"){
    # Calculate the variance based on the constrained coefficient:
    resid_variance <- sigma_hathat_c(parameter = constrained_coefs, N = n, 
                                     x=X, y=Y,
                                     no_periods = no_periods)
    # Run the Bootstrap:
    bootstrap_maxcoefs <- maxTestBoot_bootstrap(Xb = X%*%constrained_coefs, X=X, B=B,
                                                variance = resid_variance, ID = dd_data$ID,
                                                period = dd_data$period, no_placebos = length(placebo_names))
  } else {
    u_ddot <- Y - X%*%constrained_coefs
    
    # Run the Wild Bootstrap:
    bootstrap_maxcoefs <- maxTestBoot_wildbootstrap(Xb = X%*%constrained_coefs, X=X, B=B,
                                                     u_ddot = u_ddot,
                                                     ID = dd_data$ID, period = dd_data$period,
                                                     no_placebos = length(placebo_names))
  }
  
  # Find the critical value at the alpha level:
  boot_crit_value <- stats::quantile(bootstrap_maxcoefs, probs = alpha)
  
  # Reject Or Not:
  reject_H0 <- (max(abs(unconstrained_coefs[1:length(placebo_names)])) < boot_crit_value)
  
  results_list <- structure(list(placebo_coefficients = unconstrained_coefs[1:length(placebo_names)],
                                 abs_placebo_coefficients = abs(unconstrained_coefs[1:length(placebo_names)]),
                                 max_abs_coefficient =max(abs(unconstrained_coefs[1:length(placebo_names)])),
                                 bootstrap_citical_value = boot_crit_value,
                                 reject_null_hypothesis = reject_H0, 
                                 equiv_threshold = equiv_threshold,
                                 B = B, 
                                 significance_level = alpha, wild = (type=="Wild"),
                                 num_individuals = n, num_periods = no_periods, 
                                 base_period = base_period), class = "maxEquivTestBoot")
  
  return(results_list)
}


#   ---- Supporting functions for the Bootstrap Approach ----
# Calculating Double Demeaned Data:
#' @title Double Demeaning a data.frame object
#'
#' @param df The data.frame object to double-demean. It should contain a column with names "ID" and "period" representing the unit and time period, respectively, over which the data is double-demeaned.
#'
#' @return a data.frame object with the double demeaned data.
#'
#' @importFrom dplyr %>% group_by_at summarise across all_of left_join
double_demean <- function(df){
  # Finding all variables that need to be double demeaned
  placebo_names <- grep("placebo_", names(df), value=TRUE)
  X_names <- grep("X_", names(df), value=TRUE)
  names <- c(placebo_names, X_names)
  
  # Matrix:
  original_data <- df[, c("Y", names)]
  
  # Calculate the mean for each unit (ID) for each value column
  unit_means <- as.data.frame(df %>%
                                group_by_at("ID") %>%
                                summarise(across(all_of(c("Y", names)), mean, .names = "unit_means_{col}")))
  
  # Calculate the mean for each time period for each value column
  time_means <-  as.data.frame(df %>%
                                 group_by_at("period") %>%
                                 summarise(across(all_of(c("Y", names)), mean, .names = "time_means_{col}")))
  
  # Merge the calculated means back to the original dataframe 
  new_df <- df
  new_df <- new_df %>%
    left_join(unit_means, by = "ID") %>%
    left_join(time_means, by = "period")
  
  # Get the correct dimensions for W of the time means, unit means and overall means
  unit_mean_names <- grep("unit_means_", names(new_df), value=TRUE)
  time_means_names <- grep("time_means_", names(new_df), value=TRUE)
  unit_means_mat <- new_df[,unit_mean_names]
  time_means_mat <- new_df[,time_means_names]
  
  over_all_means <- colMeans(as.matrix(original_data))
  n_row <- nrow(as.matrix(original_data))
  overall_means <- as.data.frame(matrix(rep(over_all_means, times = n_row), nrow = n_row, byrow = TRUE))
  
  demeaned_data <- original_data - unit_means_mat - time_means_mat + overall_means
  ID <- df$ID
  period <- df$period
  demeaned_data <- cbind(ID, period, demeaned_data)
  colnames(demeaned_data) <- c("ID", "period", "Y", names)
  return(demeaned_data)
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


# Optimization function:
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
                                                    maxeval=2000000, xtol_rel = 1e-06),
                                        x = x, y=y, no_placebos = no_placebos, equiv_threshold=equiv_threshold)
  return(constrained_optimum)
}

# Constrained Variance:
sigma_hathat_c <- function(parameter, x, y, N, no_periods){
  Xb <- as.numeric(x%*%parameter)
  
  residuals_demeaned <- y - Xb
  
  c_sigma_hathat <- sum((residuals_demeaned^2))/((N-1)*no_periods)
  return(c_sigma_hathat)
}



# --- maxTest Error Function ---
maxTest_error <- function(type, equiv_threshold, vcov){
  
  if(length(type)!=1 && !identical(type, c("IU", "Boot", "Wild"))){
    return(list(error=TRUE, message = "type is not valid"))
  }
  
  if(length(type)==1 && !(type %in% c("IU", "Boot", "Wild"))){
    return(list(error=TRUE, message = "type is not valid"))
  }
  
  # If type = IU, vcov must be correctly specified:
  if(type == "IU" && !is.null(vcov)){
    if(!(vcov %in% c("HC", "HAC", "CL")) && !is.function(vcov)){
      return(list(error=TRUE, message = "vcov is not valid"))
    }
  }
  
  # if type != IU, delta cannot be NULL
  if(type != "IU" && is.null(equiv_threshold)){
    return(list(error=TRUE, message = "delta must be specified for types Boot and Wild."))
  }
  
  return(list(error=FALSE))
}










