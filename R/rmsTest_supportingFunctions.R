# ---- Error Function RMS ----
rmsTest_error <- function(alpha, no_lambda){
  
  if(!alpha %in% c(0.01, 0.025, 0.05, 0.1, 0.2)){
    return(list(error=TRUE, message="alpha must be one of 0.01, 0.025, 0.05, 0.1 or 0.2"))
  }
  
  if(!is.numeric(no_lambda) || no_lambda <= 0 || no_lambda != round(no_lambda)){
    return(list(error=TRUE, message="no_lambda must be a positive integer"))
  }
  
  return(list(error = FALSE))
}



# ----------- The Testing Procedure --------------------------------------------
rmsTest_func <- function(data, equiv_threshold, alpha, no_lambda, base_period){
  # Formula for plm function:
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm_formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Number of individuals in the sample
  individuals <- unique(data[,"ID"])
  N <- length(individuals)
  
  # Number of periods:
  no_periods <- length(unique(data[,"period"]))
  
  # Matrix storing the placebo coefficient RMS on 1/lambda of the data for 
  # lambda = 1,..., no.lambda:
  MS_placebo_lambda <- rep(NA, no_lambda)
  for(i in 1:no_lambda){
    lambda <- i/no_lambda
    # Drawing the subset of individuals of size lambda*N
    sub_N <- floor(lambda*N)
    subset_index <- sample(1:N, sub_N)
    subset_indiv <- individuals[subset_index]
    
    subset_data <- data[data[,"ID"] %in% subset_indiv,]
    
    # Calculating the TWFE estimators for this subset of data:
    all_coefs <- plm::plm(plm_formula, data=subset_data, effect="twoways", 
                          model="within", index=c("ID","period"))$coefficients
    
    # Calculate the mean squared error for the placebo estimates:
    placebo_beta_lambda <- all_coefs[placebo_names]
    placebo_mean_sqrd <- base::mean((placebo_beta_lambda)^2)
    # Store the mean squared placebo coefficient estimate:
    MS_placebo_lambda[i] <- placebo_mean_sqrd
    
  }
  
  betas_placebo <- placebo_beta_lambda
  
  RMS_placebo <- sqrt(MS_placebo_lambda[no_lambda])
  MS_placebo <- MS_placebo_lambda[no_lambda]
  
  # Calculating \hat{V}_n:
  diff_vec <- MS_placebo_lambda[1:(no_lambda-1)] - MS_placebo_lambda[no_lambda]
  V_n <- sqrt(mean(diff_vec^2))
  
  # Calculating the (Bootstrap) critical value of the W distribution:
  Q_W <- W_critical_value(alpha)
  
  if(!is.null(equiv_threshold)){
    # Critical value for RMS test
    MS_critical_value <- as.numeric(equiv_threshold^2 + Q_W*V_n)
    RMS_critical_value <- sqrt(MS_critical_value)
    
    # Reject or not:
    reject_H0 <- MS_placebo_lambda[no_lambda] < MS_critical_value
    
    results_list <- structure(list(placebo_coefficients = betas_placebo, rms_placebo_coefficients = RMS_placebo,
                                   rms_critical_value = RMS_critical_value,
                                   reject_null_hypothesis = reject_H0, equiv_threshold = equiv_threshold,
                                   significance_level = alpha,
                                   num_individuals = N, num_periods = no_periods, 
                                   base_period = base_period,
                                   equiv_threshold_specified = TRUE),
                              class = "rmsEquivTest")
  } else {
    min_equiv_threshold <- as.numeric(sqrt(MS_placebo-Q_W*V_n))
    results_list <- structure(list(placebo_coefficients = betas_placebo, rms_placebo_coefficients = RMS_placebo,
                                   minimum_equiv_threshold = min_equiv_threshold, significance_level = alpha,
                                   num_individuals = N, num_periods = no_periods, base_period = base_period,
                                   equiv_threshold_specified = FALSE),
                              class = "rmsEquivTest")
  }
  return(results_list)
}


# ---- Critical Value RMS ----
# Function to get critical value for a given significance level
W_critical_value <- function(significance_level) {
  crit_val_matrix <- matrix(c(
    0.010, -4.2329959,
    0.025, -2.9047318,
    0.050, -2.1431720,
    0.100, -1.4601327,
    0.200, -0.8561188,
    0.800,  0.8344549,
    0.900,  1.4928501,
    0.950,  2.2003307,
    0.975,  3.0065821,
    0.990,  4.1150156
  ), ncol = 2, byrow = TRUE)
  
  crit_val_row <- which(crit_val_matrix[,1] == significance_level)
  crit_val <- crit_val_matrix[crit_val_row, 2]
  return(crit_val)
}