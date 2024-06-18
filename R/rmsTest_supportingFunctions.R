# ---- Error Function RMS ----
rmsTest.error <- function(alpha, no.lambda){
  
  if(!alpha %in% c(0.01, 0.025, 0.05, 0.1, 0.2)){
    return(list(error=TRUE, message="alpha must be one of 0.01, 0.025, 0.05, 0.1 or 0.2"))
  }
  
  if(!is.numeric(no.lambda) || no.lambda <= 0 || no.lambda != round(no.lambda)){
    return(list(error=TRUE, message="no.lambda must be a positive integer"))
  }
  
  return(list(error = FALSE))
}



# ----------- The Testing Procedure --------------------------------------------
rmsTest.func <- function(data, delta, alpha, no.lambda, base.period){
  # Formula for plm function:
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm.formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Number of individuals in the sample
  individuals <- unique(data[,"ID"])
  N <- length(individuals)
  
  # Number of periods:
  no.periods <- length(unique(data[,"period"]))
  
  # Matrix storing the placebo coefficient RMS on 1/lambda of the data for 
  # lambda = 1,..., no.lambda:
  MS.placebo.lambda <- rep(NA, no.lambda)
  for(i in 1:no.lambda){
    lambda <- i/no.lambda
    # Drawing the subset of individuals of size lambda*N
    sub.N <- floor(lambda*N)
    subset.index <- sample(1:N, sub.N)
    subset.indiv <- individuals[subset.index]
    
    subset.data <- data[data[,"ID"] %in% subset.indiv,]
    
    # Calculating the TWFE estimators for this subset of data:
    all.coefs <- plm::plm(plm.formula, data=subset.data, effect="twoways", 
                          model="within", index=c("ID","period"))$coefficients
    
    # Calculate the mean squared error for the placebo estimates:
    placebo.beta.lambda <- all.coefs[placebo_names]
    placebo.mean.sqrd <- base::mean((placebo.beta.lambda)^2)
    # Store the mean squared placebo coefficient estimate:
    MS.placebo.lambda[i] <- placebo.mean.sqrd
    
  }
  
  betas.placebo <- placebo.beta.lambda
  
  RMS.placebo <- sqrt(MS.placebo.lambda[no.lambda])
  MS.placebo <- MS.placebo.lambda[no.lambda]
  
  # Calculating \hat{V}_n:
  diff.vec <- MS.placebo.lambda[1:(no.lambda-1)] - MS.placebo.lambda[no.lambda]
  V_n <- sqrt(mean(diff.vec^2))
  
  # Calculating the (Bootstrap) critical value of the W distribution:
  Q.W <- W_critical_value(alpha)
  
  if(!is.null(delta)){
    # Critical value for RMS test
    MS.critical.value <- as.numeric(delta^2 + Q.W*V_n)
    RMS.critical.value <- sqrt(MS.critical.value)
    
    # Reject or not:
    reject.H0 <- MS.placebo.lambda[no.lambda] < MS.critical.value
    
    results.list <- structure(list(placebo.coefs = betas.placebo, RMS = RMS.placebo,
                                   critical.value = RMS.critical.value,
                                   reject.H0 = reject.H0, delta = delta,
                                   sign.level = alpha,
                                   N = N, no.periods = no.periods, base.period = base.period,
                                   delta.specified = TRUE),
                              class = "rmsEquivTest")
  } else {
    min.delta <- as.numeric(sqrt(MS.placebo-Q.W*V_n))
    results.list <- structure(list(placebo.coefs = betas.placebo, RMS = RMS.placebo,
                                   minimum.delta = min.delta, sign.level = alpha,
                                   N = N, no.periods = no.periods, base.period = base.period,
                                   delta.specified = FALSE),
                              class = "rmsEquivTest")
  }
  return(results.list)
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