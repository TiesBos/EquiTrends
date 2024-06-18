# ----------- The Mean Test Function -------------------------------------------
meanTest.func <- function(data, delta, vcov, cluster, alpha, n, no.periods, base.period){
  # Construct the formula for the plm() function
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm.formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  #print(plm.formula)
  # Run the two-way fixed effects model:
  plm.twfe <- plm::plm(plm.formula, data=data, effect="twoways", 
                       model="within", index=c("ID","period"))
  
  # Extract the estimated coefficients:
  betas <- plm.twfe$coefficients
  betas.placebo <- c(betas[placebo_names])
  
  # Calculating the mean of the absolute placebo coefficients:
  no.placebos <- length(placebo_names)
  one.vec <- rep(1, no.placebos)
  demean.vec <- one.vec/no.placebos
  mean.placebo <- abs(as.numeric(t(demean.vec)%*%betas.placebo))
  
  # Calculate the variance-covariance matrix:
  if(is.null(vcov)){
    vcov.mat <- plm.twfe$vcov
  } else if(vcov == "HC"){
    vcov.mat <- sandwich::vcovHC(plm.twfe, type="HC1", method = "white1")
  } else if(vcov == "HAC"){
    vcov.mat <- sandwich::vcovHC(plm.twfe, type="HC3", method = "arellano")
  } else if(vcov == "CL"){
    if(is.null(cluster)){
      vcov.mat <- clubSandwich::vcovCR(plm.twfe, cluster="ID", type="CR0")
    } else {
      vcov.mat <- clubSandwich::vcovCR(plm.twfe, cluster=data[,"cluster"], type="CR0")
    }
  } else {
    vcov.mat <- vcov(plm.twfe)
    # Verifying it gives a square matrix:
    if(!is.matrix(vcov.mat) || nrow(vcov.mat) != ncol(vcov.mat)){
      stop("The vcov argument is invalid.")
    }
  }
  
  # The estimated variance of the mean placebo coefficient estimator:
  placebo.vcov <- vcov.mat[placebo_names, placebo_names]
  sum.variance <- as.numeric(t(one.vec)%*%placebo.vcov%*%one.vec)
  mean.placebo.var <- sum.variance/n
  
  # If delta is not NULL, we calculate the critical value / p-value. Else, the minimum delta 
  # for which H0 can be rejected.
  if(!is.null(delta)){
    # The Critical value at the alpha level:
    crit.val <- VGAM::qfoldnorm(alpha, mean=delta, sd=sqrt(mean.placebo.var))
    
    # The p-value
    p_value <- VGAM::pfoldnorm(mean.placebo, mean=delta, sd=sqrt(mean.placebo.var))
    
    # reject H0:
    reject.H0 <- mean.placebo < crit.val
    
    results.list <- list(placebo.coefs = betas.placebo, abs.mean.placebo = mean.placebo,
                         placebo_var = mean.placebo.var, 
                         critical.value = crit.val, p.value = p_value, 
                         reject.H0 = reject.H0,
                         delta.specified = TRUE, delta = delta,
                         sign.level = alpha,
                         no.periods = no.periods, base.period = base.period, N=n)
    
  } else {
    
    min.delta <- meanTest.optim.func(mean.placebo, sqrt(mean.placebo.var), 0.05)
    
    results.list <- list(placebo.coefs = betas.placebo, abs.mean.placebo = mean.placebo,
                         placebo_var = mean.placebo.var,
                         delta = delta,
                         delta.specified = FALSE, 
                         minimum.delta = min.delta, sign.level = alpha,
                         no.periods = no.periods, base.period = base.period, N=n)
    
  }
  return(results.list)
}

# ----------- Finding the Minimum Delta ----------------------------------------
# Function to minimize:
meanTest.obj.func <- function(coef, mean, sd, alpha){
  return(1e100*(VGAM::pfoldnorm(coef, mean, sd) - alpha)^2)
}


meanTest.optim.func <- function(coef, sd, alpha){
  obj.wrapper <- function(x) meanTest.obj.func(coef=coef, mean=x, sd=sd, alpha=alpha)
  
  
  result <- nloptr::nloptr(x0 = coef,
                           eval_f = obj.wrapper,
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