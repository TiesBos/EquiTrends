############ Supporting Functions of maxTest ###################################
# ----------- Intersection Union Approach --------------------------------------
maxTestIU <- function(data, delta, vcov, cluster, alpha, n, no.periods, base.period){
  # Construct the formula for the plm() function
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm.formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Run the two-way fixed effects model:
  IU.twfe <- plm::plm(plm.formula, data=data, effect="twoways", 
                      model="within", index=c("ID","period"))
  
  # Extract the estimated coefficients:
  betas <- IU.twfe$coefficients
  betas.placebo <- abs(c(betas[placebo_names]))
  
  # Extract the Variance-Covariance Matrix based on the user input:
  if(is.null(vcov)){
    vcov.mat <- IU.twfe$vcov
  } else if(vcov == "HC"){
    vcov.mat <- sandwich::vcovHC(IU.twfe, type="HC1", method = "white1")
  } else if(vcov == "HAC"){
    vcov.mat <- sandwich::vcovHC(IU.twfe, type="HC3", method = "arellano")
  } else if(vcov == "CL"){
    if(is.null(cluster)){
      vcov.mat <- clubSandwich::vcovCR(IU.twfe, cluster="ID", type="CR0")
    } else {
      vcov.mat <- clubSandwich::vcovCR(IU.twfe, cluster=data[,"cluster"], type="CR0")
    }
  } else {
    vcov.mat <- vcov(IU.twfe)
    if(!is.matrix(vcov.mat) || nrow(vcov.mat) != ncol(vcov.mat)){
      stop("The vcov argument is invalid")
    }
  }
  
  # Extracting the variances of the slope coefficients
  subvcov.mat <- vcov.mat[placebo_names, placebo_names]
  beta.var <- diag(subvcov.mat)
  
  # Calculating the standard errors
  beta.SE <- as.matrix(sqrt(beta.var/n))
  
  if(!is.null(delta)){
    
    # Critical Values method
    crit_values <- c(base::sapply(beta.SE, 
                                  FUN = function(x){return(VGAM::qfoldnorm(alpha, mean=delta, sd=x))}))
    reject.vec <- abs(betas.placebo) < crit_values
    
    # P-values calculated
    p_values <- c(base::mapply(FUN = function(x, y){return(VGAM::pfoldnorm(x, mean=delta, sd=y))},
                               betas.placebo, beta.SE))
    
    reject.H0 <- all(reject.vec==TRUE)
    
    results.list <- list(absolute.placebo.coefs = betas.placebo, standard.errors = beta.SE, 
                         critical.values = crit_values, p.values = p_values,
                         reject.H0 = reject.H0, delta = delta,
                         sign.level = alpha, no.periods = no.periods, 
                         base.period = base.period, N=n, coef.names = placebo_names,
                         delta.specified = TRUE)
    
    return(results.list)
  } else{
    
    # Calculating the minimum deltas
    min.deltas <- base::mapply(FUN = function(x,y){maxTestIU.optim.func(coef=x, sd=y, alpha=alpha)}, betas.placebo, beta.SE)
    # Then, the minimum delta is the maximum value over all these values:
    min.delta <- max(min.deltas)
    
    results.list <- list(absolute.placebo.coefs = betas.placebo, standard.errors = beta.SE,
                         minimum.delta = min.delta,
                         minimum.deltas = min.deltas,
                         sign.level = alpha,
                         no.periods = no.periods, base = base.period, N=n,
                         coef.names = placebo_names,
                         delta.specified = FALSE)
    
    return(results.list)
  }
}

# Functions to help find the minimum delta such that the p-value is alpha
maxTestIU.obj.func <- function(coef, mean, sd, alpha){
  return(1e20*(VGAM::pfoldnorm(coef, mean, sd) - alpha)^2)
}

maxTestIU.optim.func <- function(coef, sd, alpha){
  obj.wrapper <- function(x) maxTestIU.obj.func(coef=coef, mean=x, sd=sd, alpha=alpha)
  
  
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

# ----------- The Bootstrap Approaches -----------------------------------------
maxTestBoot <- function(data, delta, alpha, n, B, no.periods, 
                                 base.period, type){
  # Obtain the double demeaned data:
  dd.data <- double.demean(data)
  
  # find the double-demeaned independent and dependent variable:
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  X <- as.matrix(dd.data[, c(placebo_names, X_names)])
  Y <- as.matrix(dd.data$Y)
  
  # The unconstrained coefficient is:
  unconstrained.coefs <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  # its maximum is:
  max.unconstr.coef <- max(abs(unconstrained.coefs[1:length(placebo_names)]))
  
  if(max.unconstr.coef >= delta){
    constrained.coefs <- unconstrained.coefs
  } else {
    constrained.coefs <- optimization.function(x=X, y=Y, no.placebos = length(placebo_names), 
                                               delta = delta,
                                               start.val = unconstrained.coefs)$solution
  }
  
  if(type == "Boot"){
    # Calculate the variance based on the constrained coefficient:
    resid.variance <- sigma_hathat_c(parameter = constrained.coefs, N = n, 
                                     x=X, y=Y,
                                     no.periods = no.periods)
    # Run the Bootstrap:
    bootstrap.maxcoefs <- maxTestBoot_bootstrap(Xb = X%*%constrained.coefs, X=X, B=B,
                                                variance = resid.variance, ID = dd.data$ID,
                                                period = dd.data$period, no_placebos = length(placebo_names))
  } else {
    u_ddot <- Y - X%*%constrained.coefs
    
    # Run the Wild Bootstrap:
    bootstrap.maxcoefs <- maxTestBoot_wildbootstrap(Xb = X%*%constrained.coefs, X=X, B=B,
                                                     u_ddot = u_ddot,
                                                     ID = dd.data$ID, period = dd.data$period,
                                                     no_placebos = length(placebo_names))
  }
  
  # Find the critical value at the alpha level:
  boot.crit.value <- stats::quantile(bootstrap.maxcoefs, probs = alpha)
  
  # Reject Or Not:
  reject.H0 <- (max(abs(unconstrained.coefs[1:length(placebo_names)])) < boot.crit.value) 
  
  results.list <- list(absolute.placebo.coefs = abs(unconstrained.coefs[1:length(placebo_names)]),
                       max.coefs =max(abs(unconstrained.coefs[1:length(placebo_names)])),
                       critical.values = boot.crit.value,
                       reject.H0 = reject.H0, B = B, delta.specified = TRUE, 
                       delta = delta, sign.level = alpha, wild = FALSE,
                       no.periods = no.periods, base.period = base.period,
                       N=n)
  
  return(results.list)
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
double.demean <- function(df){
  # Finding all variables that need to be double demeaned
  placebo_names <- grep("placebo_", names(df), value=TRUE)
  X_names <- grep("X_", names(df), value=TRUE)
  names <- c(placebo_names, X_names)
  
  # Matrix:
  original.data <- df[, c("Y", names)]
  
  # Calculate the mean for each unit (ID) for each value column
  unit_means <- as.data.frame(df %>%
                                group_by_at("ID") %>%
                                summarise(across(all_of(c("Y", names)), mean, .names = "unit_means_{col}")))
  
  # Calculate the mean for each time period for each value column
  time_means <-  as.data.frame(df %>%
                                 group_by_at("period") %>%
                                 summarise(across(all_of(c("Y", names)), mean, .names = "time_means_{col}")))
  
  # Merge the calculated means back to the original dataframe 
  new.df <- df
  new.df <- new.df %>%
    left_join(unit_means, by = "ID") %>%
    left_join(time_means, by = "period")
  
  # Get the correct dimensions for W of the time means, unit means and overall means
  unit.mean.names <- grep("unit_means_", names(new.df), value=TRUE)
  time.means.names <- grep("time_means_", names(new.df), value=TRUE)
  unit.means.mat <- new.df[,unit.mean.names]
  time.means.mat <- new.df[,time.means.names]
  
  over.all.means <- colMeans(as.matrix(original.data))
  n.row <- nrow(as.matrix(original.data))
  overall.means <- as.data.frame(matrix(rep(over.all.means, times = n.row), nrow = n.row, byrow = TRUE))
  
  demeaned.data <- original.data - unit.means.mat - time.means.mat + overall.means
  ID <- df$ID
  period <- df$period
  demeaned.data <- cbind(ID, period, demeaned.data)
  colnames(demeaned.data) <- c("ID", "period", "Y", names)
  return(demeaned.data)
}

# The objective function of the constrained estimator:
objective.function <- function(parameter, x, y, no.placebos, delta){
  Xb <- as.numeric(x%*%parameter)
  MSE <- mean((y-Xb)^2)
  return(MSE)
}


# The constraint function for the optimization. Note that we only penalize the 
# coefficients corresponding to the placebos!
constraint.function <- function(parameter, x, y, no.placebos, delta){
  value <- -max(abs(parameter[1:no.placebos])) + delta
  return(value)
}


# Optimization function:
optimization.function <- function(x, y, no.placebos, delta, start.val){
  constrained.optimum <- nloptr::nloptr(x0 = start.val,
                                        eval_f = objective.function,
                                        eval_grad_f = NULL,
                                        lb = rep(-Inf, length(start.val)),
                                        ub = rep(Inf, length(start.val)),
                                        eval_g_ineq = constraint.function,
                                        eval_jac_g_ineq = NULL,
                                        eval_g_eq = NULL,
                                        eval_jac_g_eq = NULL,
                                        opts = list("algorithm" = "NLOPT_LN_COBYLA", 
                                                    maxeval=2000000, xtol_rel = 1e-06),
                                        x = x, y=y, no.placebos = no.placebos, delta=delta)
  return(constrained.optimum)
}

# Constrained Variance:
sigma_hathat_c <- function(parameter, x, y, N, no.periods){
  Xb <- as.numeric(x%*%parameter)
  
  residuals.demeaned <- y - Xb
  
  c.sigma.hathat <- sum((residuals.demeaned^2))/((N-1)*no.periods)
  return(c.sigma.hathat)
}



# --- maxTest Error Function ---
maxTest.error <- function(type, delta, vcov){
  
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
  if(type != "IU" && is.null(delta)){
    return(list(error=TRUE, message = "delta must be specified for types Boot and Wild."))
  }
  
  return(list(error=FALSE))
}










