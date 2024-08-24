#' @title Simulating a panel data for a binary treatment
#' 
#' @description sim.paneldata generates a panel data set with N cross-sectional units and tt time periods. The data set includes a binary treatment variable, a set of placebo variables, and a set of additional regressors. The data set can be generated under homoskedasticity or heteroskedasticity, and/or AR(1) errors.
#'
#' @param N The number of cross-sectional units in the panel-data
#' @param tt The number of time periods in the panel-data
#' @param beta The vector of coefficients for the placebo variables. Must be of size tt.
#' @param p The number of additional regressors
#' @param gamma The vector of coefficients for the additional regressors
#' @param eta The vector of fixed effects. Must be of size N.
#' @param lambda The vector of time effects. Must be of size tt.
#' @param het The heteroskedasticity parameter. Must be 0 or 1: \code{het = 1} indicates that the error terms are generated under heteroskedasticity, \code{het = 0} indicates the error terms are generated under homoscedasticity. 
#' @param phi The AR(1) parameter for the error terms. Must be in the interval [0,1). 
#' @param sd The standard deviation of the error terms. Must be a positive number.
#' @param burnins The number of burn-ins for the AR(1) process. Must be a positive integer.
#' 
#' @export
#' 
#' @return A \code{data.frame} with the following columns:
#'  \item{ID}{The cross-sectional unit identifier}
#'  \item{period}{The time period identifier}
#'  \item{Y}{The dependent variable}
#'  \item{G}{The binary treatment variable}
#'  \item{X_1, ..., X_p}{The additional regressors}
#'
#' 
#' @examples
#' sim_data <- sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p=1, 
#'                           gamma = rep(0,1), het = 1, phi = 0.5, sd = 1, 
#'                           burnins = 100)
sim_paneldata <- function(N = 500, tt = 5, beta = rep(0, tt), p=1, gamma = rep(1, p),
                          eta = rep(0, N), lambda = rep(0, tt), het = 0, 
                          phi = c(0), sd = 1, burnins = 100){
  # Check input:
  checks <- sim_check(N, tt, beta, p, gamma, eta, lambda, het, phi, sd, burnins)
  if(checks$error){stop(checks$message)}
  
  # Create the cross-sectional individual identifier
  ID=rep(1:N,each=tt)
  # Create the time period identifier:
  period=rep(1:tt,N)
  
  # create treatment dummy; only individuals with even numbers are treated
  G=ifelse(ID%%2==0,1,0)
  
  # The number of regressors:
  X <- matrix(stats::rnorm(N*tt*p), nrow = N*tt, ncol = p)
  
  # create N*tt vector of fixed effects that are constant for each i
  alpha_full <- rep(eta, each = tt)
  
  # create N*tt vector of time effects that are constant for each t
  lambda_full = rep(lambda, N)
  
  # The dummies and placebos:
  placebos <- data.frame(ID, period)
  for(dummies in 1:(tt))
  {
    placebos[, paste0("dummy_", dummies)]= ifelse(period==dummies,1,0)
    placebos[, paste0("placebo_", dummies)]=placebos[[paste0("dummy_", dummies)]]*G
  }
  
  # Creating the error-terms:
  eps = rep(NA, N*tt)
  if(phi == 0 && het == 0){
    eps=stats::rnorm(N*tt,0,sd)
  } else {
    for (i in 1:N) {
      start_idx <- (i - 1) * tt + 1
      end_idx <- i * tt
      
      # Generate AR(1) process using arima.sim()
      het.error=G[start_idx:end_idx]*het
      ar_process = stats::arima.sim(n = tt, list(ar = c(phi)), 
                                    sd = sd, n.start = burnins)
      
      # Place the block values into the result vector
      eps[start_idx:end_idx] = ar_process
    }
  }
  
  placebos <- as.matrix(placebos[paste0("placebo_", 1:tt)])
  
  
  # The dependent variable:
  Y = alpha_full + lambda_full + X %*% gamma + placebos %*% beta + eps
  
  if(p > 0){
    sim.data=data.frame(ID, period, Y, G, X)
    colnames(sim.data)=c("ID","period", "Y", "G",paste0("X_",1:p))
  } else {
    sim.data=data.frame(ID, period, Y, G)
    colnames(sim.data)=c("ID","period", "Y", "G")
  } 
  return(sim.data)
}


#' Checking input for the sim_paneldata function
#'
#' @param N The number of cross-sectional units in the panel-data
#' @param tt The number of time periods in the panel-data
#' @param beta The vector of coefficients for the placebo variables. Must be of size tt.
#' @param p The number of additional regressors
#' @param gamma The vector of coefficients for the additional regressors
#' @param eta The vector of fixed effects. Must be of size N.
#' @param lambda The vector of time effects. Must be of size tt.
#' @param het The heteroskedasticity parameter. Must be 0 or 1: \code{het = 1} indicates that the error terms are generated under heteroskedasticity, \code{het = 0} indicates the error terms are generated under homoscedasticity. 
#' @param phi The AR(1) parameter for the error terms. Must be in the interval [0,1). 
#' @param sd The standard deviation of the error terms. Must be a positive number.
#' @param burnins The number of burn-ins for the AR(1) process. Must be a positive integer.
#'
#' @return
#' A list with two elements: a logical object error indicating if an error is encountered and a message (a character string) corresponding to the error. If error is TRUE, message contains an error message. If error is FALSE, message is an empty string.
#' 
sim_check <- function(N, tt, beta, p, gamma, eta, lambda, het, phi, sd, burnins){
  if(!is.numeric(N) || N <= 0 || N != round(N)){
    return(list(error=TRUE, message="N must be a positive integer"))
  }
  if(!is.numeric(tt) || tt <= 0 || tt != round(tt)){
    return(list(error=TRUE, message="tt must be a positive integer"))
  }
  if(!is.numeric(p) || p < 0 || p != round(p)){
    return(list(error=TRUE, message="p must be a non-negative integer"))
  }
  if(!(het %in% c(0,1))){
    return(list(error=TRUE, message="het must take on a value of 0 or 1"))
  }
  if(length(phi) != 1){
    return(list(error=TRUE, message="phi must be a scalar"))
  }
  if(!(phi >= 0 && phi < 1)){
    return(list(error=TRUE, message="phi must be in the interval [0,1)"))
  }
  if(!is.numeric(sd) || sd <= 0){
    return(list(error=TRUE, message="sd must be a positive number"))
  }
  if(!is.numeric(burnins) || burnins <= 0 || burnins != round(burnins)){
    return(list(error=TRUE, message="burnins must be a positive integer"))
  }
  if(length(eta) != N){
    return(list(error=TRUE, message="eta must be a vector of length N"))
  }
  if(length(lambda) != tt){
    return(list(error=TRUE, message="lambda must be a vector of length tt"))
  }
  if(length(beta) != tt){
    return(list(error=TRUE, message="beta must be a vector of length tt"))
  }
  if(length(gamma) != p){
    return(list(error=TRUE, message="gamma must be a vector of length p"))
  }
  return(list(error=FALSE))
}


