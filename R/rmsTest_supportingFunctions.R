# ---- Error Function RMS ----
#' Additional input checks for the rmsEquivTest function
#'
#' @param alpha The significance level for the test. Must be one of 0.01, 0.025, 0.05, 0.1 or 0.2.
#' @param no_lambda see \link[EquiTrends]{rmsEquivTest}
#'
#' @return
#' A list with two elements: a logical object error indicating if an error is encountered and a message (a character string) corresponding to the error. If error is TRUE, message contains an error message. If error is FALSE, message is an empty string.
#'
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
#' An internal function of the RMS Equivalence Testing procedure
#'
#' @description This is a supporting function of the \code{rmsEquivTest} function. It calculates the placebo coefficients and the RMS of the placebo coefficients. It then calculates the critical value for the test and checks whether the null hypothesis can be rejected, according to Dette & Schumann (2024).
#'
#' @param data The data.frame object containing the data for the test. Should be of the form what is returned by the \link[EquiTrends]{EquiTrends_dataconstr} function.
#' @param equiv_threshold The equivalence threshold for the test. If NULL, the minimum equivalence threshold for which the null hypothesis can be rejected is calculated.
#' @param alpha The significance level for the test. Must be one of 0.01, 0.025, 0.05, 0.1 or 0.2.
#' @param no_lambda See \link[EquiTrends]{rmsEquivTest}.
#' @param base_period The base period for the test. Must be one of the unique periods in the data.
#' @param no_periods The number of periods in the data.
#' @param is_panel_balanced A logical value indicating whether the panel data is balanced.
#'
#' @references 
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \doi{10.1080/07350015.2024.2308121}
#'
#' @return
#' An object of class "rmsEquivTest" containing:
#' \item{\code{placebo_coefficients}}{A numeric vector of the estimated placebo coefficients,}
#' \item{\code{rms_placebo_coefs}}{the root mean squared value of the placebo coefficients,}
#' \item{\code{significance_level}}{the significance level of the test,}
#' \item{\code{num_individuals}}{the number of cross-sectional individuals in the data (n),}
#' \item{\code{num_periods}}{the number of pre-treatment periods in the data (T),}
#' \item{\code{num_observations}}{the number of observations in the data (N),}
#' \item{\code{base_period}}{the base period in the data,}
#' \item{\code{equiv_threshold_specified}}{a logical value indicating whether an equivalence threshold was specified.}
#' \item{\code{is_panel_balanced}}{a logical value indicating whether the panel data is balanced.}
#'
#' If \code{is.null(equiv_threshold)}, then additionally \code{minimum_equiv_threshold}: the minimum equivalence threshold for which the null hypothesis of non-negligible (based on the equivalence threshold) trend-differnces can be rejected. 
#' 
#' if \code{!(is.null(equiv_threshold))}, then additionally
#' \itemize{
#' \item \code{rms_critical_value}: the critical value at the alpha level,
#' \item \code{reject_null_hypothesis}: A logical value indicating whether to reject the null hypothesis,
#' \item \code{equiv_threshold}: the equivalence threshold specified.
#' }

rmsTest_func <- function(data, equiv_threshold, alpha, no_lambda, base_period, no_periods, is_panel_balanced){
  # Formula for plm function:
  placebo_names <- base::grep("placebo_",base::names(data),value=TRUE)
  X_names <- base::grep("X_", base::names(data), value=TRUE)
  plm_formula <- stats::as.formula(paste("Y~", paste(c(placebo_names, X_names), collapse=" + ")))
  
  # Number of individuals in the sample
  individuals <- unique(data[,"ID"])
  N <- length(individuals)
  
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
    placebo_names <- base::grep("placebo_",base::names(all_coefs),value=TRUE)
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
    if(MS_critical_value< 0){
      stop("The critical value is negative. Please enter a higher equivalence threshold.")
    }
    RMS_critical_value <- sqrt(MS_critical_value)
    
    # Reject or not:
    reject_H0 <- MS_placebo_lambda[no_lambda] < MS_critical_value
    
    results_list <- structure(list(placebo_coefficients = betas_placebo, rms_placebo_coefficients = RMS_placebo,
                                   rms_critical_value = RMS_critical_value,
                                   reject_null_hypothesis = reject_H0, equiv_threshold = equiv_threshold,
                                   significance_level = alpha,
                                   base_period = base_period,
                                   equiv_threshold_specified = TRUE,
                                   num_individuals = N, num_periods = no_periods,
                                   num_observations = nrow(data),
                                   is_panel_balanced = is_panel_balanced),
                              class = "rmsEquivTest")
  } else {
    min_equiv_threshold <- as.numeric(sqrt(MS_placebo-Q_W*V_n))
    results_list <- structure(list(placebo_coefficients = betas_placebo, rms_placebo_coefficients = RMS_placebo,
                                   minimum_equiv_threshold = min_equiv_threshold, significance_level = alpha,
                                   base_period = base_period,
                                   equiv_threshold_specified = FALSE,
                                   num_individuals = N, num_periods = no_periods,
                                   num_observations = nrow(data),
                                   is_panel_balanced = is_panel_balanced),
                              class = "rmsEquivTest")
  }
  return(results_list)
}


# ---- Critical Value RMS ----
# Function to get critical value for a given significance level
#' Calculating the critical value for the W distribution as construced in Dette & Schumann (2024).
#'
#' @param significance_level The significance level for the test. Must be one of 0.01, 0.025, 0.05, 0.1, 0.2, 0.8, 0.9, 0.95, 0.975, 0.99.
#' 
#' @references 
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \doi{10.1080/07350015.2024.2308121}
#'
#' @return
#' A numeric scalar with the critical value for the W distribution at the given significance level.
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