#' @title Equivalence Test for Pre-trends based on the Maximum Absolute Placebo Coefficient 
#'
#' @description This function performs an equivalence test for pre-trends based on the maximum absolute placebo coefficient from Dette & Schumann (2024). The test can be performed using the intersection-union approach (IU), a bootstrap procedure for spherical errors (Boot) and a wild bootstrap procedure (Wild).
#'
#' @param Y A numeric vector with the variable of interest. If \code{data} is supplied, \code{Y} should be a scalar indicating the column number or column-name character string that corresponds to the numeric dependent (outcome) variable in ’data’.
#' @param ID A numeric vector identifying the different cross-sectional units in the dataset. If \code{data} is supplied, \code{ID} should be a scalar indicating the column number or column-name character string that corresponds to the cross-sectional units identifier in \code{data}.
#' @param G A binary or logic vector (of the same dimension as \code{Y} and \code{ID}) indicating if the individual (e.g. as indicated by \code{ID}) receives treatment (e.g. 1 or TRUE) or not (0 or FALSE). f 'data' is supplied, \code{G} should be a scalar identifying the column number or column-name character string associated to \code{G} in \code{data}.
#' @param period A numeric vector (of the same dimension as Y) indicating time. If \code{data} is supplied, \code{period} should be a scalar indicating the column number or column-name character string that corresponds to the time identifier in \code{data}.
#' @param X  A vector, matrix, or data.frame containing the control variables. If \code{data} is supplied, \code{X} must be a vector of column numbers or column-name character strings that identifies the control variables’ columns. 
#' @param data An optional \code{data.frame} object containing the variables in Y, ID, G, T and, if supplied, X and cluster as its columns.
#' @param equiv_threshold The scalar equivalence threshold (must be positive). The default is NULL, implying that the function must look for the minimum value for which the null hypothesis of ”non-negligible differences” can still be rejected.
#' @param pretreatment_period A numeric vector identifying the pre-treatment periods that should be used for testing. \code{pretreatment_period} must be a subset of the periods included through \code{period}. The default is to use all periods that are included in \code{period}.
#' @param base_period The pre-treatment period to compare the post-treatment observation to. The default is to take the last period of the pre-treatment period.
#' @param vcov The variance-covariance matrix that needs to be used. See \emph{Details} for more details.
#' @param cluster If \code{vcov = "CL"}, a vector indicating which observations belong to the same cluster. \code{cluster} must be of the same length as Y. If \code{data} is supplied, \code{cluster} must be either the column index or column name of this vector in the data.frame/matrix. The default (\code{cluster=NULL}) assumes every unit in ID is its own cluster. Only required if \code{vcov = "CL"} and \code{type = "IU"}.
#' @param alpha Significance level of the test. The default is 0.05. Only required if \code{equiv_threshold} is not specified.
#' @param type The type of maximum test that should be performed. "IU" for the intersection-union test, "Boot" for the regular bootstrap procedure from Dette & Schumann (2023) and "Wild" for the Wild bootstrap procedure.
#' @param B If type = Boot or type = Wild, the number of bootstrap samples used. The default is 1000.
#'
#' 
#' @details The \code{vcov} parameter specifies the variance-covariance matrix to be used in the function for \code{type = "IU"}. 
#' This parameter can take two types of inputs:
#' \enumerate{
#' \item A character string specifying the type of variance-covariance matrix estimation. The options are:\itemize{
#'    \item \code{NULL}: The default variance-covariance matrix estimated by the \link[plm]{plm} function is used.
#'    \item \code{"HC"}: A heteroscedasticity-robust (HC) covariance matrix is estimated using the \code{vcovHC} function from the \code{plm} package, \link[plm]{vcovHC}, with type \code{"HC1"} and method \code{"white1"} (see White, 1980).
#'    \item \code{"HAC"}: A heteroscedasticity and autocorrelation robust (HAC) covariance matrix is estimated using the \code{vcovHC} function from the \code{plm} package, \link[plm]{vcovHC}, with type \code{"HC3"} and method \code{"arellano"} (see Arellano, 1987).
#'    \item \code{"CL"}: A cluster-robust covariance matrix is estimated using the \code{vcovCR} function from the \code{clubSandwich} package with type \code{"CR0"} (see Lian & Zegers (1986)). The cluster variable is either \code{"ID"} or a custom cluster variable provided in the \code{data} dataframe.
#' }
#'\item A function that takes an \link[plm]{plm} object as input and returns a variance-covariance matrix. 
#'    This allows for custom variance-covariance matrix estimation methods. For example, you could 
#'    use the \code{vcovHC} function from the \code{sandwich} package with a specific method and type:
#'    \preformatted{function(x) {vcovHC(x, method = "white1", type = "HC2")}}
#' }
#' If no \code{vcov} parameter is provided, the function defaults to using the variance-covariance matrix 
#' estimated by the \link[plm:plm]{plm::plm()} function.
#' 
#' Note that for the bootstrap functions (\code{type= "Boot"} and \code{type = "Wild"}), the equiv_threshold must be specified. Currently, no option exist to find the minimum equivalence threshold for these types.
#' 
#' @references
#' Arellano M (1987). “Computing Robust Standard Errors for Within-groups Estimators.” \emph{Oxford bulletin of Economics and Statistics}, 49(4), 431–434.
#' 
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \href{https://doi.org/10.1080/07350015.2024.2308121}{10.1080/07350015.2024.2308121}
#' 
#' Liang, K.-Y., & Zeger, S. L. (1986). "Longitudinal data analysis using generalized linear models." \emph{Biometrika}, 73(1), 13-22. \href{doi:10.1093/biomet/73.1.13}{doi:10.1093/biomet/73.1.13}
#' 
#' White H (1980). “A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity.” \emph{Econometrica}, 48(4), 817–838.
#'
#' @seealso \code{\link[=print.maxEquivTestBoot]{print.maxEquivTestBoot}} \code{\link[=print.maxEquivTestIU]{print.maxEquivTestIU}}
#' 
#' @return If \code{type = "IU"}, an object of class \code{maxEquivTestIU}  with \itemize{
#' \item{\code{placebo_coefficients}: A numeric vector of the estimated placebo coefficients},
#' \item{\code{abs_placebo_coefficients}: a numeric vector with the absolute values of estimated placebo coefficients,}
#' \item{\code{placebo_coefficients_se}: a numeric vector with the standard errors of the placebo coefficients,}
#' \item{\code{significance_level}: the chosen significance level of the test,}
#' \item{\code{num_individuals}: the number of cross-sectional individuals in \code{data},}
#' \item{\code{num_periods}: the number of periods in \code{data},}
#' \item{\code{base_period}: the base period in the data,}
#' \item{\code{placebo_names}: the names corresponding to the placebo coefficients,}
#' \item{\code{equiv_threshold_specified}: a logical value indicating whether an equivalence threshold was specified.}
#' \item{if \code{!(is.null(equiv_threshold))}\itemize{
#'  \item{\code{IU_critical_values}: a numeric vector with the individual critical values for each of the placebo coefficients,}
#'  \item{\code{reject_null_hypothesis}: a logical value indicating whether the null hypothesis of negligible pre-trend differences can be rejected at the specified significance level \code{alpha},}
#'  \item{\code{equiv_threshold}: the equivalence threshold employed.}
#' }}
#' \item{Additionally, if \code{is.null(equiv_threshold)}\itemize{
#' \item{\code{minimum_equiv_thresholds}: a numeric vector including for each placebo coefficient the minimum equivalence threshold for which the null hypothesis of negligible pre-trend differences can be rejected for the corresponding placebo coefficient individually,}
#' \item{\code{minimum_equiv_threshold}: a numeric scalar minimum equivalence threshold for which the null hypothesis of negligible pre-trend differences can be rejected for all placebo coefficients individually.}
#' }}}
#' 
#' if \code{type = "Boot"} or \code{type = "Wild"}, an object of class "maxEquivTestBoot"  with \itemize{
#'  \item{\code{placebo_coefficients}: a numeric vector of the estimated placebo coefficients,}
#'  \item{\code{abs_placebo_coefficients}: a numeric vector with the absolute values of estimated placebo coefficients,}
#'  \item{\code{max_abs_coefficient}: the maximum absolute estimated placebo coefficient,}
#'  \item{\code{bootstrap_critical_value}: the by bootstrap found critical value for the equivalence test based on the maximum absolute placebo coefficient,}
#'  \item{\code{reject_null_hypothesis}: a logical value indicating whether the null hypothesis of negligible pre-trend differences can be rejected at the specified significance level \code{alpha},}
#'  \item{\code{B}: the number of bootstrap samples used to find the critical value,}
#'  \item{\code{significance_level}: the chosen significance level of the test \code{alpha},}
#'  \item{\code{num_individuals}: the number of cross-sectional individuals in \code{data},}
#'  \item{\code{num_periods}: the number of periods in \code{data},}
#'  \item{\code{base_period}: the base period in the data,}
#'  \item{\code{placebo_names}: the names corresponding to the placebo coefficients,}
#'  \item{\code{equiv_threshold_specified}: a logical value indicating whether an equivalence threshold was specified.}
#' }
#' 
#' @examples
#' # generate some data
#' # Simulate data
#' sim_data <- sim_paneldata()
#' 
#' # -----------------  IU Approach -----------------
#' # Perform the test with equivalent threshold specified as 1 based on 
#' # pre-treatment periods 1:4 and homoscedastic error-terms:
#' equivalence_test <- maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                                  data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                                  base_period = 4, type = "IU")
#' print(equivalence_test)
#' 
#' # Perform the test without specifying the equivalence threshold with heteroscedastic 
#' # and autocorrelation robust variance-covariance matrix estimator:
#' equivalence_test <- maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4, vcov = "HAC", type = "IU")
#' print(equivalence_test)
#' 
#' # Perform the test without specifying the equivalence threshold with a custom
#' # variance-covariance matrix estimator:
#' vcov_func <- function(x) {plm::vcovHC(x, method = "white1", type = "HC2")}
#' 
#' equivalence_test <- maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4, vcov = vcov_func, type = "IU")
#' print(equivalence_test)
#' 
#' #-----------------  Bootstrap Approach -----------------
#' # Perform the test with equivalent threshold specified as 1 based on 
#' # pre-treatment periods 1:4 (with base period 4) with the general bootstrap procedure:
#' equivalence_test <- maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4, type = "Boot")
#' print(equivalence_test)
#' 
#' # Perform the test with equivalent threshold specified as 1 based on 
#' # pre-treatment periods 1:4 (with base period 4) with the wild bootstrap procedure:
#' equivalence_test <- maxEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4, type = "Wild")
#' print(equivalence_test)
#' 
#' @export
#' 
maxEquivTest <- function(Y, ID, G, period, X = NULL, data = NULL, equiv_threshold = NULL,  
                       pretreatment_period = NULL, base_period = NULL, 
                       vcov = NULL, cluster = NULL, alpha = 0.05, 
                       type = c("IU", "Boot", "Wild"), B = 1000){
  # If no type is specified, the type is "IU"
  if(identical(type, c("IU", "Boot", "Wild"))){
    type <- "IU"
    warning("type is not specified: Intersection-Union (IU) approach is used.")
  }
  
  # Error/Warnings:
  # Checking errors specific to this test procedure
  error_maxTest <- maxTest_error(type, equiv_threshold, vcov, B)
  if(error_maxTest$error){stop(error_maxTest$message)}
  
  
  # General error checking:
  error_test <- EquiTrends_inputcheck(Y, ID, G, period, X, data, equiv_threshold, pretreatment_period, 
                                     base_period, cluster, alpha)
  
  if(error_test$error){stop(error_test$message)}
  
  # Structuring the data:
  data_constr <- EquiTrends_dataconstr(Y, ID, G, period, X, data, pretreatment_period, base_period,
                                   cluster)
 
  # The dataframe:
  df <- data_constr$dataset

  # The base period:
  base_period <- data_constr$baseperiod
  
  # Original names of control variables:
  orig_names <- data_constr$orig_names
  
  # Panel balanced or not:
  panel_balanced <- data_constr$panel_balanced
  
  # Number of individuals
  N <- length(df[unique(df$ID), "ID"])
  # Number of Time Periods:
  no_periods <- length(df[unique(df$period), "period"])
  
  if(type=="IU"){
    test_results <- maxTestIU_func(df, equiv_threshold, vcov, cluster, alpha, N, no_periods, base_period)
  } else if (type == "Boot" || type == "Wild"){
    test_results <- maxTestBoot_func(df, equiv_threshold, alpha, N, B, no_periods, base_period, type, orig_names)
  }
  
  return(test_results)
}

#' @title Equivalence Test for Pre-trends based on the Mean Placebo Coefficient 
#'
#' @description This function performs an equivalence test for pre-trends based on the mean placebo coefficient from Dette & Schumann (2024). 
#' 
#' @param Y A numeric vector with the variable of interest. If \code{data} is supplied, \code{Y} should be a scalar indicating the column number or column-name character string that corresponds to the numeric dependent (outcome) variable in ’data’.
#' @param ID A numeric vector identifying the different cross-sectional units in the dataset. If \code{data} is supplied, \code{ID} should be a scalar indicating the column number or column-name character string that corresponds to the cross-sectional units identifier in \code{data}.
#' @param G A binary or logic vector (of the same dimension as \code{Y} and \code{ID}) indicating if the individual (e.g. as indicated by \code{ID}) receives treatment (e.g. 1 or TRUE) or not (0 or FALSE). f 'data' is supplied, \code{G} should be a scalar identifying the column number or column-name character string associated to \code{G} in \code{data}.
#' @param period A numeric vector (of the same dimension as Y) indicating time. If \code{data} is supplied, \code{period} should be a scalar indicating the column number or column-name character string that corresponds to the time identifier in \code{data}.
#' @param X  A vector, matrix, or data.frame containing the control variables. If \code{data} is supplied, \code{X} must be a vector of column numbers or column-name character strings that identifies the control variables’ columns. 
#' @param data An optional \code{data.frame} object containing the variables in Y, ID, G, T and, if supplied, X and cluster as its columns.
#' @param equiv_threshold The scalar equivalence threshold (must be positive). The default is NULL, implying that the function must look for the minimum value for which the null hypothesis of ”non-negligible differences” can still be rejected.
#' @param pretreatment_period A numeric vector identifying the pre-treatment periods that should be used for testing. \code{pretreatment_period} must be a subset of the periods included through \code{period}. The default is to use all periods that are included in \code{period}.
#' @param base_period The pre-treatment period to compare the post-treatment observation to. The default is to take the last period of the pre-treatment period.
#' @param vcov The variance-covariance matrix that needs to be used. See \emph{Details} for more details.
#' @param cluster If \code{vcov = "CL"}, a vector indicating which observations belong to the same cluster. \code{cluster} must be of the same length as Y. If \code{data} is supplied, \code{cluster} must be either the column index or column name of this vector in the data.frame/matrix. The default (\code{cluster=NULL}) assumes every unit in ID is its own cluster.
#' @param alpha Significance level of the test. The default is 0.05. Only required if \code{equiv_threshold} is not specified.
#' 
#' @seealso \code{\link[=print.meanEquivTest]{print.meanEquivTest}}
#' 
#' @details 
#' The \code{vcov} parameter specifies the variance-covariance matrix to be used in the function. 
#' This parameter can take two types of inputs:
#' \enumerate{
#' \item A character string specifying the type of variance-covariance matrix estimation. The options are:\itemize{
#'    \item \code{NULL}: The default variance-covariance matrix estimated by the \link[plm]{plm} function is used.
#'    \item \code{"HC"}: A heteroscedasticity-robust (HC) covariance matrix is estimated using the \code{vcovHC} function from the \code{plm} package, \link[plm]{vcovHC}, with type \code{"HC1"} and method \code{"white1"} (see White, 1980).
#'    \item \code{"HAC"}: A heteroscedasticity and autocorrelation robust (HAC) covariance matrix is estimated using the \code{vcovHC} function from the \code{plm} package, \link[plm]{vcovHC}, with type \code{"HC3"} and method \code{"arellano"} (see Arellano, 1987).
#'    \item \code{"CL"}: A cluster-robust covariance matrix is estimated using the \code{vcovCR} function from the \code{clubSandwich} package with type \code{"CR0"} (see Lian & Zegers (1986)). The cluster variable is either \code{"ID"} or a custom cluster variable provided in the \code{data} dataframe.
#' }
#'\item A function that takes an \link[plm]{plm} object as input and returns a variance-covariance matrix. 
#'    This allows for custom variance-covariance matrix estimation methods. For example, you could 
#'    use the \code{vcovHC} function from the \code{sandwich} package with a specific method and type:
#'    \preformatted{function(x) {vcovHC(x, method = "white1", type = "HC2")}}
#' }
#' If no \code{vcov} parameter is provided, the function defaults to using the variance-covariance matrix 
#' estimated by the \link[plm:plm]{plm::plm()} function.
#'
#' @return 
#' An object of class "meanEquivTest" containing:\itemize{
#' \item{\code{placebo_coefficients}}{A numeric vector of the estimated placebo coefficients,}
#' \item{\code{abs_mean_placebo_coefs}}{the absolute value of the mean of the placebo coefficients,}
#' \item{\code{var_mean_placebo_coef}}{the estimated variance of the mean placebo coefficient,}
#' \item{\code{significance_level}}{the significance level of the test,}
#' \item{\code{num_individuals}}{the number of cross-sectional individuals in the data,}
#' \item{\code{num_periods}}{the number of periods in the data,}
#' \item{\code{base_period}}{the base period in the data,}
#' \item{\code{equiv_threshold_specified}}{a logical value indicating whether an equivalence threshold was specified.}
#'}
#' If \code{is.null(equiv_threshold)}, then additionally \code{minimum_equiv_threshold}: the minimum equivalence threshold for which the null hypothesis of non-negligible (based on the equivalence threshold) trend-differnces can be rejected. 
#' 
#' if \code{!(is.null(equiv_threshold))}, then additionally \itemize{
#' \item \code{mean_critical_value}: the critical value at the alpha level,
#' \item \code{p_value}: the p-value of the test,
#' \item \code{reject_null_hypothesis}: A logical value indicating whether to reject the null hypothesis,
#' \item \code{equiv_threshold}: the equivalence threshold specified.
#' }
#' 
#' @examples
#' # Simulate data
#' sim_data <- sim_paneldata()
#' 
#' # Perform the test with equivalent threshold specified as 1 based on 
#' # pre-treatment periods 1:4 and homoscedastic error-terms:
#' equivalence_test <- meanEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4)
#' print(equivalence_test)
#' 
#' # Perform the test without specifying the equivalence threshold with heteroscedastic 
#' # and autocorrelation robust variance-covariance matrix estimator:
#' equivalence_test <- meanEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4, vcov = "HAC")
#' print(equivalence_test)
#' 
#' # Perform the test without specifying the equivalence threshold with a custom
#' # variance-covariance matrix estimator:
#' vcov_func <- function(x) {plm::vcovHC(x, method = "white1", type = "HC2")}
#' 
#' equivalence_test <- meanEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4, vcov = vcov_func)
#' print(equivalence_test)
#' 
#' @references
#' Arellano M (1987). “Computing Robust Standard Errors for Within-groups Estimators.” \emph{Oxford bulletin of Economics and Statistics}, 49(4), 431–434.
#' 
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \href{https://doi.org/10.1080/07350015.2024.2308121}{10.1080/07350015.2024.2308121}
#' 
#' Liang, K.-Y., & Zeger, S. L. (1986). "Longitudinal data analysis using generalized linear models." \emph{Biometrika}, 73(1), 13-22. \href{doi:10.1093/biomet/73.1.13}{doi:10.1093/biomet/73.1.13}
#' 
#' White H (1980). “A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity.” \emph{Econometrica}, 48(4), 817–838.
#' 
#'
#' 
#' @export
#'
meanEquivTest <- function(Y, ID, G, period, X = NULL, data = NULL, equiv_threshold = NULL,  
                     pretreatment_period = NULL, base_period = NULL, 
                     vcov = NULL, cluster = NULL, alpha=0.05){
  
  # Error/Warnings:
  # General error checking:
  error_test <- EquiTrends_inputcheck(Y, ID, G, period, X, data, equiv_threshold, pretreatment_period, 
                                     base_period, cluster, alpha)
  
  if(error_test$error){stop(error_test$message)}
  
  # We construct the data.frame:
  data_constr <- EquiTrends_dataconstr(Y, ID, G, period, X, data, pretreatment_period, base_period,
                                   cluster)
  df <- data_constr$data
  
  # Store the base.period:
  base_period <- data_constr$baseperiod
  
  # number of individuals in the sample:
  N <- length(df[unique(df$ID), "ID"])
  # The number of periods:
  no_periods <- length(unique(df$period))
  
  # Perform the test:
  results <- meanTest_func(df, equiv_threshold, vcov, cluster, alpha, N, no_periods, base_period)
  
  return(results)
}



#' @title Equivalence Test for Pre-trends based on the RMS Placebo Coefficient
#' 
#' @description This function performs an equivalence test for pre-trends based on the root mean squared placebo coefficient from Dette & Schumann (2024).
#'
#' @param Y A numeric vector with the variable of interest. If \code{data} is supplied, \code{Y} should be a scalar indicating the column number or column-name character string that corresponds to the numeric dependent (outcome) variable in ’data’.
#' @param ID A numeric vector identifying the different cross-sectional units in the dataset. If \code{data} is supplied, \code{ID} should be a scalar indicating the column number or column-name character string that corresponds to the cross-sectional units identifier in \code{data}.
#' @param G A binary or logic vector (of the same dimension as \code{Y} and \code{ID}) indicating if the individual (e.g. as indicated by \code{ID}) receives treatment (e.g. 1 or TRUE) or not (0 or FALSE). f 'data' is supplied, \code{G} should be a scalar identifying the column number or column-name character string associated to \code{G} in \code{data}.
#' @param period A numeric vector (of the same dimension as Y) indicating time. If \code{data} is supplied, \code{period} should be a scalar indicating the column number or column-name character string that corresponds to the time identifier in \code{data}.
#' @param X  A vector, matrix, or data.frame containing the control variables. If \code{data} is supplied, \code{X} must be a vector of column numbers or column-name character strings that identifies the control variables’ columns. 
#' @param data An optional \code{data.frame} object containing the variables in Y, ID, G, T and, if supplied, X and cluster as its columns.
#' @param equiv_threshold The scalar equivalence threshold (must be positive). The default is NULL, implying that the function must look for the minimum value for which the null hypothesis of ”non-negligible differences” can still be rejected.
#' @param pretreatment_period A numeric vector identifying the pre-treatment periods that should be used for testing. \code{pretreatment_period} must be a subset of the periods included through \code{period}. The default is to use all periods that are included in \code{period}.
#' @param base_period The pre-treatment period to compare the post-treatment observation to. The default is to take the last period of the pre-treatment period.
#' @param alpha Significance level of the test. The default is 0.05.
#' @param no_lambda Parameter specifying the number of incremental segments of the dataset over which a statistic is calculated. See \emph{Details}. The default is 5.
#'  
#' @details
#' \code{no_lambda} determines the proportions lambda/\code{no.lambda} for lambda = 1,...,\code{no_lambda} of the cross-sectional units at which the placebo coefficients are estimated. The placebo coefficients are estimated for each of these proportions and the root mean squared (RMS) of the placebo coefficients is calculated, which are then used to construct the critical value at a significance level of \code{alpha}. See Dette & Schumann (2024, s. 4.2.3.) for more details.
#' 
#' 
#' @seealso \code{\link[=print.rmsEquivTest]{print.rmsEquivTest}}
#' 
#' @return
#' An object of class "rmsEquivTest" containing:
#' \item{\code{placebo_coefficients}}{A numeric vector of the estimated placebo coefficients,}
#' \item{\code{rms_placebo_coefs}}{the root mean squared value of the placebo coefficients,}
#' \item{\code{significance_level}}{the significance level of the test,}
#' \item{\code{num_individuals}}{the number of cross-sectional individuals in the data,}
#' \item{\code{num_periods}}{the number of periods in the data,}
#' \item{\code{base_period}}{the base period in the data,}
#' \item{\code{equiv_threshold_specified}}{a logical value indicating whether an equivalence threshold was specified.}
#'
#' If \code{is.null(equiv_threshold)}, then additionally \code{minimum_equiv_threshold}: the minimum equivalence threshold for which the null hypothesis of non-negligible (based on the equivalence threshold) trend-differnces can be rejected. 
#' 
#' if \code{!(is.null(equiv_threshold))}, then additionally
#' \itemize{
#' \item \code{rms_critical_value}: the critical value at the alpha level,
#' \item \code{reject_null_hypothesis}: A logical value indicating whether to reject the null hypothesis,
#' \item \code{equiv_threshold}: the equivalence threshold specified.
#' }
#' 
#' @references
#' Dette, H., & Schumann, M. (2024). "Testing for Equivalence of Pre-Trends in Difference-in-Differences Estimation." \emph{Journal of Business & Economic Statistics}, 1–13. DOI: \href{https://doi.org/10.1080/07350015.2024.2308121}{10.1080/07350015.2024.2308121}
#' 
#' 
#' @examples
#' # Simulate data
#' sim_data <- sim_paneldata()
#' 
#' # Perform the test with equivalent threshold specified as 1 based on 
#' # pre-treatment periods 1:4 and homoscedastic error-terms:
#' equivalence_test <- rmsEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                                  data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                                  base_period = 4)
#' print(equivalence_test)
#' 
#' # Perform the test without specifying the equivalence threshold:
#' equivalence_test <- rmsEquivTest(Y = "Y", ID = "ID", G = "G", period = "period", 
#'                            data = sim_data, equiv_threshold = 1, pretreatment_period = 1:4,
#'                            base_period = 4)
#' print(equivalence_test)
#' 
#' 
#' @export
#'
#' 
rmsEquivTest <- function(Y, ID, G, period, X = NULL, data = NULL, equiv_threshold = NULL,  
                    pretreatment_period = NULL, base_period = NULL, 
                    alpha=0.05, no_lambda = 5){
  
  # rmsTest specific error checking:
  error_rmsTest <- rmsTest_error(alpha, no_lambda)
  if(error_rmsTest$error){stop(error_rmsTest$message) }
  
  # General error checking:
  error_test <- EquiTrends_inputcheck(Y, ID, G, period, X, data, equiv_threshold, pretreatment_period, 
                                     base_period, cluster = NULL, alpha)
  
  if(error_test$error){stop(error_test$message)}
  
  # Read the data:
  data_constr <- EquiTrends_dataconstr(Y, ID, G, period, X, data, pretreatment_period, 
                                       base_period, cluster=NULL)
  df <- data_constr$dataset
  base_period <- data_constr$baseperiod
  
  # Perform the test:
  test_results <- rmsTest_func(df, equiv_threshold, alpha, no_lambda, base_period)
  
  return(test_results)
}



