#' @title Equivalence Test for Pre-trends based on the Maximum Absolute Placebo Coefficient 
#'
#' @description This function performs an equivalence test for pre-trends based on the maximum absolute placebo coefficient. The test can be performed using the intersection-union approach (IU), a bootstrap procedure for spherical errors (Boot) and a wild bootstrap procedure (Wild).
#'
#' @param Y If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the numeric dependent (outcome) variable in ’data’. If 'data' is not supplied, a numeric vector with the variable of interest.
#' @param ID If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the unit numbers in ’data’. If 'data' is not supplied, a numeric vector (of the same dimension as Y) containing the unit numbers of the observations.
#' @param G If 'data' is supplied, a scalar identifying the column number or column-name character string associated to the binary or logic variable indicating if the individual receives treatment (e.g. 1 or TRUE) or not (0 or FALSE). If 'data' is not supplied, a vector (of the same dimension as Y) binary or logic indicating if the individual (e.g. the ID vector).
#' @param period If 'data' is supplied, a scalar identifying the column number or column-name character string associated with period (time) data. The time variable has to be numeric. If 'data' is not supplied, a numeric vector (of the same dimension as Y) indicating time.
#' @param X  If 'data' is supplied, a vector of column numbers or column-name character strings that identifies the control variables’ columns. If data is not supplied, a vector, matrix or data.frame containing the control variables.
#' @param data An optional data.frame object containing the variables in Y, ID, G, T and, if supplied, X and cluster as its columns.
#' @param equiv_threshold The scalar equivalence threshold (must be positive). The default is NULL, implying that the function must look for the minimum value for which the null of ”non-negligible differences” can be rejected.
#' @param pretreatment_period A numeric vector identifying the pre-treatment periods that should be used for testing. The default is to use all periods that are included in T.
#' @param base_period The pre-treatment period to compare the post-treatment observation to. The default is to take the last specified pre-treatment period.
#' @param vcov The variance-covariance matrix that needs to be used. See details for more details.
#' @param cluster If vcov = "CL", a vector indicating which observations belong to the same cluster of the same length as Y. If 'data' is supplied, 'cluster' must be either the column index or column name of this vector in the data.frame/matrix. The default (cluster=NULL) assumes every unit in ID is its own cluster.
#' @param alpha Significance level of the test. The default is 0.05.
#' @param type The type of maximum test that should be performed. "IU" for the intersection-union test, "Boot" for the regular bootstrap procedure from Dette & Schumann (2023) and "Wild" for the Wild bootstrap procedure.
#' @param B If type = Boot or type = Wild, the number of bootstrap samples used. The default is 1000.
#' 
#' @details DETAILS ON "vcov" parameter.
#'
#' @return hoi
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
  error_maxTest <- maxTest.error(type, equiv_threshold, vcov)
  if(error_maxTest$error){stop(error_maxTest$message)}
  
  # General error checking:
  error_test <- EquiTrends.inputcheck(Y, ID, G, period, X, data, equiv_threshold, pretreatment_period, 
                                     base_period, cluster, alpha)
  
  if(error_test$error){stop(error_test$message)}
  
  # Structuring the data:
  data_constr <- EquiTrends.dataconstr(Y, ID, G, period, X, data, pretreatment_period, base_period,
                                   cluster)
  # The dataframe:
  df <- data_constr$dataset
  
  # The base period:
  base_period <- data_constr$baseperiod
  
  # Number of individuals
  N <- length(df[unique(df$ID), "ID"])
  # Number of Time Periods:
  no_periods <- length(df[unique(df$period), "period"])
  
  if(type=="IU"){
    test_results <- maxTestIU(df, equiv_threshold, vcov, cluster, alpha, N, no_periods, base_period)
  } else if (type == "Boot" || type == "Wild"){
    test_results <- maxTestBoot(df, equiv_threshold, alpha, N, B, no_periods, base_period, type)
  }
  
  return(test_results)
}

#' @title Equivalence Test for Pre-trends based on the Mean Placebo Coefficient 
#'
#' @description This function performs an equivalence test for pre-trends based on the maximum absolute placebo coefficient. The test can be performed using the intersection-union approach (IU), a bootstrap procedure for spherical errors (Boot) and a wild bootstrap procedure (Wild).
#'
#' @param Y If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the numeric dependent (outcome) variable in ’data’. If 'data' is not supplied, a numeric vector with the variable of interest.
#' @param ID If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the unit numbers in ’data’. If 'data' is not supplied, a numeric vector (of the same dimension as Y) containing the unit numbers of the observations.
#' @param G If 'data' is supplied, a scalar identifying the column number or column-name character string associated to the binary or logic variable indicating if the individual receives treatment (e.g. 1 or TRUE) or not (0 or FALSE). If 'data' is not supplied, a vector (of the same dimension as Y) binary or logic indicating if the individual (e.g. the ID vector).
#' @param period If 'data' is supplied, a scalar identifying the column number or column-name character string associated with period (time) data. The time variable has to be numeric. If 'data' is not supplied, a numeric vector (of the same dimension as Y) indicating time.
#' @param X  If 'data' is supplied, a vector of column numbers or column-name character strings that identifies the control variables’ columns. If data is not supplied, a vector, matrix or data.frame containing the control variables.
#' @param data An optional data.frame object containing the variables in Y, ID, G, T and, if supplied, X and cluster as its columns.
#' @param equiv_threshold The scalar equivalence threshold (must be positive). The default is NULL, implying that the function must look for the minimum value for which the null of ”non-negligible differences” can still be rejected.
#' @param pretreatment_period A numeric vector identifying the pre-treatment periods that should be used for testing. The default is to use all periods that are included in T.
#' @param base_period The pre-treatment period to compare the post-treatment observation to. The default is to take the last specified pre-treatment period.
#' @param vcov The variance-covariance matrix that needs to be used. See details for more details.
#' @param cluster If vcov = "CL", a vector indicating which observations belong to the same cluster of the same length as Y. If 'data' is supplied, 'cluster' must be either the column index or column name of this vector in the data.frame/matrix. The default (cluster=NULL) assumes every unit in ID is its own cluster.
#' @param alpha Significance level of the test. The default is 0.05.
#'
#' @return hoi
#' @export
#'
meanEquivTest <- function(Y, ID, G, period, X = NULL, data = NULL, equiv_threshold = NULL,  
                     pretreatment_period = NULL, base_period = NULL, 
                     vcov = NULL, cluster = NULL, alpha=0.05){
  
  # Error/Warnings:
  # General error checking:
  error_test <- EquiTrends.inputcheck(Y, ID, G, period, X, data, equiv_threshold, pretreatment_period, 
                                     base_period, cluster, alpha)
  
  if(error_test$error){stop(error_test$message)}
  
  # We construct the data.frame:
  data_constr <- EquiTrends.dataconstr(Y, ID, G, period, X, data, pretreatment_period, base_period,
                                   cluster)
  df <- data_constr$data
  
  # Store the base.period:
  base.period <- data_constr$baseperiod
  
  # number of individuals in the sample:
  N <- length(df[unique(df$ID), "ID"])
  # The number of periods:
  no_periods <- length(unique(df$period))
  
  # Perform the test:
  results <- meanTest.func(df, equiv_threshold, vcov, cluster, alpha, N, no_periods, base_period)
  
  return(results)
}



#' @title Equivalence Test for Pre-trends based on the RMS Placebo Coefficient 
#' @param Y If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the numeric dependent (outcome) variable in ’data’. If 'data' is not supplied, a numeric vector with the variable of interest.
#' @param ID If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the unit numbers in ’data’. If 'data' is not supplied, a numeric vector (of the same dimension as Y) containing the unit numbers of the observations.
#' @param G If 'data' is supplied, a scalar identifying the column number or column-name character string associated to the binary or logic variable indicating if the individual receives treatment (e.g. 1 or TRUE) or not (0 or FALSE). If 'data' is not supplied, a vector (of the same dimension as Y) binary or logic indicating if the individual (e.g. the ID vector).
#' @param period If 'data' is supplied, a scalar identifying the column number or column-name character string associated with period (time) data. The time variable has to be numeric. If 'data' is not supplied, a numeric vector (of the same dimension as Y) indicating time.
#' @param X  If 'data' is supplied, a vector of column numbers or column-name character strings that identifies the control variables’ columns. If data is not supplied, a vector, matrix or data.frame containing the control variables.
#' @param data An optional data.frame object containing the variables in Y, ID, G, T and, if supplied, X and cluster as its columns.
#' @param equiv_threshold The scalar equivalence threshold (must be positive). The default is NULL, implying that the function must look for the minimum value for which the null of ”non-negligible differences” can still be rejected.
#' @param pretreatment_period A numeric vector identifying the pre-treatment periods that should be used for testing. The default is to use all periods that are included in T.
#' @param base_period The pre-treatment period to compare the post-treatment observation to. The default is to take the last specified pre-treatment period.
#' @param no_lambda The scaling factor of the uniform distribution
#' @param alpha Significance level of the test. Must be one of the following: 0.01, 0.025, 0.05, 0.1 or 0.2. The default is 0.05.
#'
#' @return hoi
#' @export
#'
#' 
rmsEquivTest <- function(Y, ID, G, period, X = NULL, data = NULL, equiv_threshold = NULL,  
                    pretreatment_period = NULL, base_period = NULL, 
                    alpha=0.05, no_lambda = 5){
  
  # rmsTest specific error checking:
  error_rmsTest <- rmsTest.error(alpha, no_lambda)
  if(error_rmsTest$error){stop(error_rmsTest$message) }
  
  # General error checking:
  error_test <- EquiTrends.inputcheck(Y, ID, G, period, X, data, equiv_threshold, pretreatment_period, 
                                     base_period, cluster = NULL, alpha)
  
  if(error_test$error){stop(error_test$message)}
  
  # Read the data:
  data_constr <- EquiTrends.dataconstr(Y, ID, G, period, X, data, pretreatment_period, 
                                       base_period, cluster=NULL)
  df <- data_constr$dataset
  base.period <- data_constr$baseperiod
  
  # Perform the test:
  test_results <- rmsTest.func(df, equiv_threshold, alpha, no_lambda, base_period)
  
  return(test_results)
}



