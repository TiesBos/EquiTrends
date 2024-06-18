#' @title Equivalence Test for Pre-trends based on the Maximum Absolute Coefficient 
#'
#' @description This function performs an equivalence test for pre-trends based on the maximum absolute placebo coefficient. The test can be performed using the intersection-union approach (IU), a bootstrap procedure for spherical errors (Boot) and a wild bootstrap procedure (Wild).
#'
#' @param Y If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the numeric dependent (outcome) variable in ’data’. If 'data' is not supplied, a numeric vector with the variable of interest.
#' @param ID If 'data' is supplied, a scalar identifying the column number or column-name character string that corresponds to the unit numbers in ’data’. If 'data' is not supplied, a numeric vector (of the same dimension as Y) containing the unit numbers of the observations.
#' @param G If 'data' is supplied, a scalar identifying the column number or column-name character string associated to the binary or logic variable indicating if the individual receives treatment (e.g. 1 or TRUE) or not (0 or FALSE). If 'data' is not supplied, a vector (of the same dimension as Y) binary or logic indicating if the individual (e.g. the ID vector).
#' @param period If 'data' is supplied, a scalar identifying the column number or column-name character string associated with period (time) data. The time variable has to be numeric. If 'data' is not supplied, a numeric vector (of the same dimension as Y) indicating time.
#' @param X  If 'data' is supplied, a vector of column numbers or column-name character strings that identifies the control variables’ columns. If data is not supplied, a vector, matrix or data.frame containing the control variables.
#' @param data An optional data.frame object containing the variables in Y, ID, G, T and, if supplied, X and cluster as its columns.
#' @param delta The scalar equivalence threshold (must be positive). The default is NULL, implying that the function must look for the minimum value for which the null of ”non-negligible differences” can still be rejected.
#' @param pretreatment.period A numeric vector identifying the pre-treatment periods that should be used for testing. The default is to use all periods that are included in T.
#' @param base.period The pre-treatment period to compare the post-treatment observation to. The default is to take the last specified pre-treatment period.
#' @param vcov The variance-covariance matrix that needs to be used. See details for more details.
#' @param cluster If vcov = "CL", a vector indicating which observations belong to the same cluster of the same length as Y. If 'data' is supplied, 'cluster' must be either the column index or column name of this vector in the data.frame/matrix. The default (cluster=NULL) assumes every unit in ID is its own cluster.
#' @param alpha Significance level of the test. The default is 0.05.
#' @param type The type of maximum test that should be performed. "IU" for the intersection-union test, "Boot" for the regular bootstrap procedure from Dette & Schumann (2023) and "Wild" for the Wild bootstrap procedure.
#' @param B If type = Boot or type = Wild, the number of bootstrap samples used. The default is 1000.
#' @param verbose A logical object indicating if test results need to be printed. The default is TRUE.
#' 
#' @importFrom Rcpp sourceCpp 
#' @importFrom RcppParallel RcppParallelLibs
#' 
#'
#' @return hoi
#' @export
#'
#' 
maxTest <- function(Y, ID, G, period, X = NULL, data = NULL, delta = NULL,  
                       pretreatment.period = NULL, base.period = NULL, 
                       vcov = NULL, cluster = NULL, alpha = 0.05, 
                       type = c("IU", "Boot", "Wild"), B = 1000, verbose=TRUE){
  # If no type is specified, the type is "IU"
  if(identical(type, c("IU", "Boot", "Wild"))){
    type <- "IU"
    warning("type is not specified: Intersection-Union (IU) approach is used.")
  }
  
  # Error/Warnings:
  # Checking errors specific to this test procedure
  error.maxTest <- maxTest.error(type, delta, vcov)
  if(error.maxTest$error){stop(error.maxTest$message)}
  
  # General error checking:
  error.test <- EquiTrends.inputcheck(Y, ID, G, period, X, data, delta, pretreatment.period, 
                                     base.period, cluster, alpha)
  
  if(error.test$error){stop(error.test$message)}
  
  # Structuring the data:
  data.constr <- EquiTrends.dataconstr(Y, ID, G, period, X, data, pretreatment.period, base.period,
                                   cluster)
  # The dataframe:
  df <- data.constr$dataset
  
  # The base period:
  base.period <- data.constr$baseperiod
  
  # Number of individuals
  N <- length(df[unique(df$ID), "ID"])
  # Number of Time Periods:
  no.periods <- length(df[unique(df$period), "period"])
  
  if(type=="IU"){
    test_results <- maxTestIU(df, delta, vcov, cluster, alpha, N, no.periods, base.period)
    class(test_results) <- "maxTestIU"
  } else if (type == "Boot" || type == "Wild"){
    test_results <- maxTestBoot(df, delta, alpha, N, B, no.periods, base.period, type)
    class(test_results) <- "maxTestBoot"
  }
  
  if(verbose){print(test_results)}
  return(invisible(test_results))
}

