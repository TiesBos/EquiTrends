# ----------- Data Construction Function ---------------------------------------
#' @title Data Construction Function for EquiTrends
#'
#' @param Y see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param ID see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param G see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param period see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param X see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param data see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param pretreatment_period see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param base_period see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param cluster see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#'
#' @return
#' A list containing the structured data.frame object used in the equivalence testing procedures, the base period for the test, a logical value indicating whether the panel is balanced and the number of periods.

EquiTrends_dataconstr <- function(Y, ID, G, period, X, data, pretreatment_period, 
                              base_period, cluster){
  orig_names <- NULL
  colnames_X <- NULL
  if(!base::is.null(data)){
    names_data <- base::colnames(data)
    Y <- data[, Y]
    ID <- data[, ID]
    G <- data[, G]
    period <- data[, period]
    if(!is.null(X)){X <- data[,X]}
    if(!is.null(cluster)){cluster <- data[,cluster]}
  }
  
  
  # Bind the input together:
  if(!is.null(X)){
    orig_names <- colnames(X)
    colnames_X <- paste0("X_", 1:ncol(as.matrix(X)))
    new_data <- data.frame(Y=Y, ID=ID, period=period, G=G, X=X)
    colnames(new_data) <- c("Y", "ID", "period", "G", colnames_X)
  } else {
    new_data <- data.frame(Y=Y, ID=ID, period=period, G=G)
    colnames(new_data) <- c("Y", "ID", "period", "G")
  }
  if(is.null(orig_names)){orig_names <- colnames_X}
  
  # Add a variable indicating the cluster if supplied
  if(!is.null(cluster)){new_data[,"cluster"] <- cluster}
  
  # Omit any rows with NAs
  original_length <- nrow(new_data)
  new_data <- stats::na.omit(new_data)
  new_length <- nrow(new_data)
  omitted_rows <- original_length - new_length
  if(omitted_rows > 0){warning(paste(omitted_rows, "rows of pre-treatment data omitted due to NAs."))}
  
  
  # Only select the pre-specified pre-selection periods
  if(!base::is.null(pretreatment_period)){
    new_data <- new_data[new_data$period %in% pretreatment_period, ]
  }
  
  # Set the base period:
  if(base::is.null(base_period)){
    base_period <- base::max(new_data$period)
  }
  
  # Create the placebo variables:
  # Collect all pre-treatment periods in time
  unique_time <- unique(new_data$period)
  # Remove the base period:
  unique_time <- unique_time[unique_time != base_period]
  # Number of placebos:
  tt <- length(unique_time)
  
  # Create placebos:
  for(l in 1:tt){
    new_data[, paste0("placebo_", unique_time[l])] <- ifelse(new_data[, "period"]==unique_time[l], 1, 0)*new_data[, "G"]
  }
  
  # Turn the ID to a 1 to N scale:
  new_data$ID <- as.integer(factor(new_data$ID , levels = unique(new_data$ID)))
  # Turn the period column to a 1:(T+1) scale:
  # period_numeric <- as.integer(factor(new_data$period, levels = unique(new_data$period)))
  # new_data$period <- period_numeric
  
  # Check if the data is balanced:
  balanced_panel_test <- is_panel_balanced(new_data)
  balanced_panel <- balanced_panel_test$balanced
  no_periods <- balanced_panel_test$no_periods
  
  # Return the final data.frame object used:
  return(list(dataset = new_data, baseperiod = base_period, orig_names = orig_names, 
              balanced_panel = balanced_panel, no_periods = no_periods))
}



# ----------- Main Error Checking Function -------------------------------------
#   This function checks the input for the test functions. It takes as input:
#' @title Input Checks Function for EquiTrends
#'
#' @param Y see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param ID see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param G see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param period see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param X see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param data see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param equiv_threshold see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param pretreatment_period see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param base_period see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param cluster see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#' @param alpha see \link[EquiTrends]{maxEquivTest}, \link[EquiTrends]{meanEquivTest} or \link[EquiTrends]{rmsEquivTest}
#'
#' @return
#' A list containing an error indicator and a message. If \code{error} is TRUE, \code{message} contains an error message. If \code{error} is FALSE, \code{message} is empty.
EquiTrends_inputcheck <- function(Y, ID, G, period, X, data, equiv_threshold, pretreatment_period, 
                                 base_period, cluster, alpha){
  
  if(is.null(data)){
    # Check if Y, ID, G, period and cluster are vectors:
    dim_Y <- dim(as.matrix(Y))[2]
    dim_ID <- dim(as.matrix(ID))[2]
    dim_G <- dim(as.matrix(G))[2]
    dim_period <- dim(as.matrix(period))[2]
    ncols_vec <- c(dim_Y, dim_ID, dim_G, dim_period)
    message_vec <- "Y, ID, period and G must be column vectors."
    if(!is.null(cluster)){
      dim_cluster <- dim(as.matrix(cluster))[2]
      ncols_vec <- c(ncols_vec, dim_cluster)
      message_vec <- "Y, ID, period, G and cluster must be column vectors."
    }
    
    if(any(ncols_vec != 1)){
      return(list(error=TRUE, message = message_vec))
    }
    
    # Check if X is not an empty matrix:
    if(!is.null(X) && dim(as.matrix(X))[2]==0){
      return(list(error=TRUE, message = "X must have at least one column."))
    }
    
    # Check if Y, ID, G, period, cluster and X are of equal length:
    # Get the lengths of the vectors
    lengths <- c(length(Y), length(G), length(ID), length(period))
    # If X is not NULL, include the number of rows of X in the lengths vector
    if (!is.null(X)) {
      lengths <- c(lengths, nrow(as.matrix(X)))
    }
    # If cluster is not NULL, include its length in the lengths vector
    if (!is.null(cluster)) {
      lengths <- c(lengths, length(cluster))
    }
    
    if(any(lengths != lengths[1])){
      return(list(error=TRUE, message="all supplied vectors and matrix must have equal length"))
    }
    
    if(any(lengths == 1)){
      warning("data possibly missing: check if `data` should be specified")
    }
    
  } else {
    # Check if the data is of the require data.frame form:
    if(!is.data.frame(data)){
      return(list(error=TRUE, message = "data must be an object of class data.frame"))
    }
    
    # For the remaining tests, we construct the variables quickly:
    Y <- data[, Y]
    ID <- data[, ID]
    period <- data[, period]
    G <- data[, G]
    if(!is.null(X)){X <- data[,X]}
    if(!is.null(cluster)){cluster <- data[,cluster]}
    
    # Check if Y, ID, G, period and cluster are vectors:
    dim_Y <- dim(as.matrix(Y))[2]
    dim_ID <- dim(as.matrix(ID))[2]
    dim_G <- dim(as.matrix(G))[2]
    dim_period <- dim(as.matrix(period))[2]
    ncols_vec <- c(dim_Y, dim_ID, dim_G, dim_period)
    message_vec <- "Y, ID, period and G must be column vectors."
    if(!is.null(cluster)){
      dim_cluster <- dim(as.matrix(cluster))[2]
      ncols_vec <- c(ncols_vec, dim_cluster)
      message_vec <- "Y, ID, period, G and cluster must be column vectors."
    }
    if(any(ncols_vec != 1)){
      return(list(error=TRUE, message = message_vec))
    }
    
    
    
  }
  
  # Check if delta is strictly positive:
  if(!is.null(equiv_threshold) && equiv_threshold < 0){
    return(list(error=TRUE, message = "equiv_threshold must be non-negative."))
  }
  
  # Check if alpha is between 0 and 1:
  if(alpha <= 0 || alpha >= 1){
    return(list(error = TRUE, message = "alpha must lie between 0 and 1."))
  }
  
  # Check if G is binary or logic:
  if(any(!(G%in%c(0, 1, TRUE, FALSE)))){
    return(list(error=TRUE, message = "Entries of G must either be logical (e.g. TRUE/FALSE) or binary (e.g. 0/1)."))
  }
  
  # Check if T is numeric:
  if (!is.numeric(period)) {
    return(list(error=TRUE, message = "period must be numeric"))
  }
  
  # time prior must be a subset of T:
  if (!is.null(pretreatment_period) && !all(pretreatment_period %in% period)) {
    return(list(error=TRUE, message="pretreatment_period must be a subset of period"))
  }
  
  if(!is.null(base_period) && length(base_period) != 1){
    return(list(error=TRUE, message = "base_period must be a scalar."))
  }
  
  # base period must lie in pretreatment_period:
  if(is.null(pretreatment_period)){pretreatment_period <- period}
  if(!is.null(base_period) && !(base_period%in%pretreatment_period)){
    return(list(error=TRUE, message = "base_period must be an element of pretreatment_period."))
  }
  
  # the pretreatment period must have at least two periods:
  if(length(unique(pretreatment_period)) < 2){
    return(list(error=TRUE, message = "pre-treatment period must have at least two unique periods."))
  }
  
  return(list(error=FALSE))
}

# ----------- Checking if the panel is balanced ---------------------------------

#' @importFrom dplyr %>%
#' @importFrom rlang .data
is_panel_balanced <- function(df) {
  # Group by the individual variable and list all unique periods for each individual
  periods_by_individual <- df %>%
    dplyr::group_by(.data$ID) %>%
    dplyr::summarise(periods = list(unique(.data$period)))
  
  # Check if all individuals have the same unique periods
  unique_periods <- periods_by_individual$periods[[1]]
  balanced <- all(sapply(periods_by_individual$periods, function(x) setequal(x, unique_periods)))
  
  if(balanced){
    no_periods <- length(unique(df$period))
  } else {
    no_periods_per_individual <- sapply(periods_by_individual$periods, length)
    no_periods <- c(min(no_periods_per_individual), max(no_periods_per_individual))
  }
  
  return(list(balanced=balanced, no_periods = no_periods))
}