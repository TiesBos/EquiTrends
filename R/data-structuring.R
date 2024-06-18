EquiTrends.dataconstr <- function(Y, ID, G, period, X, data, pretreatment.period, 
                              base.period, cluster){
  colnames.X <- NULL
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
    colnames.X <- paste0("X_", 1:ncol(as.matrix(X)))
    new.data <- data.frame(Y=Y, ID=ID, period=period, G=G, X=X)
    colnames(new.data) <- c("Y", "ID", "period", "G", colnames.X)
  } else {
    new.data <- data.frame(Y=Y, ID=ID, period=period, G=G)
    colnames(new.data) <- c("Y", "ID", "period", "G")
  }
  
  # Add a variable indicating the cluster if supplied
  if(!is.null(cluster)){new.data[,"cluster"] <- cluster}
  
  # Omit any rows with NAs
  original.length <- nrow(new.data)
  new.data <- stats::na.omit(new.data)
  new.length <- nrow(new.data)
  omitted.rows <- original.length - new.length
  if(omitted.rows > 0){warning(paste(omitted.rows, "rows of pre-treatment data omitted due to NAs."))}
  
  
  # Only select the pre-specified pre-selection periods
  if(!base::is.null(pretreatment.period)){
    new.data <- new.data[new.data$period %in% pretreatment.period, ]
  }
  
  # Set the base period:
  if(base::is.null(base.period)){
    
    base.period <- base::max(new.data$period)
  }
  
  # Create the placebo variables:
  # Collect all pre-treatment periods in time
  unique.time <- unique(new.data$period)
  # Remove the base period:
  unique.time <- unique.time[unique.time != base.period]
  # Number of placebos:
  tt <- length(unique.time)
  
  # Create placebos:
  for(l in 1:tt){
    new.data[, paste0("placebo_", unique.time[l])] <- ifelse(new.data[, "period"]==unique.time[l], 1, 0)*new.data[, "G"]
  }
  
  # Turn the ID to a 1 to N scale:
  new.data$ID <- as.integer(factor(new.data$ID , levels = unique(new.data$ID)))
  # Turn the period column to a 1:(T+1) scale:
  period.numeric <- as.integer(factor(new.data$period, levels = unique(new.data$period)))
  new.data$period <- period.numeric
  
  # Return the final data.frame object used:
  return(list(dataset = new.data, baseperiod = base.period))
}



# ----------- Main Error Checking Function -------------------------------------
#   This function checks the input for the test functions. It takes as input:
EquiTrends.inputcheck <- function(Y, ID, G, period, X, data, delta, pretreatment.period, 
                                 base.period, cluster, alpha){
  
  if(is.null(data)){
    # Check if Y, ID, G, period and cluster are vectors:
    dim.Y <- dim(as.matrix(Y))[2]
    dim.ID <- dim(as.matrix(ID))[2]
    dim.G <- dim(as.matrix(G))[2]
    dim.period <- dim(as.matrix(period))[2]
    ncols.vec <- c(dim.Y, dim.ID, dim.G, dim.period)
    message.vec <- "Y, ID, period and G must be column vectors."
    if(!is.null(cluster)){
      dim.cluster <- dim(as.matrix(cluser))[2]
      ncols.vec <- c(ncols.vec, dim.cluster)
      message.vec <- "Y, ID, period, G and cluster must be column vectors."
    }
    
    if(any(ncols.vec != 1)){
      return(list(error=TRUE, message = message.vec))
    }
    
    # Check if X is not an empty matrix:
    if(!is.null(X) && dim(as.matrix(X))[1]==0){
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
      return(list(error=TRUE, message="all supplied vectors and matrix must have the same length"))
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
    
    
  }
  
  # Check if delta is strictly positive:
  if(!is.null(delta) && delta < 0){
    return(list(error=TRUE, message = "delta must be non-negative."))
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
  if (!is.null(pretreatment.period) && !all(pretreatment.period %in% period)) {
    return(list(error=TRUE, message="pretreatment.period must be a subset of period"))
  }
  
  if(!is.null(base.period) && length(base.period) != 1){
    return(list(error=TRUE, message: "base.period must be a scalar."))
  }
  
  # base period must lie in time.prior:
  if(is.null(pretreatment.period)){pretreatment.period <- period}
  if(!is.null(base.period) && !(base.period%in%pretreatment.period)){
    return(list(error=TRUE, message = "base.period can not be found not in the pre-treatment period."))
  }
  
  # Warnings:
  
  return(list(error=FALSE))
}