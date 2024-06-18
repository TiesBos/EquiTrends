# A function to print data.frame objects:
pretty_print <- function(df) {
  # Get the column names
  col_names <- colnames(df)
  
  # Determine the maximum width for each column
  max_widths <- sapply(df, function(col) max(nchar(as.character(col))))
  max_widths <- pmax(max_widths, nchar(col_names))
  
  # Format the column names
  formatted_col_names <- mapply(function(name, width) format(name, width = width), name = col_names, width = max_widths)
  
  # Print the column names
  cat(paste(formatted_col_names, collapse = "\t"), "\n")
  
  # Check if df has only one row
  if(nrow(df) == 1){
    # Format and print the row
    formatted_row <- mapply(function(col, width) format(as.character(col), width = width), col = df[1, ], width = max_widths)
    cat(paste(formatted_row, collapse = "\t"), "\n")
  } else {
    # Format the data frame
    formatted_df <- apply(df, 2, function(col, width) format(as.character(col), width = width), width = max_widths)
    
    # Print each row
    apply(formatted_df, 1, function(row) {
      cat(paste(row, collapse = "\t"), "\n")
    })
  }
}

# ---- Summary Function for the Intersection-Union Approach --------------------
#' @title Print maxEquivTestIU objects
#'
#' @description Print method for objects of class 'maxEquivTestIU'.
#'
#' @param x An object of class 'maxEquivTestIU' containing the results of the maximum test based on the intersection-union approach.
#' @param ... Further arguments passed to or from other methods.
#' @method print maxEquivTestIU
#' @return The function prints a summary of the results of the maximum test based on the intersection-union approach.
#' @export
#'
print.maxEquivTestIU <- function(x, ...){
  cat("\n")
  width <- getOption("width")
  title <- "Dette & Schumann (2023) Equivalence Tests for Pre-trends in DiD Estimation"
  
  # Check if console width is less than title length
  if (width < nchar(title)) {
    width <- nchar(title) + 2  # Adjust width to be slightly larger than title
  }
  
  separator <- strrep(" ", floor((width - nchar(title)) / 2))
  
  # Centered title
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  cat(separator, title, "\n", sep = "")
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  
  cat("Type: Intersection Union \n")
  cat("Significance level:", x$sign.level, "\n")
  if(x$delta.specified){  
    cat("Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of", x$delta, ".\n")
    df.print <- data.frame(as.numeric(formatC(x$absolute.placebo.coefs, format = "g", digits = 4)),
                           as.numeric(formatC(x$standard.errors, format = "g", digits = 4)),
                           as.numeric(formatC(x$critical.values, format = "g", digits = 4)),
                           as.numeric(formatC(x$p.values, format = "g", digits = 4)))
    
    # Convert the formatted values back to numeric
    #df.print <- apply(df.print, MARGIN = c(1,2), FUN = function(x){as.numeric(x)})
    rownames(df.print) <- x$coef.names
    colnames(df.print) <- c("Abs. Estimate", "Std. Error", "Critical Value", "p-value")
  } else {
    cat("Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold.\n")
    cat("Minimum equivalence threshold to accept the alternative:", formatC(x$minimum.delta, format = "g", digits = 4), "\n")
    df.print <- data.frame(as.numeric(formatC(x$absolute.placebo.coefs, format = "g", digits = 4)), 
                           as.numeric(formatC(x$standard.errors, format = "g", digits = 4)), 
                           as.numeric(formatC(x$minimum.deltas, format = "g", digits = 4)))
    colnames(df.print) <- c(" Estimate", "Std. Error ", " Minimum Equivalence Threshold")
    rownames(df.print) <- x$coef.names
  }
  cat("---\n")
  pretty_print(df.print)
  cat("---\n")
  
  # Summary statistics
  cat("No. placebo coefficients estimated (T):", length(x$absolute.placebo.coefs), "\n")
  cat("No. pre-treatment periods (T+1):", x$no.periods ,"\n")
  cat("Base period:", x$base.period ,"\n")
  cat("No. Individuals (N):", x$N, "\n")
  cat("\n")
}


#' @title Print maxEquivTestBoot objects
#'
#' @param x An object of class 'maxEquivTestBoot' containing the results of the maximum test based on the bootstrap procedure.
#' @param ... Further arguments passed to or from other methods.
#' @method print maxEquivTestBoot
#' @return The function prints a summary of the results of the maximum test based on the bootstrap procedures.
#' @export
#'
print.maxEquivTestBoot <- function(x, ...){
  cat("\n")
  width <- getOption("width")
  title <- "Dette & Schumann (2023) Equivalence Tests for Pre-trends in DiD Estimation"
  
  # Check if console width is less than title length
  if (width < nchar(title)) {
    width <- nchar(title) + 2  # Adjust width to be slightly larger than title
  }
  
  separator <- strrep(" ", floor((width - nchar(title)) / 2))
  
  # Centered title
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  cat(separator, title, "\n", sep = "")
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  
  # Type of Bootstrap
  if (!x$wild) {
    cat("Type: Bootstrap for Spherical Errors  (Based on ", x$B, " bootstrap samples)\n", sep="")
  } else {
    cat("Type: Cluster Wild Bootstrap (Based on ", x$B, " bootstrap samples)\n", sep="")
  }
  
  # Additional details
  cat("Significance level:", x$sign.level, "\n")
  cat("Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of", x$delta, ".\n")
  cat("---\n")
  
  # Display results
  if (x$delta.specified) {
    output.df <- data.frame(as.numeric(formatC(x$max.coefs, format = "g", digits = 4)), 
                            as.numeric(formatC(x$critical.values, format = "g", digits = 4)),
                            x$reject.H0)
    colnames(output.df) <- c("Max. Abs. Coefficient", "Bootstrap Critical Value", "Reject H0")
    rownames(output.df) <- c("")
    pretty_print(output.df)
  }
  
  # Summary statistics
  cat("---\n")
  cat("No. placebo coefficients estimated (T):", length(x$absolute.placebo.coefs), "\n")
  cat("No. pre-treatment periods (T+1):", x$no.periods ,"\n")
  cat("Base period:", x$base.period ,"\n")
  cat("No. Individuals (N):", x$N, "\n")
  cat("\n")
}

#' @title Print meanEquivTest objects
#'
#' @param x An object of class 'meanEquivTest' containing the results of the maximum test based on the bootstrap procedure.
#' @param ... Further arguments passed to or from other methods.
#' @method print meanEquivTest
#' @return The function prints a summary of the results of the maximum test based on the bootstrap procedures.
#' @export
#'
print.meanEquivTest <- function(x, ...){
  cat("\n")
  width <- getOption("width")
  title <- "Dette & Schumann (2023) Equivalence Tests for Pre-trends in DiD Estimation"
  
  # Check if console width is less than title length
  if (width < nchar(title)) {
    width <- nchar(title) + 2  # Adjust width to be slightly larger than title
  }
  
  separator <- strrep(" ", floor((width - nchar(title)) / 2))
  
  # Centered title
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  cat(separator, title, "\n", sep = "")
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  
  cat("Type: Mean Placebo Effect \n")
  if(x$delta.specified){
    cat("Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold of", x$delta, ".\n")
    
    df.print <- data.frame(as.numeric(formatC(x$abs.mean.placebo, format = "g", digits = 4)),
                           as.numeric(formatC(sqrt(x$placebo_var), format = "g", digits = 4)), 
                           as.numeric(formatC(x$p.value, format = "g", digits = 4)))
    colnames(df.print) <- c("Abs. Mean Placebo Effect", "Std. Error", "p-value")
    rownames(df.print) <- c("")
  } else {
    cat("Significance level:", x$sign.level, "\n")
    cat("Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold.\n")
    df.print <- data.frame(as.numeric(formatC(x$abs.mean.placebo, format = "g", digits = 4)), 
                           as.numeric(formatC(sqrt(x$placebo_var), format = "g", digits = 4)), 
                           as.numeric(formatC(x$minimum.delta, format = "g", digits = 4)))
    colnames(df.print) <- c("Abs. Mean Placebo Effect", "Std. Error", "Min. Equiv. Threshold")
    rownames(df.print) <- c("")  
  }
  cat("---\n")
  pretty_print(df.print)
  cat("---\n")
  
  # Data statistics
  cat("No. placebo coefficients estimated (T):", length(x$placebo.coefs), "\n")
  cat("No. pre-treatment periods (T+1):", x$no.periods ,"\n")
  cat("Base period:", x$base.period ,"\n")
  cat("No. Individuals (N):", x$N, "\n")
  cat("\n")
}  

#' @title Print rmsEquivTest objects
#'
#' @param x An object of class 'rmsEquivTest' containing the results of the maximum test based on the bootstrap procedure.
#' @param ... Further arguments passed to or from other methods.
#' @method print rmsEquivTest
#' @return The function prints a summary of the results of the maximum test based on the bootstrap procedures.
#' @export
#'
print.rmsEquivTest <- function(x, ...){
  cat("\n")
  width <- getOption("width")
  title <- "Dette & Schumann (2023) Equivalence Tests for Pre-trends in DiD Estimation"
  
  # Check if console width is less than title length
  if (width < nchar(title)) {
    width <- nchar(title) + 2  # Adjust width to be slightly larger than title
  }
  
  separator <- strrep(" ", floor((width - nchar(title)) / 2))
  
  # Centered title
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  cat(separator, title, "\n", sep = "")
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  
  cat("Type: Mean Placebo Effect \n")
  if(x$delta.specified){
    cat("Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold of", x$delta, ".\n")
    
    df.print <- data.frame(as.numeric(formatC(x$RMS, format = "g", digits = 4)),
                           as.numeric(formatC(x$critical.value, format = "g", digits = 4)))
    colnames(df.print) <- c("RMS Placebo Effect", "Simulated Crit. Val.")
    rownames(df.print) <- c("")
  } else {
    cat("Significance level:", x$sign.level, "\n")
    cat("Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold.\n")
    df.print <- data.frame(as.numeric(formatC(x$RMS, format = "g", digits = 4)), 
                           as.numeric(formatC(x$minimum.delta, format = "g", digits = 4)))
    colnames(df.print) <- c("RMS Placebo Effect", "Min. Equiv. Threshold")
    rownames(df.print) <- c("")  
  }
  cat("---\n")
  pretty_print(df.print)
  cat("---\n")
  
  # Data statistics
  cat("No. placebo coefficients estimated (T):", length(x$placebo.coefs), "\n")
  cat("No. pre-treatment periods (T+1):", x$no.periods ,"\n")
  cat("Base period:", x$base.period ,"\n")
  cat("No. Individuals (N):", x$N, "\n")
  cat("\n")
}  


