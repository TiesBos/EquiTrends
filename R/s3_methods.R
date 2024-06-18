# ---- Summary Function for the Intersection-Union Approach --------------------
#' @title Print maxTestIU objects
#'
#' @description Print method for objects of class 'maxTestIU'.
#'
#' @param x An object of class 'maxTestIU' containing the results of the maximum test based on the intersection-union approach.
#' @param ... urther arguments passed to or from other methods.
#' @method print maxTestIU
#' @return The function prints a summary of the results of the maximum test based on the intersection-union approach.
#' @export
#'
print.maxTestIU <- function(x, ...){
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
    cat("Minimum equivalence threshold to reject H0:", x$minimum.delta, "\n")
    df.print <- data.frame(as.numeric(formatC(x$absolute.placebo.coefs, format = "g", digits = 4)), 
                           as.numeric(formatC(x$standard.errors, format = "g", digits = 4)), 
                           as.numeric(formatC(x$minimum.deltas, format = "g", digits = 4)))
    colnames(df.print) <- c(" Estimate", "Std. Error ", " Minimum Equivalence Threshold")
    rownames(df.print) <- x$coef.names
    tibble.print <- tibble::as_tibble(df.print)
  }
  cat("---\n")
  tibble::print(tibble.print, ...)
  cat("---\n")
  
  # Summary statistics
  cat("No. placebo coefficients estimated (T):", length(x$absolute.placebo.coefs), "\n")
  cat("No. pre-treatment periods (T+1):", x$no.periods ,"\n")
  cat("Base period:", x$base.period ,"\n")
  cat("No. Individuals (N):", x$N, "\n")
  cat("\n")
}


#' @title Print maxTestBoot objects
#'
#' @param x An object of class 'maxTestBoot' containing the results of the maximum test based on the bootstrap procedure.
#' @param ... Further arguments passed to or from other methods.
#' @method print maxTestBoot
#' @return The function prints a summary of the results of the maximum test based on the bootstrap procedures.
#' @export
#'
print.maxTestBoot <- function(x, ...){
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
    colnames(output.df) <- c("Max. Abs. Coefficient", "    Bootstrap Critical Value", "  Reject H0")
    rownames(output.df) <- c("")
    print.default(output.df, ...)
  }
  
  # Summary statistics
  cat("---\n")
  cat("No. placebo coefficients estimated (T):", length(x$absolute.placebo.coefs), "\n")
  cat("No. pre-treatment periods (T+1):", x$no.periods ,"\n")
  cat("Base period:", x$base.period ,"\n")
  cat("No. Individuals (N):", x$N, "\n")
  cat("\n")
}