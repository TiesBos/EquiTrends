# A Function for nice representation of small numbers:
number_rep <- function(x){
  new_x <- ifelse(abs(x) < 1e-4, format(x, scientific = TRUE), format(x, digits = 4))
  new_x <- ifelse(abs(as.numeric(new_x)) < 2e-16,  "<2e-16", new_x)
  return(new_x)
}

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
  title <- "Equivalence Tests for Pre-trends in DiD Estimation"
  
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
  if(x$equiv_threshold_specified){  
    cat("Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of", x$equiv_threshold, ".\n")
    cat("Reject null hypothesis:" , x$reject_null_hypothesis, "\n")
    cat("( Critical values are printed for the significance level:", x$significance_level, ")\n")
    df_print <- data.frame(number_rep(x$abs_placebo_coefficients),
                           number_rep(x$placebo_coefficients_se),
                           number_rep(x$IU_critical_values))
    
    rownames(df_print) <- x$placebo_coef_names
    colnames(df_print) <- c("Abs. Estimate", "Std. Error", "Critical Value")
  } else {
    cat("Significance level:", x$significance_level, "\n")
    cat("Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold.\n")
    cat("Minimum equivalence threshold to accept the alternative:", number_rep(x$minimum_equiv_threshold), "\n")
    df_print <- data.frame(number_rep(x$abs_placebo_coefficients), 
                           number_rep(x$placebo_coefficients_se), 
                           number_rep(x$minimum_equiv_thresholds))
    colnames(df_print) <- c(" Estimate", "Std. Error ", " Minimum Equivalence Threshold")
    rownames(df_print) <- x$placebo_coef_names
  }
  cat("---\n")
  pretty_print(df_print)
  cat("---\n")
  
  # Summary statistics
  cat("No. placebo coefficients estimated:", length(x$abs_placebo_coefficients), "\n")
  cat("Base period:", x$base_period ,"\n")
  cat(" \n")
  if(x$is_panel_balanced){
    cat("Balanced Panel: \n")
    cat(" + No. pre-treatment periods:", x$num_periods ,"\n")
  } else {
    cat("Unbalanced Panel: \n")
    range_periods <- paste((x$num_periods)[1], "-", (x$num_periods)[2])
    cat(" + No. pre-treatment periods (range):", range_periods ,"\n")
  }
  cat(" + No. individuals:", x$num_individuals, "\n")
  cat(" + Total no. observations:", x$num_observations,"\n")
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
  title <- "Equivalence Tests for Pre-trends in DiD Estimation"
  
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
  cat("Significance level:", x$significance_level, "\n")
  if(x$equiv_threshold_specified){
    cat("Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold of", x$equiv_threshold, ".\n")
    
    # Display results
    df_print <- data.frame(number_rep(x$max_abs_coefficient), 
                            number_rep(x$bootstrap_critical_value),
                            x$reject_null_hypothesis)
    colnames(df_print) <- c("Max. Abs. Coefficient", "Bootstrap Critical Value", "Reject H0")
    rownames(df_print) <- c("")
  } else {
    cat("Alternative hypothesis: the maximum placebo effect does not exceed the equivalence threshold.\n")
    df_print <- data.frame(number_rep(x$max_abs_coefficient), 
                           number_rep(x$minimum_equiv_threshold))
    colnames(df_print) <- c("Max. Abs. Coefficient", "Min. Equiv. Threshold")
    rownames(df_print) <- c("")
  }
  
  cat("---\n")
  pretty_print(df_print)
  cat("---\n")
  
  # Dataset statistics:
  cat("No. placebo coefficients estimated:", length(x$abs_placebo_coefficients), "\n")
  cat("Base period:", x$base_period ,"\n")
  cat(" \n")
  if(x$is_panel_balanced){
    cat("Balanced Panel:\n")
    cat(" + No. pre-treatment periods:", x$num_periods ,"\n")
  } else {
    cat("Unbalanced Panel:\n")
    range_periods <- paste((x$num_periods)[1], "-", (x$num_periods)[2])
    cat(" + No. pre-treatment periods (range):", range_periods ,"\n")
  }
  cat(" + No. individuals:", x$num_individuals, "\n")
  cat(" + Total no. observations:", x$num_observations,"\n")
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
  title <- "Equivalence Tests for Pre-trends in DiD Estimation"
  
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
  if(x$equiv_threshold_specified){
    cat("Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold of", x$equiv_threshold, ".\n")
    
    df_print <- data.frame(number_rep(x$abs_mean_placebo_coefs),
                           number_rep(sqrt(x$var_mean_placebo_coef)), 
                           number_rep(x$p_value),
                           x$reject_null_hypothesis)
    colnames(df_print) <- c("Abs. Mean Placebo Effect", "Std. Error", "p-value", "Reject H0")
    rownames(df_print) <- c("")
  } else {
    cat("Significance level:", x$significance_level, "\n")
    cat("Alternative hypothesis: the mean placebo effect does not exceed the equivalence threshold.\n")
    df_print <- data.frame(number_rep(x$abs_mean_placebo_coefs), 
                           number_rep(sqrt(x$var_mean_placebo_coef)), 
                           number_rep(x$minimum_equiv_threshold))
    colnames(df_print) <- c("Abs. Mean Placebo Effect", "Std. Error", "Min. Equiv. Threshold")
    rownames(df_print) <- c("")  
  }
  cat("---\n")
  pretty_print(df_print)
  cat("---\n")
  
  # Data statistics
  cat("No. placebo coefficients estimated:", length(x$placebo_coefficients), "\n")
  cat("Base period:", x$base_period ,"\n")
  cat(" \n")
  if(x$is_panel_balanced){
    cat("Balanced Panel: \n")
    cat(" + No. pre-treatment periods:", x$num_periods ,"\n")
  } else {
    cat("Unbalanced Panel: \n")
    range_periods <- paste((x$num_periods)[1], "-", (x$num_periods)[2])
    cat(" + No. pre-treatment periods (range):", range_periods ,"\n")
  }
  cat(" + No. individuals:", x$num_individuals, "\n")
  cat(" + Total no. observations:", x$num_observations,"\n")
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
  title <- "Equivalence Tests for Pre-trends in DiD Estimation"
  
  # Check if console width is less than title length
  if (width < nchar(title)) {
    width <- nchar(title) + 2  # Adjust width to be slightly larger than title
  }
  
  separator <- strrep(" ", floor((width - nchar(title)) / 2))
  
  # Centered title
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  cat(separator, title, "\n", sep = "")
  cat(separator, strrep("=", nchar(title)), "\n", sep = "")
  
  cat("Type: Root Mean Squared Placebo Effect \n")
  cat("Significance level:", x$significance_level, "\n")
  if(x$equiv_threshold_specified){
    cat("Alternative hypothesis: the root mean squared placebo effect does not exceed the equivalence threshold of", x$equiv_threshold, ".\n")
    
    df_print <- data.frame(number_rep(x$rms_placebo_coefficients),
                           number_rep(x$rms_critical_value),
                           x$reject_null_hypothesis)
    colnames(df_print) <- c("RMS Placebo Effect", "Simulated Crit. Val.", "Reject H0")
    rownames(df_print) <- c("")
  } else {
    cat("Alternative hypothesis: the root mean squared placebo effect does not exceed the equivalence threshold.\n")
    df_print <- data.frame(number_rep(x$rms_placebo_coefficients), 
                           number_rep(x$minimum_equiv_threshold))
    colnames(df_print) <- c("RMS Placebo Effect", "Min. Equiv. Threshold")
    rownames(df_print) <- c("")  
  }
  cat("---\n")
  pretty_print(df_print)
  cat("---\n")
  
  # Data statistics
  cat("No. placebo coefficients estimated:", length(x$placebo_coefficients), "\n")
  cat("Base period:", x$base_period ,"\n")
  cat(" \n")
  if(x$is_panel_balanced){
    cat("Balanced Panel: \n")
    cat(" + No. pre-treatment periods:", x$num_periods ,"\n")
  } else {
    cat("Panel: unbalanced\n")
    range_periods <- paste((x$num_periods)[1], "-", (x$num_periods)[2])
    cat(" + No. pre-treatment periods (range):", range_periods ,"\n")
  }
  cat(" + No. individuals:", x$num_individuals, "\n")
  cat(" + Total no. observations:", x$num_observations,"\n")
  cat("\n")
}  


