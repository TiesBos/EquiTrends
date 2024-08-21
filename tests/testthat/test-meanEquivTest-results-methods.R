
test_that("meanEquivTest return for matrix input and standard variance-covariance matrix",{
  sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
  
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, c("X_1", "X_2")]
  cluster_data <- sim_data$cluster
  
  # For the matrix input version of the function:
  # - No equivalence threshold:
  pre_treatment_period <- 1:5
  base_period <- 5
  alpha <- 0.1
  
  
  # Do the procedure by hand for comparison:
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  
  mean_var <- sum((plm_test$vcov)[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 16
  
  
  meanEquivTest_results <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                         period = period_data, X=X_data, 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period,
                                         cluster = cluster_data, 
                                         alpha = alpha)
  
  expect_equal(class(meanEquivTest_results), "meanEquivTest")
  expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(meanEquivTest_results$significance_level, alpha)
  expect_equal(meanEquivTest_results$num_individuals, 500)
  expect_equal(meanEquivTest_results$num_periods, 5)
  expect_equal(meanEquivTest_results$base_period, 5)
  expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.24250704475537498972)
  expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
  p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
  expect_equal(p_val, alpha)
  
  expect_snapshot(print(meanEquivTest_results))
  
  # Equivalence threshold specified:
  meanEquivTest_results2 <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                         period = period_data, X=X_data, 
                                         equiv_threshold = 0.2,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period,
                                         cluster = cluster_data, 
                                         alpha = alpha)
  expect_equal(class(meanEquivTest_results2), "meanEquivTest")
  expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(meanEquivTest_results2$significance_level, alpha)
  expect_equal(meanEquivTest_results2$num_individuals, 500)
  expect_equal(meanEquivTest_results2$num_periods, 5)
  expect_equal(meanEquivTest_results2$base_period, 5)
  expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
  prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
  expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
  expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
  reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
  expect_equal(reject_conclusion, FALSE)
  expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
  
  expect_snapshot(print(meanEquivTest_results2))
  # with index input:
  meanEquivTest_results3 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                         period = 3, X= c(5,6), 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period,
                                         cluster = 7, data = sim_data,
                                         alpha = alpha)
  
  expect_snapshot(print(meanEquivTest_results3))
  # Equivalence threshold specified:
  meanEquivTest_results4 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                          period = 3, X=c(5,6), 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          cluster = 7, data = sim_data,
                                          alpha = alpha)
  expect_snapshot(print(meanEquivTest_results4))
  expect_equal(meanEquivTest_results3, meanEquivTest_results)
  expect_equal(meanEquivTest_results4, meanEquivTest_results2)
  
})



test_that("meanEquivTest return for HC-type variance-covariance matrix",{
  sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
  
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, c("X_1", "X_2")]
  cluster_data <- sim_data$cluster
  
  # For the matrix input version of the function:
  # - No equivalence threshold:
  pre_treatment_period <- 1:5
  base_period <- 5
  alpha <- 0.1
  
  
  # Do the procedure by hand for comparison:
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  var_cov <- plm::vcovHC(plm_test, type="HC1", method = "white1")
  mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 16
  
  
  meanEquivTest_results <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                         period = period_data, X=X_data, 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "HC",
                                         cluster = cluster_data, 
                                         alpha = alpha)
  
  expect_equal(class(meanEquivTest_results), "meanEquivTest")
  expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(meanEquivTest_results$significance_level, alpha)
  expect_equal(meanEquivTest_results$num_individuals, 500)
  expect_equal(meanEquivTest_results$num_periods, 5)
  expect_equal(meanEquivTest_results$base_period, 5)
  expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.22687340815107762126)
  expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
  p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
  expect_equal(p_val, alpha)
  
  # Equivalence threshold specified:
  meanEquivTest_results2 <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                          period = period_data, X=X_data, 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          cluster = cluster_data, 
                                          alpha = alpha, vcov = "HC")
  expect_equal(class(meanEquivTest_results2), "meanEquivTest")
  expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(meanEquivTest_results2$significance_level, alpha)
  expect_equal(meanEquivTest_results2$num_individuals, 500)
  expect_equal(meanEquivTest_results2$num_periods, 5)
  expect_equal(meanEquivTest_results2$base_period, 5)
  expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
  prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
  expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
  expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
  reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
  expect_equal(reject_conclusion, FALSE)
  expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
  
  # with index input:
  meanEquivTest_results3 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                         period = 3, X= c(5,6), 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "HC",
                                         cluster = 7, data = sim_data,
                                         alpha = alpha)
  meanEquivTest_results4 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                          period = 3, X=c(5,6), 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          cluster = 7, data = sim_data,vcov = "HC",
                                          alpha = alpha)
  expect_equal(meanEquivTest_results3, meanEquivTest_results)
  expect_equal(meanEquivTest_results4, meanEquivTest_results2)
  
})


# test_that("meanEquivTest return for index input and HC type variance-covariance matrix",{
#   sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
#   
#   # The data in vector/matrix form:
#   Y_data <- sim_data$Y
#   ID_data <- sim_data$ID
#   G_data <- sim_data$G
#   period_data <- sim_data$period
#   X_data <- sim_data[, c("X_1", "X_2")]
#   cluster_data <- sim_data$cluster
#   
#   # For the matrix input version of the function:
#   # - No equivalence threshold:
#   pre_treatment_period <- 1:5
#   base_period <- 5
#   alpha <- 0.1
#   
#   
#   # Do the procedure by hand for comparison:
#   subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
#   test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
#   plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
#   placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
#   var_cov <- plm::vcovHC(plm_test, type="HC1", method = "white1")
#   
#   mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 500
#   
#   
#   meanEquivTest_results <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                          period = 3, X= c(5,6), 
#                                          equiv_threshold = NULL,
#                                          pretreatment_period = pre_treatment_period, 
#                                          base_period = base_period, vcov = "HC",
#                                          cluster = 7, data = sim_data,
#                                          alpha = alpha)
#   
#   expect_equal(class(meanEquivTest_results), "meanEquivTest")
#   expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
#   expect_equal(meanEquivTest_results$significance_level, alpha)
#   expect_equal(meanEquivTest_results$num_individuals, 500)
#   expect_equal(meanEquivTest_results$num_periods, 5)
#   expect_equal(meanEquivTest_results$base_period, 5)
#   expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.13751285)
#   expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
#   p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
#   expect_equal(p_val, alpha)
#   
#   # Equivalence threshold specified:
#   meanEquivTest_results2 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                           period = 3, X=c(5,6), 
#                                           equiv_threshold = 0.2,
#                                           pretreatment_period = pre_treatment_period, 
#                                           base_period = base_period,
#                                           cluster = 7, data = sim_data,vcov = "HC",
#                                           alpha = alpha)
#   expect_equal(class(meanEquivTest_results2), "meanEquivTest")
#   expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
#   expect_equal(meanEquivTest_results2$significance_level, alpha)
#   expect_equal(meanEquivTest_results2$num_individuals, 500)
#   expect_equal(meanEquivTest_results2$num_periods, 5)
#   expect_equal(meanEquivTest_results2$base_period, 5)
#   expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
#   prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
#   expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
#   expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
#   reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
#   expect_equal(reject_conclusion, TRUE)
#   expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
#   
# })

test_that("meanEquivTest return for CL-type variance-covariance matrix",{
  sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
  
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, c("X_1", "X_2")]
  cluster_data <- sim_data$cluster
  
  # For the matrix input version of the function:
  # - No equivalence threshold:
  pre_treatment_period <- 1:5
  base_period <- 5
  alpha <- 0.1
  
  
  # Do the procedure by hand for comparison:
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  var_cov <- clubSandwich::vcovCR(plm_test, cluster=cluster_data, type="CR0")
  mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 16
  
  
  meanEquivTest_results <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                         period = period_data, X=X_data, 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "CL",
                                         cluster = cluster_data, 
                                         alpha = alpha)
  
  expect_equal(class(meanEquivTest_results), "meanEquivTest")
  expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(meanEquivTest_results$significance_level, alpha)
  expect_equal(meanEquivTest_results$num_individuals, 500)
  expect_equal(meanEquivTest_results$num_periods, 5)
  expect_equal(meanEquivTest_results$base_period, 5)
  expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.16322737373597387411, tolerance = 1e-7)
  expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
  p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
  expect_equal(p_val, alpha)
  
  # Equivalence threshold specified:
  meanEquivTest_results2 <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                          period = period_data, X=X_data, 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          cluster = cluster_data, 
                                          alpha = alpha, vcov = "CL")
  expect_equal(class(meanEquivTest_results2), "meanEquivTest")
  expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(meanEquivTest_results2$significance_level, alpha)
  expect_equal(meanEquivTest_results2$num_individuals, 500)
  expect_equal(meanEquivTest_results2$num_periods, 5)
  expect_equal(meanEquivTest_results2$base_period, 5)
  expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
  prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
  expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
  expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
  reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
  expect_equal(reject_conclusion, TRUE)
  expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
  
  # For index input:
  meanEquivTest_results3 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                         period = 3, X= c(5,6), 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "CL",
                                         cluster = 7, data = sim_data,
                                         alpha = alpha)
  
  meanEquivTest_results4 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                          period = 3, X=c(5,6), 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          cluster = 7, data = sim_data,vcov = "CL",
                                          alpha = alpha)
  expect_equal(meanEquivTest_results3, meanEquivTest_results)
  expect_equal(meanEquivTest_results4, meanEquivTest_results2)
  
})


# test_that("meanEquivTest return for index input and CL type variance-covariance matrix",{
#   sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
#   
#   # The data in vector/matrix form:
#   Y_data <- sim_data$Y
#   ID_data <- sim_data$ID
#   G_data <- sim_data$G
#   period_data <- sim_data$period
#   X_data <- sim_data[, c("X_1", "X_2")]
#   cluster_data <- sim_data$cluster
#   
#   # For the matrix input version of the function:
#   # - No equivalence threshold:
#   pre_treatment_period <- 1:5
#   base_period <- 5
#   alpha <- 0.1
#   
#   
#   # Do the procedure by hand for comparison:
#   subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
#   test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
#   plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
#   placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
#   var_cov <- clubSandwich::vcovCR(plm_test, cluster=cluster_data, type="CR0")
#   
#   mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 500
#   
#   
#   meanEquivTest_results <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                          period = 3, X= c(5,6), 
#                                          equiv_threshold = NULL,
#                                          pretreatment_period = pre_treatment_period, 
#                                          base_period = base_period, vcov = "CL",
#                                          cluster = 7, data = sim_data,
#                                          alpha = alpha)
#   
#   expect_equal(class(meanEquivTest_results), "meanEquivTest")
#   expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
#   expect_equal(meanEquivTest_results$significance_level, alpha)
#   expect_equal(meanEquivTest_results$num_individuals, 500)
#   expect_equal(meanEquivTest_results$num_periods, 5)
#   expect_equal(meanEquivTest_results$base_period, 5)
#   expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.1261254, tolerance = 1e-7)
#   expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
#   p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
#   expect_equal(p_val, alpha)
#   
#   # Equivalence threshold specified:
#   meanEquivTest_results2 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                           period = 3, X=c(5,6), 
#                                           equiv_threshold = 0.2,
#                                           pretreatment_period = pre_treatment_period, 
#                                           base_period = base_period,
#                                           cluster = 7, data = sim_data,vcov = "CL",
#                                           alpha = alpha)
#   expect_equal(class(meanEquivTest_results2), "meanEquivTest")
#   expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
#   expect_equal(meanEquivTest_results2$significance_level, alpha)
#   expect_equal(meanEquivTest_results2$num_individuals, 500)
#   expect_equal(meanEquivTest_results2$num_periods, 5)
#   expect_equal(meanEquivTest_results2$base_period, 5)
#   expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
#   prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
#   expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
#   expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
#   reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
#   expect_equal(reject_conclusion, TRUE)
#   expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
#   
# })

test_that("meanEquivTest return for HAC type variance-covariance matrix",{
  sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
  
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, c("X_1", "X_2")]
  cluster_data <- sim_data$cluster
  
  # For the matrix input version of the function:
  # - No equivalence threshold:
  pre_treatment_period <- 1:5
  base_period <- 5
  alpha <- 0.1
  
  
  # Do the procedure by hand for comparison:
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  var_cov <- plm::vcovHC(plm_test, type="HC3", method = "arellano")
  
  mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 16
  
  
  meanEquivTest_results <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                         period = period_data, X= X_data, 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "HAC",
                                         cluster = cluster_data, 
                                         alpha = alpha)
  
  expect_equal(class(meanEquivTest_results), "meanEquivTest")
  expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(meanEquivTest_results$significance_level, alpha)
  expect_equal(meanEquivTest_results$num_individuals, 500)
  expect_equal(meanEquivTest_results$num_periods, 5)
  expect_equal(meanEquivTest_results$base_period, 5)
  expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.23890025676027837331, tolerance = 1e-7)
  expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
  p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
  expect_equal(p_val, alpha)
  
  # Equivalence threshold specified:
  meanEquivTest_results2 <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                          period = period_data, X= X_data, 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          cluster = cluster_data, vcov = "HAC",
                                          alpha = alpha)
  expect_equal(class(meanEquivTest_results2), "meanEquivTest")
  expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(meanEquivTest_results2$significance_level, alpha)
  expect_equal(meanEquivTest_results2$num_individuals, 500)
  expect_equal(meanEquivTest_results2$num_periods, 5)
  expect_equal(meanEquivTest_results2$base_period, 5)
  expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
  prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
  expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
  expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
  reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
  expect_equal(reject_conclusion, FALSE)
  expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
  
  meanEquivTest_results3 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                         period = 3, X= c(5,6), 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "HAC",
                                         cluster = 7, data = sim_data,
                                         alpha = alpha)
  
  meanEquivTest_results4 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                          period = 3, X=c(5,6), 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          cluster = 7, data = sim_data,vcov = "HAC",
                                          alpha = alpha)
  expect_equal(meanEquivTest_results3, meanEquivTest_results)
  expect_equal(meanEquivTest_results4, meanEquivTest_results2)
})

# test_that("meanEquivTest return for index input and HAC type variance-covariance matrix",{
#   sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
#   
#   # The data in vector/matrix form:
#   Y_data <- sim_data$Y
#   ID_data <- sim_data$ID
#   G_data <- sim_data$G
#   period_data <- sim_data$period
#   X_data <- sim_data[, c("X_1", "X_2")]
#   cluster_data <- sim_data$cluster
#   
#   # For the matrix input version of the function:
#   # - No equivalence threshold:
#   pre_treatment_period <- 1:5
#   base_period <- 5
#   alpha <- 0.1
#   
#   
#   # Do the procedure by hand for comparison:
#   subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
#   test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
#   plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
#   placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
#   var_cov <- plm::vcovHC(plm_test, type="HC3", method = "arellano")
#   
#   mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 500
#   
#   
#   meanEquivTest_results <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                          period = 3, X= c(5,6), 
#                                          equiv_threshold = NULL,
#                                          pretreatment_period = pre_treatment_period, 
#                                          base_period = base_period, vcov = "HAC",
#                                          cluster = 7, data = sim_data,
#                                          alpha = alpha)
#   
#   expect_equal(class(meanEquivTest_results), "meanEquivTest")
#   expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
#   expect_equal(meanEquivTest_results$significance_level, alpha)
#   expect_equal(meanEquivTest_results$num_individuals, 500)
#   expect_equal(meanEquivTest_results$num_periods, 5)
#   expect_equal(meanEquivTest_results$base_period, 5)
#   expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.13966959, tolerance = 1e-7)
#   expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
#   p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
#   expect_equal(p_val, alpha)
#   
#   # Equivalence threshold specified:
#   meanEquivTest_results2 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                           period = 3, X=c(5,6), 
#                                           equiv_threshold = 0.2,
#                                           pretreatment_period = pre_treatment_period, 
#                                           base_period = base_period,
#                                           cluster = 7, data = sim_data,vcov = "HAC",
#                                           alpha = alpha)
#   expect_equal(class(meanEquivTest_results2), "meanEquivTest")
#   expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
#   expect_equal(meanEquivTest_results2$significance_level, alpha)
#   expect_equal(meanEquivTest_results2$num_individuals, 500)
#   expect_equal(meanEquivTest_results2$num_periods, 5)
#   expect_equal(meanEquivTest_results2$base_period, 5)
#   expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
#   prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
#   expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
#   expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
#   reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
#   expect_equal(reject_conclusion, TRUE)
#   expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
#   
# })



test_that("meanEquivTest return for CL-type variance-covariance matrix (no cluster index)",{
  sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
  
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, c("X_1", "X_2")]
  cluster_data <- sim_data$cluster
  
  # For the matrix input version of the function:
  # - No equivalence threshold:
  pre_treatment_period <- 1:5
  base_period <- 5
  alpha <- 0.1
  
  
  # Do the procedure by hand for comparison:
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  var_cov <- clubSandwich::vcovCR(plm_test, cluster = "ID", type="CR0")
  mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 16
  
  
  meanEquivTest_results <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                         period = period_data, X=X_data, 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "CL",
                                         alpha = alpha)
  
  expect_equal(class(meanEquivTest_results), "meanEquivTest")
  expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(meanEquivTest_results$significance_level, alpha)
  expect_equal(meanEquivTest_results$num_individuals, 500)
  expect_equal(meanEquivTest_results$num_periods, 5)
  expect_equal(meanEquivTest_results$base_period, 5)
  expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.23861523182253563391, tolerance = 2e-7)
  expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
  p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
  expect_equal(p_val, alpha)
  
  # Equivalence threshold specified:
  meanEquivTest_results2 <- meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                          period = period_data, X=X_data, 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          alpha = alpha, vcov = "CL")
  expect_equal(class(meanEquivTest_results2), "meanEquivTest")
  expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(meanEquivTest_results2$significance_level, alpha)
  expect_equal(meanEquivTest_results2$num_individuals, 500)
  expect_equal(meanEquivTest_results2$num_periods, 5)
  expect_equal(meanEquivTest_results2$base_period, 5)
  expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
  expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
  expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
  expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
  expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
  expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
  prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
  expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
  expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
  reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
  expect_equal(reject_conclusion, FALSE)
  expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
  
  meanEquivTest_results3 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                         period = 3, X= c(5,6), 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period, vcov = "CL",
                                         data = sim_data,
                                         alpha = alpha)
  
  meanEquivTest_results4 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
                                          period = 3, X=c(5,6), 
                                          equiv_threshold = 0.2,
                                          pretreatment_period = pre_treatment_period, 
                                          base_period = base_period,
                                          data = sim_data,vcov = "CL",
                                          alpha = alpha)
  expect_equal(meanEquivTest_results3, meanEquivTest_results)
  expect_equal(meanEquivTest_results4, meanEquivTest_results2)
  
})


# test_that("meanEquivTest return for index input and CL type variance-covariance matrix (no cluster index)",{
#   sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
#   
#   # The data in vector/matrix form:
#   Y_data <- sim_data$Y
#   ID_data <- sim_data$ID
#   G_data <- sim_data$G
#   period_data <- sim_data$period
#   X_data <- sim_data[, c("X_1", "X_2")]
#   cluster_data <- sim_data$cluster
#   
#   # For the matrix input version of the function:
#   # - No equivalence threshold:
#   pre_treatment_period <- 1:5
#   base_period <- 5
#   alpha <- 0.1
#   
#   
#   # Do the procedure by hand for comparison:
#   subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
#   test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
#   plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
#   placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
#   var_cov <- clubSandwich::vcovCR(plm_test, cluster = "ID", type="CR0")
#   
#   mean_var <- sum(var_cov[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]) / 500
#   
#   
#   meanEquivTest_results <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                          period = 3, X= c(5,6), 
#                                          equiv_threshold = NULL,
#                                          pretreatment_period = pre_treatment_period, 
#                                          base_period = base_period, vcov = "CL",
#                                          data = sim_data,
#                                          alpha = alpha)
#   
#   expect_equal(class(meanEquivTest_results), "meanEquivTest")
#   expect_equal(meanEquivTest_results$equiv_threshold_specified, FALSE)
#   expect_equal(meanEquivTest_results$significance_level, alpha)
#   expect_equal(meanEquivTest_results$num_individuals, 500)
#   expect_equal(meanEquivTest_results$num_periods, 5)
#   expect_equal(meanEquivTest_results$base_period, 5)
#   expect_equal(meanEquivTest_results$minimum_equiv_threshold, 0.1396184, tolerance = 2e-7)
#   expect_equal(length(meanEquivTest_results$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results$placebo_coefficients, placebo_coefs)
#   p_val <- VGAM::pfoldnorm(meanEquivTest_results$abs_mean_placebo_coefs, mean=meanEquivTest_results$minimum_equiv_threshold, sd=sqrt(meanEquivTest_results$var_mean_placebo_coef))
#   expect_equal(p_val, alpha)
#   
#   # Equivalence threshold specified:
#   meanEquivTest_results2 <- meanEquivTest(Y = 1, ID= 2, G = 4, 
#                                           period = 3, X=c(5,6), 
#                                           equiv_threshold = 0.2,
#                                           pretreatment_period = pre_treatment_period, 
#                                           base_period = base_period,
#                                           data = sim_data,vcov = "CL",
#                                           alpha = alpha)
#   expect_equal(class(meanEquivTest_results2), "meanEquivTest")
#   expect_equal(meanEquivTest_results2$equiv_threshold_specified, TRUE)
#   expect_equal(meanEquivTest_results2$significance_level, alpha)
#   expect_equal(meanEquivTest_results2$num_individuals, 500)
#   expect_equal(meanEquivTest_results2$num_periods, 5)
#   expect_equal(meanEquivTest_results2$base_period, 5)
#   expect_equal(length(meanEquivTest_results2$placebo_coefficients), 4)
#   expect_equal(length(meanEquivTest_results2$abs_mean_placebo_coefs), 1)
#   expect_equal(meanEquivTest_results2$abs_mean_placebo_coefs, abs(mean(placebo_coefs)))
#   expect_equal(length(meanEquivTest_results2$var_mean_placebo_coef), 1)
#   expect_equal(meanEquivTest_results2$var_mean_placebo_coef, mean_var)
#   expect_equal(meanEquivTest_results2$placebo_coefficients, placebo_coefs)
#   prob_crit_val <- VGAM::pfoldnorm(meanEquivTest_results2$mean_critical_value, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef))
#   expect_equal(prob_crit_val, 0.1, tolerance = 1e-7)
#   expect_equal(VGAM::qfoldnorm(alpha, mean=meanEquivTest_results2$equiv_threshold, sd=sqrt(meanEquivTest_results2$var_mean_placebo_coef)), meanEquivTest_results2$mean_critical_value)
#   reject_conclusion <- meanEquivTest_results2$reject_null_hypothesis
#   expect_equal(reject_conclusion, TRUE)
#   expect_equal(meanEquivTest_results2$equiv_threshold, 0.2)
#   
# })



