test_that("results rmsEquivTest", {
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
  
  # If the equivalence threshold is not specified:
  rmsEquivTest_results <- rmsEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                         period = period_data, X=X_data, 
                                         equiv_threshold = NULL,
                                         pretreatment_period = pre_treatment_period, 
                                         base_period = base_period,
                                         alpha = alpha)
  
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  
  expect_equal(class(rmsEquivTest_results), "rmsEquivTest")
  expect_equal(rmsEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(rmsEquivTest_results$significance_level, alpha)
  expect_equal(rmsEquivTest_results$num_individuals, 500)
  expect_equal(rmsEquivTest_results$num_periods, 5)
  expect_equal(rmsEquivTest_results$base_period, 5)
  expect_equal(length(rmsEquivTest_results$placebo_coefficients), 4)
  expect_equal(length(rmsEquivTest_results$minimum_equiv_threshold), 1)
  expect_equal(rmsEquivTest_results$placebo_coefficients, placebo_coefs)
  expect_equal(rmsEquivTest_results$rms_placebo_coefficients, sqrt(mean((placebo_coefs)^2)))
  
  
  # The same but with index input:
  rmsEquivTest_results2 <- rmsEquivTest(Y = 1, ID= 2, G = 4, 
                                       period = 3, X=c(5,6), 
                                       equiv_threshold = NULL,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period, data = sim_data,
                                       alpha = alpha)
  expect_equal(class(rmsEquivTest_results2), "rmsEquivTest")
  expect_equal(rmsEquivTest_results2$equiv_threshold_specified, FALSE)
  expect_equal(rmsEquivTest_results2$significance_level, alpha)
  expect_equal(rmsEquivTest_results2$num_individuals, 500)
  expect_equal(rmsEquivTest_results2$num_periods, 5)
  expect_equal(rmsEquivTest_results2$base_period, 5)
  expect_equal(length(rmsEquivTest_results2$placebo_coefficients), 4)
  expect_equal(rmsEquivTest_results2$placebo_coefficients, placebo_coefs)
  expect_equal(rmsEquivTest_results2$rms_placebo_coefficients, sqrt(mean((placebo_coefs)^2)))
  
  
  # If the equivalence threshold is specified:
  rmsEquivTest_results3 <- rmsEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = 1,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       alpha = alpha)
  
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  
  expect_equal(class(rmsEquivTest_results3), "rmsEquivTest")
  expect_equal(rmsEquivTest_results3$equiv_threshold_specified, TRUE)
  expect_equal(rmsEquivTest_results3$significance_level, alpha)
  expect_equal(rmsEquivTest_results3$num_individuals, 500)
  expect_equal(rmsEquivTest_results3$num_periods, 5)
  expect_equal(rmsEquivTest_results3$base_period, 5)
  expect_equal(length(rmsEquivTest_results3$placebo_coefficients), 4)
  expect_equal(rmsEquivTest_results3$placebo_coefficients, placebo_coefs)
  expect_equal(rmsEquivTest_results3$rms_placebo_coefficients, sqrt(mean((placebo_coefs)^2)))
  expect_equal(rmsEquivTest_results3$equiv_threshold, 1)
  expect_equal(length(rmsEquivTest_results3$rms_critical_value), 1)
  
  
  # The same but with index input:
  rmsEquivTest_results4 <- rmsEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X=c(5,6), 
                                        equiv_threshold = 1,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period, data = sim_data,
                                        alpha = alpha)
  expect_equal(class(rmsEquivTest_results4), "rmsEquivTest")
  expect_equal(rmsEquivTest_results4$equiv_threshold_specified, TRUE)
  expect_equal(rmsEquivTest_results4$significance_level, alpha)
  expect_equal(rmsEquivTest_results4$num_individuals, 500)
  expect_equal(rmsEquivTest_results4$num_periods, 5)
  expect_equal(rmsEquivTest_results4$base_period, 5)
  expect_equal(length(rmsEquivTest_results4$placebo_coefficients), 4)
  expect_equal(rmsEquivTest_results4$placebo_coefficients, placebo_coefs)
  expect_equal(rmsEquivTest_results4$rms_placebo_coefficients, sqrt(mean((placebo_coefs)^2)))
  expect_equal(rmsEquivTest_results4$equiv_threshold, 1)
  expect_equal(length(rmsEquivTest_results4$rms_critical_value), 1)
  
  
  
  # Check the print methods:
  rmsEquivTest_results$minimum_equiv_threshold <- 0
  expect_snapshot(print(rmsEquivTest_results))
  
  rmsEquivTest_results2$minimum_equiv_threshold <- 0
  expect_snapshot(print(rmsEquivTest_results2))
  
  rmsEquivTest_results3$rms_critical_value <- 0
  expect_snapshot(print(rmsEquivTest_results3))
  
  rmsEquivTest_results4$rms_critical_value <- 0
  expect_snapshot(print(rmsEquivTest_results4))
})








