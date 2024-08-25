# Unit Tests for the output of the bootstrap version of the maxEquivTest function
test_that("Output of maxEquivTest with type = Boot",{
  skip_on_cran()
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
  
  B <- 100
  
  # For the version with data input as matrices and vectors:
  # If the equivalence threshold is specified:
  maxEquivTest_results <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = 1,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       alpha = alpha, B=B, type = "Boot")
  
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  
  expect_equal(class(maxEquivTest_results), "maxEquivTestBoot")
  expect_equal(maxEquivTest_results$equiv_threshold_specified, TRUE)
  expect_equal(maxEquivTest_results$significance_level, alpha)
  expect_equal(maxEquivTest_results$num_individuals, 500)
  expect_equal(maxEquivTest_results$num_periods, 5)
  expect_equal(maxEquivTest_results$base_period, 5)
  expect_equal(length(maxEquivTest_results$placebo_coefficients), 4)
  expect_equal(max(abs(placebo_coefs)), maxEquivTest_results$max_abs_coefficient)
  expect_equal(maxEquivTest_results$placebo_coefficients, placebo_coefs, tolerance = 1e-3)
  expect_equal(maxEquivTest_results$B, B)
  expect_equal(length(B), 1)
  expect_equal(length(maxEquivTest_results$bootstrap_critical_value), 1)
  reject_null <- max(abs(placebo_coefs)) < maxEquivTest_results$bootstrap_critical_value
  expect_equal(maxEquivTest_results$reject_null_hypothesis, reject_null)
  
  maxEquivTest_results2 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                       period = 3, X=c(5,6), 
                                       equiv_threshold = 1,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period, data = sim_data,
                                       alpha = alpha, B=B, type = "Boot")
  
  expect_equal(class(maxEquivTest_results2), "maxEquivTestBoot")
  expect_equal(maxEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(maxEquivTest_results2$significance_level, alpha)
  expect_equal(maxEquivTest_results2$num_individuals, 500)
  expect_equal(maxEquivTest_results2$num_periods, 5)
  expect_equal(maxEquivTest_results2$base_period, 5)
  expect_equal(length(maxEquivTest_results2$placebo_coefficients), 4)
  expect_equal(max(abs(placebo_coefs)), maxEquivTest_results2$max_abs_coefficient)
  expect_equal(maxEquivTest_results2$placebo_coefficients, placebo_coefs, tolerance = 1e-3)
  expect_equal(maxEquivTest_results2$B, B)
  expect_equal(length(B), 1)
  expect_equal(length(maxEquivTest_results2$bootstrap_critical_value), 1)
  reject_null <- max(abs(placebo_coefs)) < maxEquivTest_results2$bootstrap_critical_value
  expect_equal(maxEquivTest_results2$reject_null_hypothesis, reject_null)
  
  test_formula <- as.formula(Y ~ placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  maxEquivTest_results3 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, 
                                        equiv_threshold = NULL,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period, data = sim_data,
                                        alpha = alpha, B=B, type = "Boot")
  expect_equal(class(maxEquivTest_results3), "maxEquivTestBoot")
  expect_equal(maxEquivTest_results3$equiv_threshold_specified, FALSE)
  expect_equal(maxEquivTest_results3$significance_level, alpha)
  expect_equal(maxEquivTest_results3$num_individuals, 500)
  expect_equal(maxEquivTest_results3$num_periods, 5)
  expect_equal(maxEquivTest_results3$base_period, 5)
  expect_equal(length(maxEquivTest_results3$placebo_coefficients), 4)
  expect_equal(max(abs(placebo_coefs)), maxEquivTest_results3$max_abs_coefficient)
  expect_equal(maxEquivTest_results3$placebo_coefficients, placebo_coefs, tolerance = 1e-3)
  expect_equal(maxEquivTest_results3$B, B)
  expect_equal(length(B), 1)
  expect_equal(length(maxEquivTest_results3$minimum_equiv_threshold), 1)
  expect_equal((maxEquivTest_results3$minimum_equiv_threshold >= maxEquivTest_results3$max_abs_coefficient), TRUE )
  
  # Check the print functions but set the minimum equiv threshold / critical value to 0 due to randomization:
  print_maxEquivTest <- maxEquivTest_results
  print_maxEquivTest$bootstrap_critical_value <- 0
  expect_snapshot(print(print_maxEquivTest))
  
  print_maxEquivTest2 <- maxEquivTest_results2
  print_maxEquivTest2$bootstrap_critical_value <- 0
  expect_snapshot(print(print_maxEquivTest2))
  
  print_maxEquivTest3 <- maxEquivTest_results3
  print_maxEquivTest3$minimum_equiv_threshold <- 0
  expect_snapshot(print(print_maxEquivTest3))
})

