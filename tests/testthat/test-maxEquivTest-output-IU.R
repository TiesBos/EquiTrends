
test_that("maxEquivTest return for matrix input and standard variance-covariance matrix",{
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
  
  vcov_mat <- plm_test$vcov
  subcov_mat <- vcov_mat[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  beta_var <- diag(subcov_mat)
  
  # Calculating the standard errors
  beta_se <- sqrt(beta_var)
  
  maxEquivTest_results <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = NULL,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       cluster = cluster_data, 
                                       alpha = alpha,
                                       type = "IU")
  
  expect_equal(class(maxEquivTest_results), "maxEquivTestIU")
  expect_equal(maxEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(maxEquivTest_results$significance_level, alpha)
  expect_equal(maxEquivTest_results$num_individuals, 500)
  expect_equal(maxEquivTest_results$num_periods, 5)
  expect_equal(maxEquivTest_results$base_period, 5)
  expect_equal(maxEquivTest_results$minimum_equiv_threshold, 0.38987476186682140655, tolerance = 1e-6)
  expect_equal(maxEquivTest_results$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results$placebo_coefficients, placebo_coefs)
  
  
  # Equivalence threshold specified:
  maxEquivTest_results2 <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                        period = period_data, X=X_data, 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = cluster_data, 
                                        alpha = alpha,
                                        type = "IU")
  
  expect_equal(class(maxEquivTest_results2), "maxEquivTestIU")
  expect_equal(maxEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(maxEquivTest_results2$significance_level, alpha)
  expect_equal(maxEquivTest_results2$num_individuals, 500)
  expect_equal(maxEquivTest_results2$num_periods, 5)
  expect_equal(maxEquivTest_results2$base_period, 5)
  expect_equal(maxEquivTest_results2$equiv_threshold, 0.2)
  expect_equal(maxEquivTest_results2$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results2$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results2$placebo_coefficients, placebo_coefs)
  expect_equal(maxEquivTest_results2$reject_null_hypothesis, FALSE)
  
  
  
  # with index input:
  maxEquivTest_results3 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X= c(5,6), 
                                        equiv_threshold = NULL,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU")
  
  # Equivalence threshold specified:
  maxEquivTest_results4 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X=c(5,6), 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU")
  expect_equal(maxEquivTest_results3, maxEquivTest_results)
  expect_equal(maxEquivTest_results4, maxEquivTest_results2)
})

test_that("maxEquivTest return for HC-type variance-covariance matrix",{
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
  vcov_mat <- plm::vcovHC(plm_test, type="HC1", method = "white1")
 
  subcov_mat <- vcov_mat[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  beta_var <- diag(subcov_mat)
  
  # Calculating the standard errors
  beta_se <- sqrt(beta_var)
  
  maxEquivTest_results <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = NULL,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       cluster = cluster_data, 
                                       alpha = alpha,
                                       type = "IU", vcov="HC")
  
  expect_equal(class(maxEquivTest_results), "maxEquivTestIU")
  expect_equal(maxEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(maxEquivTest_results$significance_level, alpha)
  expect_equal(maxEquivTest_results$num_individuals, 500)
  expect_equal(maxEquivTest_results$num_periods, 5)
  expect_equal(maxEquivTest_results$base_period, 5)
  expect_equal(maxEquivTest_results$minimum_equiv_threshold, 0.37412870289644517552, tolerance = 1e-6)
  expect_equal(maxEquivTest_results$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results$placebo_coefficients, placebo_coefs)
  
  # Equivalence threshold specified:
  maxEquivTest_results2 <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                        period = period_data, X=X_data, 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = cluster_data, 
                                        alpha = alpha,
                                        type = "IU", vcov="HC")
  
  expect_equal(class(maxEquivTest_results2), "maxEquivTestIU")
  expect_equal(maxEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(maxEquivTest_results2$significance_level, alpha)
  expect_equal(maxEquivTest_results2$num_individuals, 500)
  expect_equal(maxEquivTest_results2$num_periods, 5)
  expect_equal(maxEquivTest_results2$base_period, 5)
  expect_equal(maxEquivTest_results2$equiv_threshold, 0.2)
  expect_equal(maxEquivTest_results2$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results2$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results2$placebo_coefficients, placebo_coefs)
  expect_equal(maxEquivTest_results2$reject_null_hypothesis, FALSE)
  
  
  # with index input:
  maxEquivTest_results3 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X= c(5,6), 
                                        equiv_threshold = NULL,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov="HC")
  # Equivalence threshold specified:
  maxEquivTest_results4 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X=c(5,6), 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov = "HC")
  expect_equal(maxEquivTest_results3, maxEquivTest_results)
  expect_equal(maxEquivTest_results4, maxEquivTest_results2)
  
})


test_that("maxEquivTest return for HAC-type variance-covariance matrix",{
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
  vcov_mat <- plm::vcovHC(plm_test, type="HC3", method = "arellano")
  
  subcov_mat <- vcov_mat[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  beta_var <- diag(subcov_mat)
  
  # Calculating the standard errors
  beta_se <- sqrt(beta_var)
  
  maxEquivTest_results <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = NULL,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       cluster = cluster_data, 
                                       alpha = alpha,
                                       type = "IU", vcov="HAC")
  
  expect_equal(class(maxEquivTest_results), "maxEquivTestIU")
  expect_equal(maxEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(maxEquivTest_results$significance_level, alpha)
  expect_equal(maxEquivTest_results$num_individuals, 500)
  expect_equal(maxEquivTest_results$num_periods, 5)
  expect_equal(maxEquivTest_results$base_period, 5)
  expect_equal(maxEquivTest_results$minimum_equiv_threshold, 0.39576374871198399807, tolerance = 1e-6)
  expect_equal(maxEquivTest_results$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results$placebo_coefficients, placebo_coefs)
  
  # Equivalence threshold specified:
  maxEquivTest_results2 <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                        period = period_data, X=X_data, 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = cluster_data, 
                                        alpha = alpha,
                                        type = "IU", vcov="HAC")
  
  expect_equal(class(maxEquivTest_results2), "maxEquivTestIU")
  expect_equal(maxEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(maxEquivTest_results2$significance_level, alpha)
  expect_equal(maxEquivTest_results2$num_individuals, 500)
  expect_equal(maxEquivTest_results2$num_periods, 5)
  expect_equal(maxEquivTest_results2$base_period, 5)
  expect_equal(maxEquivTest_results2$equiv_threshold, 0.2)
  expect_equal(maxEquivTest_results2$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results2$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results2$placebo_coefficients, placebo_coefs)
  expect_equal(maxEquivTest_results2$reject_null_hypothesis, FALSE)
  
  
  # with index input:
  maxEquivTest_results3 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X= c(5,6), 
                                        equiv_threshold = NULL,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov="HAC")
  # Equivalence threshold specified:
  maxEquivTest_results4 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X=c(5,6), 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov = "HAC")
  expect_equal(maxEquivTest_results3, maxEquivTest_results)
  expect_equal(maxEquivTest_results4, maxEquivTest_results2)
  
})


test_that("maxEquivTest return for CL-type variance-covariance matrix without cluster vector",{
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
  vcov_mat <- clubSandwich::vcovCR(plm_test, cluster="ID", type="CR0")
  
  subcov_mat <- vcov_mat[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  beta_var <- diag(subcov_mat)
  
  # Calculating the standard errors
  beta_se <- sqrt(beta_var)
  
  maxEquivTest_results <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = NULL,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       alpha = alpha,
                                       type = "IU", vcov="CL")
  
  expect_equal(class(maxEquivTest_results), "maxEquivTestIU")
  expect_equal(maxEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(maxEquivTest_results$significance_level, alpha)
  expect_equal(maxEquivTest_results$num_individuals, 500)
  expect_equal(maxEquivTest_results$num_periods, 5)
  expect_equal(maxEquivTest_results$base_period, 5)
  expect_equal(maxEquivTest_results$minimum_equiv_threshold, 0.39537894825602459825, tolerance = 1e-6)
  expect_equal(maxEquivTest_results$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results$placebo_coefficients, placebo_coefs)
  
  # Equivalence threshold specified:
  maxEquivTest_results2 <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                        period = period_data, X=X_data, 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        alpha = alpha,
                                        type = "IU", vcov="CL")
  
  expect_equal(class(maxEquivTest_results2), "maxEquivTestIU")
  expect_equal(maxEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(maxEquivTest_results2$significance_level, alpha)
  expect_equal(maxEquivTest_results2$num_individuals, 500)
  expect_equal(maxEquivTest_results2$num_periods, 5)
  expect_equal(maxEquivTest_results2$base_period, 5)
  expect_equal(maxEquivTest_results2$equiv_threshold, 0.2)
  expect_equal(maxEquivTest_results2$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results2$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results2$placebo_coefficients, placebo_coefs)
  expect_equal(maxEquivTest_results2$reject_null_hypothesis, FALSE)
  
  
  # with index input:
  maxEquivTest_results3 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X= c(5,6), 
                                        equiv_threshold = NULL,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov="CL")
  # Equivalence threshold specified:
  maxEquivTest_results4 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X=c(5,6), 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov = "CL")
  expect_equal(maxEquivTest_results3, maxEquivTest_results)
  expect_equal(maxEquivTest_results4, maxEquivTest_results2)
  
})

test_that("maxEquivTest return for CL-type variance-covariance matrix with cluster vector",{
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
  vcov_mat <- clubSandwich::vcovCR(plm_test, cluster= cluster_data, type="CR0")
  
  subcov_mat <- vcov_mat[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), drop = FALSE]
  beta_var <- diag(subcov_mat)
  
  # Calculating the standard errors
  beta_se <- sqrt(beta_var)
  
  maxEquivTest_results <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = NULL,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       cluster = cluster_data, 
                                       alpha = alpha,
                                       type = "IU", vcov="CL")
  
  expect_equal(class(maxEquivTest_results), "maxEquivTestIU")
  expect_equal(maxEquivTest_results$equiv_threshold_specified, FALSE)
  expect_equal(maxEquivTest_results$significance_level, alpha)
  expect_equal(maxEquivTest_results$num_individuals, 500)
  expect_equal(maxEquivTest_results$num_periods, 5)
  expect_equal(maxEquivTest_results$base_period, 5)
  expect_equal(maxEquivTest_results$minimum_equiv_threshold, 0.36172344249370752545, tolerance = 1e-6)
  expect_equal(maxEquivTest_results$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results$placebo_coefficients, placebo_coefs)
  
  # Equivalence threshold specified:
  maxEquivTest_results2 <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                        period = period_data, X=X_data, 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = cluster_data, 
                                        alpha = alpha,
                                        type = "IU", vcov="CL")
  
  expect_equal(class(maxEquivTest_results2), "maxEquivTestIU")
  expect_equal(maxEquivTest_results2$equiv_threshold_specified, TRUE)
  expect_equal(maxEquivTest_results2$significance_level, alpha)
  expect_equal(maxEquivTest_results2$num_individuals, 500)
  expect_equal(maxEquivTest_results2$num_periods, 5)
  expect_equal(maxEquivTest_results2$base_period, 5)
  expect_equal(maxEquivTest_results2$equiv_threshold, 0.2)
  expect_equal(maxEquivTest_results2$placebo_coefficients_se, beta_se, tolerance = 1e-6)
  expect_equal(length(maxEquivTest_results2$placebo_coefficients), 4)
  expect_equal(maxEquivTest_results2$placebo_coefficients, placebo_coefs)
  expect_equal(maxEquivTest_results2$reject_null_hypothesis, FALSE)
  
  
  # with index input:
  maxEquivTest_results3 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X= c(5,6), 
                                        equiv_threshold = NULL,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov="CL")
  # Equivalence threshold specified:
  maxEquivTest_results4 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X=c(5,6), 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU", vcov = "CL")
  expect_equal(maxEquivTest_results3, maxEquivTest_results)
  expect_equal(maxEquivTest_results4, maxEquivTest_results2)
  
})


# Print functions:
test_that("maxEquivTest print for standard variance-covariance matrix",{
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
  alpha <- 0.05
  
  
  # Do the procedure by hand for comparison:
  subdata <- sim_data[,c("ID", "period", "Y", "placebo_1", "placebo_2", "placebo_3", "placebo_4", "X_1", "X_2")]
  test_formula <- as.formula(Y ~ X_1 + X_2 + placebo_1 + placebo_2 + placebo_3 + placebo_4)
  plm_test <- plm::plm(test_formula, data=subdata, effect="twoways", model="within", index=c("ID","period"))
  placebo_coefs <- plm_test$coefficients[c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  
  vcov_mat <- plm_test$vcov
  subcov_mat <- vcov_mat[c("placebo_1", "placebo_2", "placebo_3", "placebo_4"), c("placebo_1", "placebo_2", "placebo_3", "placebo_4")]
  beta_var <- diag(subcov_mat)
  
  # Calculating the standard errors
  beta_se <- sqrt(beta_var)
  
  maxEquivTest_results <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                       period = period_data, X=X_data, 
                                       equiv_threshold = NULL,
                                       pretreatment_period = pre_treatment_period, 
                                       base_period = base_period,
                                       cluster = cluster_data, 
                                       alpha = alpha,
                                       type = "IU")
  
  expect_snapshot(print(maxEquivTest_results))
  
  # Equivalence threshold specified:
  maxEquivTest_results2 <- maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, 
                                        period = period_data, X=X_data, 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = cluster_data, 
                                        alpha = alpha,
                                        type = "IU")
  
  expect_snapshot(print(maxEquivTest_results2))
  
  # with index input:
  maxEquivTest_results3 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X= c(5,6), 
                                        equiv_threshold = NULL,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU")
  expect_snapshot(print(maxEquivTest_results3))
  
  # Equivalence threshold specified:
  maxEquivTest_results4 <- maxEquivTest(Y = 1, ID= 2, G = 4, 
                                        period = 3, X=c(5,6), 
                                        equiv_threshold = 0.2,
                                        pretreatment_period = pre_treatment_period, 
                                        base_period = base_period,
                                        cluster = 7, data = sim_data,
                                        alpha = alpha,
                                        type = "IU")

  expect_snapshot(print(maxEquivTest_results4))
  
})






