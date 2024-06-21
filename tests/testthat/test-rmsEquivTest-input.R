test_that("Input data rmsEquivTest",{
  # Data
  sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, paste0("X_", 1:2)]
  
  # Test if the function does not give an error if we use the data input and if we do not use the data option:
  expect_no_error(rmsEquivTest(Y = 1, ID=2, G = 4, period = 3, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                               base_period = 5))
  expect_no_error(rmsEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = 1, pretreatment_period = 1:5,
                               base_period = 5))
  
  # expect an error if both data AND Y, ID, G, period, X 
  expect_error(rmsEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, 
                            equiv_threshold = 1, pretreatment_period = 1:5,
                            base_period = 5, data = sim_data))
  
  # Expect error if data is supplied and a vector:
  expect_error(rmsEquivTest(Y= Y_data, ID= 1, G = 4, period = 2, X= c(5,6), equiv_threshold = 100, pretreatment_period = 1:5,
                            base_period = 5, data = sim_data))
  expect_error(rmsEquivTest(Y = 1, ID=ID_data, G = 4, period = 3, X=c(5,6), data = sim_data, equiv_threshold = 100, pretreatment_period = 1:5,
                            base_period = 5))
  expect_error(rmsEquivTest(Y = 1, ID=2, G = G_data, period = 3, X=c(5,6), data = sim_data, equiv_threshold = 100, pretreatment_period = 1:5,
                            base_period = 5))
  expect_error(rmsEquivTest(Y = 1, ID=2, G = 4, period = period_data, X=c(5,6), data = sim_data, equiv_threshold = 100, pretreatment_period = 1:5,
                            base_period = 5))
  expect_error(rmsEquivTest(Y = 1, ID=2, G = 4, period = 3, X=X_matrix, data = sim_data, equiv_threshold = 100, pretreatment_period = 1:5,
                            base_period = 5))
  
  # Expect error if Y, ID, G, period, X are not vectors/matrices of equal length:
  expect_error(rmsEquivTest(Y= Y_data[1:100], ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = 100, pretreatment_period = 1:5,
                            base_period = 5))
  expect_error(rmsEquivTest(Y= Y_data, ID= ID_data[1:100], G = G_data, period = period_data, X=X_data, equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= Y_data, ID= ID_data, G = G_data[1:100], period = period_data, X=X_data, equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data[1:100], X=X_data, equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data[1:100,], equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if Y, ID, G, period are multi-column matrices:
  expect_error(rmsEquivTest(Y= c(1,3), ID= 3, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= 1, ID= c(2,2), G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = c(4,4), period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = c(3,3), X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if G is not a logical vector:
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 1, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if period is not numeric:
  string_period <- as.character(period_data)
  expect_error(rmsEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = string_period, X=X_data, 
                            equiv_threshold = 100, pretreatment_period = 1:5,
                            base_period = 5))
})

test_that("Remaining data inputs rmsEquivTest",{
  # Data
  sim_data <- readRDS(test_path("fixtures", "test_data.rds"))
  cluster_data <- ifelse(sim_data$ID %% 3 ==0, 1, 2)
  sim_data$cluster <- cluster_data
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, paste0("X_", 1:2)]
  
  # expect an error if pre-treatment period is not a subset of period:
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6),  data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:8, base_period = 5))
  
  # expect an error if the base period is not a scalar or not in pre_treatment period:
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = c(4,5)))
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 6))
  
  # Expect error if equiv threshold is negative or not an numeric value or not a scalar:
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = -1, pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = "a", pretreatment_period = 1:5, base_period = 5))
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = c(1,2), pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if alpha is not between in 0.01, 0.025, 0.05, 0.1 or 0.2:
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5, alpha = 0.15))
  # Expect error if alpha is not numeric:
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6),  data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5, alpha = "a"))
  
  # Expect error if no_lambda is non-integer or non-numeric:
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5, no_lambda = 1.5))
  expect_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5, no_lambda = "a"))
  
  expect_no_error(rmsEquivTest(Y= 1, ID= 2, G = 4, period = 3, X=c(5,6), data = sim_data,
                            equiv_threshold = 100, pretreatment_period = 1:5, base_period = 5, no_lambda = 10))
})


