
# Test for input validation
test_that("sim_paneldata input", {
  #N is not a positive integer
  expect_error(sim_paneldata(N = -1, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 1), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = 1, burnins = 100))
  
  #tt is not a positive integer
  expect_error(sim_paneldata(N = 500, tt = 0, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = 1, burnins = 100))
  
  # beta length does not match tt
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 4), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = 1, burnins = 100))
  
  # phi not in interval [0,1)
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 0, 
                             phi = 1, sd = 1, burnins = 100))
  
  # all inputs correct
  expect_no_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1),
                                eta = rep(0, 500), lambda = rep(0, 5), het = 0, 
                                phi = 0.5, sd = 1, burnins = 100))
  
  # p is not a non-negative integer
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = -1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = 1, burnins = 100))
  
  # het is not 0 or 1
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 2, phi = 0, 
                             sd = 1, burnins = 100))
  
  # sd is not a positive number
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = -1, burnins = 100))
  
  # burnins is not a positive integer
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = 1, burnins = -1))
  
  # eta length does not match N
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 499), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = 1, burnins = 100))
  
  # lambda length does not match tt
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), 
                             eta = rep(0, 500), lambda = rep(0, 4), het = 0, phi = 0, 
                             sd = 1, burnins = 100))
  
  # gamma length does not match p
  expect_error(sim_paneldata(N = 500, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 2), 
                             eta = rep(0, 500), lambda = rep(0, 5), het = 0, phi = 0, 
                             sd = 1, burnins = 100))
})

# Test for function output structure
test_that("sim_paneldata output", {
  N <- 500
  tt <- 5
  p <- 1
  sim_data <- sim_paneldata(N = N, tt = tt, beta = rep(0, tt), p = p, gamma = rep(1, p), het = 1, phi = 0.5, sd = 1, burnins = 100)
  
  # Check if output is a data frame
  expect_s3_class(sim_data, "data.frame")
  
  # Check the number of rows and columns
  expect_equal(nrow(sim_data), N * tt)
  expect_equal(ncol(sim_data), 4 + p)
  
  # Check the column names
  expected_colnames <- c("ID", "period", "Y", "G", paste0("X_", 1:p))
  expect_equal(colnames(sim_data), expected_colnames)
  
  # Check column types
  expect_type(sim_data$ID, "integer")
  expect_type(sim_data$period, "integer")
  expect_type(sim_data$Y, "double")
  
  for(i in 1:p) {
    expect_type(sim_data[[paste0("X_", i)]], "double")
  }
})

# Additional tests

# Boundary Testing: Minimum values
test_that("sim_paneldata works with minimum values", {
  sim_data <- sim_paneldata(N = 1, tt = 1, beta = rep(0, 1), p = 1, gamma = rep(1, 1), het = 0, phi = 0, sd = 1, burnins = 1)
  expect_equal(nrow(sim_data), 1)
  expect_equal(ncol(sim_data), 5)
  expect_s3_class(sim_data, "data.frame")
})

# Randomness Check: Ensure reproducibility with set seed
test_that("sim_paneldata produces consistent output with set seed", {
  set.seed(123)
  sim_data1 <- sim_paneldata(N = 10, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), het = 0, phi = 0, sd = 1, burnins = 100)
  
  set.seed(123)
  sim_data2 <- sim_paneldata(N = 10, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), het = 0, phi = 0, sd = 1, burnins = 100)
  
  expect_equal(sim_data1, sim_data2)
})

# Heteroskedasticity and AR(1) Process Check
test_that("sim_paneldata handles heteroskedasticity and AR(1) process correctly", {
  sim_data <- sim_paneldata(N = 10, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), het = 1, phi = 0.5, sd = 1, burnins = 100)
  expect_s3_class(sim_data, "data.frame")
  expect_equal(nrow(sim_data), 10 * 5)
})

# Zero Regressors
test_that("sim_paneldata works with zero additional regressors", {
  sim_data <- sim_paneldata(N = 10, tt = 5, beta = rep(0, 5), p = 0, gamma = numeric(0), het = 0, phi = 0, sd = 1, burnins = 100)
  expect_equal(ncol(sim_data), 4)
  expect_s3_class(sim_data, "data.frame")
})

# Extreme Values
test_that("sim_paneldata works with extreme values", {
  sim_data <- sim_paneldata(N = 1000, tt = 100, beta = rep(0, 100), p = 10, gamma = rep(1, 10), het = 0, phi = 0, sd = 10, burnins = 100)
  expect_s3_class(sim_data, "data.frame")
  expect_equal(nrow(sim_data), 1000 * 100)
  expect_equal(ncol(sim_data), 4 + 10)
})

# Edge Cases for phi
test_that("sim_paneldata handles edge cases for phi", {
  sim_data <- sim_paneldata(N = 10, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), het = 0, phi = 0, sd = 1, burnins = 100)
  expect_s3_class(sim_data, "data.frame")
  
  sim_data <- sim_paneldata(N = 10, tt = 5, beta = rep(0, 5), p = 1, gamma = rep(1, 1), het = 0, phi = 0.999, sd = 1, burnins = 100)
  expect_s3_class(sim_data, "data.frame")
})

