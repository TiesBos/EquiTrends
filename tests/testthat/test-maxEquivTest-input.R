test_that("vcov input maxEquivTest, type = IU", {
  # Data
  sim_data <- sim_paneldata(N = 100, tt = 10, p = 2, gamma = c(1, 1), beta = rep(0, 10), phi = 0, het = 0, sd = 1, burnins = 50)
  
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, paste0("X_", 1:2)]
  
  # For the IU procedure:
    # no error should be obtained for the different variance-covariance:
    # - No vcov argument:  
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5, type = "IU"))
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                  base_period = 5, type = "IU"))
  expect_no_error(maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5, type = "IU"))
  expect_no_error(maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                  base_period = 5, type = "IU"))
    # - vcov argument: "HC:
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = "HC"))
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = "HC"))
  
  # - vcov argument: "HAC":
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = "HAC"))
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = "HAC"))

  # - vcov argument: "CL":
  #   + No cluster:
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = "HAC"))
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = "HAC"))
  #   + Cluster:
  cluster_data <- ifelse(sim_data$ID %% 3 ==0, 1, 2)
  expect_no_error(maxEquivTest(Y = Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = "CL", cluster = cluster_data))
  
  # - custom vcovs:
  vcov_func1 <- function(x) {plm::vcovHC(x, method = "white1", type = "HC2")}
  expect_no_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5, type = "IU", vcov = vcov_func1))
  
  # invalid vcov argument:
  vcov_func2 <- function(x){x}
  expect_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
               base_period = 5, type = "IU", vcov = vcov_func2))
  
  expect_error(maxEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
               base_period = 5, type = "IU", vcov = "INV"))
  
})
