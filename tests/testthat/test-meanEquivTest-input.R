test_that("Input data meanEquivTest",{
  # Data
  sim_data <- sim_paneldata(N = 100, tt = 5, p = 2, gamma = c(1, 1), beta = rep(0, 5), phi = 0, het = 0, sd = 1, burnins = 50)
  cluster_data <- ifelse(sim_data$ID %% 3 ==0, 1, 2)
  sim_data$cluster <- cluster_data
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, paste0("X_", 1:2)]
  
  # Test if the function does not give an error if we use the data input and if we do not use the data option:
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5))
  expect_no_error(meanEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, cluster = cluster_data, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5))
  
  # expect an error if both data AND Y, ID, G, period, X 
  expect_error(meanEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, 
                             cluster = cluster_data, equiv_threshold = 1, pretreatment_period = 1:5,
                             base_period = 5, data = sim_data))
  
  # Expect error if data is supplied and a vector:
  expect_error(meanEquivTest(Y= Y_data, ID= 1, G = 4, period = 2, X= c(5,6), cluster = 7, equiv_threshold = 1, pretreatment_period = 1:5,
                  base_period = 5, data = sim_data))
  expect_error(meanEquivTest(Y = 3, ID=ID_data, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  expect_error(meanEquivTest(Y = 3, ID=1, G = G_data, period = 2, X=c(5,6), cluster =7, data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  expect_error(meanEquivTest(Y = 3, ID=1, G = 4, period = period_data, X=c(5,6), cluster =7, data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  expect_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=X_matrix, cluster =7, data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  expect_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), cluster =cluster_data, data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  
  # Expect error if Y, ID, G, period, X are not vectors/matrices of equal length:
  expect_error(meanEquivTest(Y= Y_data[1:100], ID= ID_data, G = G_data, period = period_data, X=X_data, cluster = cluster_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  expect_error(meanEquivTest(Y= Y_data, ID= ID_data[1:100], G = G_data, period = period_data, X=X_data, cluster = 
                               cluster_data, equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= Y_data, ID= ID_data, G = G_data[1:100], period = period_data, X=X_data, cluster = 
                               cluster_data, equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data[1:100], X=X_data, cluster = 
                               cluster_data, equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data[1:100,], cluster = 
                               cluster_data, equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, cluster = 
                               cluster_data[1:100], equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if Y, ID, G, period and cluster are multi-column matrices:
  expect_error(meanEquivTest(Y= c(1,3), ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= 3, ID= c(1,1), G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= 3, ID= 1, G = c(4,4), period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = c(2,2), X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =c(7,7), data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if G is not a logical vector:
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 1, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if period is not numeric:
  string_period <- as.character(period_data)
  expect_error(meanEquivTest(Y= Y_data, ID= ID_data, G = G_data, period = string_period, X=X_data, 
                             cluster = cluster_data, equiv_threshold = 1, pretreatment_period = 1:5,
                             base_period = 5))
}
          
          
)

test_that("Remaining data inputs meanEquivTest",{
  # Data
  sim_data <- sim_paneldata(N = 100, tt = 5, p = 2, gamma = c(1, 1), beta = rep(0, 5), phi = 0, het = 0, sd = 1, burnins = 50)
  cluster_data <- ifelse(sim_data$ID %% 3 ==0, 1, 2)
  sim_data$cluster <- cluster_data
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, paste0("X_", 1:2)]
  
  # expect an error if pre-treatment period is not a subset of period:
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:8, base_period = 5))
  
  # expect an error if the base period is not a scalar or not in pre_treatment period:
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = c(4,5)))
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 6))
  
  # Expect error if equiv threshold is negative or not an numeric value or not a scalar:
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = -1, pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = "a", pretreatment_period = 1:5, base_period = 5))
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = c(1,2), pretreatment_period = 1:5, base_period = 5))
  
  # Expect error if alpha is not between 0 and 1:
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5, alpha = 0))
  # Expect error if alpha is not numeric:
  expect_error(meanEquivTest(Y= 3, ID= 1, G = 4, period = 2, X=c(5,6), cluster =7, data = sim_data,
                             equiv_threshold = 1, pretreatment_period = 1:5, base_period = 5, alpha = "a"))
})


test_that("vcov input meanEquivTest", {
  # Data
  sim_data <- sim_paneldata(N = 100, tt = 5, p = 2, gamma = c(1, 1), beta = rep(0, 5), phi = 0, het = 0, sd = 1, burnins = 50)
  
  # The data in vector/matrix form:
  Y_data <- sim_data$Y
  ID_data <- sim_data$ID
  G_data <- sim_data$G
  period_data <- sim_data$period
  X_data <- sim_data[, paste0("X_", 1:2)]
  
  # For the procedure:
  # no error should be obtained for the different variance-covariance:
  # - No vcov argument:  
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                                base_period = 5))
  expect_no_error(meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5))
  expect_no_error(meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                                base_period = 5))
  # - vcov argument: "HC:
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5, vcov = "HC"))
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                                base_period = 5, vcov = "HC"))
  
  # - vcov argument: "HAC":
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5, vcov = "HAC"))
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                                base_period = 5, vcov = "HAC"))
  
  # - vcov argument: "CL":
  #   + No cluster:
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5, vcov = "HAC"))
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                                base_period = 5, vcov = "HAC"))
  #   + Cluster:
  cluster_data <- ifelse(sim_data$ID %% 3 ==0, 1, 2)
  expect_no_error(meanEquivTest(Y = Y_data, ID= ID_data, G = G_data, period = period_data, X=X_data, equiv_threshold = NULL, pretreatment_period = 1:5,
                                base_period = 5, vcov = "CL", cluster = cluster_data))
  
  # - custom vcovs:
  vcov_func1 <- function(x) {plm::vcovHC(x, method = "white1", type = "HC2")}
  expect_no_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                                base_period = 5, vcov = vcov_func1))
  
  # invalid vcov argument:
  vcov_func2 <- function(x){x}
  expect_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                             base_period = 5, vcov = vcov_func2))
  
  expect_error(meanEquivTest(Y = 3, ID=1, G = 4, period = 2, X=c(5,6), data = sim_data, equiv_threshold = 1, pretreatment_period = 1:5,
                             base_period = 5, vcov = "INV"))
  
})



