#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]

// Function that calculates the mean for each group:

// [[Rcpp::export]]
arma::vec groupMeans(arma::vec x, arma::vec group) {
  arma::vec unique_groups = arma::unique(group);
  arma::vec means(unique_groups.n_elem);
  
  for(unsigned int i = 0; i < unique_groups.n_elem; i++) {
    means[i] = arma::mean(x.elem(find(group == unique_groups[i])));
  }
  
  return means;
}

struct BootstrapWorker : public RcppParallel::Worker {
  // input data
  const arma::vec& Xb;
  const arma::mat& X;
  const arma::vec& ID;
  const arma::vec& period;
  const double variance;
  const int no_placebos;
  const int B;
  const arma::mat& error_terms;
  
  // output data
  arma::vec& max_abs_bootstrap_coef;
  
  // constructor
  BootstrapWorker(const arma::vec& Xb, const arma::mat& X, const arma::vec& ID, const arma::vec& period, const double variance, const int no_placebos, const int B, const arma::mat& error_terms, arma::vec& max_abs_bootstrap_coef)
    : Xb(Xb), X(X), ID(ID), period(period), variance(variance), no_placebos(no_placebos), B(B), error_terms(error_terms), max_abs_bootstrap_coef(max_abs_bootstrap_coef) {}
  
  // worker function
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t b = begin; b < end; b++) {
      arma::vec boot_u_ddot = error_terms.col(b);
      arma::vec boot_y_ddot = Xb + boot_u_ddot;
      
      arma::vec boot_Y_ID_mean = groupMeans(boot_y_ddot, ID);
      arma::vec boot_Y_period_mean = groupMeans(boot_y_ddot, period);
      double boot_Y_mean_total = arma::mean(boot_y_ddot);
      
      arma::uvec ID_indeces = arma::conv_to<arma::uvec>::from(ID-1);    
      arma::uvec period_indeces = arma::conv_to<arma::uvec>::from(period-1);
      arma::vec double_demeaned_boot_Y = boot_y_ddot - boot_Y_ID_mean.elem(ID_indeces) - boot_Y_period_mean.elem(period_indeces) + boot_Y_mean_total;
      
      arma::colvec boot_coef = solve(X.t() * X, X.t() * double_demeaned_boot_Y);
      arma::colvec placebo_boot_coefs = boot_coef.head(no_placebos);
      max_abs_bootstrap_coef[b] = max(abs(placebo_boot_coefs));
      //    Rcpp::Rcout << "Processing iteration " << b << ": " << max(abs(placebo_boot_coefs)) << std::endl;
    }
  }
};

// [[Rcpp::export]]
arma::vec maxTestBoot_bootstrap(const arma::vec& Xb, const arma::mat& X,
                                 const int B,
                                 const double variance,
                                 const arma::vec& ID,
                                 const arma::vec& period,
                                 const int no_placebos){
  
  //Creating a matrix with the error terms for each bootstrap iteration:
  arma::mat error_terms = arma::randn<arma::mat>(X.n_rows, B) * sqrt(variance);
  
  arma::vec max_abs_bootstrap_coef(B);
  
  BootstrapWorker worker(Xb, X, ID, period, variance, no_placebos, B, error_terms, max_abs_bootstrap_coef);
  RcppParallel::parallelFor(0, B, worker);
  
  return worker.max_abs_bootstrap_coef;
}



// The Wild Bootstrap:
struct WildBootstrapWorker : public RcppParallel::Worker {
  // source data
  const arma::vec& Xb;
  const arma::mat& X;
  const arma::vec& u_ddot;
  const arma::vec& ID;
  const arma::vec& period;
  const int no_placebos;
  const arma::vec& unique_ID;
  const int N;
  
  // destination data
  arma::vec& max_abs_bootstrap_coef;
  
  
  // initialize with source and destination
  WildBootstrapWorker(const arma::vec& Xb, const arma::mat& X, const arma::vec& u_ddot, const arma::vec& ID, const arma::vec& period, const int no_placebos, const arma::vec& unique_ID, const int N, arma::vec& max_abs_bootstrap_coef) 
    : Xb(Xb), X(X), u_ddot(u_ddot), ID(ID), period(period), no_placebos(no_placebos), unique_ID(unique_ID), N(N), max_abs_bootstrap_coef(max_abs_bootstrap_coef) {}
  
  // take the square root of the range of elements requested
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t b = begin; b < end; b++) {
      // Obtain the rademacher variable:
      arma::vec R_vars =  arma::randi<arma::vec>(N, arma::distr_param(0,1));
      R_vars.transform([](double val){return val == 0 ? -1 : 1; });
      arma::uvec ID_indeces = arma::conv_to<arma::uvec>::from(ID-1);
      arma::vec R_vars_full = R_vars.elem(ID_indeces);
      
      // Draw the error-term, boot.u.ddot, from a random normal distribution with mean 0 and variance "variance":
      arma::vec boot_u_ddot = u_ddot % R_vars_full;
      // Calculate the bootstrapped dependent variable, boot.y.ddot:
      arma::vec boot_y_ddot = Xb + boot_u_ddot;
      
      // We double-demean the dependent variable
      arma::vec boot_Y_ID_mean = groupMeans(boot_y_ddot, ID);
      arma::vec boot_Y_period_mean = groupMeans(boot_y_ddot, period);
      double boot_Y_mean_total = mean(boot_y_ddot);
      
      arma::colvec double_demeaned_boot_Y(X.n_rows);
      arma::uvec period_indeces = arma::conv_to<arma::uvec>::from(period-1);
      double_demeaned_boot_Y = boot_y_ddot - boot_Y_ID_mean.elem(ID_indeces) - boot_Y_period_mean.elem(period_indeces) + boot_Y_mean_total;
      
      // for (int i = 0; i < X.n_rows; i++) {
      //   double_demeaned_boot_Y[i] = boot_y_ddot[i] - boot_Y_ID_mean[ID[i] - 1] - boot_Y_period_mean[period[i] - 1] + boot_Y_mean_total_vec;
      // }
      
      // Calculate the bootstrapped coefficients:
      arma::colvec boot_coef = solve(X.t() * X, X.t() * double_demeaned_boot_Y);
      // obtaining the first no.placebo elements of the vector
      // i.e. those corresponding to the placebo variables:
      arma::colvec placebo_boot_coefs = boot_coef.head(no_placebos);
      max_abs_bootstrap_coef[b] = max(abs(placebo_boot_coefs));
    }
  }
};

// [[Rcpp::export]]
arma::vec maxTestBoot_wildbootstrap(const arma::vec& Xb, const arma::mat& X,
                                    int B,
                                    const arma::vec& u_ddot,
                                    const arma::vec& ID,
                                    const arma::vec& period,
                                    int no_placebos){
  arma::vec unique_ID = arma::unique(ID);
  int N = unique_ID.n_elem;
  // Matrix to store the bootstrap coefficients in:
  arma::vec max_abs_bootstrap_coef(B);
  
  // create the worker
  WildBootstrapWorker worker(Xb, X, u_ddot, ID, period, no_placebos, unique_ID, N, max_abs_bootstrap_coef);
  
  // call it with parallelFor
  parallelFor(0, B, worker);
  
  return max_abs_bootstrap_coef;
}


