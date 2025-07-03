///////////////////////////////////////////////////////////////////////////
// MCMC sampler C++ script for Bayesian spatio-temporal SAE model        //
// Author: Elliot S. Shannon                                             //
// Email: shann125@msu.edu                                               //
// Affiliation: Michigan State University                                //
// Date Created: July 02, 2025                                           //
// Date Modified: July 02, 2025                                          //
///////////////////////////////////////////////////////////////////////////

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

///////////////////////////////////////////////////////////
// Function to simulate multivariate normal distribution //
///////////////////////////////////////////////////////////

arma::vec rmvn(const arma::vec& mu, const arma::mat& Sigma) {
  arma::mat L = arma::chol(Sigma, "lower");
  arma::vec z = arma::randn<arma::vec>(mu.n_elem);
  return mu + L * z;
}

////////////////////////////
// Inverse-logit function //
////////////////////////////

double logit_inv(double z, double a, double b) {
  return b - (b-a)/(1+exp(z));
}

////////////////////
// Logit function //
////////////////////

double logit(double theta, double a, double b) {
  return log((theta-a)/(b-theta));
}

// [[Rcpp::export]]
void sampler(const int& n_samples,
             const arma::vec& keep, // length n_keep
             const int& batch_length,
             const int& J, // number of areal units (counties)
             const int& T, // number of discrete time steps (years)
             const int& P, // number of covariates associated with beta_t
             const arma::vec& y, // length N, nested as plots within counties within years
             const arma::sp_mat& X, // dimension N x T*P
             const arma::sp_mat& X_tilde, // dimension N x J
             const arma::sp_mat& X_mu, // dimension J*T x T*P
             const arma::sp_mat& X_tilde_mu, // dimension J*T x J
             const arma::sp_mat& A, // dimension N x J*T
             const arma::sp_mat& B, // dimension N x T
             const arma::vec& times, // length N
             const arma::vec& times_indx, // length T 
             const arma::vec& times_count, // length T
             const arma::vec& counties, // length N
             const arma::vec& years, // length N
             const arma::vec& years_indx, // length T
             const arma::vec& years_count, // length T
             const arma::mat& W, // dimension J x J
             const arma::vec& D, // length J
             const arma::mat& S_1, // dimension J x J
             const arma::mat& S_2, // dimension J x J
             const arma::vec& lambda, // length J
             const arma::vec& mu_0, // length P
             const arma::mat& Sigma_0, // dimension P x P
             const arma::mat& H_beta, // dimension P x P
             const double& h_beta, // hyperparameter
             const double& sigma_sq_a, // hyperparameter
             const double& sigma_sq_b, // hyperparameter
             const double& tau_sq_eta_a, // hyperparameter
             const double& tau_sq_eta_b, // hyperparameter
             const double& tau_sq_w_a, // hyperparameter
             const double& tau_sq_w_b, // hyperparameter
             const double& rho_eta_c, // hyperparameter
             const double& rho_eta_d, // hyperparameter
             const double& rho_eta_tuning, // tuning parameter
             const double& rho_w_c, // hyperparameter
             const double& rho_w_d, // hyperparameter
             const double& rho_w_tuning, // tuning parameter
             arma::vec& mu_s, // length J*T
             arma::vec& beta_s, // length T*P
             arma::vec& eta_s, //length J
             arma::vec& u_s, // length J*T
             arma::vec& sigma_sq_s, // length T
             arma::mat& Sigma_beta_s, // dimension P x P
             arma::vec& beta_0_s, // length P
             double& tau_sq_eta_s,
             double& rho_eta_s,
             arma::vec& tau_sq_w_s, // length T
             double& rho_w_s,
             arma::mat& mu_samples, // dimension J*T x n_keep
             arma::mat& beta_samples, // dimension T*P x n_keep
             arma::mat& eta_samples, // dimension J x n_keep
             arma::mat& u_samples, // dimension J*T x n_keep
             arma::mat& sigma_sq_samples, // dimension T x n_keep
             arma::cube& Sigma_beta_samples, // dimension P x P x n_keep
             arma::mat& beta_0_samples, // dimension P x n_keep
             arma::vec& tau_sq_eta_samples, // length n_keep
             arma::vec& rho_eta_samples, // length n_keep
             arma::mat& tau_sq_w_samples, // dimension T x n_keep
             arma::vec& rho_w_samples // length n_keep
) {
  
  std::cout << "Prepping data..." << std::endl;
  std::cout << "----------------------------------" << std::endl;
  
  ///////////////////////
  // Declare variables //
  ///////////////////////
  
  int t; // year
  int s; // sample
  int N = y.size(); // number of response data
  int s_keep = 0; // track number of samples kept
  int batch_iter = 0; // track batch size
  int batch_accept_w = 0; // track acceptance rate for rho_w
  int batch_accept_eta = 0; // track acceptance rate for rho_eta
  arma::vec N_sigma_sq_inv(N); // Large covariance matrix values
  arma::sp_mat Sigma_inv(N, N); // Large covariance matrix
  arma::mat V_beta_0(P, P); // Covariance matrix for updating beta_0
  arma::vec v_beta_0(P); // Vv is mean vector for updating beta_0
  arma::mat V_beta(P, P); // Covariance matrix for updating beta_t
  arma::vec v_beta(P); // Vv is mean vector for updating beta_t
  arma::mat V_eta(J, J); // V is covariance matrix for updating eta
  arma::vec v_eta(J); // Vv is mean vector for updating eta
  arma::mat Q_u_inv(J, J); // Inverse of Q_u
  arma::mat Q_eta_inv(J, J); // Inverse of Q_eta
  arma::mat V_u(J, J); // V is covariance matrix for updating u_t
  arma::vec v_u(J); // Vv is mean vector for updating u_t
  double tau_sq_w_a_tilde; // shape parameter for updating tau_sq_w
  double tau_sq_w_b_tilde; // rate parameter for updating tau_sq_w
  double sigma_sq_a_tilde; // shape parameter for updating sigma_sq
  double sigma_sq_b_tilde; // rate parameter for updating sigma_sq
  double tau_sq_eta_a_tilde; // shape parameter for updating tau_sq_eta
  double tau_sq_eta_b_tilde; // rate parameter for updating tau_sq_eta
  double h_beta_tilde; // degrees of freedom for updating Sigma_beta
  arma::mat H_beta_tilde(P, P); // symmetric PD matrix for updating Sigma_beta
  double rho_eta_cand; // candidate value for rho_eta
  double rho_eta_current_ltd; // log-target-density of current rho_eta 
  double rho_eta_cand_ltd; // log-target-density of candidate rho_eta
  double rho_w_cand; // candidate rho_w
  double rho_w_current_ltd; // log-target-density of current rho_w
  double rho_w_cand_ltd; // log-target-density of candidate rho_w
  double unif; // uniform random variable
  arma::mat rho_eta_current_Q_inv(J, J); // Inverse of matrix Q with current rho_eta
  arma::mat rho_eta_cand_Q_inv(J, J); // Inverse of matrix Q with candidate rho_eta
  arma::mat rho_w_current_Q_inv(J, J); // Inverse of matrix Q with current rho_w
  arma::mat rho_w_cand_Q_inv(J, J); // Inverse of matrix Q with candidate rho_w
  
  ////////////////////////////////////////////
  // Build large diagonal covariance matrix //
  ////////////////////////////////////////////
  
  N_sigma_sq_inv = B * (1.0 / sigma_sq_s);
  Sigma_inv.diag() = N_sigma_sq_inv;
  
  ////////////////////////////
  // Pring starting message //
  ////////////////////////////
  
  std::cout << "Starting sampler..." << std::endl;
  std::cout << "----------------------------------" << std::endl;
  
  auto start = std::chrono::high_resolution_clock::now(); // start timer
  
  //////////////////////////
  // Loop through sampler //
  //////////////////////////
  
  for (s = 0; s < n_samples; s++) {
    
    ///////////////
    // Update mu //
    ///////////////
    
    mu_s = X_mu * beta_s + X_tilde_mu * eta_s + u_s;
    
    ///////////////////
    // Update beta_0 //
    ///////////////////
    
    V_beta_0 = arma::inv(arma::inv(Sigma_0) + arma::inv(Sigma_beta_s));
    v_beta_0 = arma::inv(Sigma_0) * mu_0 + arma::inv(Sigma_beta_s) * beta_s.rows(0, P-1);
    beta_0_s = rmvn(V_beta_0 * v_beta_0, V_beta_0);
    
    /////////////////////////////////////
    // Update parameters in time order //
    /////////////////////////////////////
    
    ///////////
    // t = 1 //
    ///////////
    
    ///////////////////
    // Update beta_1 //
    ///////////////////
    
    V_beta = arma::inv(X.submat(0, 0, times_indx(1) - 1, 1).t() * Sigma_inv.submat(0, 0, times_indx(1) - 1, times_indx(1) - 1) * X.submat(0, 0, times_indx(1) - 1, 1) + 2*arma::inv(Sigma_beta_s));
    v_beta = X.submat(0, 0, times_indx(1) - 1, 1).t() * Sigma_inv.submat(0, 0, times_indx(1) - 1, times_indx(1) - 1) * (y.rows(0, times_indx(1) - 1) - X_tilde.rows(0, times_indx(1) - 1) * eta_s - A.rows(0, times_indx(1) - 1) * u_s) + arma::inv(Sigma_beta_s) * (beta_s.rows(2, 3) + beta_0_s);
    beta_s.rows(0, 1) = rmvn(V_beta * v_beta, V_beta);
    
    ////////////////
    // Update u_1 //
    ////////////////
    
    Q_u_inv = S_1 - rho_w_s * S_2; // fast trick, see computing notes
    V_u = arma::inv_sympd(A.submat(0, 0, times_indx(1) - 1, J - 1).t() * Sigma_inv.submat(0, 0, times_indx(1) - 1, times_indx(1) - 1) * A.submat(0, 0, times_indx(1) - 1, J - 1) + (1/tau_sq_w_s(0) * Q_u_inv) + (1/tau_sq_w_s(1) * Q_u_inv));
    v_u = A.submat(0, 0, times_indx(1) - 1, J - 1).t() * Sigma_inv.submat(0, 0, times_indx(1) - 1, times_indx(1) - 1) * (y.rows(0, times_indx(1) - 1) - X.submat(0, 0, times_indx(1) - 1, 1) * beta_s.rows(0, 1) - X_tilde.rows(0, times_indx(1) - 1) * eta_s) + (1/tau_sq_w_s(1) * Q_u_inv) * u_s.rows(J, 2*J - 1);
    u_s.rows(0, J-1) = rmvn(V_u * v_u, V_u);
    
    ///////////////////////
    // Update sigma_sq_1 //
    ///////////////////////
    
    sigma_sq_a_tilde = sigma_sq_a + 0.5 * years_count(0);
    sigma_sq_b_tilde = sigma_sq_b + 0.5 * arma::dot(y.rows(years_indx(0), years_indx(0) + years_count(0) - 1) - X.submat(years_indx(0), 0, years_indx(0) + years_count(0) - 1, 1) * beta_s.rows(0,1) - X_tilde.rows(years_indx(0), years_indx(0) + years_count(0) - 1) * eta_s - A.rows(years_indx(0), years_indx(0) + years_count(0) - 1) * u_s, y.rows(years_indx(0), years_indx(0) + years_count(0) - 1) - X.submat(years_indx(0), 0, years_indx(0) + years_count(0) - 1, 1) * beta_s.rows(0,1) - X_tilde.rows(years_indx(0), years_indx(0) + years_count(0) - 1) * eta_s - A.rows(years_indx(0), years_indx(0) + years_count(0) - 1) * u_s);
    sigma_sq_s.row(0) = 1/arma::randg(1, arma::distr_param(sigma_sq_a_tilde, 1/sigma_sq_b_tilde));
    
    ///////////////////////
    // Update tau_sq_w_1 //
    ///////////////////////
    
    tau_sq_w_a_tilde = tau_sq_w_a + J/2;
    tau_sq_w_b_tilde = tau_sq_w_b + 0.5 * arma::dot(u_s.rows(0, J - 1).t() * Q_u_inv, u_s.rows(0, J - 1));
    tau_sq_w_s.row(0) = 1/arma::randg(1, arma::distr_param(tau_sq_w_a_tilde, 1/tau_sq_w_b_tilde));
    
    //////////////////////
    // t = 2, ... T - 1 //
    //////////////////////
    
    for (t = 1; t < (T-1); t++) {
      
      ///////////////////
      // Update beta_t //
      ///////////////////
      
      V_beta = arma::inv(X.submat(times_indx(t), 2*t, times_indx(t+1) - 1, 2*t + 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), times_indx(t+1) - 1, times_indx(t+1) - 1) * X.submat(times_indx(t), 2*t, times_indx(t+1) - 1, 2*t + 1) + 2*arma::inv(Sigma_beta_s));
      v_beta = X.submat(times_indx(t), 2*t, times_indx(t+1) - 1, 2*t + 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), times_indx(t+1) - 1, times_indx(t+1) - 1) * (y.rows(times_indx(t), times_indx(t+1) - 1) - X_tilde.rows(times_indx(t), times_indx(t+1)-1) * eta_s - A.rows(times_indx(t), times_indx(t+1) - 1) * u_s) + arma::inv(Sigma_beta_s) * (beta_s.rows(2*(t+1), 2*(t+1)+1) + beta_s.rows(2*(t-1), 2*(t-1)+1));
      beta_s.rows(2*t, 2*t + 1) = rmvn(V_beta * v_beta, V_beta);
      
      ////////////////
      // Update u_t //
      ////////////////
      
      V_u = arma::inv_sympd(A.submat(times_indx(t), t*J, times_indx(t+1) - 1, (t+1)*J - 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), times_indx(t+1) - 1, times_indx(t+1) - 1) * A.submat(times_indx(t), t*J, times_indx(t+1) - 1, (t+1)*J - 1) + (1/tau_sq_w_s(t) * Q_u_inv) + (1/tau_sq_w_s(t+1) * Q_u_inv));
      v_u = A.submat(times_indx(t), t*J, times_indx(t+1) - 1, (t+1)*J - 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), times_indx(t+1) - 1, times_indx(t+1) - 1) * (y.rows(times_indx(t), times_indx(t+1) - 1) - X.submat(times_indx(t), 2*t, times_indx(t+1) - 1, 2*t + 1) * beta_s.rows(2*t, 2*t + 1) - X_tilde.rows(times_indx(t), times_indx(t+1) - 1) * eta_s) + (1/tau_sq_w_s(t) * Q_u_inv) * u_s.rows((t-1)*J, t*J - 1) + (1/tau_sq_w_s(t+1) * Q_u_inv) * u_s.rows((t+1)*J, (t+2)*J - 1);
      u_s.rows(t*J, (t+1)*J-1) = rmvn(V_u * v_u, V_u);
      
      ///////////////////////
      // Update sigma_sq_t //
      ///////////////////////
      
      sigma_sq_a_tilde = sigma_sq_a + 0.5 * years_count(t);
      sigma_sq_b_tilde = sigma_sq_b + 0.5 * arma::dot(y.rows(years_indx(t), years_indx(t) + years_count(t) - 1) - X.submat(years_indx(t), 2*t, years_indx(t) + years_count(t) - 1, 2*t + 1) * beta_s.rows(2*t, 2*t + 1) - X_tilde.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * eta_s - A.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * u_s, y.rows(years_indx(t), years_indx(t) + years_count(t) - 1) - X.submat(years_indx(t), 2*t, years_indx(t) + years_count(t) - 1, 2*t + 1) * beta_s.rows(2*t,2*t + 1) - X_tilde.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * eta_s - A.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * u_s);
      sigma_sq_s.row(t) = 1/arma::randg(1, arma::distr_param(sigma_sq_a_tilde, 1/sigma_sq_b_tilde));
      
      ///////////////////////
      // Update tau_sq_w_t //
      ///////////////////////
      
      tau_sq_w_b_tilde = tau_sq_w_b + 0.5 * arma::dot((u_s.rows(t*J, (t+1)*J - 1) - u_s.rows((t-1)*J, t*J - 1)).t() * Q_u_inv, (u_s.rows(t*J, (t+1)*J - 1) - u_s.rows((t-1)*J, t*J - 1)));
      tau_sq_w_s.row(t) = 1/arma::randg(1, arma::distr_param(tau_sq_w_a_tilde, 1/tau_sq_w_b_tilde));
    }
    
    ///////////
    // t = T //
    ///////////
    
    t = T-1;
    
    ///////////////////
    // Update beta_T //
    ///////////////////
    
    V_beta = arma::inv(X.submat(times_indx(t), 2*t, N - 1, 2*t + 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), N - 1, N - 1) * X.submat(times_indx(t), 2*t, N - 1, 2*t + 1) + arma::inv(Sigma_beta_s));
    v_beta = X.submat(times_indx(t), 2*t, N - 1, 2*t + 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), N - 1, N - 1) * (y.rows(times_indx(t), N - 1) - X_tilde.rows(times_indx(t), N - 1) * eta_s - A.rows(times_indx(t), N - 1) * u_s) + arma::inv(Sigma_beta_s) * beta_s.rows(2*(t-1), 2*(t-1)+1);
    beta_s.rows(2*t, 2*t + 1) = rmvn(V_beta * v_beta, V_beta);
    
    ////////////////
    // Update u_T //
    ////////////////
    
    V_u = arma::inv_sympd(A.submat(times_indx(t), t*J, N - 1, T*J - 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), N - 1, N - 1) * A.submat(times_indx(t), t*J, N - 1, T*J - 1) + (1/tau_sq_w_s(t) * Q_u_inv));
    v_u = A.submat(times_indx(t), t*J, N - 1, T*J - 1).t() * Sigma_inv.submat(times_indx(t), times_indx(t), N - 1, N - 1) * (y.rows(times_indx(t), N - 1) - X.submat(times_indx(t), 2*t, N - 1, 2*t + 1) * beta_s.rows(2*t, 2*t + 1) - X_tilde.rows(times_indx(t), N - 1) * eta_s) + (1/tau_sq_w_s(t) * Q_u_inv) * u_s.rows((t-1)*J, (t)*J - 1);
    u_s.rows(t*J, T*J - 1) = rmvn(V_u * v_u, V_u);
    
    ///////////////////////
    // Update sigma_sq_T //
    ///////////////////////
    
    sigma_sq_a_tilde = sigma_sq_a + 0.5 * years_count(t);
    sigma_sq_b_tilde = sigma_sq_b + 0.5 * arma::dot(y.rows(years_indx(t), years_indx(t) + years_count(t) - 1) - X.submat(years_indx(t), 2*t, years_indx(t) + years_count(t) - 1, 2*t + 1) * beta_s.rows(2*t, 2*t + 1) - X_tilde.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * eta_s - A.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * u_s, y.rows(years_indx(t), years_indx(t) + years_count(t) - 1) - X.submat(years_indx(t), 2*t, years_indx(t) + years_count(t) - 1, 2*t + 1) * beta_s.rows(2*t,2*t + 1) - X_tilde.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * eta_s - A.rows(years_indx(t), years_indx(t) + years_count(t) - 1) * u_s);
    sigma_sq_s.row(t) = 1/arma::randg(1, arma::distr_param(sigma_sq_a_tilde, 1/sigma_sq_b_tilde));
    
    ///////////////////////
    // Update tau_sq_w_T //
    ///////////////////////
    
    tau_sq_w_b_tilde = tau_sq_w_b + 0.5 * arma::dot((u_s.rows(t*J, J*T - 1) - u_s.rows((t-1)*J, t*J - 1)).t() * Q_u_inv, (u_s.rows(t*J, J*T - 1) - u_s.rows((t-1)*J, t*J - 1)));
    tau_sq_w_s.row(t) = 1/arma::randg(1, arma::distr_param(tau_sq_w_a_tilde, 1/tau_sq_w_b_tilde));
    
    ////////////////
    // Update eta //
    ////////////////
    
    V_eta = arma::inv_sympd(X_tilde.t() * Sigma_inv * X_tilde + (1/tau_sq_eta_s * (S_1 - rho_eta_s * S_2)));
    v_eta = X_tilde.t() * Sigma_inv * (y - X * beta_s - A * u_s);
    eta_s = rmvn(V_eta * v_eta, V_eta);
    
    //////////////////////////////////////////// 
    // Build large diagonal covariance matrix //
    ////////////////////////////////////////////
    
    N_sigma_sq_inv = B * (1.0 / sigma_sq_s);
    Sigma_inv.diag() = N_sigma_sq_inv;
    
    ///////////////////////
    // Update tau_sq_eta //
    ///////////////////////
    
    Q_eta_inv = S_1 - rho_eta_s * S_2;
    tau_sq_eta_a_tilde = tau_sq_eta_a + 0.5 * J;
    tau_sq_eta_b_tilde = tau_sq_eta_b + 0.5 * arma::dot(eta_s.t() * Q_eta_inv, eta_s);
    tau_sq_eta_s = 1/arma::randg(1, arma::distr_param(tau_sq_eta_a_tilde, 1/tau_sq_eta_b_tilde))(0);
    
    ///////////////////////
    // Update Sigma_beta //
    ///////////////////////
    
    h_beta_tilde = h_beta + T;
    H_beta_tilde = H_beta + (beta_s.rows(0, 1) - beta_0_s) * (beta_s.rows(0, 1) - beta_0_s).t();
    for (t = 1; t < T; t++) {
      H_beta_tilde += (beta_s.rows(2*t, 2*t + 1) - beta_s.rows(2*t - 2, 2*t - 1)) * (beta_s.rows(2*t, 2*t + 1) - beta_s.rows(2*t - 2, 2*t - 1)).t();
    }
    Sigma_beta_s = arma::iwishrnd(H_beta_tilde, h_beta_tilde);
    
    ////////////////////
    // Update rho_eta //
    ////////////////////
    
    rho_eta_cand = logit_inv(arma::as_scalar(randn(1, arma::distr_param(logit(rho_eta_s, rho_eta_c, rho_eta_d), sqrt(rho_eta_tuning)))), rho_eta_c, rho_eta_d);
    
    rho_eta_current_Q_inv = S_1 - rho_eta_s * S_2;
    rho_eta_cand_Q_inv = S_1 - rho_eta_cand * S_2;
    
    rho_eta_current_ltd = log(rho_eta_s - rho_eta_c) + log(rho_eta_d - rho_eta_s);
    rho_eta_current_ltd += -0.5 * (-J * log(1/tau_sq_eta_s) - sum(log(D % (1 - rho_eta_s * lambda))));
    rho_eta_current_ltd += -0.5 * 1/tau_sq_eta_s * arma::dot(eta_s.t() * rho_eta_current_Q_inv, eta_s);
    rho_eta_cand_ltd = log(rho_eta_cand - rho_eta_c) + log(rho_eta_d - rho_eta_cand);
    rho_eta_cand_ltd += -0.5 * (-J * log(1/tau_sq_eta_s) - sum(log(D % (1 - rho_eta_cand * lambda))));
    rho_eta_cand_ltd += -0.5 * 1/tau_sq_eta_s * arma::dot(eta_s.t() * rho_eta_cand_Q_inv, eta_s);
    
    unif = arma::randu();
    
    if (unif < exp(rho_eta_cand_ltd - rho_eta_current_ltd)) {
      rho_eta_s = rho_eta_cand;
      batch_accept_eta += 1;
    }
    
    //////////////////
    // Update rho_w //
    //////////////////
    
    rho_w_cand = logit_inv(arma::as_scalar(randn(1, arma::distr_param(logit(rho_w_s, rho_w_c, rho_w_d), sqrt(rho_w_tuning)))), rho_w_c, rho_w_d);
    
    rho_w_current_Q_inv = S_1 - rho_w_s * S_2;
    rho_w_cand_Q_inv = S_1 - rho_w_cand * S_2;
    
    rho_w_current_ltd = log(rho_w_s - rho_w_c) + log(rho_w_d - rho_w_s);
    rho_w_cand_ltd = log(rho_w_cand - rho_w_c) + log(rho_w_d - rho_w_cand);
    
    rho_w_current_ltd += -0.5 * (-J * log(1/tau_sq_w_s(0)) - sum(log(D % (1 - rho_w_s * lambda))));
    rho_w_current_ltd += -0.5 * 1/tau_sq_w_s(0) * arma::dot(u_s.rows(0, J - 1).t() * rho_w_current_Q_inv, u_s.rows(0, J - 1));
    rho_w_cand_ltd += -0.5 * (-J * log(1/tau_sq_w_s(0)) - sum(log(D % (1 - rho_w_cand * lambda))));
    rho_w_cand_ltd += -0.5 * 1/tau_sq_w_s(0) * arma::dot(u_s.rows(0, J - 1).t() * rho_w_cand_Q_inv, u_s.rows(0, J - 1));
    
    for (t = 1; t < T; t++) {
      rho_w_current_ltd += -0.5 * (-J * log(1/tau_sq_w_s(t)) - sum(log(D % (1 - rho_w_s * lambda))));
      rho_w_current_ltd += -0.5 * 1/tau_sq_w_s(t) * arma::dot((u_s.rows(t*J, (t+1)*J - 1) - u_s.rows((t-1)*J, (t)*J - 1)).t() * rho_w_current_Q_inv, (u_s.rows(t*J, (t+1)*J - 1) - u_s.rows((t-1)*J, (t)*J - 1)));
      rho_w_cand_ltd += -0.5 * (-J * log(1/tau_sq_w_s(t)) - sum(log(D % (1 - rho_w_cand * lambda))));
      rho_w_cand_ltd += -0.5 * 1/tau_sq_w_s(t) * arma::dot((u_s.rows(t*J, (t+1)*J - 1) - u_s.rows((t-1)*J, (t)*J - 1)).t() * rho_w_cand_Q_inv, (u_s.rows(t*J, (t+1)*J - 1) - u_s.rows((t-1)*J, (t)*J - 1)));
    }
    
    unif = arma::randu();
    
    if (unif < exp(rho_w_cand_ltd - rho_w_current_ltd)) {
      rho_w_s = rho_w_cand;
      batch_accept_w += 1;
    }
    
    ///////////////////
    // Store samples //
    ///////////////////
    
    if (std::count(keep.begin(), keep.end(), s) > 0) {
      mu_samples.col(s_keep) = mu_s;
      beta_0_samples.col(s_keep) = beta_0_s;
      beta_samples.col(s_keep) = beta_s;
      eta_samples.col(s_keep) = eta_s;
      u_samples.col(s_keep) = u_s;
      sigma_sq_samples.col(s_keep) = sigma_sq_s;
      tau_sq_w_samples.col(s_keep) = tau_sq_w_s;
      tau_sq_eta_samples(s_keep) = tau_sq_eta_s;
      Sigma_beta_samples.slice(s_keep) = Sigma_beta_s;
      rho_eta_samples(s_keep) = rho_eta_s;
      rho_w_samples(s_keep) = rho_w_s;
      s_keep += 1;
    }
    
    //////////////////////////////
    // Report batch information //
    //////////////////////////////
    
    batch_iter += 1;
    if (batch_iter == batch_length) {
      std::cout << std::fixed << std::setprecision(2) << "            Complete: " << 100*(s+1)/n_samples << "%" << std::endl;
      std::cout << std::fixed << std::setprecision(2) << "Batch accept   rho_w: " << 100*batch_accept_w/batch_length << std::endl;
      std::cout << std::fixed << std::setprecision(2) << "Batch accept rho_eta: " << 100*batch_accept_eta/batch_length << std::endl;
      auto stop = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::ratio<3600>> elapsed_time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<3600>>>(stop - start);
      std::cout << std::fixed << std::setprecision(3) << "                Time: " << elapsed_time.count() << " hours" << std::endl;
      std::cout << "----------------------------------" << std::endl;
      
      batch_accept_w = 0;
      batch_accept_eta = 0;
      batch_iter = 0;
    }
  }
}