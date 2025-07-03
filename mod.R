################################################################
### R script for running Bayesian spatio-temporal SAE model  ###
### Author: Elliot S. Shannon                                ###
### Email: shann125@msu.edu                                  ###
### Affiliation: Michigan State University                   ###
### Date Created: July 02, 2025                              ###
### Date Modified: July 03, 2025                             ###
### R version: 4.3.1 (2023-06-16)                            ###
################################################################

rm(list = ls())

######################
### Load Libraries ###
######################

library(sf) # sf_1.0-19
library(dplyr) # dplyr_1.1.4
library(Rcpp) # Rcpp_1.0.13-1
library(RcppArmadillo) # RcppArmadillo_14.2.2-1
library(Matrix) # Matrix_1.6-1

set.seed(123)

############################
### Compile the C++ code ###
############################

Rcpp::sourceCpp("sampler.cpp")

########################
### Load in the data ###
########################

data <- readRDS("data.rds")

########################
### Prepare the data ###
########################

plots <- data$plots %>% arrange(t)
counties <- data$counties
W <- data$W
tcc <- data$x_full # tree canopy cover
y <- plots$y # response
x <- plots$x # covariate
J <- nrow(counties) # number of areal units (counties)
T <- max(plots$t) # number of discrete time steps (years)
P <- 2 # number of covariates (including intercept)
N <- length(y)

############################################
### CAR precision matrix pre-computation ###
############################################

D <- diag(rowSums(W))
D.sqrt.inv <- diag(1 / sqrt(diag(D)))
eig <- eigen(D.sqrt.inv %*% W %*% D.sqrt.inv)
DP <- sqrt(D) %*% eig$vectors
lambda <- eig$values

S_1 <- matrix(0, nrow = J, ncol = J)
S_2 <- matrix(0, nrow = J, ncol = J)
for (j in 1:J) {
  S_1 <- S_1 + (DP[, j] %*% t(DP[, j]))
  S_2 <- S_2 + lambda[j] * (DP[, j] %*% t(DP[, j]))
}

#####################################################################
### N x (P*T) design matrix of covariates with temporally-varying ###
### impact on response y                                          ###
#####################################################################

X <- sparseMatrix(
  c(1:N, 1:N),
  c(plots$t * 2 - 1, plots$t * 2),
  x = c(rep(1, times = N), x)
)

####################################################################
### N x (Q*J) design matrix of covariates with spatially-varying ###
### impact on response y (we assume Q == 1)                      ###
####################################################################

X_tilde <- sparseMatrix(1:N, plots$j, x = x)

#########################################################################
### (J*T) x (P*T) design matrix of covariates with temporally-varying ###
### impact on latent mean mu                                          ###
#########################################################################

X_mu <- sparseMatrix(
  c(1:(J * T), 1:(J * T)),
  c(tcc$t * 2 - 1, tcc$t * 2),
  x = c(rep(1, times = J * T), tcc$x)
)

##########################################################################
### (J*T) x (Q*J) design matrix of covariates  with spatially-varying  ###
### impact on latent mean mu                                           ###
##########################################################################

X_tilde_mu <- sparseMatrix(1:(J * T), tcc$j, x = tcc$x)

################################
### N x (J*T) mapping matrix ###
################################

A <- sparseMatrix(
  c(1:N),
  J * (plots$t - 1) + plots$j,
  x = rep(1, times = N)
)

############################
### N x T mapping matrix ###
############################

B <- sparseMatrix(
  c(1:N),
  plots$t,
  x = rep(1, times = N)
)

#########################
### Priors and tuning ###
#########################

mu_0 <- rep(0, times = P)
Sigma_0 <- diag(100, nrow = P)

H_beta <- diag(100, nrow = P)
h_beta <- 10

sigma_sq_a <- 2
sigma_sq_b <- 100

tau_sq_eta_a <- 2
tau_sq_eta_b <- 100

tau_sq_w_a <- 2
tau_sq_w_b <- 100

rho_eta_c <- 0
rho_eta_d <- 1

rho_w_c <- 0
rho_w_d <- 1

rho_eta_tuning <- 10
rho_w_tuning <- 2

#######################
### Starting values ###
#######################

mu.s <- rep(0, times = J * T)
beta.s <- rep(0, times = T * P)
eta.s <- rep(0, times = J)
u.s <- rep(0, times = J * T)
sigma_sq.s <- rep(100, times = T)
Sigma_beta.s <- diag(100, nrow = P)
beta_0.s <- rep(0, times = P)
tau_sq_eta.s <- 100
rho_eta.s <- 0.5
tau_sq_w.s <- rep(10, times = T)
rho_w.s <- 0.5

#####################
### Store samples ###
#####################

n_samples <- 50000 # number of samples to take
n_burn <- 25000 # number of burn-in samples
iter <- 25 # thinning interval
keep <- seq(from = n_burn, to = n_samples, by = iter)
n_keep <- length(keep)
batch_length <- 500 # batch length to report on console

mu.samples <- matrix(-999, nrow = J * T, ncol = n_keep)
beta.samples <- matrix(-999, nrow = T * P, ncol = n_keep)
eta.samples <- matrix(-999, nrow = J, ncol = n_keep)
u.samples <- matrix(-999, nrow = J * T, ncol = n_keep)
sigma_sq.samples <- matrix(-999, nrow = T, ncol = n_keep)
Sigma_beta.samples <- array(-999, dim = c(P, P, n_keep))
beta_0.samples <- matrix(-999, nrow = P, ncol = n_keep)
tau_sq_eta.samples <- rep(-999, times = n_keep)
rho_eta.samples <- rep(-999, times = n_keep)
tau_sq_w.samples <- matrix(-999, nrow = T, ncol = n_keep)
rho_w.samples <- rep(-999, times = n_keep)

#######################
### Run the sampler ###
#######################

sampler(
  n_samples = n_samples,
  keep = keep - 1,
  batch_length = batch_length,
  J = J,
  T = T,
  P = P,
  y = y,
  X = X,
  X_tilde = X_tilde,
  X_mu = X_mu,
  X_tilde_mu = X_tilde_mu,
  A = A,
  B = B,
  times = plots$t,
  times_indx = match(unique(plots$t), plots$t) - 1,
  times_count = unname(table(plots$t)),
  counties = plots$j,
  years = plots$t,
  years_indx = (plots %>%
    group_by(t) %>%
    count() %>%
    pull(n) %>%
    cumsum() %>%
    append(0) %>%
    sort())[1:T],
  years_count = plots %>%
    group_by(t) %>%
    count() %>%
    pull(n),
  W = W,
  D = rowSums(W),
  S_1 = S_1,
  S_2 = S_2,
  lambda = lambda,
  mu_0 = mu_0,
  Sigma_0 = Sigma_0,
  H_beta = H_beta,
  h_beta = h_beta,
  sigma_sq_a = sigma_sq_a,
  sigma_sq_b = sigma_sq_b,
  tau_sq_eta_a = tau_sq_eta_a,
  tau_sq_eta_b = tau_sq_eta_b,
  tau_sq_w_a = tau_sq_w_a,
  tau_sq_w_b = tau_sq_w_b,
  rho_eta_c = rho_eta_c,
  rho_eta_d = rho_eta_d,
  rho_eta_tuning = rho_eta_tuning,
  rho_w_c = rho_w_c,
  rho_w_d = rho_w_d,
  rho_w_tuning = rho_w_tuning,
  mu_s = mu.s,
  beta_s = beta.s,
  eta_s = eta.s,
  u_s = u.s,
  sigma_sq_s = sigma_sq.s,
  Sigma_beta_s = Sigma_beta.s,
  beta_0_s = beta_0.s,
  tau_sq_eta_s = tau_sq_eta.s,
  rho_eta_s = rho_eta.s,
  tau_sq_w_s = tau_sq_w.s,
  rho_w_s = rho_w.s,
  mu_samples = mu.samples,
  beta_samples = beta.samples,
  eta_samples = eta.samples,
  u_samples = u.samples,
  sigma_sq_samples = sigma_sq.samples,
  Sigma_beta_samples = Sigma_beta.samples,
  beta_0_samples = beta_0.samples,
  tau_sq_eta_samples = tau_sq_eta.samples,
  rho_eta_samples = rho_eta.samples,
  tau_sq_w_samples = tau_sq_w.samples,
  rho_w_samples = rho_w.samples
)

#######################
### Collect samples ###
#######################

samples <- list(
  "mu.samples" = mu.samples,
  "beta.samples" = beta.samples,
  "eta.samples" = eta.samples,
  "u.samples" = u.samples,
  "sigma_sq.samples" = sigma_sq.samples,
  "Sigma_beta.samples" = Sigma_beta.samples,
  "beta_0.samples" = beta_0.samples,
  "tau_sq_eta.samples" = tau_sq_eta.samples,
  "rho_eta.samples" = rho_eta.samples,
  "tau_sq_w.samples" = tau_sq_w.samples,
  "rho_w.samples" = rho_w.samples
)

####################
### Save samples ###
####################

saveRDS(samples, "samples.rds")
