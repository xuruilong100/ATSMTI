library(R6)
library(dplyr)
library(FKF)

source('SimulateMethod.R')
source('AffineModel.R')
source('KalmanFilterAffineModel.R')

# param for Vasicek ----

L <- 10
dt <- 1/12
measure_error <- 0.001

# 1 FACTOR
# factor_num <- 1
# kappa     <- c(0.06)
# theta     <- c(0.05)
# sigma     <- c(0.02)
# lambda    <- c(-0.20)

# 2 FACTOR
# factor_num <- 2
# kappa     <- c(0.06, 0.7)
# theta     <- c(0.05, 0.01)
# sigma     <- c(0.02, 0.05)
# lambda    <- c(-0.20, -0.50)

# 3 FACTOR
factor_num <- 3
kappa     <- c(0.06, 0.3, 0.7)
theta     <- c(0.01, 0.02, 0.04)
sigma     <- c(0.02, 0.05, 0.03)
lambda    <- c(-0.20, -0.50, -0.15)

tau  <- c(1/12 ,1/4, 1/2, 1, 2, 3, 4, 5, 7, 10, 12, 15, 20, 30)

transMat <- VasicekTransMat$new()
euler <- Euler$new()
milstein <- Milstein$new()

model <- Vasicek$new(
    kappa, theta, sigma, lambda, milstein)

set.seed(123)

idx <- seq(0, L - dt, dt)
zeros <- model$zeroRate(10, dt, tau, measure_error)

plot(x = idx, y = zeros[,1], type = 'l', ylim = c(0,0.2))
for (i in 2:length(tau))
{
    lines(x = idx, y = zeros[, i], col = palette(rainbow(length(tau)))[i])
}
lines(x = idx, y = zeros[,1], type = 'l', lwd = 2)

plot(x = tau, y = zeros[100,], type = 'o', lwd = 2)

# Test Kalman Filter

kappa_init <- runif(factor_num, min = 0.0, max = 0.8)
theta_init <- runif(factor_num, min = 0.0, max = 0.1)
sigma_init <- runif(factor_num, min = 0.0, max = 0.1)
lambda_init <- runif(factor_num, min = -1.0, max = 0.0)
measurement_err_init <- runif(length(tau), min = 0.0, max = 0.1)

initial_param <- c(
    kappa_init,
    theta_init,
    sigma_init,
    lambda_init,
    measurement_err_init)

optim_controls <- list(
    trace = 100, eval.max = 3000, iter.max = 6000, rel.tol = 1e-5)
upper_bound <- c(
    rep(c(1.0, 1.0, 1.0, 1.0), each = factor_num),
    rep(0.1, length(tau)))
lower_bound <- c(
    rep(c(0.0001, 0.0001, 0.0001, -1.0 ), each = factor_num),
    rep(0.0001, length(tau)))

helper <- VasicekHelper$new(factor_num)

vkf <- KalmanFilterVasicek$new(
    helper, dt, tau)

rst <- vkf$kalmanFilterEstimate(
    zeros, initial_param, optim_controls,
    lower_bound, upper_bound)

rst$convergence # rst$convergence == 0 means successful convergence.
rst$par[1:(4*factor_num)] %>% matrix(ncol = factor_num, byrow = TRUE)
rst$par[-1:(-4*factor_num)]
