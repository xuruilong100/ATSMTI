# Affine Class and derivatives ----

AffineModel <- R6Class(
    classname = 'AffineModel',
    private = list(
        kappa_ = NA,
        theta_ = NA,
        sigma_ = NA,
        lambda_ = NA,
        factor_num_ = NA,
        method_ = NA),
    public = list(
        initialize = function(kappa, theta,
                              sigma, lambda,
                              method)
        {
            private$kappa_ = kappa
            private$theta_ = theta
            private$sigma_ = sigma
            private$lambda_ = lambda
            private$factor_num_ = length(kappa)
            private$method_ = method
        },
        kappa = function(){return(private$kappa_)},
        theta = function(){return(private$theta_)},
        sigma = function(){return(private$sigma_)},
        lambda = function(){return(private$lambda_)},
        factorNum = function(){return(private$factor_num_)},
        drift = function(x){return(NA)},
        diffusion = function(x){return(NA)},
        diffusionDerivative = function(x){return(NA)},
        init = function(x){return(NA)},
        A = function(tau){return(NA)},
        B = function(tau){return(NA)},
        factor = function(terminal, dt)
        {
            f <- private$method_$simulate(
                self, terminal, dt)
            return(f)
        },
        shortRate = function(terminal, dt)
        {
            n <- terminal / dt
            sr <- self$factor(terminal, dt)
            if (private$factor_num_ == 1)
            {
                return(sr)
            }
            else
            {
                sr %<>% apply(1, sum) %>% matrix(nrow = n)
                return(sr)
            }
        },
        zeroRate = function(terminal, dt, tau,
                            measurement_error = NULL)
        {
            n <- terminal / dt

            f <- self$factor(
                terminal, dt)

            m = length(tau)

            A <- -self$A(tau) / tau
            H <- self$B(tau) / matrix(tau, nrow = m, ncol = private$factor_num_)

            if (is.null(measurement_error))
            {
                nu <- 0
            }
            else
            {
                nu <- rnorm(
                    m*n, mean = 0, sd = measurement_error) %>%
                    matrix(ncol = n)
            }

            z = matrix(A, nrow = m, ncol = n) + H %*% t(f) + nu

            return(t(z))
        })
)

Vasicek <- R6Class(
    classname = 'Vasicek',
    inherit = AffineModel,
    private = list(),
    public = list(
        A = function(tau)
        {
            B <- self$B(tau)
            gamma = private$kappa_^2*(private$theta_ - (private$sigma_*private$lambda_)/private$kappa_) - private$sigma_^2/2

            A <- matrix(
                0.0,
                nrow = length(tau),
                ncol = 1)
            for (i in 1:ncol(B))
            {
                A <- A + gamma[i] * (B[,i] - tau)/(private$kappa_[i]^2) - private$sigma_[i]^2*B[,i]^2/(4*private$kappa_[i])
            }
            return(A)
        },
        B = function(tau)
        {
            B <- matrix(
                nrow = length(tau),
                ncol = private$factor_num_)
            for (i in 1:ncol(B))
            {
                B[,i] <- 1/private$kappa_[i]*(1 - exp(-private$kappa_[i] * tau))
            }
            return(B)
        },
        drift = function(x)
        {
            return(private$kappa_*(private$theta_ - x))
        },
        diffusion = function(x)
        {
            return(private$sigma_)
        },
        diffusionDerivative = function(x)
        {
            return(0)
        },
        init = function(x)
        {
            return(private$theta_)
        })
)

CIR <- R6Class(
    classname = 'CIR',
    inherit = AffineModel,
    private = list(),
    public = list(
        A = function(tau)
        {
            kappa <- private$kappa_
            theta <- private$theta_
            sigma <- private$sigma_
            lambda <- private$lambda_

            gamma <- sqrt((kappa + lambda)^2 + 2*sigma^2)

            A <- matrix(
                0.0,
                nrow = length(tau),
                ncol = 1)

            for (i in 1:private$factor_num_)
            {
                numerator <- 2*gamma[i]*exp((gamma[i] + kappa[i] + lambda[i])*tau/2)
                denominator <- (gamma[i] + kappa[i] + lambda[i])*(exp(gamma[i]*tau) - 1) + 2*gamma[i]
                A <- A + log(numerator / denominator) * (2*kappa[i]*theta[i]/sigma[i]^2)
            }
            return(A)
        },
        B = function(tau)
        {
            kappa <- private$kappa_
            sigma <- private$sigma_
            lambda <- private$lambda_

            B <- matrix(
                nrow = length(tau),
                ncol = private$factor_num_)

            gamma <- sqrt((kappa + lambda)^2 + 2*sigma^2)

            for (i in 1:ncol(B))
            {
                numerator <- 2*(exp(gamma[i]*tau) - 1)
                denominator <- (gamma[i] + kappa[i] + lambda[i])*(exp(gamma[i]*tau) - 1) + 2*gamma[i]
                B[,i] <- numerator / denominator
            }
            return(B)
        },
        drift = function(x)
        {
            return(private$kappa_*(private$theta_ - x))
        },
        diffusion = function(x)
        {
            return(private$sigma_ * sqrt(x))
        },
        diffusionDerivative = function(x)
        {
            return(0.5 * private$sigma_ / sqrt(x))
        },
        init = function(x)
        {
            return(private$theta_)
        })
)

Helper <- R6Class(
    classname = 'Helper',
    private = list(
        model_ = NA,
        factor_num_ = NA),
    public = list(
        initialize = function(m)
        {
            private$factor_num_ <- m
        },
        setParam = function(kappa, theta, sigma, lambda) {},
        A = function(tau)
        {
            return(private$model_$A(tau))
        },
        B = function(tau)
        {
            return(private$model_$B(tau))
        },
        factorNum = function()
        {
            return(private$factor_num_)
        }
    )
)

VasicekHelper <- R6Class(
    classname = 'VasicekHelper',
    inherit = Helper,
    private = list(),
    public = list(
        initialize = function(m)
        {
            super$initialize(m)
            private$model_ <- Vasicek$new(
                NA, NA, NA, NA, NA)
        },
        setParam = function(kappa, theta, sigma, lambda)
        {
            private$model_ <- Vasicek$new(
                kappa, theta, sigma, lambda, NA)
        })
)

CIRHelper <- R6Class(
    classname = 'CIRHelper',
    inherit = Helper,
    private = list(),
    public = list(
        initialize = function(m)
        {
            super$initialize(m)
            private$model_ <- CIR$new(
                NA, NA, NA, NA, NA)
        },
        setParam = function(kappa, theta, sigma, lambda)
        {
            private$model_ <- CIR$new(
                kappa, theta, sigma, lambda, NA)
        })
)
