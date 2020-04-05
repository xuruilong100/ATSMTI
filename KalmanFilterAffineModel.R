library(optimx)
# KalmanFilterAffineModel Class and derivatives ----

KalmanFilterAffineModel <- R6Class(
    classname = 'KalmanFilterAffineModel',
    private = list(
        factor_num_ = NA,
        helper_ = NA,
        dt_ = NA,
        tau_ = NA,
        # m be the dimension of the state variable,
        # d be the dimension of the observations,
        # n the number of observations.
        transInterceptMat = function(kappa, theta, sigma, lambda)
        {
            m <- private$factor_num_
            dt <- private$dt_

            # intercept of state transition equation (dt: m x 1)
            C <- matrix(
                theta*(1 - exp(-kappa*dt)),
                nrow = m, ncol = 1)

            return(C)
        },
        transMat = function(kappa, theta, sigma, lambda)
        {
            m <- private$factor_num_
            dt <- private$dt_
            # factor of transition equation (Tt: m x m x 1)
            F_ <- diag(exp(-kappa*dt), m, m)

            return(F_)
        },
        measureInterceptMat = function(kappa, theta, sigma, lambda)
        {
            d <- length(private$tau_)
            tau <- private$tau_
            private$helper_$setParam(
                kappa, theta, sigma, lambda)

            tauMat <- matrix(
                data = tau, nrow = d, ncol = 1, byrow = FALSE)

            # intercept of measurement equation (ct: d x 1)
            A <- -private$helper_$A(tau) / tauMat

            return(A)
        },
        measureMat = function(kappa, theta, sigma, lambda)
        {
            m <- private$factor_num_
            d <- length(private$tau_)
            tau <- private$tau_
            private$helper_$setParam(
                kappa, theta, sigma, lambda)

            tauMat <- matrix(
                data = tau, nrow = d, ncol = m, byrow = FALSE)

            # factor of measurement equation (Zt: d x m x 1)
            H <- private$helper_$B(tau) / tauMat

            return(H)
        },
        nagtiveLoglike = function(zeros){return(NA)}),
    public = list(
        initialize = function(helper, dt, tau)
        {
            private$helper_ <- helper
            private$dt_ <- dt
            private$tau_ <- tau
            private$factor_num_ <- helper$factorNum()
        },
        kalmanFilterEstimate = function(zeros, initialParam,
                                        optimControls,
                                        lowerBound, upperBound)
        {
            loglike <- private$nagtiveLoglike(zeros)

            fitted_model <- nlminb(
                initialParam,
                loglike,
                control = optimControls,
                lower = lowerBound,
                upper = upperBound)

            return(fitted_model)
        }
    )
)

KalmanFilterVasicek <- R6Class(
    classname = 'KalmanFilterVasicek',
    inherit = KalmanFilterAffineModel,
    private = list(
        transVar = function(kappa, theta, sigma, lambda)
        {
            # m be the dimension of the state variable,
            # d be the dimension of the observations,
            # n the number of observations.

            m <- private$factor_num_
            dt <- private$dt_
            # variance of innovations of transition (HHt: m x m x 1)
            Q <- diag(sigma^2/(2*kappa)*(1 - exp(-2*kappa*dt)), m, m)

            return(Q)
        },
        nagtiveLoglike = function(zeros)
        {
            func <- function(x)
            {
                dt <- private$dt_
                tau <- private$tau_

                # m be the dimension of the state variable,
                # d be the dimension of the observations,
                # n the number of observations.

                m <- private$factor_num_
                d <- length(private$tau_) # same as ncol(zeros)
                n <- nrow(zeros)

                kappa <- x[1:m]
                theta <- x[(m + 1):(2*m)]
                sigma <- x[(2*m + 1):(3*m)]
                lambda <- x[(3*m + 1):(4*m)]
                measurement_error <- x[(4*m + 1):length(x)]

                # initial state variable (a0: m x 1)
                r_init <- as.vector(theta)

                # variance of state variable (P0: m x m)
                P_init <- diag(sigma^2/(2*kappa), m, m) # unconditional variance of state variable

                # intercept of state transition equation (dt: m x 1)
                C <- private$transInterceptMat(kappa, theta, sigma, lambda)
                #cat('C:', dim(C), '\n')

                # factor of transition equation (Tt: m x m x 1)
                F_ <- private$transMat(kappa, theta, sigma, lambda)
                #cat('F:', dim(F_), '\n')

                # factor of measurement equation (Zt: d x m x 1)
                H <- private$measureMat(kappa, theta, sigma, lambda)
                #cat('H:', dim(H), '\n')

                # intercept of measurement equation (ct: d x 1)
                A <- private$measureInterceptMat(kappa, theta, sigma, lambda)
                #cat('A:', dim(A), '\n')

                # variance of innovations of transition (HHt: m x m x 1)
                Q <- private$transVar(kappa, theta, sigma, lambda)
                #cat('Q:', dim(Q), '\n')

                # variance of measurement error (GGt: d x d x 1)
                R <- diag(measurement_error^2, nrow = d, ncol = d)
                #cat('R:', dim(R), '\n')

                filtered_process <- fkf(
                    a0 = r_init,
                    P0 = P_init,
                    dt = C,
                    ct = A,
                    Tt = F_,
                    Zt = H,
                    HHt = Q,
                    GGt = R,
                    yt = t(zeros))

                return(-filtered_process$logLik)
            }

            return(func)
        }),
    public = list()
)

KalmanFilterCIR <- R6Class(
    classname = 'KalmanFilterCIR',
    inherit = KalmanFilterAffineModel,
    private = list(
        transVar = function(kappa, theta, sigma, lambda, y)
        {
            # m be the dimension of the state variable,
            # d be the dimension of the observations,
            # n the number of observations.

            m <- private$factor_num_
            dt <- private$dt_
            # variance of innovations of transition (HHt: m x m x 1)
            Q <- diag(
                theta*sigma^2/(2*kappa)*(1 - exp(-kappa*dt))^2 + sigma^2/kappa*(exp(-kappa*dt) - exp(-2*kappa*dt))*y,
                m, m)

            return(Q)
        },
        nagtiveLoglike = function(zeros)
        {
            func <- function(x)
            {
                dt <- private$dt_
                tau <- private$tau_

                # m be the dimension of the state variable,
                # d be the dimension of the observations,
                # n the number of observations.

                m <- private$factor_num_
                d <- length(private$tau_) # same as ncol(zeros)
                n <- nrow(zeros)

                kappa <- x[1:m]
                theta <- x[(m + 1):(2*m)]
                sigma <- x[(2*m + 1):(3*m)]
                lambda <- x[(3*m + 1):(4*m)]
                measurement_error <- x[(4*m + 1):length(x)]

                loglike <- 0.0

                A <- private$measureInterceptMat(kappa, theta, sigma, lambda)
                H <- private$measureMat(kappa, theta, sigma, lambda)
                I <- diag(1, nrow = m, ncol = m)
                C <- private$transInterceptMat(kappa, theta, sigma, lambda)
                F_ <- private$transMat(kappa, theta, sigma, lambda)

                for (i in 1:n)
                {
                    if (i == 1)
                    {
                        Ey <- matrix(theta, nrow = m, ncol = 1)
                        vary <- diag(sigma^2*theta/(2*kappa), nrow = m, ncol = m)
                    }

                    Q <- private$transVar(kappa, theta, sigma, lambda, as.vector(Ey))

                    Ez <- A + H %*% Ey
                    varz <- H %*% vary %*% t(H) + diag(measurement_error^2, nrow = d, ncol = d)

                    zeta <- matrix(zeros[i,], ncol = 1) - Ez

                    loglike <- loglike + log(det(varz)) + t(zeta) %*% solve(varz) %*% zeta

                    # update part

                    K <- vary %*% t(H) %*% solve(varz)
                    Ey_filter <- Ey + K %*% zeta

                    Ey_filter[Ey_filter < 0, ] <- 0

                    vary_filter <- (I - K %*% H) %*% vary

                    Ey <- C + F_ %*% Ey_filter
                    vary <- vary - F_ %*% vary_filter %*% t(F_) + Q
                }

                return(loglike)
            }
            return(func)
        }),
    public = list()
)
