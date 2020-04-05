# SimulateMethod Class and derivatives ----

SimulateMethod <- R6Class(
    classname = 'SimulateMethod',
    private = list(),
    public = list(
        simulate = function(model, terminal, dt)
        {
            return(NA)
        })
)

VasicekTransMat <- R6Class(
    classname = 'VasicekTransMat',
    inherit = SimulateMethod,
    private = list(),
    public = list(
        initialize = function(){},
        simulate = function(model, terminal, dt)
        {
            kappa <- model$kappa()
            theta <- model$theta()
            sigma <- model$sigma()
            factor_num <- model$factorNum()

            n <- terminal / dt
            r <- matrix(
                theta,
                nrow = n,
                ncol = factor_num,
                byrow = TRUE)

            for (i in 2:n)
            {
                C = theta*(1 - exp(-kappa * dt))
                F_ = exp(-kappa*dt)
                Q = sqrt(sigma^2/(2*kappa)*(1 - exp(-2*kappa*dt)))
                e = rnorm(factor_num, mean = 0, sd = Q)
                r[i,] = C + F_ * r[i - 1,] + e
            }

            return(r)
        })
)

CIRTransMat <- R6Class(
    classname = 'CIRTransMat',
    inherit = SimulateMethod,
    private = list(),
    public = list(
        initialize = function(){},
        simulate = function(model, terminal, dt)
        {
            kappa <- model$kappa()
            theta <- model$theta()
            sigma <- model$sigma()
            factor_num <- model$factorNum()

            n <- terminal / dt
            r <- matrix(
                theta,
                nrow = n,
                ncol = factor_num,
                byrow = TRUE)

            for (i in 2:n)
            {
                C = theta*(1 - exp(-kappa * dt))
                F_ = exp(-kappa*dt)
                Q = theta*sigma^2/(2*kappa)*(1 - exp(-kappa * dt))^2 + sigma^2/kappa * (exp(-kappa * dt) - exp(-2*kappa * dt))*r[i - 1,]
                e = rnorm(factor_num, mean = 0, sd = sqrt(Q))
                r[i,] = C + F_ * r[i - 1,] + e
            }

            return(r)
        })
)

Euler <- R6Class(
    classname = 'Euler',
    inherit = SimulateMethod,
    private = list(),
    public = list(
        initialize = function(){},
        simulate = function(model, terminal, dt)
        {
            init <- model$init()
            factor_num <- model$factorNum()

            n <- terminal / dt
            r <- matrix(
                init,
                nrow = n,
                ncol = factor_num,
                byrow = TRUE)

            for (i in 2:n)
            {
                dw = sqrt(dt) * rnorm(factor_num)
                r[i,] = r[i - 1,] + model$drift(r[i - 1,]) * dt  + model$diffusion(r[i - 1,]) * dw
            }

            return(r)
        })
)

Milstein <- R6Class(
    classname = 'Milstein',
    inherit = SimulateMethod,
    private = list(),
    public = list(
        initialize = function(){},
        simulate = function(model, terminal, dt)
        {
            init <- model$init()
            factor_num <- model$factorNum()

            n <- terminal / dt
            r <- matrix(
                init,
                nrow = n,
                ncol = factor_num,
                byrow = TRUE)

            for (i in 2:n)
            {
                dw = sqrt(dt) * rnorm(factor_num)
                r[i,] = r[i - 1,] +
                    model$drift(r[i - 1,]) * dt  +
                    model$diffusion(r[i - 1,]) * dw +
                    0.5*model$diffusionDerivative(r[i - 1,])*model$diffusion(r[i - 1,])*(dw^2 - dt)
            }

            return(r)
        })
)
