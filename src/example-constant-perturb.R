set.seed(1)
source("example-distance-est-core.R")

EndemicEquilSIR <- function(beta=(R0 * (mu + gamma)), eta=17/5e4,
                            gamma=365/22, mu=1/50, p=0,  R0=17, verbose=FALSE) {
  ## Computes the endemic equilibrium of an SIR model with immigration, Eq. 6.
  ##
  ## Args:
  ##   beta: numeric. The transmission rate.
  ##   eta: numeric. The rate of infection from outside.
  ##   gamma: numeric. The recovery rate.
  ##   mu: numeric. The birth rate.
  ##   p: numeric. The vaccination uptake.
  ##   R0: numeric. The basic reproduction number.
  ##
  ## Returns:
  ##   A list with numeric elements S, I, and R, coresponding to the
  ##   equilibrium fractions of the population in the
  ##   susceptible, infected, and removed states.
  stopifnot(c(beta, eta, gamma, p, R0) >= 0)
  stopifnot(p <= 1)
  a <- - beta * (gamma + mu)
  b <- beta * mu * (1 - p) - (gamma + mu) * (eta + mu)
  c <- mu * (1 - p) * eta
  eq <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
  i.star <- ifelse(p == 1, 0, eq)
  s.star <- ifelse(p == 1, 0, mu * (1 - p)/ (beta * i.star + eta + mu))
  if (verbose) {
    ds.star <- mu *(1 - p) - beta * s.star * i.star - eta * s.star - mu * s.star
    di.star <- beta * s.star * i.star + eta * s.star - (gamma + mu) * i.star
    cat('dS = ', ds.star, '\n')
    cat('dI = ', di.star, '\n')
  }
  return(list(S=s.star, I=i.star, R=1 - i.star - s.star))
}



get_drift_mat <- function(R0 = 17, N_0 = 2e6, eta = 2e-8, tstep = 1 / 52,
                             tlength = 1000, burninyrs = 10, dnoise = 0,
                             fnoise = 0, savemem = FALSE){

  params <- c(b_0 = NA, d = 0.02, sigma_env_d = dnoise, sigma_env_f = fnoise,
              gamma = 365 / 22, eta = eta, mu = N_0 * 0.02, N_0 = N_0, p = 0,
              dbdt = 0)
  params["b_0"] <- with(as.list(params), (gamma + d) * R0)

  equil <- with(as.list(params), EndemicEquilSIR(beta = b_0, eta = eta,
                                           gamma = gamma, mu = d, p = p))

  params["X1_0"] <- equil$S
  params["X2_0"] <- equil$I
  params["X3_0"] <- equil$R

  drift <- with(as.list(params),
                get_drift(beta = b_0, eta = eta, d = d, gamma = gamma, S = X1_0,
                          I = X2_0, sigma_env_f = sigma_env_f,
                          sigma_env_d = sigma_env_d))
  drift
}

repnum_seq <- 2^seq(-4, 5, length.out = 100)

F <- lapply(repnum_seq, get_drift_mat)

D <- diag(2)

calc_Sigma <- function(Fmat) {
    T <- eigen(Fmat)$vectors
    Lambda <- eigen(Fmat)$values
    
    Dt <- solve(T) %*% D %*% t(solve(T))
    denom <- outer(Lambda, Lambda, "+")
    Sigmat <- -Dt / denom
    T %*% Sigmat %*% t(T)
    }

Sigma <- lapply(F, calc_Sigma)

S11 <- Re(sapply(Sigma, "[", 1, 1))
S22 <- Re(sapply(Sigma, "[", 2, 2))
gvar <- sapply(Sigma, function(x) det(Re(x)))

plot(repnum_seq, S11, log = 'xy', type = 'l', xlab = expression(R[0]), ylab = "Var(S)")
png('varI-constant-perturb.png')
plot(repnum_seq, S22, log = 'xy', type = 'l', xlab = expression(R[0]), ylab = "Var(I)")
dev.off()
plot(repnum_seq, gvar, log = 'xy', type = 'l')


P <- function(t, p.s, p.0) {
  ## Computes the vaccination level as function of time.
  val <- p.s + p.0 * t
  val <- ifelse(val > 1, 1, val)
  ifelse(val < 0, 0, val)
}

RunSDESIR <- function(beta=(R0 * (gamma + mu)), eta=2e-8,
                      gamma=(365 / 22), imports=100,  mu=(1 / 50), p.s=0,
                      p.0=0, pop.size=5e6, R0=17,
                      n.intervals=1e4, start.time=0, stop.time=1000) {
  ## Solves for the trajectory of the stochastic fast-slow SIS
  ## model corresponding to Eq. 7 using the Euler method.
  ##
  ## Args:
  ##   beta: numeric. The transmission rate.
  ##   eta: numeric. The rate of infection from outside.
  ##   gamma: numeric. The recovery rate.
  ##   mu: numeric. The birth rate.
  ##   imports: numeric. The expected number of imported cases.
  ##   p.s: numeric. The vaccination uptake at time 0.
  ##   p.0: numeric. The rate of change in vaccination uptake per unit time.
  ##   pop.size: numeric. The population size.
  ##   R0: numeric. The basic reproduction number (with no vaccination).
  ##   init.vars: numeric of initial number of susceptible and infected
  ##     individuals with names 'X' and 'Y'.
  ##   n.intervals: integer. Number of intervals in returned time series.
  ##   start.time: numeric. Time at which to begin solving.
  ##   stop.time: numeric. Time at which to stop solving.
  ##
  ## Returns:
  ##   numeric matrix. The sampling times and numbers of susceptible,
  ##   infected, and removed individuals.

  R0 <- beta / (gamma + mu) # In case user provides beta instead of R0
  stopifnot(c(beta, eta, gamma, imports, p.s, pop.size, R0) >= 0)
  stopifnot(p.s <= 1)
  stopifnot(stop.time > start.time)

  equil <- EndemicEquilSIR(beta = beta, eta = eta, gamma = gamma, mu = mu, p = 0)

  init.vars <- c(X = equil$S * pop.size, Y = equil$I * pop.size)
    
  sqrtM22 <- function(M){
    ## Computes the square root of 2 by 2 matrix
    A <- M[1, 1]
    B <- M[1, 2]
    C <- M[2, 1]
    D <- M[2, 2]
    tau <- A + D
    delta <- A * D - B * C
    s <- sqrt(delta)
    t <- sqrt(tau + 2 * s)
    stopifnot(t != 0)
    rbind(c(A + s, B), c(C, D + s)) / t
  }

  X <- numeric(n.intervals + 1)
  Y <- numeric(n.intervals + 1)
  X[1] <- init.vars['X']
  Y[1] <- init.vars['Y']
  time.range <- stop.time - start.time
  dt <- time.range / n.intervals
  delta.w1 <- rnorm(n.intervals, 0, sqrt(dt))
  delta.w2 <- rnorm(n.intervals, 0, sqrt(dt))
  time <- seq(start.time, stop.time, dt)
  p.t <- P(time, p.s, p.0)
  KeepInBounds <- function(x) {
    if (x > pop.size) {
      x <- pop.size
    } else if (x < 0) {
      x <- 0
    }
    x
  }
  for (i in 1:n.intervals) {
    rate.trans <- (beta * Y[i] / pop.size + eta ) * X[i]
    EdX <- -rate.trans + mu * (pop.size * (1 - p.t[i]) - X[i])
    EdY <- rate.trans - (gamma + mu) * Y[i]
    B <- sqrtM22(D)
    stopifnot(isTRUE(all.equal(B %*% B, D, check.attributes = FALSE)))

    X[i + 1] <- X[i] + EdX * dt + B[1, 1] * delta.w1[i] + B[1, 2] * delta.w2[i]
    Y[i + 1] <- Y[i] + EdY * dt + B[2, 1] * delta.w1[i] + B[2, 2] * delta.w2[i]
    X[i + 1] <- KeepInBounds(X[i + 1])
    Y[i + 1] <- KeepInBounds(Y[i + 1])
    Z <- pop.size - X[i + 1] - Y[i + 1]
    if (Z < 0){
      ## keep Z positive
      X[i + 1] <- X[i + 1] * (pop.size * (1 - 1e-6) / (pop.size - Z))
      Y[i + 1] <- Y[i + 1] * (pop.size * (1 - 1e-6) / (pop.size - Z))
    }
  }
  cbind(time=time, X=X, Y=Y, Z=(pop.size - X - Y))
}

get_numerical_var <- function(repnum){
    var(RunSDESIR(R0=repnum, n.intervals = 5e4)[-seq(1,1e4),-1])
}

sigma_num <- lapply(repnum_seq, get_numerical_var)

S11n <- Re(sapply(sigma_num, "[", 1, 1))
S22n <- Re(sapply(sigma_num, "[", 2, 2))
gvarn <- sapply(sigma_num, function(x) det(Re(x[seq(1,2), seq(1,2)])))

plot(S11n, S11, log='xy')
plot(S22n, S22, log='xy')
plot(gvarn, gvar, log = 'xy')

plot(repnum_seq, S11n, log='xy')
lines(repnum_seq, S11)

plot(repnum_seq, S22n, log='xy')
lines(repnum_seq, S22)

plot(repnum_seq, gvarn, log='xy')
lines(repnum_seq, gvar)
