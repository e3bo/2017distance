#!/usr/bin/Rscript

source("example-distance-est-core.R")

set.seed(2)

linear_sir_stepfun <- pomp::Csnippet("
  double dw1 = rnorm(0, sqrt(dt));
  double dw2 = rnorm(0, sqrt(dt));
  double Z1new = f11 * Z1 * dt + f12 * Z2 * dt + b11 * dw1 + b12 * dw2 + Z1;
  double Z2new = f21 * Z1 * dt + f22 * Z2 * dt + b12 * dw1 + b22 * dw2 + Z2;
  cases += Z2 * gamma * dt;
  Z1 = Z1new;
  Z2 = Z2new;
")

sir_initializer <- function(params, t0, ...) {
  statenames <- c("Z1", "Z2", "cases")
  x0 <- c(params[c("Z1_0", "Z2_0")], 0)
  names(x0) <- statenames
  x0
}

create_sir_sim <- function(times = seq(0, 9),
                                t0 = min(times),
                                params, delta.t = 0.01){
  statenames <- c("Z1", "Z2", "cases")
  rp_names <- c("f11", "f12", "f21", "f22", "b11", "b12", "b22", "gamma")
  ivp_names <- c("Z1_0", "Z2_0")
  paramnames <- c(rp_names, ivp_names)
  data <- data.frame(time = times, reports = NA)
  simp <- pomp::pomp(data = data, times = "time", t0 = t0, params = params,
                     rprocess = pomp::euler.sim(step.fun = linear_sir_stepfun,
                     delta.t = delta.t),
                     #zeronames = "cases",
                     obsnames = "reports",
                     statenames = statenames,
                     paramnames = paramnames,
                     initializer = sir_initializer)
  simp
}

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
                             fnoise = 0, p = 0){

  params <- c(b_0 = NA, d = 0.02, sigma_env_d = dnoise, sigma_env_f = fnoise,
              gamma = 365 / 22, eta = eta, mu = N_0 * 0.02, N_0 = N_0, p = p,
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

pvac_seq <- seq(0, 1, by = 0.01)

Fl <- lapply(pvac_seq, function(x) get_drift_mat(p = x))

Dl <- list(diag(2))

mat2pars <- function(F, D, opars = c(gamma = 1, Z1_0 = 10, Z2_0 = 0)){
  B <- sqrtM22(D)
  mpars <- c(f11 = F[1,1], f12 = F[1,2], f21 = F[2,1], f22 = F[2,2], b11 = B[1,1], b12 = B[1,2], b22 = B[2,2])
  c(mpars, opars)
}

#Fl <- list(matrix(c(-1, 0, 0, -.3), 2, 2))
pars <- Map(mat2pars, Fl, Dl)

sim <- create_sir_sim(params = pars[[1]], delta.t = 1e-3)

tseq <- seq(0, 10, length.out=100)
out <- simulate(sim, states = TRUE, times = tseq, nsim = 1e4)

sigman_seq <- apply(out, 3, function(x) var(t(x)))
ind <- seq(1, ncol(sigman_seq))
sigma_num_list <- lapply(ind, function(x) matrix(sigman_seq[,x], 3, 3))

calc_Sigma <- function(Fmat, Dmat, t = Inf) {
    T <- eigen(Fmat)$vectors
    Lambda <- eigen(Fmat)$values

    Dt <- solve(T) %*% Dmat %*% t(solve(T))
    denom <- outer(Lambda, Lambda, "+")
    Sigmat <- -Dt / denom + diag(exp(Lambda * t)) %*% Dt %*% diag(exp(Lambda * t)) / denom
    Re(T %*% Sigmat %*% t(T))
    }


tfine <- seq(0, max(tseq), len = 100)
Sigmal <- Map(calc_Sigma, Fl[1], Dl[1], t = tfine)

S11 <- sapply(Sigmal, "[", 1, 1)
S12 <- sapply(Sigmal, "[", 1, 2)
S22 <- sapply(Sigmal, "[", 2, 2)


par(mfrow = c(3, 1))
plot(tfine, S11, type = 'l')
points(tseq, sigman_seq[1,])
plot(tfine, S12, type = 'l')
points(tseq, sigman_seq[2,])
plot(tfine, S22, type = 'l')
points(tseq, sigman_seq[5,])


## check solution for mean of cases

plot(tseq, colMeans(out["cases", ,]))
plot(tseq[-1], diff(colMeans(out["cases", ,])))
plot(tseq, colMeans(out["Z2", ,]))

ExpectedVals <- function(Z1_0, Z2_0, t, F){
  Lambda <- eigen(F)$values
  T <- eigen(F)$vectors
  Re(T %*% diag(exp(Lambda * t)) %*% solve(T) %*% c(Z1_0, Z2_0))
}

EZ2 <- sapply(tseq, function(x) ExpectedVals(10, 0, t = x, Fl[[1]]))[2,]
lines(tseq, EZ2)

ExpectedCases <- function(Z1_0, Z2_0, t, F){
  Lambda <- eigen(F)$values
  T <- eigen(F)$vectors
  Re(T %*% diag((exp(Lambda * t) - 1) / Lambda) %*% solve(T) %*% c(Z1_0, Z2_0))[2]
}

par(mfrow = c(1,1))
Ecases <- sapply(tseq, function(x) ExpectedCases(10, 0, t = x, Fl[[1]]))
plot(tseq, Ecases, type = 'l')
points(tseq, colMeans(out["cases", ,]))
