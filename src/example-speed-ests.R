#!/usr/bin/Rscript

source("../src/example-distance-est-core.R")

get_drift_eq_at_t <- function(params, tpoint){
    beta <- params["b_0"] + params["dbdt"] * tpoint
    equil <- with(as.list(params), get_equil(beta = beta, eta = eta,
                                             gamma = gamma, d = d, p = p,
                                             sigma_env_d = sigma_env_d,
                                             sigma_env_f = sigma_env_f))
    params["X1_0"] <- equil$sstar
    params["X2_0"] <- equil$istar
    params["X3_0"] <- 1 - equil$sstar - equil$istar
    drift <- with(as.list(params),
                  get_drift(beta = beta, eta = eta, d = d, gamma = gamma, S = X1_0,
                            I = X2_0, sigma_env_f = sigma_env_f,
                            sigma_env_d = sigma_env_d))
    lambda <- eigen(drift)$values
    list(equil = equil, drift = drift, lambda = lambda)
}

get_speed_est <- function(R0_0 = 0, N_0 = 1e7, eta = 400e-7, tstep = 1 / 52,
                          tlength = 520 * 4, burninyrs = 10, dnoise = 0,
                          fnoise = 0, dbdt = 0, lag.max = 52){

  params <- c(b_0 = NA, d = 0.02, sigma_env_d = dnoise, sigma_env_f = fnoise,
              gamma = 365 / 22, eta = eta, mu = N_0 * 0.02, N_0 = N_0, p = 0,
              dbdt = dbdt)
  params["b_0"] <- with(as.list(params), (gamma + d) * R0_0)

  times <- burninyrs + seq(from = 0, length.out = tlength, by = tstep)
  del <- lapply(times, get_drift_eq_at_t, params = params)
  del0 <- get_drift_eq_at_t(params = params, tpoint = 0)

  params["X1_0"] <- del0$equil$sstar
  params["X2_0"] <- del0$equil$istar
  params["X3_0"] <- 1 - del0$equil$sstar - del0$equil$istar

  sim <- create_sir_sim(t0 = 0, params = params, times = times)

  out <- pomp::simulate(sim, states = TRUE, params = params)
  outI <- ts(out["X2", 1, ], deltat = tstep, start = burninyrs)

  cut <- (max(time(outI)) - min(time(outI))) / 2  + min(time(outI))
  w1 <- window(outI, end = cut)
  w2 <- window(outI, start = cut)

  acestw1 <- acf(w1, lag.max = 1040 - 30, plot = FALSE)$acf[, 1, 1]
  acestw2 <- acf(w2, lag.max = 1040 - 30, plot = FALSE)$acf[, 1, 1]

  if(any(c(is.na(acestw1), is.na(acestw2)))) {
    fitw1 <- fitw2 <- list(coef=c(omega = NA, gamma = NA, a = NA))
  } else {
    fitw1 <- get_fit(acestw1, tstep = tstep)
    fitw2 <- get_fit(acestw2, tstep = tstep)
  }
  ret <- list(estsw1 = fitw1$coef, estsw2 = fitw2$coef,
              del = del, del0 = del0, acestw1 = acestw1, acestw2 = acestw2,
              fitw1 = fitw1$fit, fitw2 = fitw2$fit,
              simts = out, altfitw1 = fitw1$altfit,
              altfitw2 = fitw2$altfit, w1 = w1, w2 = w2, outI = outI)
  ret
}

simulate_ests <- function(R0_0 = 0, dbdt = 8 / 50){
    nreplicates <- 1000
    tmpf <- function(...) get_speed_est(dbdt = dbdt, R0_0 = R0_0)
    parallel::mclapply(integer(nreplicates), tmpf)
}

system.time(se <- simulate_ests())

save.image("example-speed-ests.RData")
