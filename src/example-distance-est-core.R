#!/usr/bin/Rscript

set.seed(2)
options(mc.cores = 20)

sir_stepfun <- pomp::Csnippet("
  double dN01, dN03;
  double dN10, dN12;
  double dN20, dN23;
  double dN30;

  double mu10, mu12;
  double mu20, mu23;
  double mu30;

  double d_perturb, mu12_perturb;

  double rate[5], trans[5], dW;

  dW = rgammawn(sigma_env_d, dt);
  d_perturb = d * dW / dt;

  mu10 = d_perturb;
  mu20 = d_perturb;
  mu12 = b * X2 / N_0 + eta;
  mu23 = gamma;
  mu30 = d_perturb;

  dW = rgammawn(sigma_env_f, dt);
  mu12_perturb = mu12 * dW / dt;

  rate[0] = mu10;
  rate[1] = mu12_perturb;
  rate[2] = mu20;
  rate[3] = mu23;
  rate[4] = mu30;
  reulermultinom(2, X1, &rate[0], dt, &trans[0]);
  reulermultinom(2, X2, &rate[2], dt, &trans[2]);
  reulermultinom(1, X3, &rate[4], dt, &trans[4]);

  dN01 = rpois(mu * (1 - p) * dt);
  dN10 = trans[0];
  dN12 = trans[1];
  dN20 = trans[2];
  dN23 = trans[3];
  dN30 = trans[4];
  dN03 = rpois(mu * p * dt);

  X1 = X1 - dN12 - dN10 + dN01;
  X2 = X2 + dN12 - dN20        - dN23;
  X3 = X3        - dN30 + dN03 + dN23;
  cases += dN23;

  b += dbdt * dt;
  if (b < 0) b = 0;
")

sir_vectorfield <- pomp::Csnippet("
  double lambda;
  double tau;
  double dmean;
  double lambdamean;

  lambda = b * X2 / N_0 + eta;
  if (sigma_env_f > 0) {
    tau = pow(sigma_env_f, 2);
    lambdamean = log(1 + tau * lambda) / tau;
  } else {
    lambdamean = lambda;
  }
  if (sigma_env_d > 0) {
    tau = pow(sigma_env_d, 2);
    dmean = log(1 + tau * d) / tau;
  } else {
    dmean = d;
  }

  DX1 =  -lambdamean * X1 - dmean * X1              + mu * (1 - p);
  DX2 =   lambdamean * X1 - dmean * X2 - gamma * X2;
  DX3 =                   - dmean * X3 + gamma * X2 + mu * p;
  Dcases = gamma * X2;
  if (b > 0) {
    Db = dbdt;
  } else {
    Db = 0;
  }
")

sir_initializer <- function(params, t0, ...) {
  statenames <- c("X1", "X2", "X3", "cases", "b")
  x0 <- params[c("X1_0", "X2_0", "X3_0")]
  x0 <- c(round(x0 / sum(x0) * params["N_0"]), 0, params["b_0"])
  names(x0) <- statenames
  x0
}

create_sir_sim <- function(times = seq(0, 9),
                                t0 = min(times),
                                params, delta.t = 0.0001){
  statenames <- c("X1", "X2", "X3", "cases", "b")
  rp_names <- c("d", "sigma_env_d", "sigma_env_f", "gamma", "eta", "mu",
                "p", "dbdt")
  ivp_names <- c("X1_0", "X2_0", "X3_0", "N_0", "b_0")
  paramnames <- c(rp_names, ivp_names)
  data <- data.frame(time = times, reports = NA)
  simp <- pomp::pomp(data = data, times = "time", t0 = t0, params = params,
                     rprocess = pomp::euler.sim(step.fun = sir_stepfun,
                     delta.t = delta.t),
                     skeleton = vectorfield(sir_vectorfield),
                     zeronames = "cases",
                     obsnames = "reports",
                     statenames = statenames,
                     paramnames = paramnames,
                     initializer = sir_initializer)
  simp
}


##

get_l1ac <- function(x, freq){
  ts <- ts(x, freq=freq)
  foo <- acf(ts, lag.max = frequency(ts), plot = FALSE)
  ind <- which(foo$lag[,1,1] == 1.0)
  foo$acf[ind]
}

get_equil <- function(beta, d, gamma, eta, sigma_env_f, sigma_env_d, p){
  if (sigma_env_f > 0) {
      tau_f <- sigma_env_f ^ 2
      f <- function(x) 1 / tau_f * log (1 + tau_f * (beta * x + eta))
  } else {
      f <- function(x) beta * x + eta
  }
  if (sigma_env_d > 0) {
      tau_d <- sigma_env_d ^ 2
      dmean <- 1 / tau_d * log (1 + tau_d * d)
  } else {
      dmean <- d
  }
  obj <- function(x) dmean * (1 - p) - x * (dmean + gamma) -  (dmean + gamma) * dmean * x / f(x)
  ans <- uniroot(f = obj, interval = c(0, 1), tol = 1e-8)
  istar <- ans$root
  sstar <- istar * (dmean + gamma) / f(istar)
  didt <- sstar * f(istar)  -  istar * (dmean + gamma)
  dsdt <- -sstar * f(istar) - dmean * sstar + dmean
  list(ans = ans, istar = istar, sstar = sstar, didt = didt, dsdt = dsdt)
}

get_drift <- function(beta, eta, d, gamma, mu, S, I, sigma_env_f,
                      sigma_env_d){
    ## S and I are densities, not counts
    if (sigma_env_f > 0) {
        tau <- sigma_env_f ^ 2
        lambda <- log(1 + tau * (beta * I + eta)) / tau
        dlambda_dI <- beta / (1 + tau * (beta * I + eta))
    } else {
        lambda <- beta * I + eta
        dlambda_dI <- beta
    }
    if (sigma_env_d > 0) {
        tau <- sigma_env_d ^ 2
        dmean <- log(1 + tau * d) / tau
    } else {
        dmean <- d
    }
    ret <- rbind(c(-lambda - dmean, -dlambda_dI * S),
                 c(lambda, dlambda_dI * S - (gamma + dmean)))
  dimnames(ret) <- list(c('S', 'I'), c('S', 'I'))
  ret
}

get_diffu <- function(beta, eta, d, gamma, mu, S, I, p, N_0, sigma_env_f, sigma_env_d){
  ## S and I are densities, not counts
  lambda <- beta * I + eta
  if (sigma_env_f > 0) {
      tau <- sigma_env_f ^ 2
      prob <- log(1 + tau * lambda) / tau
      m2_trans <- S * prob + S / tau * (S * N_0 - 1) * (2 * log(1 + lambda * tau) - log(1 + 2 * lambda * tau))
  } else {
      m2_trans <- S * lambda
  }
  if (sigma_env_d > 0) {
      tau <- sigma_env_d ^ 2
      dmean <- log(1 + tau * d) / tau
      m2_death_S <- S * dmean + S / tau * (S * N_0 - 1) * (2 * log(1 + d * tau) - log(1 + 2 * d * tau))
      m2_death_I <- I * dmean + I / tau * (I * N_0 - 1) * (2 * log(1 + d * tau) - log(1 + 2 * d * tau))
      m2_death_SI <- S / tau * I * N_0 * (2 * log(1 + d * tau) - log(1 + 2 * d * tau))
  } else {
      m2_death_S <- S * d
      m2_death_I <- I * d
      m2_death_SI <- 0
  }
  ret <- rbind(c(mu * (1 - p) / N_0 + m2_trans + m2_death_S,
                 -m2_trans + m2_death_SI),
               c(-m2_trans + m2_death_SI,
                 m2_trans + gamma * I + m2_death_I))
  dimnames(ret) <- list(c('S', 'I'), c('S', 'I'))
  ret
}

get_fit <- function(y, tstep, est_K = FALSE) {
  x <- seq_along(y) * tstep
  start <- list()
  start$gamma <- unname(coef(lm(log(I(abs(y)))~x))["x"])
  spec <- spectrum(y, plot = FALSE)
  start$omega <- spec$freq[which.max(spec$spec)] / tstep
  start$a <- 0
  fit_osc <- try(minpack.lm::nlsLM(
      y~sqrt(1 + a^2) * exp(x * gamma) * sin(x * omega + atan2(1, a)),
      start = start,
      control = minpack.lm::nls.lm.control(maxiter = 1000)))
  if (est_K) {
      fit_decay <- try(minpack.lm::nlsLM(y~K * exp(x * gamma),
                    start = list(gamma = start$gamma, K = y[1]),
                    control = minpack.lm::nls.lm.control(maxiter = 1000)))
  } else {
      K <- y[1]
      fit_decay <- try(minpack.lm::nlsLM(y~K * exp(x * gamma),
                    start = list(gamma = start$gamma),
                    control = minpack.lm::nls.lm.control(maxiter = 1000)))
  }
  if (inherits(fit_osc , "try-error")) {
      e_osc <- Inf
  } else {
      e_osc <- fit_osc$m$resid()
  }
  if (inherits(fit_decay, "try-error")) {
      e_decay <- Inf
  } else {
      e_decay <- fit_decay$m$resid()
  }
  nll <- function(resids) {
      n <- length(resids)
      log(sum(resids ^ 2)) * n
  }
  aic <- c(constant = nll(y), fit_decay = nll(e_decay) + 2 * (1 + est_K),
           fit_osc = nll(e_osc) + 2 * 3)
  fits <- list(constant = "contant_y=0", fit_decay = fit_decay, fit_osc = fit_osc)

  names(aic)[which.min(aic)]
  coefests <- try(coef(fits[[which.min(aic)]])[c("omega", "gamma", "a")])
  if (inherits(coefests, "try-error")){
      coefests <- c(NA, NA, NA)
  }
  names(coefests) <- c("omega", "gamma", "a")
  c(list(coef = coefests), fits)
}

##

extract_eigen <- function(fit, tstep){
    pars <- coef(fit)
    if (length(pars) == 0) {
        ret <- c(NA, NA, NA)
    } else if ("ar2" %in% names(pars)) {
        if (pars["ar1"] ^ 2 + pars["ar2"] * 4 >= 0) {
            ret <- c(NA, NA, NA)
        } else {
            R <- sqrt(-pars["ar2"])
            gamma <- log(R)  / tstep
            omega <- acos(pars["ar1"]  / (2 * R))
            ret <- c(omega, gamma, NA)
        }
    } else {
        ret <- c(0, log(pars["ar1"]) / tstep, 0)
    }
    names(ret) <- c("omega", "gamma", "a")
    ret
}

get_arma_fit <- function(devs, q = FALSE, tstep){
    orders <- list(c(0, 0, 0), c(1, 0, q), c(2, 0, q))
    tmpf <- function(ord) {
        arima(x = devs, order = ord, include.mean = FALSE, method = "ML")
    }
    fits <- lapply(orders, tmpf)
    scores <- sapply(fits, AIC)
    coefests <- extract_eigen(fits[[which.min(scores)]], tstep)
    list(fits, scores, coefests)
}

get_distance_est <- function(R0 = 17, N_0 = 2e6, eta = 2e-4, tstep = 1 / 52,
                             tlength = 1000, burninyrs = 10, dnoise = 0,
                             fnoise = 0, savemem = FALSE, prob_rep = 1){

  params <- c(b_0 = NA, d = 0.02, sigma_env_d = dnoise, sigma_env_f = fnoise,
              gamma = 365 / 22, eta = eta, mu = N_0 * 0.02, N_0 = N_0, p = 0,
              dbdt = 0)
  params["b_0"] <- with(as.list(params), (gamma + d) * R0)

  equil <- with(as.list(params), get_equil(beta = b_0, eta = eta,
                                           gamma = gamma, d = d, p = p,
                                           sigma_env_d = sigma_env_d,
                                           sigma_env_f = sigma_env_f))
  params["X1_0"] <- equil$sstar
  params["X2_0"] <- equil$istar
  params["X3_0"] <- 1 - equil$sstar - equil$istar

  drift <- with(as.list(params),
                get_drift(beta = b_0, eta = eta, d = d, gamma = gamma, S = X1_0,
                          I = X2_0, sigma_env_f = sigma_env_f,
                          sigma_env_d = sigma_env_d))
  lambda <- eigen(drift)$values

  times <- burninyrs + seq(from = 0, length.out = tlength, by = tstep)
  sim <- create_sir_sim(t0 = 0, params = params, times = times)

  out <- pomp::simulate(sim, states = TRUE, params = params)
  outI <- ts(out["X2", 1, ], deltat = tstep)
  outS <- ts(out["X1", 1, ], deltat = tstep)
  outR <- ts(out["X3", 1, ], deltat = tstep)
  outC <- ts(out["cases", 1, -1], deltat = tstep)
  outC <- sapply(outC, function(x) rbinom(n = 1, size = x, prob = prob_rep))

  lag.max <- dim(out)[3] - 30
  acesti <- acf(outI, lag.max = lag.max, plot = FALSE)$acf[, 1, 1]
  acests <- acf(outS, lag.max = lag.max, plot = FALSE)$acf[, 1, 1]
  acestc <- acf(outC, lag.max = lag.max - 1, plot = FALSE)$acf[-1, 1, 1]

  if(any(is.na(acesti)) | any(is.na(acests))) {
    fitc <- fiti <- fits <- list(coef=c(omega = NA, gamma = NA, a = NA))
  } else {
    fiti <- get_fit(acesti, tstep = tstep, est_K = FALSE)
    fits <- get_fit(acests, tstep = tstep, est_K = FALSE)
    fitc <- get_fit(acestc, tstep = tstep, est_K = TRUE)
  }

    if (FALSE) {
    armafiti <- try(get_arma_fit(devs = outI - mean(outI), tstep = tstep))
    if (!inherits(armafiti, "try-error")){
        arma_ests <- armafiti[[3]]
    } else {
        arma_ests <- c(omega = NA, gamma = NA, a = NA)
    }
    }
    #armafits <- try(get_arma_fit(devs = outS - mean(outS), tstep = tstep))
    
    ret <- list(estsi = fiti$coef, estss = fits$coef, estsc = fitc$coef,
                lambda = lambda, acesti = acesti, acests = acests,
                acestc = acestc, fiti = fiti$fit, fits = fits$fit,
                fitc = fitc$fit, simts = out, altfiti = fiti$altfit,
                altfits = fits$altfit, altffitc = fitc$altfit)
    if (savemem){
        ret$fits <- ret$fiti <- ret$fitc <- ret$simts <- ret$altfiti <- NA
        ret$altfits <- ret$altfitc <- NA
    }
    ret
}

simulate_ests <- function(repnum, N_0, eta, tstep, dnoise, fnoise, prob_rep = 1){
  replicate(100, get_distance_est(tlength = 520 * 2, R0 = repnum, N_0 = N_0,
                                  eta = eta, tstep = tstep, dnoise = dnoise,
                                  fnoise = fnoise, prob_rep = prob_rep),
            simplify = FALSE)
}

makedf <- function(x, repnum, N_0, tstep, dnoise, fnoise, prob_rep = 1){
  pars <- t(sapply(x, "[[", "lambda"))

  estsi <- t(sapply(x, "[[", "estsi"))
  res <- data.frame(estsi, var = "I")

  estss <- t(sapply(x, "[[", "estss"))
  res2 <- data.frame(estss, var = "S")
  names(res2) <- names(res)
  res <- rbind(res, res2)

  estsc <- t(sapply(x, "[[", "estsc"))
  res2 <- data.frame(estsc, var = "C")
  names(res2) <- names(res)
  res <- rbind(res, res2)

  res$omega[is.na(res$omega)] <- 0
  res$a[is.na(res$a)] <- 0

  data.frame(res, lambda1 = pars[1, 1], lambda2 = pars[1, 2], repnum = repnum,
             N_0 = N_0, tstep = tstep, dnoise = dnoise, fnoise = fnoise,
             prob_rep = prob_rep)
}

calc_acf <- function(R0 = 17, N_0 = 2e6, eta = 2e-4, tstep = 1 / 52,
                     tstop = 1000, burninyrs = 10, dnoise = 0, fnoise = 0){

  params <- c(b_0 = NA, d = 0.02, sigma_env_d = dnoise, sigma_env_f = fnoise,
              gamma = 365 / 22, eta = eta, mu = N_0 * 0.02, N_0 = N_0, p = 0,
              dbdt = 0)

  params["b_0"] <- with(as.list(params), (gamma + d) * R0)

  equil <- with(as.list(params), get_equil(beta = b_0, eta = eta,
                                           gamma = gamma, d = d, p = p,
                                           sigma_env_d = sigma_env_d,
                                           sigma_env_f = sigma_env_f))
  params["X1_0"] <- equil$sstar
  params["X2_0"] <- equil$istar
  params["X3_0"] <- 1 - equil$sstar - equil$istar

  drift <- with(as.list(params),
                get_drift(beta = b_0, eta = eta, d = d, gamma = gamma, S = X1_0,
                          I = X2_0, sigma_env_f = sigma_env_f,
                          sigma_env_d = sigma_env_d))
  diffu <- with(as.list(params),
                get_diffu(beta = b_0, eta = eta, d = d, gamma = gamma,
                          mu = mu, S = X1_0, I = X2_0, p = p, N_0 = N_0,
                          sigma_env_f = sigma_env_f, sigma_env_d = sigma_env_d))

  Lambda <- eigen(drift)$values
  T <- eigen(drift)$vectors
  Dt <- solve(T) %*% diffu %*% solve(t(T))
  factor <- 1 / outer(Lambda, Lambda, "+")
  Sigmat <- -Dt * factor
  Sigma <- T %*% Sigmat %*% t(T)
  acf <- function(x) T %*% diag(exp(Lambda * x)) %*% Sigmat %*% t(T)

  samp_times <- seq(from = 0, to = tstop, by = tstep)
  times <- burninyrs + samp_times
  sim <- create_sir_sim(t0 = 0, params = params, times = times)

  out <- pomp::simulate(sim, states = TRUE)
  outI <- ts(out["X2", 1, ], deltat = tstep)

  lmax <- length(outI) - 30
  lags <- head(samp_times, n = lmax + 1)
  numeric <- stats::acf(outI, lag.max = lmax, plot = FALSE)$acf[, 1, 1]
  tmpf <- function(x) acf(x)[2, 2] / Sigma[2, 2]
  analyt <- sapply(lags, tmpf)
  data.frame(lag = lags, numeric = numeric, analyt = analyt, R0 = R0,
             N_0 = N_0, eta = eta, tstep = tstep, dnoise = dnoise,
             fnoise = fnoise)
}

palette <- c("X, susceptibles"="#999999", "Y, infecteds"="#E69F00", "C, case reports"="#56B4E9")
palette2 <- c("#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
