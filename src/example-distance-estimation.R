#!/usr/bin/Rscript

set.seed(2)

sir_stepfun <- pomp::Csnippet("
  double dN01, dN03;
  double dN10, dN12;
  double dN20, dN23;
  double dN30;

  double mu10, mu12;
  double mu20, mu23;
  double mu30;

  double d_perturb;

  double rate[5], trans[5], dW;

  dW = rgammawn(sigma_env, dt);
  d_perturb = d * dW / dt;

  mu10 = d_perturb;
  mu20 = d_perturb;
  mu12 = b * X2 / (X1 + X2 + X3) + eta;
  mu23 = gamma;
  mu30 = d_perturb;

  rate[0] = mu10;
  rate[1] = mu12;
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
")

sir_vectorfield <- pomp::Csnippet("
  double lambda;

  lambda = b * X2 / (X1 + X2 + X3) + eta;

  DX1 =  -lambda * X1 - d * X1                     + mu * (1 - p);
  DX2 =   lambda * X1 - d * X2      - gamma * X2;
  DX3 =               - d * X3      + gamma * X2   + mu * p;
")

sir_initializer <- function(params, t0, ...) {
    statenames <- c("X1", "X2", "X3", "cases")
    x0 <- params[c("X1_0", "X2_0", "X3_0")]
    x0 <- c(round(x0 / sum(x0) * params["N_0"]), 0)
    names(x0) <- statenames
    x0
}

create_sir_sim <- function(times = seq(0, 9),
                                t0 = min(times),
                                params, delta.t = 0.0001){
  statenames <- c("X1", "X2", "X3", "cases")
  rp_names <- c("b", "d", "sigma_env", "gamma", "eta", "mu", "p")
  ivp_names <- c("X1_0", "X2_0", "X3_0", "N_0")
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

get_drift <- function(beta, eta, d, gamma, mu, S, I){
  ret <- rbind(c(-beta * I - eta - d, -beta * S),
               c(beta * I + eta, beta * S - (gamma + d)))
  dimnames(ret) <- list(c('S', 'I'), c('S', 'I'))
  ret
}

get_diffu <- function(beta, eta, d, gamma, mu, S, I){
  rate_trans <- (beta * I + eta ) * S
  ret <- rbind(c(rate_trans + d  + d * S, -rate_trans),
               c(-rate_trans, rate_trans + (gamma + d) * I))
  dimnames(ret) <- list(c('S', 'I'), c('S', 'I'))
  ret
}

get_diffu_ext <- function(d, sigma, S, I){
  ret <- rbind(c(S ^ 2, S * I),
               c(I * S, I ^ 2))
  dimnames(ret) <- list(c('S', 'I'), c('S', 'I'))
  ret * (sigma * d) ^ 2
}

get_fit <- function(y, tstep) {
  x <- seq_along(y) * tstep
  start <- list()
  start$gamma <- unname(coef(lm(log(I(abs(y + 1)))~x))["x"])
  spec <- spectrum(y, plot = FALSE)
  start$omega <- spec$freq[which.max(spec$spec)]
  start$a <- 0
  fit <- minpack.lm::nlsLM(
      y~sqrt(1 + a^2) * exp(x * gamma) * sin(x * omega + atan2(1, a)),
      start = start,
      control = minpack.lm::nls.lm.control(maxiter = 1000))
  fit1 <- minpack.lm::nlsLM(
      y~exp(x * gamma),
      start = list(gamma = start$gamma),
      control = minpack.lm::nls.lm.control(maxiter = 1000))
  e <- fit$m$resid()
  e1 <- fit1$m$resid()
  if (sum(e^2) + 2 * 3 <  sum(e1^2) + 2){
    bestfit <- fit
    worstfit <- fit1
  } else {
    bestfit <- fit1
    worstfit <- fit
  }
  coefests <- coef(bestfit)[c("omega", "gamma", "a")]
  names(coefests) <- c("omega", "gamma", "a")
  list(coef = coefests, fit = bestfit, altfit = worstfit)
}

##

get_distance_est <- function(R0 = 17, N_0 = 2e6, eta = 2e-4, tstep = 1 / 52,
                             tlength = 1000, burninyrs = 10, savemem = TRUE){

  params <- c(b = 365, d = 0.02, sigma_env = 0, gamma = 365 / 22, eta = eta,
              mu = N_0 * 0.02, N_0 = N_0, p = 0)
  params["b"] <- with(as.list(params), (gamma + d) * R0)

  equil <- with(as.list(params), EndemicEquilSIR(beta = b, eta = eta,
                                                 gamma = gamma, mu = d, p = p))
  params["X1_0"] <- equil$S
  params["X2_0"] <- equil$I
  params["X3_0"] <- 1 - equil$S - equil$I

  drift <- with(as.list(params),
                get_drift(beta = b, eta = eta, d = d, gamma = gamma, S = X1_0,
                          I = X2_0))
  lambda <- eigen(drift)$values

  freq <- 100
  times <- burninyrs + seq(from = 0, length.out = tlength, by = tstep)
  sim <- create_sir_sim(t0 = 0, params = params, times = times)

  out <- pomp::simulate(sim, states = TRUE, params = params)
  outI <- ts(out["X2", 1, ], deltat = tstep)
  outS <- ts(out["X1", 1, ], deltat = tstep)
  outR <- ts(out["X3", 1, ], deltat = tstep)
  outC <- ts(out["cases", 1, -1], deltat = tstep)

  lag.max <- dim(out)[3] - 30
  acesti <- acf(outI, lag.max = lag.max, plot = FALSE)$acf[, 1, 1]
  acests <- acf(outS, lag.max = lag.max, plot = FALSE)$acf[, 1, 1]
  acestc <- acf(outC, lag.max = lag.max - 1, plot = FALSE)$acf[, 1, 1]

  if(any(is.na(acesti)) | any(is.na(acests))) {
    fitc <- fiti <- fits <- list(coef=c(omega = NA, gamma = NA, a = NA))
  } else {
    fiti <- get_fit(acesti, tstep = tstep)
    fits <- get_fit(acests, tstep = tstep)
    fitc <- get_fit(acestc, tstep = tstep)
  }
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

repnum_seq <- c(0.1, 0.5, 0.9, 2, 4, 8, 16, 32)
des <- data.frame(repnum = repnum_seq, N_0 = 1e7)
des$eta <- 400 / des$N_0

repnum_seq <- c(0.1, 2, 16)
popsize_seq <- c(2e5, 2e4, 2e3, 2e2)

repnum_seq <- c(0.1, 2, 16)
popsize_seq <- c(1e6, 1e5, 1e4, 1e3, 1e2)
des_next <- expand.grid(repnum = repnum_seq, N_0 = popsize_seq)
des_next$eta <- 400 / des_next$N_0

des <- rbind(des, des_next)
des$tstep <- 1 / 52

tstep_seq <- c(1 / 365, 1 / 52, 4 / 52, 1)
des_next <- expand.grid(repnum = repnum_seq, tstep = tstep_seq)
des_next$N_0 <- 1e7
des_next$eta <- 400 / des_next$N_0

des <- rbind(des, des_next)

simulate_ests <- function(repnum, N_0, eta, tstep){
  replicate(100, get_distance_est(tlength = 520 * 2, R0 = repnum, N_0 = N_0,
                                 eta = eta, tstep = tstep), simplify = FALSE)
}
estl <- mapply(simulate_ests, repnum = des$repnum, N_0 = des$N_0,
               eta = des$eta, tstep = des$tstep, SIMPLIFY = FALSE)

makedf <- function(x, repnum, N_0, tstep){
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
             N_0 = N_0, tstep = tstep)
}
pdf <- with(des, mapply(makedf, x = estl, repnum = repnum, N_0 = N_0,
                        tstep = tstep, SIMPLIFY = FALSE))
estpdf <- do.call(rbind, pdf)

library(ggplot2)

sub <- estpdf[estpdf$N_0 < 1e7 & estpdf$var != "C", ]
sub$R0 <- sub$repnum
g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = N_0, color = var))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5)
g <- g + geom_point(data = sub, shape=4, color=1, aes(y=Mod(lambda1)))
g <- g + scale_x_log10()
g <- g + ylab("Distance to threshold\n")
g <- g + xlab("\nPopulation size")
g <- g + guides(color=guide_legend(title="Variable\nused for\nestimate"))
g <- g + theme_minimal()
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5))
g <- g + facet_wrap(~R0, labeller = "label_both")
ggsave("distance-vs-popsize.pdf", plot=g, width=8, height=5)

is_weekly <- abs(estpdf$tstep - 1 / 52) < 1 / 52
sub <- estpdf[estpdf$N_0 > 1e6 & estpdf$var != "C" & is_weekly, ]
g <- ggplot(data=sub, aes(x=abs(omega), y=abs(gamma),
                color=as.factor(repnum)))
g <- g + geom_point()
g <- g + geom_point(data=sub, shape=2,
                    aes(x=abs(Im(lambda1)), y=abs(Re(lambda1)),
                        color=as.factor(repnum)))
g <- g + theme_minimal()
ggsave("gamma-omega-scatter.pdf", plot=g)

g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = repnum, color = var))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5)
g <- g + geom_point(data=sub, shape=4, color=1, aes(y=Mod(lambda1)))
g <- g + scale_x_log10(breaks=c(0.1, 0.5, 0.9, 2, 4, 8, 16, 32))
g <- g + ylab("Distance to threshold\n")
g <- g + xlab("\nBasic reproduction number")
g <- g + guides(color=guide_legend(title="Variable\nused for\nestimate"))
g <- g + theme_minimal()
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5))
ggsave("distance-to-threshold-scatter.pdf", plot=g, width=4, height=5)

sub <- estpdf[estpdf$N_0 > 1e6 & estpdf$repnum %in% c(0.1, 2, 16), ]
sub$R0 <- sub$repnum
g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = 1 / tstep,
                color = var))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5)
g <- g + geom_point(data = sub, shape=4, color=1, aes(y=Mod(lambda1)))
g <- g + scale_x_log10() + scale_y_log10()
g <- g + ylab("Distance to threshold\n")
g <- g + xlab("\nObservations per year")
g <- g + guides(color=guide_legend(title="Variable\nused for\nestimate"))
g <- g + theme_minimal()
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5))
g <- g + facet_wrap(~R0, labeller = "label_both")
ggsave("distance-vs-sampling-freq.pdf", plot=g, width=8, height=5)

save.image("example-distance-estimation.RData")
