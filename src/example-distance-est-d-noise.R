#!/usr/bin/Rscript
source("../src/example-distance-est-core.R")

repnum_seq <- c(0.1, 2, 16)
dnoise_seq <- c(0.01, 0.1, 1)
N_0_seq <- c(1e7, 1e9)
des <- expand.grid(repnum = repnum_seq, dnoise = dnoise_seq, N_0 = N_0_seq)

des$eta <- 1 / 3 / sqrt(1e7)
des$tstep <- 1 / 52
des$fnoise <- 0

estl <- parallel::mcMap(simulate_ests, repnum = des$repnum, N_0 = des$N_0,
                        eta = des$eta, tstep = des$tstep, dnoise = des$dnoise,
                        fnoise = des$fnoise)
pdf <- with(des, mapply(makedf, x = estl, repnum = repnum, N_0 = N_0,
                        tstep = tstep, dnoise = dnoise, fnoise = fnoise,
                        SIMPLIFY = FALSE))
estpdf <- do.call(rbind, pdf)
estpdf$N_0 <- as.factor(estpdf$N_0)

makevardf <- function(x, repnum, N_0, tstep, dnoise){
  pars <- t(sapply(x, "[[", "lambda"))

  getvar <- function(x) var(x$simts["X2", 1, ])
  var <- sapply(x, getvar)

  data.frame(variance = var, var = "I", lambda1 = pars[1, 1],
             lambda2 = pars[1, 2], repnum = repnum, N_0 = N_0, tstep = tstep,
             dnoise = dnoise)
}

vdf <- with(des, mapply(makevardf, x = estl, repnum = repnum, N_0 = N_0,
                        tstep = tstep, dnoise = dnoise, SIMPLIFY = FALSE))
vardf <- do.call(rbind, vdf)


sub <- des[des$repnum < .11,]
acfdf <- with(sub, parallel::mcMap(calc_acf, R0 = repnum, N_0 = N_0, eta = eta,
                                   tstep = tstep, dnoise = dnoise))
acfdf <- do.call(rbind, acfdf)
acfdf$N_0 <- as.factor(acfdf$N_0)
acfdf$dnoise <- as.factor(acfdf$dnoise)

acfdf20 <- with(sub, parallel::mcMap(calc_acf, R0 = repnum, N_0 = N_0,
                                     eta = eta, tstep = tstep, dnoise = dnoise,
                                     tstop = 20))
acfdf20 <- do.call(rbind, acfdf20)
acfdf20$N_0 <- as.factor(acfdf20$N_0)
acfdf20$dnoise <- as.factor(acfdf20$dnoise)

save.image("example-distance-est-dnoise.RData")

