#!/usr/bin/Rscript
source("../src/example-distance-est-core.R")

repnum_seq <- c(0.1, 2, 16)
fnoise_seq <- c(0.01, 0.1, 1)
N_0_seq <- c(1e7, 1e9)
des <- expand.grid(repnum = repnum_seq, fnoise = fnoise_seq, N_0 = N_0_seq)
des$dnoise <- 0
des$eta <- 400 / des$N_0
des$tstep <- 1 / 52

estl <- parallel::mcMap(simulate_ests, repnum = des$repnum, N_0 = des$N_0,
                        eta = des$eta, tstep = des$tstep, dnoise = des$dnoise,
                        fnoise = des$fnoise)
pdf <- with(des, Map(makedf, x = estl, repnum = repnum, N_0 = N_0,
                        tstep = tstep, dnoise = dnoise, fnoise = fnoise))
estpdf <- do.call(rbind, pdf)

###

sub <- des[des$repnum < .11,]
acfdf <- with(sub, parallel::mcMap(calc_acf, R0 = repnum, N_0 = N_0, eta = eta,
                                   tstep = tstep, dnoise = dnoise,
                                   fnoise = fnoise))
acfdf <- do.call(rbind, acfdf)

acfdf20 <- with(sub, parallel::mcMap(calc_acf, R0 = repnum, N_0 = N_0, eta = eta,
                                     tstep = tstep, dnoise = dnoise, fnoise = fnoise,
                                     tstop = 20))
acfdf20 <- do.call(rbind, acfdf20)

save.image("example-distance-est-fnoise.RData")
