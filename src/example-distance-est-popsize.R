#!/usr/bin/Rscript
source("../src/example-distance-est-core.R")

repnum_seq <- c(0.1, 2, 16)
popsize_seq <- c(1e6, 1e5, 1e4, 1e3, 1e2)
des <- expand.grid(repnum = repnum_seq, N_0 = popsize_seq)
des$eta <- 400 / des$N_0
des$tstep <- 1 / 52
des$fnoise <- 0
des$dnoise <- 0

estl <- parallel::mcMap(simulate_ests, repnum = des$repnum, N_0 = des$N_0,
                        eta = des$eta, tstep = des$tstep, fnoise = des$fnoise,
                        dnoise = des$dnoise)
pdf <- with(des, Map(makedf, x = estl, repnum = repnum, N_0 = N_0,
                     tstep = tstep, dnoise = dnoise, fnoise = fnoise))
estpdf <- do.call(rbind, pdf)
save.image("example-distance-est-popsize.RData")
