#!/usr/bin/Rscript
source("../src/example-distance-est-core.R")

repnum_seq <- c(0.9, 2, 16)
tstep_seq <- c(1 / 365, 1 / 52, 4 / 52, 1)
des <- expand.grid(repnum = repnum_seq, tstep = tstep_seq)
des$N_0 <- 1e7
des$eta <- 400 / des$N_0
des$dnoise <- 0
des$fnoise <- 0

estl <- parallel::mcMap(simulate_ests, repnum = des$repnum, N_0 = des$N_0,
                        eta = des$eta, tstep = des$tstep, dnoise = des$dnoise,
                        fnoise = des$fnoise)
pdf <- with(des, Map(makedf, x = estl, repnum = repnum, N_0 = N_0,
                     tstep = tstep, dnoise = dnoise, fnoise = fnoise))
estpdf <- do.call(rbind, pdf)
save.image("example-distance-est-sampling.RData")
