#!/usr/bin/Rscript
source("../src/example-distance-est-core.R")

repnum_seq <- c(0.1, 0.5, 0.9, 2, 4, 8, 16, 32)
des <- expand.grid(repnum = repnum_seq, N_0 = 1e7)
des$eta <- 1 / 3 / sqrt(des$N_0)
des$tstep <- 1 / 52
des$fnoise <- 0
des$dnoise <- 0

estl <- parallel::mcMap(simulate_ests, repnum = des$repnum, N_0 = des$N_0,
                        eta = des$eta, tstep = des$tstep, fnoise = des$fnoise,
                        dnoise = des$dnoise)
pdf <- with(des, Map(makedf, x = estl, repnum = repnum, N_0 = N_0,
                        tstep = tstep, fnoise = fnoise, dnoise = dnoise))
estpdf <- do.call(rbind, pdf)
save.image("example-distance-est-R0.RData")
