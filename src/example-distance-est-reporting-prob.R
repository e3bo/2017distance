#!/usr/bin/Rscript
source("../src/example-distance-est-core.R")

prob_rep_seq <- c(0.01, 0.25, 0.5, 0.75, 1)
repnum_seq <- c(0.5, 0.9)

tstep_seq <- c(1 / 52)
des <- expand.grid(prob_rep = prob_rep_seq, repnum = repnum_seq)
des$repnum <- 0.9
des$tstep <- 1 / 52
des$N_0 <- 1e7
des$eta <- 400 / des$N_0
des$dnoise <- 0
des$fnoise <- 0

estl <- parallel::mcMap(simulate_ests, repnum = des$repnum, N_0 = des$N_0,
                        eta = des$eta, tstep = des$tstep, dnoise = des$dnoise,
                        fnoise = des$fnoise, prob_rep = des$prob_rep)
pdf <- with(des, Map(makedf, x = estl, repnum = repnum, N_0 = N_0,
                     tstep = tstep, dnoise = dnoise, fnoise = fnoise,
                     prob_rep = prob_rep))
estpdf <- do.call(rbind, pdf)
save.image("example-distance-est-reporting-prob.RData")
