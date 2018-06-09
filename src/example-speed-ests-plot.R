#!/usr/bin/Rscript

library(ggplot2)
library(tikzDevice)
load("example-speed-ests.RData")

e1 <- sapply(se, "[[", "estsw1")
e2 <- sapply(se, "[[", "estsw2")

imaginary <- !is.na(e1["omega", ]) | !is.na(e2["omega", ])
if (any(imaginary)){
    print("Warning: imaginary parts were estimated in the following replicates:")
    print(which(imaginary))
    print("These are the estimates:")
    print(cbind(e1[c("gamma", "omega"), imaginary], e2[c("gamma", "omega"), imaginary]))
}

dist2 <- ifelse(is.na(e2["omega", ]), sqrt(e2["gamma", ]^2), sqrt(e2["gamma", ]^2 - e2["omega", ]^2))
dist1 <- ifelse(is.na(e1["omega", ]), sqrt(e1["gamma", ]^2), sqrt(e1["gamma", ]^2 - e1["omega", ]^2))
speed_ests <- dist2 - dist1

target_dist <- sapply(se[[1]]$del, "[[", "lambda")[1, ]
itime <- time(se[[1]]$outI)
target_dist <- -ts(target_dist, start = min(itime), end = max(itime), deltat = itime[2] - itime[1])
target_change <- target_dist[round(itime * 52) == 40 * 52] - target_dist[round(itime * 52) == 20 * 52]

tikz("change-in-distance.tex", width = 6.5, height = 6, standAlone = TRUE)
par(mfrow = c(3, 1), mar = c(4, 5, 1, 1), cex = 1)
palette2 <- c("#009E73", "#F0E442")
plot(target_dist, xlab = "Time (y)", ylab = "Distance to\nthreshold", frame = FALSE, yaxt = "n")
axis(2, c(10, 15))
plot(se[[1]]$outI, xlab = "Time (y)", ylab = "Infecteds, Y", type = 'n', frame = FALSE)
lines(se[[1]]$w1, col = palette2[1], lwd = 2)
lines(se[[1]]$w2, col = palette2[2], lwd = 2)
legend("topleft", legend = c("Window 1", "Window 2"), lwd = 2, col = palette2,
       ncol = 2, box.lty = 0)
hist(speed_ests, xlab = "Estimated change\nin distance", main = "")
abline(v=target_change, col= "grey", lwd = 2)
dev.off()
