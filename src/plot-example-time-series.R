#!/usr/bin/Rscript

library(ggplot2)
library(tikzDevice)
library(reshape2)
load("example-distance-est-R0.RData")

selected_repnums <- c(0.1, 0.9, 2, 4)

get_df <- function(repnum){
  ind <- which(abs(repnum - repnum_seq) < .Machine$double.eps)
  res <- data.frame(t(estl[[ind]][[1]]$simts[c("X1", "X2"), 1, ]))
  res$repnum <- repnum
  res$time <- seq(0, by = des$tstep[ind], length.out = nrow(res))
  res$X1 <- res$X1 - mean(res$X1)
  res$X2 <- res$X2 - mean(res$X2)
  colnames(res)[1:2] <- c("Susceptibles, $X$",
                          "Infecteds, $Y$")
  res
}

dfl <- lapply(selected_repnums, get_df)
df <- do.call(rbind, dfl)

mdf <- melt(df, id.vars = c("repnum", "time"))
mdf$R0 <- factor(paste0("$R_0 = ", mdf$repnum, "$"),
                 levels = paste0("$R_0 = ", sort(unique(mdf$repnum)), "$"))

g <- ggplot(data = mdf, aes(x = time, y = value))
g <- g + geom_line()
g <- g + facet_grid(variable ~ R0, scales = "free_y")
g <- g + theme_minimal()
g <- g + labs(x = "Time (y)", y = "Deviation from mean")
g <- g + theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0)))
g <- g + theme(axis.title.x = element_text(margin = margin(10, 0, 0, 0)))

tikz("example-time-series.tex", width = 6.5, height = 4, standAlone = TRUE)
print(g)
dev.off()
