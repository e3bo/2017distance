#!/usr/bin/Rscript

library(ggplot2)
library(tikzDevice)
load("example-distance-est-reporting-prob.RData")

sub <- estpdf[estpdf$var == "C", ]

sub$R0 <- factor(paste0("$R_0 = ", sub$repnum, "$"),
                 levels = paste0("$R_0 = ", sort(unique(sub$repnum)), "$"))

g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = prob_rep))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5, color = palette[1])
g <- g + geom_point(data = sub, shape=4, color=1, aes(y=Mod(lambda1)))
g <- g + scale_y_log10()
g <- g + ylab("Distance to threshold\n")
g <- g + xlab("\nReporting probability")
g <- g + theme_minimal()
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5))
g <- g + facet_wrap(~R0)

tikz("distance-vs-reporting-prob.tex", width = 6.5, height = 5, standAlone = TRUE)
print(g)
dev.off()
