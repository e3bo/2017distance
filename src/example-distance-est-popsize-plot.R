#!/usr/bin/Rscript

library(ggplot2)
library(tikzDevice)

load("example-distance-est-popsize.RData")

sub <- estpdf[estpdf$var != "C", ]

sub$R0 <- factor(paste0("$R_0 = ", sub$repnum, "$"),
                 levels = paste0("$R_0 = ", sort(unique(sub$repnum)), "$"))

sub$var_plot <- factor(sub$var, levels=c("S", "I"), labels=c("X", "Y"))

names(palette) <- c("X", "Y", "C")
g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = N_0, color = var_plot))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5)
g <- g + geom_point(data = sub, shape=4, color=1, aes(y=Mod(lambda1)))
breaks <- c(1e2, 1e3, 1e5, 1e6)
labs <- paste0(paste0("$", paste0("10^", log10(breaks))), "$")
g <- g + scale_x_continuous(breaks = breaks, labels = labs, trans = "log10")
g <- g + ylab("Distance to threshold\n")
g <- g + xlab("\nPopulation size, $N_0$")
g <- g + guides(color=guide_legend(title="Variable"))
g <- g + theme_minimal()
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5))
g <- g + facet_wrap(~R0)
g <- g + scale_color_manual(values = palette,
                            labels = list("X" = "Susceptibles, $X$",
                                          "Y" = "Infecteds, $Y$"))

tikz("distance-vs-popsize.tex", width = 6.5, height = 5, standAlone = TRUE)
print(g)
dev.off()
