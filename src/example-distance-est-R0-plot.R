#!/usr/bin/Rscript

library(ggplot2)
library(tikzDevice)

load("example-distance-est-R0.RData")

sub <- estpdf[estpdf$var != "C", ]
g <- ggplot(data= sub, aes(x = abs(omega), y = abs(gamma),
                color = as.factor(repnum)))
g <- g + geom_point()
g <- g + geom_point(data = sub, shape=2,
                    aes(x = abs(Im(lambda1)), y = abs(Re(lambda1)),
                        color = as.factor(repnum)))
g <- g + theme_minimal()
ggsave("gamma-omega-scatter.pdf", plot = g)

sub$var_plot <- factor(sub$var, levels=c("S", "I", "C"), labels=c("X", "Y", "C"))
names(palette) <- c("X", "Y", "C")

g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = repnum, color = var_plot))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5)
g <- g + geom_point(data = sub, shape = 4, color = 1, aes(y = Mod(lambda1)))
g <- g + scale_x_log10(breaks=c(0.1, 0.5, 0.9, 2, 4, 8, 16, 32))
g <- g + ylab("Distance to threshold\n")
g <- g + xlab("\nBasic reproduction number, $R_0$")
g <- g + guides(color=guide_legend(title="Variable"))
g <- g + theme_minimal()
g <- g + scale_color_manual(values = palette)
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5))
g <- g + theme(legend.position = "top")
g <- g + scale_color_manual(values = palette,
                            labels = list("X" = "Susceptibles, $X$",
                                          "Y" = "Infecteds, $Y$",
                                          "C" = "Case reports, $C$"))


tikz("distance-to-threshold-scatter.tex", width = 3.25, height = 3.25, standAlone = TRUE)
print(g)
dev.off()
