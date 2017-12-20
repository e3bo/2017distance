#!/usr/bin/Rscript

library(ggplot2)
library(tikzDevice)
load("example-distance-est-fnoise.RData")

sub <- estpdf[estpdf$var != "C", ]
sub$R0 <- factor(paste0("$R_0 = ", sub$repnum, "$"),
                 levels = paste0("$R_0 = ", sort(unique(sub$repnum)), "$"))

names(palette) <- c("X", "Y", "C")
sub$var_plot <- factor(sub$var, levels=c("S", "I"), labels=c("X", "Y"))

g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = fnoise,
                color = var_plot, shape = as.factor(N_0)))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5)
g <- g + geom_point(data = sub, shape = 4, color = 1, aes(y = Mod(lambda1)))
g <- g + scale_y_log10() + scale_x_log10()
g <- g + ylab("Distance to threshold\n")
g <- g + xlab("\nStandard deviation in force of infection, $\\tau_f^{1/2}$")
g <- g + guides(color = guide_legend(title = "Variable"))
g <- g + guides(shape = guide_legend(title = "Population size, $N_0$"))
g <- g + theme_minimal()
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + facet_wrap(~R0)
g <- g + scale_color_manual(values = palette,
                            labels = list("X" = "Susceptibles, $X$",
                                          "Y" = "Infecteds, $Y$"))
g <- g + scale_shape_discrete(labels = list("1e+07"="$10^7$",
                                            "1e+09"="$10^9$"))

tikz("distance-vs-fnoise.tex", width = 6.5, height = 5, standAlone = TRUE)
print(g)
dev.off()


g <- ggplot(data = acfdf[acfdf$lag <= 1, ],
            aes(x = lag, y = analyt, color = as.factor(fnoise)))
g <- g + geom_line()
g <- g + geom_point(aes(y = numeric))
g <- g + ylab("Autocorrelation\n")
g <- g + xlab("\nLag (years)")
g <- g + guides(color = guide_legend(title = "Standard dev.\nof F.O.I."))
g <- g + theme_minimal()
g <- g + theme(legend.position = "top")
ggsave("autocor-by-fnoise.pdf", plot = g, width = 4, height = 5)

g <- ggplot(data = acfdf20[acfdf20$lag <= 1, ],
            aes(x = lag, y = analyt, color = as.factor(fnoise)))
g <- g + geom_line()
g <- g + geom_point(aes(y = numeric))
g <- g + ylab("Autocorrelation\n")
g <- g + xlab("\nLag (years)")
g <- g + guides(color = guide_legend(title = "Standard dev.\nof F.O.I."))
g <- g + theme_minimal()
g <- g + theme(legend.position = "top")
ggsave("autocor-by-fnoise-20yr.pdf", plot = g, width = 4, height = 5)
