#!/usr/bin/Rscript

library(ggplot2)
library(tikzDevice)
load("example-distance-est-dnoise.RData")

sub <- estpdf[estpdf$var != "C", ]
sub$R0 <- factor(paste0("$R_0 = ", sub$repnum, "$"),
                 levels = paste0("$R_0 = ", sort(unique(sub$repnum)), "$"))

names(palette) <- c("X", "Y", "C")
sub$var_plot <- factor(sub$var, levels=c("S", "I"), labels=c("X", "Y"))
g <- ggplot(data = sub,
            aes(y = sqrt((omega) ^ 2 + (gamma) ^ 2), x = dnoise,
                color = var_plot, shape = N_0))
g <- g + geom_jitter(width = 0.05, height = 0, alpha = 0.5)
g <- g + geom_point(data = sub, shape = 4, color = 1, aes(y = Mod(lambda1)))
g <- g + scale_y_log10() + scale_x_log10()
g <- g + ylab("Distance to threshold")
g <- g + xlab(label = "Standard deviation in death rate, $\\tau_d^{1/2}$")
g <- g + guides(color = guide_legend(title = "Variable"))
g <- g + guides(shape = guide_legend(title = "Population size, $N_0$"))
g <- g + theme_minimal()
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + theme(axis.title.y = element_text(margin = margin(0, 10, 0, 0)))
g <- g + theme(axis.title.x = element_text(margin = margin(10, 0, 0, 0)))
g <- g + facet_wrap(~R0)
g <- g + scale_color_manual(values = palette,
                            labels = list("X" = "Susceptibles, $X$",
                                          "Y" = "Infecteds, $Y$"))
g <- g + scale_shape_discrete(labels = list("1e+07"="$10^7$",
                                            "1e+09"="$10^9$"))

tikz("distance-vs-dnoise.tex", width = 6.5, height = 5, standAlone = TRUE)
print(g)
dev.off()

pdf <- acfdf[acfdf$lag <= 1, ]

pdf$N0 <- factor(pdf$N_0, levels=c("1e+07", "1e+09"),
                 labels = c("$N_0 = 10^7$", "$N_0 = 10^9$"))

g <- ggplot(data = pdf,
            aes(x = lag, y = analyt, color = dnoise))
g <- g + geom_line()
g <- g + geom_point(aes(y = numeric))
g <- g + ylab("Autocorrelation")
g <- g + xlab("Lag (years)")
title <- "Standard deviation of death rate, $\\tau_d^{1/2}$"
g <- g + guides(color = guide_legend(title = title, title.position = "top"))
g <- g + theme_minimal()
g <- g + facet_grid(~N0)
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + theme(axis.title.y = element_text(margin = margin(0,10,0,0)))
g <- g + theme(axis.title.x = element_text(margin = margin(10,0,0,0)))
g <- g + scale_color_manual(values = palette2)

tikz("autocor-by-dnoise.tex", width = 3.25, height = 4, standAlone = TRUE)
print(g)
dev.off()

pdf <- acfdf20[acfdf20$lag <= 1, ]

pdf$N0 <- factor(pdf$N_0, levels=c("1e+07", "1e+09"),
                 labels = c("$N_0 = 10^7$", "$N_0 = 10^9$"))

g <- ggplot(data = pdf,
            aes(x = lag, y = analyt, color = dnoise))
g <- g + geom_line()
g <- g + geom_point(aes(y = numeric))
g <- g + ylab("Autocorrelation")
g <- g + xlab("Lag (years)")
title <- "Standard deviation of death rate, $\\tau_d^{1/2}$"
g <- g + guides(color = guide_legend(title = title, title.position = "top"))
g <- g + theme_minimal()
g <- g + facet_grid(~N0)
g <- g + theme(legend.position = "top")
g <- g + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + theme(axis.title.y = element_text(margin = margin(0,10,0,0)))
g <- g + theme(axis.title.x = element_text(margin = margin(10,0,0,0)))
g <- g + scale_color_manual(values = palette2)

tikz("autocor-by-dnoise-20yr.tex", width = 3.25, height = 4, standAlone = TRUE)
print(g)
dev.off()
