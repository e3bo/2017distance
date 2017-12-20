all: \
change-in-distance.pdf \
distance-vs-fnoise.pdf \
distance-vs-dnoise.pdf \
distance-to-threshold-scatter.pdf \
autocor-by-dnoise.pdf \
autocor-by-dnoise-20yr.pdf \
distance-vs-sampling-freq.pdf \
distance-vs-popsize.pdf \
example-time-series.pdf \
distance-to-threshold-scatter.pdf \
noise-map-adjusted.pdf

%.pdf: %.tex
	/usr/bin/pdflatex -halt-on-error $(<F) >$(<F).pdflatex.stdout 2>$(<F).pdflatex.stderr

noise-map-adjusted.tex: make-noise-map-figure.R adjust-text.R
	./$< >noise-map.stdout 2>noise-map.stderr && ./$(word 2,$^) >adjust.stdout 2>adjust.stderr

change-in-distance.tex: example-speed-ests-plot.R example-speed-ests.RData
	./$< >speed-ests-plot.stdout 2>speed-ests-plot.stderr

example-speed-ests.RData: example-speed-ests.R example-distance-est-core.R
	./$< >speed-ests.stdout 2>speed-ests.stderr

distance-vs-fnoise.tex: example-distance-est-foi-noise-plot.R example-distance-est-fnoise.RData
	./$< >fnoise-plot.stdout 2>fnoise-plot.stderr

example-distance-est-fnoise.RData: example-distance-est-foi-noise.R example-distance-est-core.R
	./$< >est-fnoise.stdout 2>est-fnoise.stderr

autocor-by-dnoise-20yr.tex autocor-by-dnoise.tex: distance-vs-dnoise.tex

distance-vs-dnoise.tex: example-distance-est-d-noise-plot.R example-distance-est-dnoise.RData
	./$< >dnoise-plot.stdout 2>dnoise-plot.stderr

example-distance-est-dnoise.RData: example-distance-est-d-noise.R example-distance-est-core.R
	./$< >est-dnoise.stdout 2>est-dnoise.stderr

distance-vs-sampling-freq.tex: example-distance-est-sampling-plot.R example-distance-est-sampling.RData
	./$< >sampling-plot.stdout 2>sampling-plot.stderr

example-distance-est-sampling.RData: example-distance-est-sampling.R example-distance-est-core.R
	./$< >est-sampling.stdout 2>est-sampling.stderr

distance-vs-popsize.tex: example-distance-est-popsize-plot.R example-distance-est-popsize.RData
	./$< >est-popsize-plot.stdout 2>est-popsize-plot.stderr

example-distance-est-popsize.RData: example-distance-est-popsize.R example-distance-est-core.R
	./$< >est-popsize.stdout 2>est-popsize.stderr

example-time-series.tex: plot-example-time-series.R example-distance-est-R0.RData
	./$< >ex-ts-plot.stdout 2>ex-ts-plot.stderr

distance-to-threshold-scatter.tex: example-distance-est-R0-plot.R example-distance-est-R0.RData
	./$< >est-R0-plot.stdout 2>est-R0-plot.stderr

example-distance-est-R0.RData: example-distance-est-R0.R example-distance-est-core.R
	./$< >est-R0.stdout 2>est-R0.stderr

clean:
	rm *.stdout *.stderr
