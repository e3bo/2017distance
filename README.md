These files are available to reproduce the results for the manuscript
"Estimating the distance to an epidemic threshold." The most important files are the R scripts that
contain the code for all of the calculations. The Makefile shows
the correct order to run the scripts in to generate all of the
results. The packrat directory contains a packrat archive of all of
the R packages used by the scripts. The Docker file gives the
instructions for building a docker image in which the R scripts may be
run. A Docker image built with this script is available on Docker
Hub. So, on a system with Docker installed, running the bash script
``run-scripts`` from within the directory containing this README will
pull down the image and run the Makefile inside of it to reproduce all
of the results. The results will be placed in a time-stamped
subdirectory. Reproducing the results takes less the 60 minutes on
workstation with 20 2.1GHz processors and 188 G of memory.
