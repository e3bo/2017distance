FROM r-base
MAINTAINER Eamon O'Dea <[last name without apostrophe]35@gmail.com>

RUN apt-get update && apt-get install -y -q --no-install-recommends \
  preview-latex-style \
  texlive-base \
  texlive-fonts-recommended \
  texlive-pictures

RUN install2.r --error packrat
RUN mkdir -p /home/docker/
COPY ./packrat /home/docker/packrat
COPY ./.Rprofile /home/docker/
RUN /usr/bin/Rscript -e "setwd(\"/home/docker\"); source(\".Rprofile\"); packrat::restore()"
RUN chown -R docker:docker /home/docker/





