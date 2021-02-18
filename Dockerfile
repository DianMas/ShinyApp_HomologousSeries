FROM rocker/shiny:latest

MAINTAINER Steffen Neumann <sneumann@ipb-halle.de>

LABEL Description="Shiny GUI for nontarget."

##
## Install system libraries
##
RUN apt-get -y update && apt-get -y install \
  libxml2-dev libssl-dev

##
## Install R dependencies
##
ADD install.R /tmp
RUN R -e "source('/tmp/install.R')"

WORKDIR /srv/shiny-server

RUN rm -rf * 
ADD app.R .
ADD homol*.csv .


