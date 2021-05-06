############################################################
# Dockerfile to build SVClone container image
# Based on Ubuntu
############################################################

FROM rocker/rstudio

MAINTAINER Belal_C

RUN mkdir /opt/svclone

COPY ./assign_remaining.R /opt/svclone

WORKDIR /opt/svclone

RUN chmod +x /opt/svclone/assign_remaining.R

ENV USER=ubuntu