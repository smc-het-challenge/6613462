############################################################
# Dockerfile to build SVClone container image
# Based on Ubuntu
############################################################

FROM ubuntu:14.04

MAINTAINER Belal_C

###USER LOGNAME normally set by login shell; we must still use nss wrapper to solve uid issue!

RUN apt-get update
RUN apt-get install -y python python-pip python-dev python-setuptools build-essential
RUN apt-get install -y python-numpy python-matplotlib python-pandas
RUN apt-get install -y ipython ipython-notebook python-ipdb
RUN apt-get install -y gfortran fort77 libgsl0ldbl cmake make gcc sudo

RUN pip install --upgrade pip
RUN pip install --upgrade virtualenv
RUN pip install --upgrade wheel

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get install -y zlib1g-dev

RUN easy_install https://pypi.python.org/packages/3.4/p/pymc/pymc-2.3.5-py3.4-linux-x86_64.egg

RUN pip install cython
RUN pip install pysam

RUN pip install pyvcf

ADD svclone /opt/svclone

WORKDIR /opt/svclone
RUN python setup.py install

RUN chmod -R 755 /usr/local/lib/python2.7/dist-packages/SVclone-0.1.1-py2.7.egg 

ENV USER=ubuntu