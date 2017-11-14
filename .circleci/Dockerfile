FROM ubuntu:xenial

RUN apt-get -qq update && apt-get upgrade --no-install-recommends -y
RUN apt-get install --no-install-recommends -y git wget bzip2 ca-certificates gfortran ssh tar gzip

RUN wget -q "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O $HOME/miniconda.sh && bash $HOME/miniconda.sh -b -f -p $HOME/miniconda

USER root
ENV PATH="${PATH}:/root/miniconda/bin"
