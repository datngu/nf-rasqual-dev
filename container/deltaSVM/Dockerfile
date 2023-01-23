FROM nfcore/base:2.1

LABEL authors="Dat T Nguyen - ndat<at>utexas.edu" \
      description="Docker image containing all requirements for running gkmsvm" 

ENV TZ='Europe/Oslo'
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt-get update && apt-get install -y \
   build-essential \
   libboost-dev \
   libgsl0-dev \
   gfortran \
   zlib1g-dev \
   libcurl4-openssl-dev \
   libblas-dev \
   libboost-all-dev \
   liblapack-dev \
   libbz2-dev \
   unzip
 

ADD environment.yml /
RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/gkmsvm/bin:$PATH

ADD bsgenome /bsgenome
ADD bin/genNullSeqs_tfr.R /usr/bin/
ADD bin/trf_dat2bed.py /usr/bin/
