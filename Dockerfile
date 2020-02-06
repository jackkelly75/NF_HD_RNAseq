############################################################
# Dockerfile for NF_HD_RNAseq
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu:18.04

# File Author / Maintainer
MAINTAINER Jack Kelly <jackkelly75@gmail.com>

# salmon binary will be installed in /home/salmon/bin/salmon


WORKDIR /home


### salmon

RUN apt-get update
RUN apt remove -y libcurl4
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
	apt-utils \
	git \
	gcc \
	make \
	g++ \
	libboost-all-dev \
	liblzma-dev \
	libbz2-dev \
	ca-certificates \
	zlib1g-dev \
	libcurl4-openssl-dev \
	curl \
	unzip \
	autoconf \
	apt-transport-https \
	ca-certificates \
	gnupg \
	software-properties-common \
	wget

RUN apt-get clean

RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add -

RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'

RUN apt-get update

RUN apt-key --keyring /etc/apt/trusted.gpg del C1F34CDD40CD72DA

RUN apt-get install kitware-archive-keyring

RUN apt-get install -y cmake

RUN curl -k -L https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz -o salmon-v1.1.0.tar.gz

RUN tar xzf salmon-v1.1.0.tar.gz

ENV PATH /home/salmon-latest_linux_x86_64/bin:${PATH}
ENV LD_LIBRARY_PATH "/usr/local/lib:${LD_LIBRARY_PATH}"

RUN echo "export PATH=$PATH" > /etc/environment
RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > /etc/environment


RUN DEBIAN_FRONTEND=noninteractive apt-get install -y r-base pandoc vim

RUN Rscript -e 'install.packages(c("rmarkdown", "pheatmap"), repos="https://cran.ma.imperial.ac.uk/")'

##"DMwR", "sjPlot", "stringr"
#RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite(c("DESeq2", "biomaRt", "EnsDb.Hsapiens.v86","IHW", "tximport"))'

# compile and install FastqPuri
RUN cd /home && git clone https://github.com/jengelmann/FastqPuri
RUN cd /home/FastqPuri && cmake -H. -Bbuild/ -DRSCRIPT=/usr/bin/Rscript
RUN cd /home/FastqPuri/build && make && make install







