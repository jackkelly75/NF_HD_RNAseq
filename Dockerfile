############################################################
# Dockerfile for NF_HD_RNAseq
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM ubuntu:18.04

# File Author / Maintainer
MAINTAINER Jack Kelly <jackkelly75@gmail.com>

# salmon binary will be installed in /home/salmon/bin/salmon


ENV PACKAGES git build-essential libtbb-dev libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev libcurl4-openssl-dev curl unzip autoconf \
    apt-transport-https ca-certificates gnupg software-properties-common \
    libz-dev wget apt-utils \
    clang vim libxml2-dev libssl-dev

ENV SALMON_VERSION 1.1.0

WORKDIR /home


RUN apt-get update && \
    apt remove -y libcurl4 && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean


#install pandoc
RUN wget https://github.com/jgm/pandoc/releases/download/2.7.3/pandoc-2.7.3-1-amd64.deb
RUN wget https://github.com/jgm/pandoc/releases/download/2.7.3/pandoc-2.7.3-linux.tar.gz
RUN dpkg -i pandoc-2.7.3-1-amd64.deb
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends texlive 



RUN apt install -y default-jre
RUN apt install -y openjdk-11-jre-headless
RUN apt install -y openjdk-8-jre-headless
#RUN apt install -y openjdk-9-jre-headless
RUN apt install -y default-jre

RUN  apt install -y --no-install-recommends


RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add -

RUN apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
RUN apt-get update

RUN apt-key --keyring /etc/apt/trusted.gpg del C1F34CDD40CD72DA

RUN apt-get install kitware-archive-keyring

RUN apt-get update

#RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends cmake

RUN wget https://github.com/Kitware/CMake/releases/download/v3.15.5/cmake-3.15.5-Linux-x86_64.sh \
      -q -O /tmp/cmake-install.sh \
      && chmod +x /tmp/cmake-install.sh \
      && mkdir /usr/bin/cmake \
      && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
      && rm /tmp/cmake-install.sh
ENV PATH="/usr/bin/cmake/bin:${PATH}"



#fastqpuri install
# Install the packages necessary to add a new repository over HTTPS:
RUN apt install -y apt-transport-https software-properties-common
#Enable the CRAN repository and add the CRAN GPG key to your system
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends gfortran r-base libpng-dev libjpeg-dev liblapack-dev libblas-dev

RUN Rscript -e 'install.packages(c("rmarkdown", "pheatmap", "DMwR", "stringr", "nlme"), repos="https://cran.ma.imperial.ac.uk/")'
RUN Rscript -e 'update.packages(ask = FALSE)'
RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' -e 'BiocManager::install(c("RCurl", "DESeq2", "openssl", "biomaRt", "EnsDb.Hsapiens.v86", "IHW", "tximport"))'


RUN curl -k -L https://github.com/jengelmann/FastqPuri/archive/v1.0.6.tar.gz -o FastqPuri-1.0.6.tar.gz && \
   tar xzvf FastqPuri-1.0.6.tar.gz && \
    cd FastqPuri-1.0.6 && \
    mkdir build && \
    cd build && \
    cmake .. -DRSCRIPT=/usr/bin/Rscript -DCMAKE_INSTALL_PREFIX=/usr/local && make && make install


#salmon install
RUN curl -k -L https://github.com/COMBINE-lab/salmon/releases/download/v1.2.0/salmon-1.2.0_linux_x86_64.tar.gz -o salmon-1.2.0_linux_x86_64.tar.gz && \
    tar xzvf salmon-1.2.0_linux_x86_64.tar.gz

ENV PATH /home/salmon-latest_linux_x86_64/bin:${PATH}
ENV LD_LIBRARY_PATH "/usr/local/lib:${LD_LIBRARY_PATH}"

RUN echo "export PATH=$PATH" > /etc/environment
RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > /etc/environment
