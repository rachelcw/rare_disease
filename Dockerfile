# Choose a base image
FROM ubuntu:18.04

LABEL maintainer="Racheli Cohen" 

# Install dependencies
RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
    build-essential \
        cmake \
        curl \
        git \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        libhdf5-serial-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        python \
        python-pip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

#python modules
RUN pip install --upgrade pip setuptools
RUN pip install 'pystan==2.17.1.0' 'numpy==1.15.4' 

# Clone the SPOT repository
RUN git clone --branch v1.0.2 https://github.com/BennyStrobes/SPOT.git /opt/spot

# Set the working directory
WORKDIR /opt/spot

# Set environment variables
ENV PYTHONPATH=/opt/spot/src
# ENV SPOT_CONFIG=/opt/spot/config.yaml 
