# USAGE INSTSRUCTIONS
# 
# Build image (run this in the main directory (where .git folder is)):
# docker build -t nemocrys/dolfinx-openfoam -f ./docker/Dockerfile .
#
# Run image:
# docker run -it --rm nemocrys/dolfinx-openfoam bash
#
# Run image and map current directory to home/workdir on Linux:
# docker run -it --rm -v $PWD:/home/workdir -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/dolfinx-openfoam bash
# and on Windows:
# docker run -it --rm -v ${PWD}:/home/workdir nemocrys/dolfinx-openfoam bash

# FROM ubuntu:20.04
FROM dolfinx/dolfinx:v0.2.0a0
SHELL ["/bin/bash", "-c"]

# Install additional stuff / just update apt-get
RUN apt-get update

# Install python packages
RUN pip3 install \
        meshio==4.4.*\
        numpy==1.20.*\
        matplotlib==3.4.*\
        scipy==1.6.*\
        pandas==1.2.*\
        pyyaml==5.4.*

# install dependencies
RUN apt-get -qq update && apt-get -qq install \
   software-properties-common \
   wget apt-transport-https
#    && \
#    rm -rf /var/lib/apt/lists/*

# disables cache (I think); this was in the DOckerfile.openfoam-adapter; not sure if needed
# ARG CACHEBUST

# install openFOAM
# see https://github.com/precice/openfoam-adapter/wiki/Notes-on-OpenFOAM#get-openfoam
# Add the signing key, add the repository, add ca-certificates (because...), update, install:
RUN wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | bash && \
        apt-get install -y ca-certificates
        
RUN apt-get update && \ 
        apt-get install -y openfoam2106-dev && \
        source /usr/lib/openfoam/openfoam2106/etc/bashrc

# Install preCICE
RUN wget https://github.com/precice/precice/releases/download/v2.3.0/libprecice2_2.3.0_focal.deb && \
        apt install -y ./libprecice2_2.3.0_focal.deb

# Install OpenFoam Adapter
RUN source /usr/lib/openfoam/openfoam2106/etc/bashrc && \
        wget https://github.com/precice/openfoam-adapter/releases/download/v1.0.0/openfoam-adapter_v1.0.0_OpenFOAMv1812-v2106.tar.gz && \
        tar -xzf openfoam-adapter_v1.0.0_OpenFOAMv1812-v2106.tar.gz && \
        cd openfoam-adapter_v1.0.0_OpenFOAMv1812-v2012 && \
        ./Allwmake

# Install pyprecice
RUN pip3 install pyprecice==2.3.*

RUN mkdir /home/workdir
WORKDIR /home/workdir

# modify user id and group
# see https://techflare.blog/permission-problems-in-bind-mount-in-docker-volume/
RUN apt-get -y install gosu 
COPY  ./docker/entrypoint.sh /
ENTRYPOINT [ "/entrypoint.sh" ]
