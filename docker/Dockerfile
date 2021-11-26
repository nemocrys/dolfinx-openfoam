# USAGE INSTSRUCTIONS
# 
# Build image (run this in the main directory (where .git folder is)):
# docker build -t nemocrys/dolfinx-openfoam -f ./docker/Dockerfile .
#
# Run image:
# docker run -it --rm nemocrys/dolfinx-openfoam:vX.X.X bash
#
# Run image and map current directory to home/workdir on Linux:
# docker run -it --rm -v $PWD:/home/workdir -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/dolfinx-openfoam bash
# and on Windows:
# docker run -it --rm -v ${PWD}:/home/workdir nemocrys/dolfinx-openfoam bash

# FROM ubuntu:20.04
FROM dolfinx/dolfinx:v0.2.0a0

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

# install openFOAM

# Install preCICE



RUN mkdir /home/workdir
WORKDIR /home/workdir

# modify user id and group
# see https://techflare.blog/permission-problems-in-bind-mount-in-docker-volume/
COPY  ./docker/entrypoint.sh /
RUN apt-get -y install gosu 
ENTRYPOINT [ "/entrypoint.sh" ]
