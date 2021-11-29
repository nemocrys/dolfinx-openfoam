# dolfinx-openfoam

Coupling of DOLFINx to OpenFOAM using the preCICE library.

## Docker

Build image with:

```shell
docker build -t nemocrys/dolfinx-openfoam -f ./docker/Dockerfile .
```

Execute container on Windows:

```shell
docker run -it --rm -v ${PWD}:/home/workdir nemocrys/dolfinx-openfoam bash
```

Execute container on Linux:

```shell
docker run -it --rm -v $PWD:/home/workdir -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/dolfinx-openfoam bash
```

## Acknowledgements

This code is based on the preCICE tutorial [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate).
