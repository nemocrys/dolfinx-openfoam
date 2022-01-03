# dolfinx-openfoam

Coupling of DOLFINx to OpenFOAM using the preCICE library.

## Docker

The fenicsx-adapter repository has to be in the same directory with this repository.

Build image with:

```shell
docker build -t nemocrys/dolfinx-openfoam -f ./docker/Dockerfile .
```

Execute container on Windows:

```shell
docker run -it --rm -v ${PWD}:/home/workdir -v ${PWD}/../fenicsx-adapter:/home/fenicsx-adapter nemocrys/dolfinx-openfoam bash
```

Execute container on Linux:

```shell
docker run -it --rm -v $PWD:/home/workdir -v $PWD/../fenicsx-adapter:/home/fenicsx-adapter -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/dolfinx-openfoam bash
```

## Acknowledgements

This project is based on the preCICE tutorial [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate and on the [dofinx dockerfile](https://github.com/FEniCS/dolfinx/blob/main/docker/Dockerfile).
