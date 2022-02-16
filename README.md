# dolfinx-openfoam

Coupling of DOLFINx to OpenFOAM using the [preCICE library](https://precice.org/).

The project is developed and maintained by the [**Model experiments group**](https://www.ikz-berlin.de/en/research/materials-science/section-fundamental-description#c486) at the Leibniz Institute for Crystal Growth (IKZ).

### Referencing
If you use this code in your research, please cite our article (available with open access):

> A. Enders-Seidlitz, J. Pal, and K. Dadzis, Model experiments for Czochralski crystal growth processes using inductive and resistive heating *IOP Conference Series: Materials Science and Engineering*, 1223 (2022) 012003. https://doi.org/10.1088/1757-899X/1223/1/012003.


## Docker configuration

The [fenicsx-adapter repository](https://github.com/precice/fenicsx-adapter) has to be in the same directory with this repository (don't forget to select the right branch). Then, if you execute it with the command bellow, the adapter will be mounted into the container and installed in editable mode.

```
Parent_directory/
├── dolfinx-openfoam
│   └── ...
├── fenicsx-adapter
│   └── ...
...
```

Build image with (execute in the `dolfinx-openfoam` directory):

```shell
docker build -t nemocrys/dolfinx-openfoam -f ./docker/Dockerfile .
```

Or [get it from dockerhub](https://hub.docker.com/r/nemocrys/dolfinx-openfoam):

```shell
docker pull nemocrys/dolfinx-openfoam
```

Execute container (in the `dolfinx-openfoam` directory) on Windows:

```shell
docker run -it --rm -v ${PWD}:/home/workdir -v ${PWD}/../fenicsx-adapter:/home/fenicsx-adapter nemocrys/dolfinx-openfoam bash
```

Execute container on Linux:

```shell
docker run -it --rm -v $PWD:/home/workdir -v $PWD/../fenicsx-adapter:/home/fenicsx-adapter -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/dolfinx-openfoam bash
```


## Flow over heated plate

The coupling is tested with the preCICE tutorial case [flow over heated plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate). Execute it with:

```shell
cd solid-fenics
python3 solid.py &> run.log &
cd ../fluid-openfoam
./run.sh &> run.log &
```


## Acknowledgements

This project is based on the preCICE tutorial [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate) and on the [dofinx dockerfile](https://github.com/FEniCS/dolfinx/blob/main/docker/Dockerfile).


[This project](https://www.researchgate.net/project/NEMOCRYS-Next-Generation-Multiphysical-Models-for-Crystal-Growth-Processes) has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 851768).

<img src="https://raw.githubusercontent.com/nemocrys/pyelmer/master/EU-ERC.png">
