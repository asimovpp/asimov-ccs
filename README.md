# ASiMoV-CCS

ASiMoV-CCS is a CFD and combustion code designed to scale to large numbers of cores. It follows a "separation of concerns' design that separates interfaces from implementations, and physics from parallelisation. There is a distinction between "user" and "back-end" code, i.e. case specific code (such as the setup of a particular test) is the user code, while the core functionality provided by ASiMoV-CCS is the back-end code.

ASiMoV-CCS is implemented in a modular fashion by separating the interface declarations contained within modules from their implementation in submodules. As a result, it is possible to implement multiple physics models and parallelisation strategies by writing separate submodules, each providing a distinct solution.

## Prerequisites

- `MPI`
- `PETSc`
- `makedepf90` - version 2.9.0 required, source can be obtained at https://salsa.debian.org/science-team/makedepf90
- `adios2` - with `hdf5` support, https://adios2.readthedocs.io/
- `fortran-yaml-cpp` - https://github.com/Nicholaswogan/fortran-yaml-cpp
- `python` - with the `pyyaml` module (and optionally the `lit` module to run tests)
- `ParHIP` - https://github.com/KaHIP/KaHIP



## Building

Set the following environment variables:

- `PETSC_DIR` to point to the PETSc install directory 
- `FYAML` to point to the root of your fortran-yaml-cpp build directory
- `ADIOS2` to point to the ADIOS2 install directory


With the prerequisites in place, ASiMoV-CCS can be built from the root directory with
```
make CMP=<compiler> all
```
where `<compiler>` is one of: `gnu`, `intel` or `cray`. `CMP` sets the compilation environment according to files in `build_tools/archs/Makefile.<compiler>`, which can be customised according to your environment. 

Tests can be run with
```
make CMP=<compiler> tests
```


## Configuring

The executable `ccs_app` that is built/linked is configured by `config.yaml`. The Fortran program is specified by `main:` where the available options currently reside in `case_setup` and are `poisson`, `scalar_advection` and `ldc`.

The build system automatically links the required modules, but submodules need to be specified in the configuration file. `base:` specifies a collection of basic submodules that work together; available options are `mpi` and `mpi_petsc`. Further customisation is available via the `options:` settings for cases where multiple implementations of the same functionality exist (none at the time of writing). 
All possible configuration settings can be found in `build_tools/config_mapping.yaml`. 


## Running
The generated application can be run as a normal MPI application, for example
```
mpirun -n 4 ./ccs_app
```

`ccs_app` accepts a number of runtime command line arguments, see `ccs_app --ccs_help` for details. 
If built with PETSc, the normal PETSc command line arguments can be passed to `ccs_app` as well.

### Running the Lid Driven Cavity case
Change into the directory when the Lid Driven Cavity configuration file (`LidDrivenCavity_config.yaml`) resides, i.e.
```
cd case_setup/LidDrivenCavity
```
You can change the values in the configuration file to customise your setup. For example, by default the number of iterations is set to `10`, but you might want to change this to something like `1000`. 

The command line option `--ccs_m` allows you to set the size of the mesh (the default is `50x50`). You can run the case as follows:
```
mpirun -n 4 ../../ccs_app --ccs_m 129 --ccs_case LidDrivenCavity
```
