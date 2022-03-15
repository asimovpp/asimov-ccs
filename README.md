# ASiMoV-CCS

## Prerequisites

- `MPI`
- `PETSc`
- `makedepf90` - version 2.9.0 reqired, source can be obtained at https://salsa.debian.org/science-team/makedepf90
- `adios2` - with `hdf5` support, https://adios2.readthedocs.io/
- `fortran-yaml-cpp` - https://github.com/Nicholaswogan/fortran-yaml-cpp
- `python` - with the `pyyaml` module (and optionally the `lit` module to run tests)



## Building

Set the following environment variables:

- `PETSC_DIR` to point to the PETSc install directory 
- `FYAML` to point to the root of your fortran-yaml-cpp build directory
- `ADIOS2` to point to the ADIOS2 install directory


With the prerequisites in place, ASiMoV-CCS can be built from the `src` directory with
```
make CMP=<compiler> all
```
where `<compiler>` is one of: `gnu`, `intel` or `cray`. `CMP` sets the compilation environment according to files in `src/build_tools/archs/Makefile.<compiler>`, which can be customised according to your environment. 

Tests can be run with
```
make CMP=<compiler> tests
```


## Configuring

The executable `ccs_app` that is built/linked is configured by `src/config.yaml`. The Fortran program is specified by `main:` where the available options currently reside in `src/case_setup` and are `poisson`, `tgv` and `scalar_advection`.

The build system automatically links the required modules, but submodules need to be specified in the configuration file. `base:` specifies a collection of basic submodules that work together; available options are `mpi` and `mpi_petsc`. Further customisation is available via the `options:` settings for cases where multiple implementations of the same functionality exist (none at the time of writing). 
All possible configuration settings can be found in `src/build_tools/config_mapping.yaml`. 


## Running
The generated application can be run as a normal MPI application, for example
```
mpirun -n 4 ccs_app
```

`ccs_app` accepts a number of runtime command line arguments, see `ccs_app --ccs_help` for details. 
If built with PETSc, the normal PETSc command line arguments can be passed to `ccs_app` as well.

