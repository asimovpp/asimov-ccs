# ASiMoV-CCS

ASiMoV-CCS is a CFD and combustion code designed to scale to large numbers of cores. It follows a "separation of concerns' design that separates interfaces from implementations, and physics from parallelisation. There is a distinction between "user" and "back-end" code, i.e. case specific code (such as the setup of a particular test) is the user code, while the core functionality provided by ASiMoV-CCS is the back-end code.

ASiMoV-CCS is implemented in a modular fashion by separating the interface declarations contained within modules from their implementation in submodules. As a result, it is possible to implement multiple physics models and parallelisation strategies by writing separate submodules, each providing a distinct solution.

## Prerequisites

- `MPI`
- `PETSc` - version 3.17 or newer required
- `makedepf90` - version 2.9.0 required, source can be obtained at https://salsa.debian.org/science-team/makedepf90
- `adios2` - with `hdf5` support, https://adios2.readthedocs.io/
- `fortran-yaml-c` - https://github.com/Nicholaswogan/fortran-yaml-c
- `python` - with the `pyyaml` module (and optionally the `lit` module to run tests)
- `ParHIP` - https://github.com/KaHIP/KaHIP



## Building

Set the following environment variables:

- `PETSC_DIR` to point to the PETSc install directory 
- `FYAMLC` to point to the root of your fortran-yaml-c build directory
- `ADIOS2` to point to the ADIOS2 install directory
- `PARHIP` to point to the root of the ParHIP install directory


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

The executable `ccs_app` that is built/linked is configured by `config.yaml`. The Fortran program is specified by `main:` where the available options currently reside in `case_setup` and are `poisson`, `scalar_advection`, `ldc`, `tgv` and `tgv2d` (for 3-D and 2-D Taylor-Green Vortex, respectively).

The build system automatically links the required modules, but submodules need to be specified in the configuration file. `base:` specifies a collection of basic submodules that work together; available options are `mpi` and `mpi_petsc`. Further customisation is available via the `options:` settings for cases where multiple implementations of the same functionality exist (none at the time of writing). 
All possible configuration settings can be found in `build_tools/config_mapping.yaml`. 


## Running
The generated application can be run as a normal MPI application, for example
```
mpirun -n 4 ./ccs_app
```

`ccs_app` accepts a number of runtime command line arguments, see `ccs_app --ccs_help` for details. 
If built with PETSc, the normal PETSc command line arguments can be passed to `ccs_app` as well, in particular if you observe stagnating residuals being printed to `stdout` you might want to try tighter tolerances for the linear solvers by passing `-ksp_atol` and `-ksp_rtol` at runtime, e.g.
```
mpirun -n ${NP} /path/to/ccs_app <ccs_options> -ksp_atol 1.0e-16 -ksp_rtol 1.0e-50
```
should give behaviour similar to a direct solver.

### Running the Taylor Green Vortex case
Change into the directory where the Taylor Green Vortex configuration file (`TaylorGreenVortex2D_config.yaml`) resides, i.e.
```
cd src/case_setup/TaylorGreenVortex
```
You can change the values in the configuration file to customise your setup. 

The command line option `--ccs_m` allows you to set the size of the mesh (the default is `50x50`). You can run the case as follows:
```
mpirun -n 4 ../../../ccs_app --ccs_m 100 --ccs_case TaylorGreenVortex2D
```

To run the 3D TGV case, you need to change `config.yaml` to state `main: tgv` instead of `main:tgv2d`, rebuild and use the run line:
```
mpirun -n 4 ../../../ccs_app --ccs_m 32 --ccs_case TaylorGreenVortex3D
```
The default 3D problem size is `16x16x16`.

In both cases time evolution of the kinetic energy and enstrophy are written to files `tgv-ke.log` and `tgv-ens.log` (`tgv2d-*.log` in the 2-D case).
In the 2-D case there also exists an analytical solution which is used to compute the RMS error in the velocity solution, note that this requires setting the viscosity `mu` in `tgv2d.f90` to match the value of `diffusion_factor` in `fv_common.f90`.
Using this error log the convergence rate of the central and upwind schemes can be confirmed by systematic refinement of the grid, noting that for the central scheme the Peclet number must be less than 2
```
Pe = rho U dx / mu < 2
```
noting that here `rho` and `U` are both 1, the grid spacing `dx` is equal to the domain length (PI) divided by the number of cells in each direction (default 50).

To change the simulated time of either case, edit the appropriate program file and either set the timestep directly or the `CFL` number.

### Running the Lid Driven Cavity case

To run the LDC case, you need to change `config.yaml` to state `main: ldc` and rebuild.

Change into the directory where the Lid Driven Cavity configuration file (`LidDrivenCavity_config.yaml`) resides, i.e.
```
cd src/case_setup/LidDrivenCavity
```
You can change the values in the configuration file to customise your setup. For example, by default the number of iterations is set to `10`, but you might want to change this to something like `1000`. 

The command line option `--ccs_m` allows you to set the size of the mesh (the default is `50x50`). You can run the case as follows:
```
mpirun -n 4 ../../../ccs_app --ccs_m 129 --ccs_case LidDrivenCavity
```

The results can be plotted against Ghia's reference data by running the `plot-structured.py` scrip from the results directory
```
python ../../../scripts/plot-structured.py
```
this requires installation of `h5py`, which can be installed via `pip install h5py`.

## Further documentation

### Code documentation
CCS uses FORD for code documentation. You can install FORD using `pip install ford`. Documentation can then be generated using `make ford`. This will generate HTML documentation in the `doc` directory with entry point `index.html`.

See https://github.com/Fortran-FOSS-Programmers/ford and https://forddocs.readthedocs.io/en/latest/ for more information on FORD.

### Developer documentation
A developer and style guide can be generated using `make dev_guide`. This requires `latex`. The output can be found in `dev_guide/ccs_dev_guide.pdf`.

Documentation describing the build system, testing framework and linting in CCS can be found [here](build_tools/build_system_readme.md).
