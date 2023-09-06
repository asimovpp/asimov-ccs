# Changelog

## v0.4

- Added reading of external mesh files
- Updated partitioning, including option to use ParMETIS
- Added local cells reordering
- Added support for scalar transport
- Added gamma and linear upwind schemes
- Replaced globally-sized arrays with shared memory to address memory capacity issues
- Added the Backwards Facing Step testcase
- Added verification tests for space and time discretisation and Poiseuille flow test case
- Improved user-friendliness of configuration and run


## v0.3

- Implemented 2D and 3D Taylor Green Vortex (TGV) use cases
- Added first and second order timestepping
- Added residual, kinetic energy and enstrophy calculation and logging
- Added convergence testing
- Split mesh into topology and geometry objects
- Added partitioning
- Implemented boundary conditions: periodic, dirichlet, neumann, extrapolate and wall
- Added outputting of solution in a ParaView readable format


## v0.2

- Lid Driven Cavity case
- SIMPLE algorithm
- ADIOS2 and YAML I/O
- LIT test suite
- Build system improvements
- Developer guide
- Improvements to back end support for matrices, vectors and fields
- Accessors for the mesh


## v0.1.1

- Various code tidying


## v0.1

- Poisson solver


