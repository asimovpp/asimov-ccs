# Changelog

## v0.5

### Added
- Added support for arbitrary cell types and improved discretisation on distorted, unstructured and non-Cartesian meshes.
- Added high-level wall, inflow and outflow boundary conditions, named mesh boundaries and utilities for locating boundary faces.
- Added non-time-accurate restart and steady-to-unsteady restart support.
- Added user-defined fixed and linearised source terms for momentum and scalar-transport equations.
- Added variable-density and variable-viscosity fields, including runtime configuration of reference properties.
- Added the Sandia flow case.
- Added Z-order space-filling-curve and no-reordering options.
- Added runtime reporting of partition quality and matrix bandwidth, plus configurable ParHIP imbalance and partitioning mode.
- Added per-variable configuration of linear solvers, preconditioners, relaxation factors, gradient relaxation, residual targets and L2/L-infinity residual norms.
- Added maximum and average CFL reporting.
- Added per-timestep solution files containing timestep and simulation-time metadata, and XDMF metadata for generated meshes.
- Added a common profiling interface with timer, Caliper and LIKWID backends, Caliper annotations and CSV timing output.
- Added graceful SIGTERM handling.
- Added support for building CCS as the static library libccs.a.

### Changed
- Introduced a generic core configuration, mesh/flow initialisation and solver driver, substantially simplifying case implementations.
- Reworked finite-volume assembly around testable advection, diffusion and transient kernels, equation objects and equation payloads.
- Made field creation and field properties dynamic, removing static field identifiers and output lists.
- Improved SIMPLE pressure–velocity coupling, pressure-Poisson consistency and pressure-correction source-imbalance handling.
- Improved gradient calculation, under-relaxation and non-orthogonal, eccentricity, pressure-gradient and mass-flux corrections.
- Reduced partitioning memory requirements by deriving distribution data locally and removing the global partition array.
- Added a general parallel data shuffle and substantially optimised face connectivity, adjacency construction, field interpolation and gradient updates.
- Added read-only views of field and old-time data and OpenMP parallelism in thread-safe kernels, gradients, residual normalisation and source initialisation.
- Changed PETSc reordering to the default; RCM remains available through the external RCM-f90 library.
- Made ParHIP and ParMETIS individually optional at build time, while continuing to require at least one partitioner.
- Select output filename extensions from the configured ADIOS2 engine instead of hard-coding them.
- Abstracted application logging and expanded runtime configuration summaries.

### Fixed
- Corrected Gamma and linear-upwind advection interpolation and coefficient calculations.
- Corrected deferred diffusion corrections and diffusion boundary-condition handling.
- Corrected momentum and scalar residual calculations.
- Corrected pressure-correction behaviour on non-orthogonal meshes.
- Fixed non-orthogonal and eccentricity corrections degradation on distorted meshes.
- Prevented out-of-bounds interpolation and improved face-interpolation-factor calculation.
- Added closed-mesh validation and diagnostics when reading meshes.
- Added validation and clearer errors for missing or unreadable ADIOS2 fields and files.
- Fixed scalar-transport, Backward-Facing Step, Taylor–Green Vortex and restart-related regressions.
- Fixed field and mesh deallocation, shared-memory allocation and other memory leaks.
- Added explicit destruction paths for compilers that do not invoke finalisers reliably on polymorphic objects (notably Cray).
- Improved runtime viscosity and reference-value parsing.

### Build, Compatibility and Quality
- Added PETSc 3.24 support and dropped support for PETSc versions earlier than 3.23.
- Updated for the ADIOS2 2.9 interface.
- Added single-precision support consistent with the PETSc build.
- Added an LLVM build configuration and improved GNU, Intel/ifx, Cray and macOS compatibility.
- Added ReFrame-based HPC regression testing for ARCHER2 and Cirrus.
- Expanded unit, parallel mesh, restart, kernel, boundary-condition, distorted-mesh, unstructured-discretisation, verification and performance testing.
- Added theory and user documentation covering case structure, finite-volume kernels, transient terms and source terms.
- Improved linting, formatting, CI dependency alignment and issue-reporting infrastructure.
- Release builds against external dependencies tag [2026.09](https://github.com/asimovpp/ccs-dependencies/releases/tag/2026.09).

### Removed
- Removed the superseded Scalar Advection example case.
- Removed bundled RCM implementation code in favour of the external RCM-f90 dependency.
- Removed obsolete field wrappers, static field IDs, output-list state and unused mesh connectivity structures.

---

## Acknowledgements

Contributions from the ASiMoV development team and community members.

For detailed documentation, see the [Developer Guide](dev_guide/ccs_dev_guide.tex) and [Theory Documentation](theory/).

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


