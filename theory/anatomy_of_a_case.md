---
title: Anatomy of a case
---

A case in `ASiMoV-CCS` is given by the combination of the case source code and its corresponding
configuration file.
Each case is, in itself, a program; the case program defines the subroutines `get_init_flow`,
`get_init_mass_flux`, `postproc` and `eval_sources`.
These subroutines are used in setting the case-specific initial flow field, performing case-specific
in-situ analysis of the solution and computing source terms; each of these is described in more
detail below.

# Initialisation

The case program begins by performing initialisation, starting with creating the parallel
environment through a call to `initialise_parallel_environment`.
This is followed by reading the case configuration by calling `get_config`, returning the runtime
options for the case in a `ccs_options` object.
These runtime options are used to determine how to create the shared environment by splitting the
basic parallel environment.

After configuring the parallelism for the case and obtaining the runtime optins the remainder of the
initialisation concerns setting up the mesh and setting the intial values of the flow.
The call to `initialise_mesh` will either: read a mesh from a file or build a 2-D or 3-D square/cube
mesh using the internal mesh generator, based on the runtime options; in either case the returned
mesh is partitioned across the parallel environment.
The local mesh (plus any halo cells) determines the extent of the various flow fields which are
created by `initialise_fields`, creating both user-defined, and default-required, fields.
Following field creation, their initial values are set by `initialise_flow`; if the user chooses a
restart these values are read from a solution file, otherwise the `get_init_flow` and
`get_init_mass_flux` subroutines which are passed as arguments to `initialise_flow` are used to
compute the initial values.
The first of these, `get_init_flow`, is used to initialise the cell-centred fields and takes a cell
locator and field name as input arguments, returning the initial field value at a point as its
return value.
The `get_init_mass_flux` is similar, but acts on face values, specifically the mass-flux field.

# Running the simulation

The core of the program is the simulation itself.
The `run_solver` subroutine, as the name suggests, runs the solver, including the outer time loop
when the case is unsteady.
Besides the solver flow fields, parallel environment and runtime options, this accepts the `eval_sources` and `postproc` subroutines as arguments.
The `eval_sources` subroutine computes the linear and fixed source vectors that are added to the equation systems, returning the integrated values, _i.e._ \(S(i) = S_P V_P\).
The `postproc` subroutine is called at the end of the timestep, performing any analysis required by
the case.

# Finalising the case

After the solver completes the case performs cleanup: `dealloc_fluid_fields` frees the memory used
by the flow, `nullify_mesh_object` frees the mesh pointer and `cleanup_parallel_environment`
finalises the parallel environment.
Before finalising the parallel environment, calling `timer_print_all` will output the timings
reported by the instrumentation of the case and `ASiMoV-CCS`.
