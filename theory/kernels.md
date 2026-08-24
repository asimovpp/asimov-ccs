---
title: Implementation of discretisation kernels
---

The proposal is to replace the current subroutines with `kernel` objects, representing the different
discretisation operations that can be applied at a point, for instance `advection_kernel`,
`diffusion_kernel`, `transient_kernel`, etc. 
These can then be combined into `equation` objects which build discrete representation of physics at
a cell, e.g. `momentum_equation`, `scalar_equation`, `poisson_equation`. 
Finally the equations all fit within the framework of an `equation_system` which handles evaluating
the equation at a cell and storing its coefficients into a row of the linear system.
In simple summary:
- `kernel` == discretisation
- `equation` == physics
- `equation_system` == storage

We want the kernels to be testable, and we want the application of equations to be high performance.
This means working with simple data structures. 
Depending on their purpose, different kernels require different data: `advection` and `diffusion`
kernels will require the field values (and gradients) straddling a face, `transient` kernels require
one or more `old` values depending on the scheme, etc. 
The `equation` object should therefore build a `payload` - a simple 1D array of the values required
by the kernels, and the offsets into this so it can retrieve the data and hand this off to kernels
as required. 
In summary:
- `kernel`s should be queriable for the `width` of data they require (e.g. stencil `width` = 1 for
  advection/diffusion -> 1 neighbour depth, stencil `width` = 2 for transient -> 2 old values)
- `equation`s should query their kernels for this data to define the `payload`
- `equation`s gather the `flow` data into the `payload` for each equation application
- `equation`s extract the `payload` data required by each `kernel` invocation
in this way the kernels remain ignorant of the physics they are used for (e.g. the `gradient` kernel
doesn't know about the pressure gradient source term in the momentum equations).

Design
======

Kernel
------

A `kernel` should fit the prototype
```
type :: kernel
contains
  procedure :: eval_coeffs   ! Returns the implicit contributions to the discretisation
  procedure :: eval_explicit ! Returns the explicit contributions to the discretisation
  procedure :: get_width     ! Returns the stencil width
  procedure :: get_order     ! Returns the theoretical order of the discretisation
end type
```

Equation
--------

An `equation` should fit the prototype
```
type :: equation
  procedure :: apply  ! Applies the equation at a cell, evaluating the coefficients/RHS
  procedure :: init   ! Defines the `payload`
  procedure :: gather ! Builds the `payload`
```

Equation System
---------------
The `equation_system` could be a subroutine
```
subroutine equation_system(flow, eq, M, b)

  call eq%init() ! This might occur externally, as part of `equation` declaration
  do i = 1, n
    call eq%gather(flow, i)          ! Gathers data for the `payload`
    call eq%apply(i, coeffs, rhs) ! Evaluates the equation coefficients/rhs
    
    ! Push coeffs/rhs -> row i of M/b
  end do

end subroutine
```

Equation object contain payload object to allow for efficient iteration over the matrix rows
