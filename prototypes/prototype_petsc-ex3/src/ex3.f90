!!! -*- mode: F90 -*-
!!! vim: set syntax=fortran:
!!!
!!!        FILE: ex3.f90
!!! DESCRIPTION: Port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code - this is to help
!!!              determine how to interface our design with PETSc.
!!!

program ex3

#include <petsc/finclude/petscksp.h>
  use petscksp
  
  implicit none

  integer :: ierr
  
  !! Initialise program
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  if (ierr /= 0) then
     print *, "Unable to initialise PETSc"
     stop
  end if
  
  !! Create stiffness matrix
  !! Assemble matrix
  !! Create right-hand-side and solution vectors
  !! Assemble right-hand-side vector
  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  !! Create linear solver & set options
  !! Solve linear system
  !! Check solution

  !! Clean up
  call PetscFinalize(ierr)
  
end program ex3
