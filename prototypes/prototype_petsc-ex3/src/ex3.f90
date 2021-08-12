!> @brief Program file PETSc ex3
!>
!> @details Port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code - this is to help
!>          determine how to interface our design with PETSc.

program ex3

  !! External uses
#include <petsc/finclude/petsc.h>
  use petsc, only : PetscInitialize, PetscFinalize, PETSC_NULL_CHARACTER

  !! ASiMoV-CCS uses
  use accsvec, only : vector, vector_init_data, create_vector
  
  implicit none

  class(vector), allocatable :: u
  type(vector_init_data) :: vec_sizes

  integer, parameter :: m = 100 ! XXX: temporary - this should be read from input
  
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
  vec_sizes%nloc = -1
  vec_sizes%nglob = (m+1)**2
  call create_vector(vec_sizes, u)
  !! Assemble right-hand-side vector
  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  !! Create linear solver & set options
  !! Solve linear system
  !! Check solution

  !! Clean up
  call PetscFinalize(ierr)
  
end program ex3
