!> @brief Program file PETSc ex3
!>
!> @details Port of PETSc ksp/tutorial/ex3.c to ASiMoV-CCS style code - this is to help
!>          determine how to interface our design with PETSc.

program ex3

  !! External uses
#include <petsc/finclude/petsc.h>
  use petsc, only : PetscInitialize, PetscFinalize, PETSC_NULL_CHARACTER

  !! ASiMoV-CCS uses
  use accs_kinds, only : accs_real, accs_int, accs_err
  use accs_types, only : vector_init_data, vector
  use accsvec, only : create_vector
  use accs_utils, only : accs_free

  implicit none

  class(vector), allocatable :: u, b
  type(vector_init_data) :: vec_sizes

  integer(accs_int), parameter :: m = 100 ! XXX: temporary - this should be read from input
  
  integer(accs_err) :: ierr
  integer(accs_int) :: istart, iend
  real(accs_real) :: h
  
  !! Initialise program
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
  if (ierr /= 0) then
     print *, "Unable to initialise PETSc"
     stop
  end if
  istart = 0; iend = 0
  
  !! Create stiffness matrix
  !! Assemble matrix
  !! Create right-hand-side and solution vectors
  vec_sizes%nloc = -1
  vec_sizes%nglob = (m+1)**2
  call create_vector(vec_sizes, u)
  call create_vector(vec_sizes, b)

  !! Assemble right-hand-side vector
  call assemble_rhs(b)
  
  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  !! Create linear solver & set options
  !! Solve linear system
  !! Check solution

  !! Clean up
  call accs_free(u)
  call accs_free(b)
  
  call PetscFinalize(ierr)

contains

  subroutine assemble_rhs(b)

    type(vector), intent(inout) :: b

    integer :: i
    real(accs_real) :: x, y

    do i = istart, iend
       x = h * modulo(i, m)
       y = h * (i / m)
    end do
    
  end subroutine assemble_rhs
  
end program ex3
