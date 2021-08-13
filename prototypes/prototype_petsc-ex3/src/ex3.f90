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
  use accsvec, only : create_vector, axpy
  use accs_utils, only : accs_free, update

  implicit none

  class(vector), allocatable :: u, b, ustar
  type(vector_init_data) :: vec_sizes

  integer(accs_int), parameter :: m = 100 ! XXX: temporary - this should be read from input
  
  integer(accs_err) :: ierr
  integer(accs_int) :: istart, iend
  real(accs_real) :: h
  integer(accs_int), parameter :: npe = 4 ! Points per element
  
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
  call update(u) ! Performs the parallel assembly
  
  !! Evaluate right-hand-side vector
  call eval_rhs(b)
  call update(b) ! Performs the parallel assembly
  
  !! Modify matrix and right-hand-side vector to apply Dirichlet boundary conditions
  !! Create linear solver & set options
  !! Solve linear system
  !! Check solution
  call create_vector(vec_sizes, ustar)
  call update(ustar) ! Performs the parallel assembly
  call axpy(-1.0_accs_real, u, ustar)
  
  !! Clean up
  call accs_free(u)
  call accs_free(b)
  call accs_free(ustar)
  
  call PetscFinalize(ierr)

contains

  subroutine eval_rhs(b)

    use accs_constants, only : add_mode
    use accs_types, only : vector_values
    use accs_utils, only : set_values
    
    type(vector), intent(inout) :: b

    integer(accs_int) :: i
    real(accs_real) :: x, y

    type(vector_values) :: val_dat

    val_dat%mode = add_mode
    allocate(val_dat%idx(npe))
    allocate(val_dat%val(npe))
    
    do i = istart, iend
       x = h * modulo(i, m)
       y = h * (i / m)

       call element_indices(i, val_dat%idx)
       call eval_element_rhs(x, y, h**2, val_dat%val)
       call set_values(val_dat, b)
    end do

    deallocate(val_dat%idx)
    deallocate(val_dat%val)
    
  end subroutine eval_rhs

  pure subroutine element_indices (i, idx)

    integer(accs_int), intent(in) :: i
    integer(accs_int), dimension(npe), intent(out) :: idx

    idx(1) = (m + 1) * (i / m) + modulo(i, m)
    idx(2) = idx(1) + 1
    idx(3) = idx(2) + (m + 1)
    idx(4) = idx(3) - 1
    
  end subroutine element_indices

  pure subroutine eval_element_rhs (x, y, H, r)
    !> @brief Apply forcing function
    
    real(accs_real), intent(in) :: x, y, H
    real(accs_real), dimension(npe), intent(out) :: r

    r(:) = 0.0
    
  end subroutine eval_element_rhs
  
end program ex3
