!> @brief Module file kinds.mod
!>
!> @details Defines kinds for use in primitive variable definitions 
!! in ASiMoV-CCS, e.g. integer(kind=accs_int) :: i

module kinds

  use iso_fortran_env

#ifdef ACCS_PETSC
#include <petsc/finclude/petscsys.h>
#endif

  implicit none

#ifdef ACCS_PETSC
  PetscReal x
  PetscInt i
  PetscErrorCode ierr
#else
  real(kind=real32) :: x
  integer(kind=int32) :: i
  integer :: ierr
#endif
  
  integer, public, parameter :: accs_real = kind(x)   !> Real kind to be used in ASiMoV-CCS
  integer, public, parameter :: accs_int = kind(i)    !> Integer kind to be used in ASiMoV-CCS
  integer, public, parameter :: accs_err = kind(ierr) !> Error kind to be used in ASiMoV-CCS

  private
  
end module kinds
