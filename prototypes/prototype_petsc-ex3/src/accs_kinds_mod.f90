!> @brief Module file accs_kinds.mod
!>
!> @details Defines kinds for use in primitive variable definitions in ASiMoV-CCS, e.g.
!>          integer(kind=accs_int) :: i

module accs_kinds

  use iso_fortran_env
#ifdef ACCS_PETSC
#include <petsc/finclude/petsc.h>
  use petsc
#endif

  implicit none

  private

#ifdef ACCS_PETSC
  PetscReal x
  PetscInt i
  PetscErrorCode ierr
#else
  real(kind=real32) :: x
  integer(kind=int32) :: i
  integer :: ierr
#endif
  
  integer, public, parameter :: accs_real = kind(x)
  integer, public, parameter :: accs_int = kind(i)
  integer, public, parameter :: accs_err = kind(ierr)
  
end module accs_kinds
