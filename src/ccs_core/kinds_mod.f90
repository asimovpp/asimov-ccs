!v Module file kinds.mod
!
!  Defines kinds for use in primitive variable definitions
!  in ASiMoV-CCS, e.g. integer(kind=ccs_int) :: i

module kinds

  use iso_fortran_env
  use mpi

#ifdef ACCS_PETSC
#include <petsc/finclude/petscsys.h>
#include <petscconf.h>
  use petscsys
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

  integer, public, parameter :: ccs_real = kind(x)       !< Real kind to be used in ASiMoV-CCS
  integer, public, parameter :: ccs_int = kind(i)        !< Integer kind to be used in ASiMoV-CCS
  integer, public, parameter :: ccs_long = kind(1_int64) !< Long integer kind to be used in ASiMoV-CCS
  integer, public, parameter :: ccs_err = kind(ierr)     !< Error kind to be used in ASiMoV-CCS

#if defined(PETSC_USE_REAL_SINGLE)
  integer, public, parameter :: CCS_MPI_PRECISION = MPI_REAL
  character(len=6), public, parameter :: CCS_PRECISION_STR = 'single'
#else
  integer, public, parameter :: CCS_MPI_PRECISION = MPI_DOUBLE_PRECISION
  character(len=6), public, parameter :: CCS_PRECISION_STR = 'double'
#endif

  private

end module kinds
