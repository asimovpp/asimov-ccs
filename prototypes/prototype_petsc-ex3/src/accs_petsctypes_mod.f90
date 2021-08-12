!> @brief Module file accs_petsctypes.mod
!>
!> @details Provides petsc-extended types.

module accs_petsctypes

#include <petsc/finclude/petscvec.h>
  use petscvec, only : tVec

  use accs_types, only : vector
  
  implicit none

  private
  
  type, public, extends(vector) :: vector_petsc
     type(tVec) :: v
  end type vector_petsc
  
end module accs_petsctypes
