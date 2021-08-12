!> @brief Module file accs_types.mod
!>
!> @details Provides concrete types and bases of extensible types.

module accs_types

  implicit none

  private

  type, public :: vector
  end type vector

  type, public :: vector_init_data
     integer :: nglob, nloc
  end type vector_init_data
  
end module accs_types
