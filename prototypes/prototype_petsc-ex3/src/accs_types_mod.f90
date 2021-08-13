!> @brief Module file accs_types.mod
!>
!> @details Provides concrete types and bases of extensible types.

module accs_types

  use accs_kinds, only : accs_int, accs_real
  
  implicit none

  private

  type, public :: vector
  end type vector

  type, public :: vector_init_data
     integer(accs_int) :: nglob, nloc
  end type vector_init_data

  type, public :: vector_values
     integer(accs_int), dimension(:), allocatable :: idx
     real(accs_real), dimension(:), allocatable :: val
     integer(accs_int) :: mode
  end type vector_values
  
end module accs_types
