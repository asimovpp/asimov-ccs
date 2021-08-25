!> @brief Module file accs_types.mod
!>
!> @details Provides concrete types and bases of extensible types.

module accs_types

  use accs_kinds, only : accs_int, accs_real
  
  implicit none

  private

  type, public :: vector
  end type vector
  type, public :: matrix
  end type matrix

  type, public :: vector_init_data
     integer(accs_int) :: nglob, nloc
     integer :: comm
  end type vector_init_data

  type, public :: vector_values
     integer(accs_int), dimension(:), allocatable :: idx
     real(accs_real), dimension(:), allocatable :: val
     integer(accs_int) :: mode
  end type vector_values

  type, public :: matrix_init_data
     integer(accs_int) :: rloc, cloc
     integer(accs_int) :: rglob, cglob
     integer :: comm
  end type matrix_init_data

  type, public :: matrix_values
     integer(accs_int), dimension(:), allocatable :: rglob
     integer(accs_int), dimension(:), allocatable :: cglob
     real(accs_real), dimension(:), allocatable :: val
     integer(accs_int) :: mode
  end type matrix_values
  
end module accs_types
