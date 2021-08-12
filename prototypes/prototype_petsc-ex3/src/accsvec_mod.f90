!> @brief Module file accsvec.mod
!>
!> @details An interface to operations on vector objects (creation, destruction, setting and
!>          getting, ...)

module accsvec

  use accs_types, only : vector, vector_init_data
  
  implicit none

  private

  interface
     module subroutine create_vector(vec_dat, v)
       !> @brief Creates a vector given the local or global size.
       !>
       !> @param[in] vector_init_data vec_dat - Data structure containing the global and local sizes
       !>                                       of the vector, -1 is interpreted as unset. If both
       !>                                       are set the local size is used.
       !> @param[out] vector v - The vector returned allocated, but (potentially) uninitialised.
       type(vector_init_data), intent(in) :: vec_dat
       class(vector), allocatable, intent(out) :: v
     end subroutine
  end interface

  public :: create_vector
  
end module accsvec
