!> @brief Module file accsvec.mod
!>
!> @details An interface to operations on vector objects (creation, destruction, setting and
!>          getting, ...)

module accsvec

  use accs_kinds, only : accs_real, accs_int
  use accs_types, only : vector, vector_init_data
  
  implicit none

  private

  public :: create_vector, set_vector_values, update_vector, begin_update_vector, end_update_vector, axpy, norm

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

     module subroutine set_vector_values(val_dat, v)

       class(*), intent(in) :: val_dat
       class(vector), intent(inout) :: v

     end subroutine

     module subroutine update_vector(v)
       class(vector), intent(inout) :: v
     end subroutine
     module subroutine begin_update_vector(v)
       class(vector), intent(inout) :: v
     end subroutine
     module subroutine end_update_vector(v)
       class(vector), intent(inout) :: v
     end subroutine

     module subroutine axpy(alpha, x, y)
       real(accs_real), intent(in) :: alpha
       class(vector), intent(in) :: x
       class(vector), intent(inout) :: y
     end subroutine

     module real(accs_real) function norm(v, norm_type)
       class(vector), intent(in) :: v
       integer(accs_int), intent(in) :: norm_type
     end function
     
  end interface

contains
  
end module accsvec
