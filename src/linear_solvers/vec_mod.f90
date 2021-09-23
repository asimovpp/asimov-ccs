!> @brief Module file vec.mod
!
!> @details An interface to operations on vector objects (creation, destruction, setting and
!!          getting, ...)

module vec

  use kinds, only : accs_real, accs_int
  use types, only : vector, vector_init_data
  
  implicit none

  private

  public :: create_vector
  public :: set_vector_values
  public :: update_vector
  public :: begin_update_vector
  public :: end_update_vector
  public :: axpy
  public :: norm

  interface
     
     !> @brief Interface to create a new vector object.
     !
     !> @param[in] vector_init_data vec_dat - Data structure containing the global and local sizes
     !!                                       of the vector, -1 is interpreted as unset. If both
     !!                                       are set the local size is used.
     !> @param[out] vector v - The vector returned allocated, but (potentially) uninitialised.
     module subroutine create_vector(vec_dat, v)
       type(vector_init_data), intent(in) :: vec_dat
       class(vector), allocatable, intent(out) :: v
     end subroutine

     !> @brief Interface to set values in a vector.
     !
     !> @param[in]     val_dat - contains the values, their indices and the mode to use for setting
     !!                          them.
     !> @param[in/out] v       - the vector.
     module subroutine set_vector_values(val_dat, v)
       class(*), intent(in) :: val_dat
       class(vector), intent(inout) :: v
     end subroutine

     !> @brief Interface to perform a parallel update of a vector.
     !
     !> @param[in/out] v - the vector
     module subroutine update_vector(v)
       class(vector), intent(inout) :: v
     end subroutine
     !> @brief Interface to begin a parallel update of a vector.
     !
     !> @param[in/out] v - the vector
     !
     !> @details Begins the parallel update to allow overlapping comms and compute.
     module subroutine begin_update_vector(v)
       class(vector), intent(inout) :: v
     end subroutine
     !> @brief Interface to end a parallel update of a vector.
     !
     !> @param[in/out] v - the vector
     !
     !> @details Ends the parallel update to allow overlapping comms and compute.
     module subroutine end_update_vector(v)
       class(vector), intent(inout) :: v
     end subroutine

     !> @brief Interface to perform the AXPY vector operation.
     !
     !> @param[in]     alpha - a scalar value
     !> @param[in]     x     - an input vector
     !> @param[in/out] y     - vector serving as input, overwritten with result
     !
     !> @details Performs the AXPY operation
     !!          y[i] = a * x[i] + y[i]
     module subroutine axpy(alpha, x, y)
       real(accs_real), intent(in) :: alpha
       class(vector), intent(in) :: x
       class(vector), intent(inout) :: y
     end subroutine

     !> @brief Interface to compute the norm of a vector
     !
     !> @param[in]  v         - the vector
     !> @param[in]  norm_type - which norm to compute? Currently supported is the 2 norm:
     !!                         norm_type=2.
     !> @param[out] n         - the computed norm returned as the result of the function
     !!                         call.
     module function norm(v, norm_type) result(n)
       class(vector), intent(in) :: v
       integer(accs_int), intent(in) :: norm_type
       real(accs_real) :: n
     end function
     
  end interface

contains
  
end module vec
