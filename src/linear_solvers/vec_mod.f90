!> @brief Module file vec.mod
!
!> @details An interface to operations on vector objects (creation, destruction, setting and
!!          getting, ...)

module vec

  use kinds, only : accs_real, accs_int
  use types, only : vector, vector_init_data, vector_values
  use parallel_types, only: parallel_environment
  
  implicit none

  private

  public :: create_vector
  public :: set_vector_values
  public :: update_vector
  public :: begin_update_vector
  public :: end_update_vector
  public :: vec_axpy
  public :: vec_norm
  public :: pack_one_vector_element
  public :: initialise_vector
  public :: set_global_vector_size
  public :: set_local_vector_size
  public :: get_vector_data
  public :: reset_vector_data
  public :: vec_reciprocal

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
    !> @param[in,out] v       - the vector.
    module subroutine set_vector_values(val_dat, v)
      class(*), intent(in) :: val_dat
      class(vector), intent(inout) :: v
    end subroutine

    !> @brief Interface to perform a parallel update of a vector.
    !
    !> @param[in,out] v - the vector
    module subroutine update_vector(v)
      class(vector), intent(inout) :: v
    end subroutine

    !> @brief Interface to begin a parallel update of a vector.
    !
    !> @details Begins the parallel update to allow overlapping comms and compute.
    !
    !> @param[in,out] v - the vector
    module subroutine begin_update_vector(v)
      class(vector), intent(inout) :: v
    end subroutine
    
    !> @brief Interface to end a parallel update of a vector.
    !
    !> @details Ends the parallel update to allow overlapping comms and compute.
    !
    !> @param[in,out] v - the vector
    module subroutine end_update_vector(v)
      class(vector), intent(inout) :: v
    end subroutine end_update_vector

    !> @brief Interface to store one vector element and its index for later setting.
    !
    !> @param[in/out] val_dat - object storing the elements, their indices and mode to use when
    !!                          setting them.
    !> @param[in]     ent     - which entry in the index/element arrays to set?
    !> @oaram[in]     idx     - vector element index
    !> @param[in]     val     - vector element value
    !
    !> @details Stores a vector element and associated index for later setting, ensuring they are
    !!          set appropriately for the backend.
    module subroutine pack_one_vector_element(val_dat, ent, idx, val)
      type(vector_values), intent(inout) :: val_dat
      integer(accs_int), intent(in) :: ent
      integer(accs_int), intent(in) :: idx
      real(accs_real), intent(in) :: val
    end subroutine pack_one_vector_element

    !> @brief Interface to perform the AXPY vector operation.
    !
    !> @details Performs the AXPY operation
    !!          y[i] = a * x[i] + y[i]
    !
    !> @param[in]     alpha - a scalar value
    !> @param[in]     x     - an input vector
    !> @param[in,out] y     - vector serving as input, overwritten with result
    module subroutine vec_axpy(alpha, x, y)
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
    module function vec_norm(v, norm_type) result(n)
      class(vector), intent(in) :: v
      integer(accs_int), intent(in) :: norm_type
      real(accs_real) :: n
    end function
    
    !> @brief Constructor for default vector values
    !
    !> param[in/out] vector_descriptor - the initialised vector values
    module subroutine initialise_vector(vector_descriptor)
      type(vector_init_data), intent(inout) :: vector_descriptor
    end subroutine initialise_vector

    !> @brief Setter for global vector size
    !
    !> param[in/out] vector_descriptor - the vector data object
    !> param[in] size                  - the global vector size
    !> param[in] par_env               - the parallel environment 
    !!                                   where the vector resides
    module subroutine set_global_vector_size(vector_descriptor, size, par_env)
      type(vector_init_data), intent(inout) :: vector_descriptor
      integer(accs_int), intent(in) :: size
      class(parallel_environment), allocatable, target, intent(in) :: par_env
    end subroutine

    !> @brief Setter for local vector size
    !
    !> param[in/out] vvector_descriptor - the vector data object
    !> param[in] size                   - the local vector size
    !> param[in] par_env                - the parallel environment 
    !!                                    where the vector resides
    module subroutine set_local_vector_size(vector_descriptor, size, par_env)
      type(vector_init_data), intent(inout) :: vector_descriptor
      integer(accs_int), intent(in) :: size
      class(parallel_environment), allocatable, target, intent(in) :: par_env
    end subroutine

    !> @brief Gets the data in a given vector
    !
    !> @param[in] vec   - the vector to get data from
    !> @param[in] array - an array to store the data in
    module subroutine get_vector_data(vec, array)
      class(vector), intent(in) :: vec
      real(accs_real), dimension(:), intent(out) :: array
    end subroutine get_vector_data

    !> @brief Resets the vector data if required for further processing
    !
    !> @param[in] vec   - the vector to reset
    !> @param[in] array - an array storing the data
    module subroutine reset_vector_data(vec, array)
      class(vector), intent(in) :: vec
      real(accs_real), dimension(:), intent(in) :: array
    end subroutine reset_vector_data

    !> @brief Replace each component of a vector by its reciprocal
    !
    !> @param[in]  vec - the vector
    !> @param[out] vec - the vector reciprocal
    module subroutine vec_reciprocal(vec)
      class(vector), intent(inout) :: vec
    end subroutine vec_reciprocal

    end interface
  
end module vec
