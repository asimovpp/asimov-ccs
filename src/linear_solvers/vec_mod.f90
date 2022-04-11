!> @brief Module file vec.mod
!
!> @details An interface to operations on vector objects (creation, destruction, setting and
!!          getting, ...)

module vec

  use kinds, only : accs_real, accs_int
  use types, only : mesh, vector, vector_init_data, vector_values
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
  public :: set_vector_size
  public :: get_vector_data
  public :: restore_vector_data
  public :: set_vector_location
  public :: vec_reciprocal
  public :: zero_vector
  public :: mult_vec_vec
  public :: scale_vec
  ! public :: vec_view
  
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
    module subroutine pack_one_vector_element(ent, idx, val, val_dat)
      integer(accs_int), intent(in) :: ent
      integer(accs_int), intent(in) :: idx
      real(accs_real), intent(in) :: val
      type(vector_values), intent(inout) :: val_dat
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

    !> @brief Setter for vector size
    !
    !> param[in]     par_env           - the parallel environment 
    !!                                   where the vector resides
    !> param[in]     geometry          - the mesh - contains the
    !!                                   information to set the
    !!                                   vector size
    !> param[in/out] vector_descriptor - the vector data object
    module subroutine set_vector_size(par_env, geometry, vector_descriptor)
      class(parallel_environment), allocatable, target, intent(in) :: par_env
      class(mesh), target, intent(in) :: geometry
      type(vector_init_data), intent(inout) :: vector_descriptor
    end subroutine set_vector_size

    !> @brief Gets the data in a given vector
    !
    !> @param[in] vec   - the vector to get data from
    !> @param[in] array - an array to store the data in
    module subroutine get_vector_data(vec, array)
      class(vector), intent(in) :: vec
      real(accs_real), dimension(:), pointer, intent(out) :: array
    end subroutine get_vector_data

    !> @brief Resets the vector data if required for further processing
    !
    !> @param[in] vec   - the vector to reset
    !> @param[in] array - the array containing the data to restore
    module subroutine restore_vector_data(vec, array)
      class(vector), intent(in) :: vec
      real(accs_real), dimension(:), pointer, intent(in) :: array
    end subroutine restore_vector_data

    !> @brief Set vector values to be located at either cell-centre or face
    !
    module subroutine set_vector_location(loc, vector_descriptor)
      integer(accs_int), intent(in) :: loc
      type(vector_init_data), intent(inout) :: vector_descriptor
    end subroutine set_vector_location

    !> @brief Replace each component of a vector by its reciprocal
    !
    !> @param[in]  vec - the vector
    !> @param[out] vec - the vector reciprocal
    module subroutine vec_reciprocal(vec)
      class(vector), intent(inout) :: vec
    end subroutine vec_reciprocal

    module subroutine zero_vector(vec)
      class(vector), intent(inout) :: vec
    end subroutine zero_vector

    module subroutine mult_vec_vec(a, b)
      class(vector), intent(in) :: a
      class(vector), intent(inout) :: b
    end subroutine mult_vec_vec

    module subroutine scale_vec(alpha, v)
      real(accs_real), intent(in) :: alpha
      class(vector), intent(inout) :: v
    end subroutine scale_vec
    
    ! module subroutine vec_view(vec_dat, vec)
    !   type(vector_init_data), intent(in) :: vec_dat
    !   class(vector), intent(in) :: vec
    ! end subroutine vec_view

  end interface
  
end module vec
