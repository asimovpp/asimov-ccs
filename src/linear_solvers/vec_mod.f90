!v Module file vec.mod
!
!  An interface to operations on vector objects (creation, destruction, setting and getting, ...)

module vec

  use kinds, only: ccs_real, ccs_int
  use types, only: ccs_mesh, ccs_vector, vector_spec, vector_values
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: create_vector
  public :: set_vector_values
  public :: clear_vector_values_entries
  public :: set_vector_values_entry
  public :: create_vector_values
  public :: set_vector_values_mode
  public :: set_vector_values_row
  public :: update_vector
  public :: begin_update_vector
  public :: end_update_vector
  public :: vec_axpy
  public :: vec_aypx
  public :: vec_norm
  public :: initialise_vector
  public :: set_vector_size
  public :: get_vector_data, get_vector_data_readonly
  public :: restore_vector_data, restore_vector_data_readonly
  public :: set_vector_location
  public :: vec_reciprocal
  public :: zero_vector
  public :: mult_vec_vec
  public :: scale_vec
  public :: get_natural_data_vec
  public :: get_global_data_vec
  public :: reorder_data_vec

  interface

    !> Interface to create a new vector object.
    module subroutine create_vector(vec_properties, v, name)
      type(vector_spec), intent(in) :: vec_properties  !< Data structure containing the global and local sizes
      !< of the vector, -1 is interpreted as unset. If both
      !< are set the local size is used.
      class(ccs_vector), allocatable, intent(out) :: v !< The vector returned allocated,
      !< but (potentially) uninitialised.
      character(len=*), optional, intent(in) :: name !< Name of the vector object
    end subroutine

    !> Interface to set values in a vector.
    module subroutine set_vector_values(val_dat, v)
      class(*), intent(in) :: val_dat       !< contains the values, their indices and the mode to use for setting them.
      class(ccs_vector), intent(inout) :: v !< the vector.
    end subroutine set_vector_values

    pure module subroutine clear_vector_values_entries(val_dat)
      type(vector_values), intent(inout) :: val_dat
    end subroutine clear_vector_values_entries

    pure module subroutine set_vector_values_entry(val, val_dat)
      real(ccs_real), intent(in) :: val
      type(vector_values), intent(inout) :: val_dat
    end subroutine set_vector_values_entry

    !> Interface to create a vector values object.
    pure module subroutine create_vector_values(nrows, val_dat)
      integer(ccs_int), intent(in) :: nrows       !< how many rows will be set?
      type(vector_values), intent(out) :: val_dat !< the vector values object
    end subroutine create_vector_values

    pure module subroutine set_vector_values_mode(mode, val_dat)
      integer(ccs_int), intent(in) :: mode
      type(vector_values), intent(inout) :: val_dat
    end subroutine set_vector_values_mode

    !v Interface to set the row currently being worked on by vector values.
    !
    !  Sets the current row in the vector value object, the implementation of this is
    !  backend-dependent as it should immediately convert to the correct indexing
    !  (whether that's 0, 1 or X-based) as used by the backend.
    pure module subroutine set_vector_values_row(row, val_dat)
      integer(ccs_int), intent(in) :: row           !< the row
      type(vector_values), intent(inout) :: val_dat !< the vector values object
    end subroutine set_vector_values_row

    !> Interface to perform a parallel update of a vector.
    module subroutine update_vector(v)
      class(ccs_vector), intent(inout) :: v !< the vector
    end subroutine

    !v Interface to begin a parallel update of a vector.
    !
    !  Begins the parallel update to allow overlapping comms and compute.
    module subroutine begin_update_vector(v)
      class(ccs_vector), intent(inout) :: v !< the vector
    end subroutine

    !v Interface to end a parallel update of a vector.
    !
    !  Ends the parallel update to allow overlapping comms and compute.
    module subroutine end_update_vector(v)
      class(ccs_vector), intent(inout) :: v !< the vector
    end subroutine end_update_vector

    !v Interface to perform the AXPY vector operation.
    !
    !  Performs the AXPY operation
    !
    !        y[i] = a * x[i] + y[i]
    module subroutine vec_axpy(alpha, x, y)
      real(ccs_real), intent(in) :: alpha   !< a scalar value
      class(ccs_vector), intent(in) :: x    !< an input vector
      class(ccs_vector), intent(inout) :: y !< vector serving as input, overwritten with result
    end subroutine

    !v Interface to perform the AYPX vector operation.
    !
    !  Performs the AYPX operation
    !
    !        y[i] = x[i] + beta * y[i]
    module subroutine vec_aypx(x, beta, y)
      class(ccs_vector), intent(in) :: x    !< an input vector
      real(ccs_real), intent(in) :: beta    !< a scalar value
      class(ccs_vector), intent(inout) :: y !< vector serving as input, overwritten with result
    end subroutine

    !> Interface to compute the norm of a vector
    module function vec_norm(v, norm_type) result(n)
      class(ccs_vector), intent(in) :: v !< the vector
      integer(ccs_int), intent(in) :: norm_type !< which norm to compute?
      !< Currently supported is the 2 norm: norm_type=2.
      real(ccs_real) :: n !< the computed norm returned as the result of the function call.
    end function

    !> Constructor for default vector values
    pure module subroutine initialise_vector(vec_properties)
      type(vector_spec), intent(inout) :: vec_properties !< the initialised vector values
    end subroutine initialise_vector

    !> Setter for vector size
    module subroutine set_vector_size(par_env, mesh, vec_properties)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< the parallel environment
      !< where the vector resides
      class(ccs_mesh), target, intent(in) :: mesh !< the mesh - contains the information to set the vector size
      type(vector_spec), intent(inout) :: vec_properties !< the vector data object
    end subroutine set_vector_size

    !> Gets the data in a given vector with possibility to overwrite the data.
    module subroutine get_vector_data(vec, array)
      class(ccs_vector), intent(inout) :: vec                        !< the vector to get data from
      real(ccs_real), dimension(:), pointer, intent(out) :: array !< an array to store the data in
    end subroutine get_vector_data

    !> Resets the vector data if required for further processing
    module subroutine restore_vector_data(vec, array)
      class(ccs_vector), intent(inout) :: vec                       !< the vector to reset
      real(ccs_real), dimension(:), pointer, intent(in) :: array !< the array containing the data to restore
    end subroutine restore_vector_data

    !> Gets the data in a given vector with readonly access.
    module subroutine get_vector_data_readonly(vec, array)
      class(ccs_vector), intent(inout) :: vec                        !< the vector to get data from
      real(ccs_real), dimension(:), pointer, intent(out) :: array !< an array to store the data in
    end subroutine get_vector_data_readonly

    !> Restores the vector data with readonly access.
    module subroutine restore_vector_data_readonly(vec, array)
      class(ccs_vector), intent(inout) :: vec                       !< the vector to reset
      real(ccs_real), dimension(:), pointer, intent(in) :: array !< the array containing the data to restore
    end subroutine restore_vector_data_readonly
    
    !> Set vector values to be located at either cell-centre or face
    pure module subroutine set_vector_location(loc, vec_properties)
      integer(ccs_int), intent(in) :: loc
      type(vector_spec), intent(inout) :: vec_properties
    end subroutine set_vector_location

    !> Replace each component of a vector by its reciprocal
    module subroutine vec_reciprocal(vec)
      class(ccs_vector), intent(inout) :: vec !< the vector, both input and output reciprocal
    end subroutine vec_reciprocal

    module subroutine zero_vector(vec)
      class(ccs_vector), intent(inout) :: vec
    end subroutine zero_vector

    module subroutine mult_vec_vec(a, b)
      class(ccs_vector), intent(in) :: a
      class(ccs_vector), intent(inout) :: b
    end subroutine mult_vec_vec

    module subroutine scale_vec(alpha, v)
      real(ccs_real), intent(in) :: alpha
      class(ccs_vector), intent(inout) :: v
    end subroutine scale_vec

    !> Interface to return the vector data in natural ordering
    module subroutine get_natural_data_vec(par_env, mesh, v, data)
      class(parallel_environment), intent(in) :: par_env
      type(ccs_mesh), intent(in) :: mesh
      class(ccs_vector), intent(inout) :: v
      real(ccs_real), dimension(:), allocatable, intent(out) :: data !< The returned vector data in
      !< natural ordering. Note the use
      !< of allocatable + intent(out),
      !< this ensures it will be
      !< de/reallocated by this subroutine.
    end subroutine get_natural_data_vec

    !> Interface to return the vector data in global ordering
    module subroutine get_global_data_vec(par_env, mesh, v, data)
      class(parallel_environment), intent(in) :: par_env
      type(ccs_mesh), intent(in) :: mesh
      class(ccs_vector), intent(inout) :: v
      real(ccs_real), dimension(:), allocatable, intent(out) :: data !< The returned vector data in
      !< global ordering. Note the use
      !< of allocatable + intent(out),
      !< this ensures it will be
      !< de/reallocated by this subroutine.
    end subroutine get_global_data_vec

    module subroutine reorder_data_vec(par_env, data_from, idx_to, data_to)
      class(parallel_environment), intent(in) :: par_env
      real(ccs_real), dimension(:), intent(in) :: data_from
      integer(ccs_int), dimension(:), intent(in) :: idx_to
      real(ccs_real), dimension(:), intent(out) :: data_to
    end subroutine reorder_data_vec
    
  end interface

end module vec
