!> @brief Module file utils.mod
!
!> @details Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!!          and call type-specific implementations of the interface in other modules.

module utils

  use iso_c_binding

  use vec, only : set_vector_values, update_vector, begin_update_vector, end_update_vector, &
                  initialise_vector, set_vector_size,         &
                  set_vector_values_mode, set_vector_values_row, set_vector_values_entry, &
                  clear_vector_values_entries, &
                  mult_vec_vec, scale_vec, zero_vector
  use mat, only : set_matrix_values, update_matrix, begin_update_matrix, end_update_matrix, &
                  initialise_matrix, finalise_matrix, set_matrix_size, &
                  set_matrix_values_mode, set_matrix_values_row, set_matrix_values_entry, &
                  clear_matrix_values_entries, pack_one_matrix_coefficient, zero_matrix
  use solver, only: initialise_equation_system
  
  implicit none

  private

  public :: set_values
  public :: set_entry
  public :: clear_entries
  public :: begin_update
  public :: end_update
  public :: update
  public :: finalise
  public :: pack_entries
  public :: initialise
  public :: set_size
  public :: mult
  public :: scale
  public :: zero
  public :: set_mode
  public :: set_row

  !> @brief Generic interface to set values on an object.
  interface set_values
     module procedure set_vector_values
     module procedure set_matrix_values
  end interface set_values

  interface set_entry
    module procedure set_vector_values_entry
    module procedure set_matrix_values_entry
  end interface set_entry
 
  interface set_mode
    module procedure set_vector_values_mode
    module procedure set_matrix_values_mode
  end interface set_mode

  interface set_row
    module procedure set_vector_values_row
    module procedure set_matrix_values_row
  end interface set_row
  
  interface finalise
    module procedure finalise_matrix
  end interface finalise

  interface clear_entries
    module procedure clear_vector_values_entries
    module procedure clear_matrix_values_entries
  end interface clear_entries

  !> @brief Generic interface to perform parallel update of an object.
  interface update
     module procedure update_vector
     module procedure update_matrix
  end interface update

  !> @brief Generic interface to begin parallel update of an object.
  !
  !> @details This is to allow overlapping comms and compute.
  interface begin_update
     module procedure begin_update_vector
     module procedure begin_update_matrix
  end interface begin_update

  !> @brief Generic interface to end parallel update of an object.
  !
  !> @details This is to allow overlapping comms and compute.
  interface end_update
    module procedure end_update_vector
    module procedure end_update_matrix
  end interface end_update

  !> @brief Generic interface to initialse vectors, matrices and linear systems
  interface initialise
    module procedure initialise_vector
    module procedure initialise_matrix
    module procedure initialise_equation_system
  end interface initialise

  !> @brief Generic interface to set vector and matrix sizes
  interface set_size
    module procedure set_vector_size
    module procedure set_matrix_size
  end interface set_size

  !> @brief Generic interface to pack entries (elements, coefficients) into a computational object.
  !
  !> @details Stores the entries and elements in an object for later setting, this ensures the
  !!          storage and values of indices in particular are set appropriately for each backend.
  interface pack_entries
    module procedure pack_one_matrix_coefficient
  end interface pack_entries

  !> @brief Generic interface to perform multiplications
  interface mult
    module procedure mult_vec_vec
  end interface mult

  !> @brief Generic interface to scale an object by a scalar constant
  interface scale
    module procedure scale_vec
  end interface scale

  !> @brief Generic interface to zero an object
  interface zero
    module procedure zero_vector
    module procedure zero_matrix
  end interface zero
  
contains

end module utils
