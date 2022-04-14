!> @brief Module file utils.mod
!
!> @details Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!!          and call type-specific implementations of the interface in other modules.

module utils

  use iso_c_binding

  use vec, only : set_vector_values, update_vector, begin_update_vector, end_update_vector, &
       initialise_vector, set_vector_size, pack_one_vector_element, &
       mult_vec_vec, zero_vector
  use mat, only : set_matrix_values, update_matrix, begin_update_matrix, end_update_matrix, &
                  initialise_matrix, finalise_matrix, set_matrix_size, &
                  pack_one_matrix_coefficient, zero_matrix
  use solver, only : initialise_equation_system
  use kinds, only : ccs_int, ccs_real

  
  implicit none

  private

  public :: set_values
  public :: begin_update
  public :: end_update
  public :: update
  public :: finalise
  public :: pack_entries
  public :: initialise
  public :: set_size
  public :: mult
  public :: zero
  public :: str
  public :: debug_print

  !> @brief Generic interface to set values on an object.
  interface set_values
     module procedure set_vector_values
     module procedure set_matrix_values
  end interface set_values

  interface finalise
    module procedure finalise_matrix
  end interface finalise

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
    module procedure pack_one_vector_element
    module procedure pack_one_matrix_coefficient
  end interface pack_entries

  !> @brief Generic interface to perform multiplications
  interface mult
    module procedure mult_vec_vec
  end interface mult

  !> @brief Generic interface to zero an object
  interface zero
    module procedure zero_vector
    module procedure zero_matrix
  end interface zero
  
  !> @brief Generic interface to converting numbers to strings
  interface str
    module procedure int2str
    module procedure real2str
  end interface str

  !> @brief Generic interface to debug printer
  interface debug_print
    module procedure dprint_print
    module procedure noop
  end interface debug_print
  
contains

  !> @brief Print a message, along with with its location
  subroutine dprint_print(msg, filename, line)
    character(*), intent(in) :: msg
    character(*), intent(in) :: filename
    integer, intent(in) :: line
  
    print *, trim(filename), "(", trim(adjustl(int2str(line))), ") : ", msg
  end subroutine
  
  !> @brief No-op routine, does nothign
  subroutine noop()
  end subroutine
  
  !> @brief Convert integer to string
  function int2str(in_int, format_str) result(out_string)
    integer(ccs_int), intent (in)       :: in_int
    character(*), optional, intent(in)  :: format_str
    character(32)                       :: tmp_string
    character(:), allocatable           :: out_string
    
    if (present(format_str)) then
      write(tmp_string, format_str) in_int
    else
      write(tmp_string, *) in_int
    end if
    out_string = trim(adjustl(tmp_string))
  end function
  
  !> @brief Convert real to string
  function real2str(in_real, format_str) result(out_string)
    real(ccs_real), intent (in)         :: in_real
    character(*), optional, intent(in)  :: format_str
    character(32)                       :: tmp_string
    character(:), allocatable           :: out_string
    
    if (present(format_str)) then
      write(tmp_string, format_str) in_real
    else
      write(tmp_string, *) in_real
    end if
    out_string = trim(adjustl(tmp_string))
  end function

end module utils
