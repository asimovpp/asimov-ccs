!> @brief Module file utils.mod
!
!> @details Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!!          and call type-specific implementations of the interface in other modules.

module utils

  use iso_c_binding

  use vec, only : set_vector_values, update_vector, begin_update_vector, end_update_vector, &
       initialise_vector, set_vector_size, pack_one_vector_element, &
       mult_vec_vec, scale_vec, zero_vector
  use mat, only : set_matrix_values, update_matrix, begin_update_matrix, end_update_matrix, &
                  initialise_matrix, finalise_matrix, set_matrix_size, &
                  pack_one_matrix_coefficient, zero_matrix
  use solver, only: initialise_linear_system
  
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
  public :: ccs_init
  public :: ccs_finalise
  public :: scale
  public :: zero

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
    module procedure initialise_linear_system
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

  !> @brief ASiMoV-CCS initialisation routine.
  !
  !> @details Currently implemented by calling PetscInitialize. Symbol must be bound using
  !!          ISO_C_BINDING to tie into pFUnit initialiser.
  subroutine ccs_init() bind (C, name="ccs_init_")

    use petsc, only : PetscInitialize, PETSC_NULL_CHARACTER
    use kinds, only : ccs_err
    
    integer(ccs_err) :: ierr
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    
    if (ierr /= 0) then
       print *, "Unable to initialise PETSc"
       stop
    end if
    
  end subroutine ccs_init

  !> @brief ASiMoV-CCS finalisation routine.
  !
  !> @details Currently implemented by calling PetscFinalize. Symbol must be bound using
  !!          ISO_C_BINDING to tie into pFUnit finaliser.
  subroutine ccs_finalise() bind (C, name="ccs_finalise_")

    use petsc, only : PetscFinalize
    use kinds, only : ccs_err

    integer(ccs_err) :: ierr

    call PetscFinalize(ierr)

  end subroutine ccs_finalise
  
end module utils
