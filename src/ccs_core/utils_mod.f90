!> @brief Module file utils.mod
!>
!> @details Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!>          and call type-specific implementations of the interface in other modules.

module utils

  use iso_c_binding

  use vec, only : set_vector_values, update_vector, begin_update_vector, end_update_vector
  use mat, only : set_matrix_values, update_matrix, begin_update_matrix, end_update_matrix
  use types, only : vector, matrix
  
  implicit none

  private

  public :: set_values, begin_update, end_update , update
  public :: accs_init, accs_finalise

  !> @brief Generic interface to set values on an object.
  interface set_values
     module procedure set_vector_values
     module procedure set_matrix_values
  end interface set_values

  !> @brief Generic interface to perform parallel update of an object.
  interface update
     module procedure update_vector
     module procedure update_matrix
  end interface update

  !> @brief Generic interface to begin parallel update of an object.
  !>
  !> @details This is to allow overlapping comms and compute.
  interface begin_update
     module procedure begin_update_vector
     module procedure begin_update_matrix
  end interface begin_update

  !> @brief Generic interface to end parallel update of an object.
  !>
  !> @details This is to allow overlapping comms and compute.
  interface end_update
     module procedure end_update_vector
     module procedure end_update_matrix
  end interface end_update
  
contains

  !> @brief ASiMoV-CCS initialisation routine.
  !>
  !> @details Currently implemented by calling PetscInitialize. Symbol must be bound using
  !>          ISO_C_BINDING to tie into pFUnit initialiser.
  subroutine accs_init() bind (C, name="accs_init_")

    use petsc, only : PetscInitialize, PETSC_NULL_CHARACTER
    use kinds, only : accs_err
    
    integer(accs_err) :: ierr
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr /= 0) then
       print *, "Unable to initialise PETSc"
       stop
    end if
    
  end subroutine accs_init

  !> @brief ASiMoV-CCS finalisation routine.
  !>
  !> @details Currently implemented by calling PetscFinalize. Symbol must be bound using
  !>          ISO_C_BINDING to tie into pFUnit finaliser.
  subroutine accs_finalise() bind (C, name="accs_finalise_")

    use petsc, only : PetscFinalize
    use accs_kinds, only : accs_err

    integer(accs_err) :: ierr

    call PetscFinalize(ierr)

  end subroutine accs_finalise
  
end module utils
