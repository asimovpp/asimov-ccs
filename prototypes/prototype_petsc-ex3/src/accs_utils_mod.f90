!> @brief Module file accs_utils.mod
!>
!> @details Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!>          and call type-specific implementations of the interface in other modules.

module accs_utils

  use iso_c_binding
  
  use accs_types, only : vector
  
  implicit none

  private

  public :: set_values, begin_update, end_update, update
  public :: accs_init

contains

  subroutine set_values(val_dat, obj)
    !> @brief Sets values in an object
    !>
    !> @details Given an object and a struct of values to place in that object, call the appropriate
    !> setter.
    
    use accsvec, only : set_vector_values

    class(*), intent(in) :: val_dat
    class(*), intent(inout) :: obj
    select type (obj)
    class is (vector)
       call set_vector_values(val_dat, obj)
    end select
    
  end subroutine set_values

  subroutine update(obj)
    !> @brief Combine begin and end update for a parallel object
    !>
    !> @details Just tidies up code that would otherwise require explicit begin/end steps.

    class(*), intent(inout) :: obj

    call begin_update(obj)
    call end_update(obj)

  end subroutine update
  
  subroutine begin_update(obj)
    !> @brief Begin updating values for a parallel object
    !>
    !> @details Separates the start of updating a parallel object from finalisation, allows
    !>          overlapping comms and computation.

    use accsvec, only : begin_update_vector
    
    class(*), intent(inout) :: obj

    select type (obj)
    class is (vector)
       call begin_update_vector(obj)
    end select
    
  end subroutine begin_update

  subroutine end_update(obj)
    !> @brief End updating values for a parallel object
    !>
    !> @details Separates the start of updating a parallel object from finalisation, allows
    !>          overlapping comms and computation.

    use accsvec, only : end_update_vector
    
    class(*), intent(inout) :: obj

    select type (obj)
    class is (vector)
       call end_update_vector(obj)
    end select
    
  end subroutine end_update

  subroutine accs_init() bind (C, name="accs_init_")

    use petsc, only : PetscInitialize, PETSC_NULL_CHARACTER
    use accs_kinds, only : accs_err
    
    integer(accs_err) :: ierr
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    if (ierr /= 0) then
       print *, "Unable to initialise PETSc"
       stop
    end if
    
  end subroutine accs_init
  
end module accs_utils
