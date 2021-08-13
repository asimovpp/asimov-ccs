!> @brief Module file accs_utils.mod
!>
!> @details Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!>          and call type-specific implementations of the interface in other modules.

module accs_utils

  use accs_types, only : vector
  
  implicit none

  private

  public :: accs_free, set_values

contains
  
  subroutine accs_free(obj)
    !> @brief Frees/destroys an object
    !>
    !> @details Given some object, call the appropriate destructor.

    use accsvec, only : free_vector
    
    class(*), intent(inout) :: obj

    select type (obj)
    class is (vector)
       call free_vector(obj)
    end select
    
  end subroutine accs_free

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
  
end module accs_utils
