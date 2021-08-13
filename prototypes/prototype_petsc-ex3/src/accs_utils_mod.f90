!> @brief Module file accs_utils.mod
!>
!> @details Provides utility functions for ASiMoV-CCS, these should be polymorphic on their input
!>          and call type-specific implementations of the interface in other modules.

module accs_utils

  use accs_types, only : vector
  use accsvec, only : free_vector
  
  implicit none

  private

  public :: accs_free

contains
  
  subroutine accs_free(obj)
    !> @brief Frees/destroys an object
    !>
    !> @details Given some object, call the appropriate destructor.
    
    class(*), intent(inout) :: obj

    select type (obj)
    class is (vector)
       call free_vector(obj)
    end select
    
  end subroutine accs_free
  
end module accs_utils
