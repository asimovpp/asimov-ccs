module accs_utils

  use accs_types, only : vector
  use accsvec, only : free_vector
  
  implicit none

  private

  public :: accs_free

contains
  
  subroutine accs_free(obj)
    class(*), intent(inout) :: obj

    select type (obj)
    class is (vector)
       call free_vector(obj)
    end select
    
  end subroutine
  
end module accs_utils
