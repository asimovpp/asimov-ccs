!> @brief Module file parallel.mod
!>
!> @details Module that defines the parallel interace for ASiMoV-CCS

module parallel

  implicit none

  private

  interface
    !> @brief Create the parallel environment
    module subroutine setup_parallel_environment(comm, rank, size)
      integer, intent(out) :: comm
      integer, intent(out) :: rank
      integer, intent(out) :: size
      
    end subroutine

    !> @brief Cleanup the parallel environment
    module subroutine cleanup_parallel_environment()
    end subroutine

    end interface

  public :: setup_parallel_environment
  public :: cleanup_parallel_environment
    
end module parallel