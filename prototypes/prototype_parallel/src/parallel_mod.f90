!> @brief Module file parallel.mod
!>
!> @details Module that defines the parallel interace for ASiMoV-CCS

module parallel

  implicit none

  private

  interface
    !> @brief Create the parallel environment
    module subroutine setup_parallel_environment(comm, rank, numprocs)
      integer, intent(out) :: comm
      integer, intent(out) :: rank
      integer, intent(out) :: numprocs
      
    end subroutine

    !> @brief Cleanup the parallel environment
    module subroutine cleanup_parallel_environment()
    end subroutine

    !> @brief Synchronise the parallel environment
    module subroutine sync(comm)
      integer, intent(in) :: comm
    end subroutine


  end interface

  public :: setup_parallel_environment
  public :: cleanup_parallel_environment
  public :: sync
    
end module parallel