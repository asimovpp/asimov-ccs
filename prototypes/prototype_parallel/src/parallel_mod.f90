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
    module subroutine cleanup_parallel_environment(comm)
      integer, intent(in) :: comm
    end subroutine

    !> @brief Synchronise the parallel environment
    module subroutine sync(comm)
      integer, intent(in) :: comm
    end subroutine

    !> @brief Timer for parallel environment
    module subroutine timer(tick)
      double precision, intent(inout) :: tick
    end subroutine

    !> @brief Global sum of integer scalars
    module subroutine global_sum_integer(x, sum, comm)
      integer, intent(in) :: x
      integer, intent(out) :: sum
      integer, intent(in) :: comm
    end subroutine

    !> @brief Global sum of double precision real scalars
    module subroutine global_sum_double(x, sum, comm)
      double precision, intent(in) :: x
      double precision, intent(out) :: sum
      integer, intent(in) :: comm
    end subroutine

    !> @brief Error handling for parallel environment
    module subroutine error_handling(error_code, comm)
      integer, intent(in) :: error_code
      integer, intent(in) :: comm
    end subroutine

    end interface

    interface global_sum
       module procedure global_sum_integer
       module procedure global_sum_double
    end interface global_sum
    
    public :: setup_parallel_environment
    public :: cleanup_parallel_environment
    public :: sync
    public :: timer
    public :: global_sum
    
end module parallel
