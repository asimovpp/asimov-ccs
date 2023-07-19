!v Submodule file parallel_env_mpi.smod
!
!  Implementation of the parallel environment using MPI
!
!  @build mpi

submodule(parallel) parallel_env_mpi
#include "ccs_macros.inc"

  use utils, only: exit_print
  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

  !> Create the MPI parallel environment
  module subroutine initialise_parallel_environment(par_env)

    class(parallel_environment), allocatable, intent(out) :: par_env !< parallel_environment_mpi

    integer :: ierr ! Error code

    allocate (parallel_environment_mpi :: par_env)

    select type (par_env)

    type is (parallel_environment_mpi)
      call mpi_init(ierr)
      call error_handling(ierr, "mpi", par_env)

      par_env%comm = MPI_COMM_WORLD
      call set_mpi_parameters(par_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  !> Cleanup the MPI parallel environment
  module subroutine cleanup_parallel_environment(par_env)

    class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi

    integer :: ierr ! Error code

    select type (par_env)

    type is (parallel_environment_mpi)
      call mpi_finalize(ierr)
      call error_handling(ierr, "mpi", par_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

end submodule parallel_env_mpi
