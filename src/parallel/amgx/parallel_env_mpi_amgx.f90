!v Submodule file parallel_env_mpi_amgx.smod
!
!  Implementation of the parallel environment using MPI and AMGX
!
!  @build mpi amgx
submodule(parallel) parallel_env_mpi_amgx
#include "ccs_macros.inc"

  use utils, only: exit_print
  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

contains

  !> Create the MPI and AMGX parallel environments
  module subroutine initialise_parallel_environment(par_env)

    class(parallel_environment), allocatable, intent(out) :: par_env !< parallel_environment_mpi

    integer :: ierr ! Error code

    allocate (parallel_environment_mpi :: par_env)

    select type (par_env)

    type is (parallel_environment_mpi)
      call mpi_init(ierr)
      call error_handling(ierr, "mpi", par_env)

      par_env%comm = MPI_COMM_WORLD

      call initialise_amgx(par_env)

      call mpi_comm_rank(par_env%comm, par_env%proc_id, ierr)
      call error_handling(ierr, "mpi", par_env)

      call mpi_comm_size(par_env%comm, par_env%num_procs, ierr)
      call error_handling(ierr, "mpi", par_env)

      call par_env%set_rop()

      par_env%root = 0

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  !> Cleanup the AMGX and MPI parallel environments
  module subroutine cleanup_parallel_environment(par_env)

    class(parallel_environment), intent(in) :: par_env !< parallel_environment_mpi

    integer :: ierr ! Error code

    select type (par_env)

    type is (parallel_environment_mpi)
      call finalise_amgx(par_env)
      call mpi_finalize(ierr)
      call error_handling(ierr, "mpi", par_env)

    class default
      call error_abort("Unsupported parallel environment")

    end select

  end subroutine

  !> Initalise AMGX
  subroutine initialise_amgx(par_env)

    type(parallel_environment_mpi), intent(in) :: par_env !< parallel_environment_mpi

    integer :: ierr ! Error code

    ! TODO: type for config, and test
    call AMGX_initialize(ierr)
    call AMGX_config_create(config, "communicator=MPI")
    call AMGX_resources_create(resources, config, par_env%comm, 1, (/ 0 /) )

    if (ierr /= 0) then
      call error_handling(ierr, "amgx", par_env)
    end if

  end subroutine

  !> Finalise AMGX
  subroutine finalise_amgx(par_env)

    type(parallel_environment_mpi), intent(in) :: par_env !< parallel_environment_mpi

    integer :: ierr ! Error code

    call AMGX_finalize(ierr)

    if (ierr /= 0) then
      call error_handling(ierr, "amgx", par_env)
    end if

  end subroutine

end submodule parallel_env_mpi_amgx
