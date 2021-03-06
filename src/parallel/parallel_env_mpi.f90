!>  Submodule file parallel_env_mpi.smod
!
!>  @build mpi
!
!>  Implementation of the parallel environment using MPI

submodule (parallel) parallel_env_mpi

  use mpi
  use parallel_types_mpi, only: parallel_environment_mpi

  implicit none

  contains

  !>  Create the MPI parallel environment
  module subroutine initialise_parallel_environment(par_env)

    class(parallel_environment), allocatable, intent(out) :: par_env

    integer :: ierr !< Error code

    allocate(parallel_environment_mpi :: par_env)

    select type (par_env)

      type is (parallel_environment_mpi)   
        call mpi_init(ierr)
        call error_handling(ierr, "mpi", par_env)

        par_env%comm = MPI_COMM_WORLD

        call mpi_comm_rank(par_env%comm, par_env%proc_id, ierr)
        call error_handling(ierr, "mpi", par_env)

        call mpi_comm_size(par_env%comm, par_env%num_procs, ierr)
        call error_handling(ierr, "mpi", par_env)

        call par_env%set_rop()
      
        par_env%root=0
    
      class default
        write(*,*) "Unsupported parallel environment"
        stop 1
    
    end select

  end subroutine

  !>  Cleanup the MPI parallel environment
  !
  !> @param[in] parallel_environment_mpi par_env
  module subroutine cleanup_parallel_environment(par_env)

    class(parallel_environment), intent(in) :: par_env
    
    integer :: ierr !< Error code

    select type (par_env)

      type is (parallel_environment_mpi)   
        call mpi_finalize(ierr)
        call error_handling(ierr, "mpi", par_env)
    
      class default
        write(*,*) "Unsupported parallel environment"
        stop 1
    
    end select

    end subroutine

  end submodule parallel_env_mpi
