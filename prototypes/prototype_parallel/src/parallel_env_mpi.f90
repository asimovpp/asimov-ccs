!> @brief Submodule file parallel_env_mpi.smod
!>
!> @details Implementation of the parallel environment using MPI

submodule (parallel) parallel_env_mpi

  use mpi_f08

  implicit none

  contains

  !> @brief Create the MPI parallel environment
  !>
  !> @param[out] parallel_environment_mpi par_env
  module subroutine initialise_parallel_environment(par_env)

    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    class(parallel_environment), allocatable, intent(out) :: par_env
    allocate(parallel_environment_mpi :: par_env)

    select type (par_env)

    type is (parallel_environment_mpi)   
      par_env%comm = MPI_COMM_WORLD
      call mpi_init(ierr)
      call error_handling(ierr, par_env)

      call mpi_comm_rank(par_env%comm, par_env%proc_id, ierr)
      call error_handling(ierr, par_env)

      call mpi_comm_size(par_env%comm, par_env%num_procs, ierr)
      call error_handling(ierr, par_env)

      ! call set_reduction_operators(par_env%set(MPI_OP))
      call par_env%set_rop()
      
      par_env%root=0
    
    class default
      write(*,*) "Unsupported parallel environment"
    
    end select

  end subroutine

  !> @brief Cleanup the MPI parallel environment
  !>
  !> @param[in] parallel_environment_mpi par_env
  module subroutine cleanup_parallel_environment(par_env)

    class(parallel_environment), intent(in) :: par_env
    integer :: ierr, length, tmp_ierr
    character(len = MPI_MAX_ERROR_STRING) :: error_message

    select type (par_env)

    type is (parallel_environment_mpi)   
      call mpi_finalize(ierr)
      call error_handling(ierr, par_env)
    
    class default
      write(*,*) "Unsupported parallel environment"
    
    end select

    end subroutine

  end submodule parallel_env_mpi