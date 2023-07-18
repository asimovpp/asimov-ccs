program split_communicator
#include "ccs_macros.inc"
  
  use mpi
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: initialise_parallel_environment, create_new_par_env, cleanup_parallel_environment
  use utils, only: debug_print, str, exit_print
  
  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: par_env_shared
  logical :: split_flag
  call initialise_parallel_environment(par_env)

  split_flag = .true.
  call create_new_par_env(par_env, MPI_COMM_TYPE_SHARED, split_flag, par_env_shared)

  select type (par_env)
  type is (parallel_environment_mpi)
    select type (par_env_shared)
    type is (parallel_environment_mpi)
      print *, "global rank ", par_env%proc_id, " shared rank ", par_env_shared%proc_id, " size ", par_env_shared%num_procs
    class default
      call error_abort("Unsupported parallel environment")
    end select
  class default
    call error_abort("Unsupported parallel environment")
  end select

  !call cleanup_parallel_environment(par_env_L3)
  call cleanup_parallel_environment(par_env)

end program split_communicator