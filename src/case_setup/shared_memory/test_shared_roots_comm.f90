program split_communicator
#include "ccs_macros.inc"
  
  use mpi
  use constants, only: ccs_split_type_shared, ccs_split_type_low_high
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: initialise_parallel_environment, create_new_par_env, &
                      cleanup_parallel_environment, create_shared_roots_comm, &
                      is_valid
  use utils, only: debug_print, str, exit_print
  
  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  class(parallel_environment), allocatable, target :: roots_env
  logical :: split_flag
  integer :: colour

  call initialise_parallel_environment(par_env)

  select type (par_env)
  type is (parallel_environment_mpi)
    split_flag = .false.
    call create_new_par_env(par_env, ccs_split_type_low_high, split_flag, shared_env)
    call create_shared_roots_comm(par_env, shared_env, roots_env)
    
    select type (roots_env)
    type is (parallel_environment_mpi)
      print *, "parallel env mpi"
    class default
      print *, "something else"
    end select

    select type (shared_env)
    type is (parallel_environment_mpi)
      select type (roots_env)
      type is (parallel_environment_mpi)
        if (is_valid(roots_env)) then
          print *, "global rank ", par_env%proc_id, " shared rank ", shared_env%proc_id, " roots rank ", roots_env%proc_id, " size ", roots_env%num_procs
        end if

      class default
        call error_abort("Unsupported parallel environment")
      end select

    class default
      call error_abort("Unsupported parallel environment")
    end select

  class default
    call error_abort("Unsupported parallel environment")
  end select

  call cleanup_parallel_environment(par_env)

end program split_communicator