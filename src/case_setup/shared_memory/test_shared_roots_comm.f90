program split_communicator
#include "ccs_macros.inc"
  
  use mpi
  use constants, only: ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: initialise_parallel_environment, create_new_par_env, create_new_par_env_wrapper, &
                      cleanup_parallel_environment, create_shared_roots_comm, &
                      is_valid, is_root
  use utils, only: debug_print, str, exit_print
  
  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: shared_env
  class(parallel_environment), allocatable, target :: roots_env
  logical :: split_flag
  integer :: colour

  call initialise_parallel_environment(par_env)

  !select type (par_env)
  !type is (parallel_environment_mpi)
  !  split_flag = .false.
  !  call create_new_par_env(par_env, ccs_split_type_low_high, split_flag, shared_env)
  !  call create_shared_roots_comm(par_env, shared_env, roots_env)
  !  
  !  select type (roots_env)
  !  type is (parallel_environment_mpi)
  !    print *, "parallel env mpi"
  !    print *, " roots rank ", roots_env%proc_id, " size ", roots_env%num_procs
  !  type is (parallel_environment)
  !    print *, "regular parallel env"
  !  class default
  !    print *, "something else"
  !  end select

  !  select type (shared_env)
  !  type is (parallel_environment_mpi)
  !    select type (roots_env)
  !    type is (parallel_environment_mpi)
  !      if (is_valid(roots_env)) then
  !        print *, "global rank ", par_env%proc_id, " shared rank ", shared_env%proc_id, " roots rank ", roots_env%proc_id, " size ", roots_env%num_procs
  !      end if

  !    class default
  !      call error_abort("Unsupported parallel environment")
  !    end select

  !  class default
  !    call error_abort("Unsupported parallel environment")
  !  end select

  !class default
  !  call error_abort("Unsupported parallel environment")
  !end select
  
  split_flag = .false.
  !call create_new_par_env(par_env, ccs_split_type_low_high, split_flag, shared_env)
  call create_new_par_env_wrapper(par_env, ccs_split_type_low_high, split_flag, shared_env)
  call create_shared_roots_comm(par_env, shared_env, roots_env)
  !if (is_root(shared_env)) then
  !  colour = 1
  !else 
  !  colour = ccs_split_undefined
  !end if

  !!split_flag = .false.

  !call create_new_par_env(par_env, colour, split_flag, roots_env)
  
  select type (roots_env)
  type is (parallel_environment_mpi)
        call dprint("parallel env mpi")
  type is (parallel_environment)
        call dprint("parallel env")
  class default
        call dprint("something else")
  end select
  
  select type (shared_env)
  type is (parallel_environment_mpi)
        call dprint("parallel env mpi")
  type is (parallel_environment)
        call dprint("parallel env")
  class default
        call dprint("something else")
  end select
  
  select type (shared_env)
      type is (parallel_environment_mpi)
      select type (roots_env)
      type is (parallel_environment_mpi)
          if (is_valid(roots_env)) then
            !print *, "global rank ", par_env%proc_id, " shared rank ", shared_env%proc_id, " roots rank ", roots_env%proc_id, " size ", roots_env%num_procs
            call dprint("global rank " // str(par_env%proc_id) // " shared rank " // str(shared_env%proc_id) // " roots rank " // str(roots_env%proc_id)  // " size " // str(roots_env%num_procs))
          end if
        class default
          call error_abort("Unsupported parallel environment")
        end select
  class default
      call error_abort("Unsupported parallel environment")
  end select

  call cleanup_parallel_environment(par_env)

end program split_communicator