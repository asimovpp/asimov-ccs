program test_shared_roots_comm
#include "ccs_macros.inc"
  
  use mpi
  use constants, only: ccs_split_type_shared, ccs_split_type_low_high, ccs_split_undefined
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: initialise_parallel_environment, create_new_par_env, &
                      cleanup_parallel_environment, create_shared_roots_comm, &
                      is_valid, is_root
  use utils, only: debug_print, str, exit_print
  
  use testing_lib

  implicit none

  class(parallel_environment), allocatable, target :: test_roots_env, test_shared_env
  logical :: use_mpi_splitting

  call init()

  select type (par_env)
  type is (parallel_environment_mpi)
    use_mpi_splitting = .false.
    call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, test_shared_env)
    call create_shared_roots_comm(par_env, test_shared_env, test_roots_env)
    
    select type (test_shared_env)
    type is (parallel_environment_mpi)
      select type (test_roots_env)
      type is (parallel_environment_mpi)
        if (is_valid(test_roots_env)) then
          call dprint("shared rank " // str(test_shared_env%proc_id) // " roots rank " // str(test_roots_env%proc_id) // " size " // str(roots_env%num_procs))
        end if

      class default
        call stop_test("Unsupported parallel environment")
      end select

    class default
      call stop_test("Unsupported parallel environment")
    end select

  class default
    call stop_test("Unsupported parallel environment")
  end select

  call fin()

end program test_shared_roots_comm
