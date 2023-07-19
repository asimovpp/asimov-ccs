program test_communicator_validity
#include "ccs_macros.inc"
  
  use mpi
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: initialise_parallel_environment, create_new_par_env, cleanup_parallel_environment, is_valid
  use utils, only: debug_print, str, exit_print
  
  implicit none

  class(parallel_environment), allocatable, target :: par_env
  class(parallel_environment), allocatable, target :: par_env_shared
  class(parallel_environment), allocatable, target :: par_env_uninitialised
  logical :: split_flag
  integer :: colour
  call initialise_parallel_environment(par_env)

  select type (par_env)
  type is (parallel_environment_mpi)
    if (modulo(par_env%proc_id, 2) == 1) then
      colour = 1
    else
      colour = -1
    end if
    split_flag = .false.
    call create_new_par_env(par_env, colour, split_flag, par_env_shared)

    select type (par_env_shared)
    type is (parallel_environment_mpi)
      print *, "colour ", colour, " check if par_env is valid ", is_valid(par_env), " shared ", is_valid(par_env_shared)
    class default
      call error_abort("Unsupported parallel environment")
    end select
  class default
    call error_abort("Unsupported parallel environment")
  end select
      
  call cleanup_parallel_environment(par_env)

end program test_communicator_validity