program split_communicator
#include "ccs_macros.inc"
  
  use mpi
  use constants, only: ccs_split_type_shared, ccs_split_type_low_high
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: initialise_parallel_environment, create_new_par_env, cleanup_parallel_environment
  use utils, only: debug_print, str, exit_print
  
  use testing_lib

  implicit none

  class(parallel_environment), allocatable, target :: par_env_shared
  logical :: split_flag
  integer :: colour

  call init()

  select type (par_env)
  type is (parallel_environment_mpi)
    ! Testing user split functionality with MPI_COMM_NULL communicator for even ranks.
    if (modulo(par_env%proc_id, 2) == 1) then
      colour = 1
    else
      colour = -1
    end if

    split_flag = .false.
    call create_new_par_env(par_env, colour, split_flag, par_env_shared)

    select type (par_env_shared)
    type is (parallel_environment_mpi)
      print *, "global rank ", par_env%proc_id, " shared rank ", par_env_shared%proc_id, " size ", par_env_shared%num_procs
    class default
      call stop_test("Unsupported parallel environment")
    end select
  class default
    call stop_test("Unsupported parallel environment")
  end select

  call fin()

end program split_communicator
