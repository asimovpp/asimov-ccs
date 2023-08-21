program test_communicator_validity
#include "ccs_macros.inc"
  
  use mpi
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: initialise_parallel_environment, create_new_par_env, cleanup_parallel_environment, is_valid
  use utils, only: debug_print, str, exit_print
 
  use testing_lib

  implicit none

  class(parallel_environment), allocatable, target :: par_env_shared
  logical :: split_flag
  integer :: colour

  call init()

  ! XXX: Are these select types necessary?
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
      if (colour /= -1) then
        ! Expect to be part of the shared environment
        if (.not. is_valid(par_env_shared)) then
          call stop_test("Shared environment invalid on process that was assigned to it!")
        end if
      else
        ! Expect not to be part of the shared environment
        if (is_valid(par_env_shared)) then
          call stop_test("Shared environment active on process that was excluded from it!")
        end if
      end if
    class default
      call stop_test("Unsupported parallel environment")
    end select
  class default
    call stop_test("Unsupported parallel environment")
  end select
 
  call fin()

end program test_communicator_validity
