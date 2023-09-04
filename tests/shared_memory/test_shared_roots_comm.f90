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

  use_mpi_splitting = .false.
  call create_new_par_env(par_env, ccs_split_type_low_high, use_mpi_splitting, test_shared_env)
  call create_shared_roots_comm(par_env, test_shared_env, test_roots_env)
  
  ! Only root of the shared environment(s) should be part of the roots env
  if (is_root(test_shared_env)) then
    if (.not. is_valid(test_roots_env)) then
      call stop_test("Root of shared environment not included in roots environment")
    end if
  else
    if (is_valid(test_roots_env)) then
      call stop_test("Non-root of shared environment included in roots environment")
    end if
  end if

  call fin()

end program test_shared_roots_comm
