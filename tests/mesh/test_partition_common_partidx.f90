!!! Test the computation of the global indices for our local partition.
!
!! Given a new partitioning
!!   part = [7, 0, 1, 3, 3, 0, ...]
!! based on our current naive (i.e. sliced) partitioning, determine the global index sets for each
!! processor and assemble them as a CSR-style data structure.
program test_partition_common_partidx

  use testing_lib, only: stop_test

  use kinds, only: ccs_int, ccs_long

  use partitioning, only: compute_global_indices_partition

  implicit none

  integer, parameter :: nlocal = 6
  integer, parameter :: nproc = 5
  integer(ccs_long), dimension(nlocal), parameter :: partition = [ &
       2, 4, 0, 0, 1, 2 &
       ]
  integer(ccs_int), dimension(nproc), parameter :: proc_ctr = [ &
       2, 1, 2, 0, 1 &
       ]
  integer(ccs_long), dimension(nproc + 1), parameter :: vtxdist = [ &
       1, 3, 4, 6, 6, 7 &
       ]
  integer(ccs_long), parameter :: global_idx_start = 42
  integer(ccs_long), dimension(nlocal), parameter :: global_indices_exp = [ &
       3, 4, 5, 1, 6, 2 &
       ] + global_idx_start

  integer(ccs_long), dimension(:), allocatable :: global_indices 

  global_indices = compute_global_indices_partition(partition, proc_ctr, vtxdist, global_idx_start)

  if (any(global_indices /= global_indices_exp)) then
     call stop_test("Partitioned global indices do not match expectation")
  end if
  
end program test_partition_common_partidx
