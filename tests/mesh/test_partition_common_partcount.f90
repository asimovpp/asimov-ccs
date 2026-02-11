!!! Tests computing the processor count for a partition

program test_partition_common_partcount

  use testing_lib, only: stop_test

  use kinds, only: ccs_int, ccs_long

  use partitioning, only: partition_count

  implicit none

  integer, parameter :: nlocal = 6
  integer, parameter :: nproc = 5
  integer(ccs_long), dimension(nlocal), parameter :: partition = [ &
       2, 4, 0, 0, 1, 2 &
       ]
  integer(ccs_int), dimension(nproc), parameter :: proc_ctr_exp = [ &
       2, 1, 2, 0, 1 &
       ]
  integer(ccs_int), dimension(:), allocatable :: proc_ctr

  proc_ctr = partition_count(nproc, partition)

  if (size(proc_ctr) /= nproc) then
     call stop_test("Process counts array misallocated")
  end if

  if (any(proc_ctr /= proc_ctr_exp)) then
     call stop_test("Process counts do not match expectation")
  end if

end program test_partition_common_partcount
