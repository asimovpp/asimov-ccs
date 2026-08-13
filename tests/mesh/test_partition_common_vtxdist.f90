!> Tests computation of the graph vertex distribution
program test_partition_common_vtxdist

  use testing_lib, only: stop_test

  use kinds, only: ccs_int, ccs_long

  use partitioning, only: proccnt_to_vtxdist

  implicit none

  integer, parameter :: nproc = 5
  integer(ccs_int), dimension(nproc), parameter :: proc_ctr = [ &
       2, 1, 2, 0, 1 &
       ]
  integer(ccs_long), dimension(nproc + 1), parameter :: vtxdist_expect = &
       [1, 3, 4, 6, 6, 7]
  integer(ccs_long), dimension(:), allocatable :: vtxdist

  vtxdist = proccnt_to_vtxdist(proc_ctr)

  if (size(vtxdist) /= (nproc + 1)) then
     call stop_test("vtxdist size is wrong!")
  end if
  if (any(vtxdist /= vtxdist_expect)) then
     call stop_test("vtxdist does not match expection")
  end if
  
end program test_partition_common_vtxdist
