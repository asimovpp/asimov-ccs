!!!v Test program for data reordering.

program test_reorder_data

  use testing_lib
  
  use utils, only: reorder_data
  
  implicit none

  integer, parameter :: nglobal = 10
  integer :: nlocal, istart
  integer :: proc_id
  
  integer, dimension(:), allocatable :: to ! Destination address
  real(ccs_real), dimension(:), allocatable :: src_data ! Data I hold
  real(ccs_real), dimension(:), allocatable :: data ! Data I recieve

  logical :: passing

  integer :: i
  
  call init()
  
  proc_id = par_env%proc_id
  if(proc_id == 0) then
    nlocal = 2
    istart = 1
    to = [ 7, 1 ]
  else if (proc_id == 1) then
    nlocal = 3
    istart = 3
    to = [ 2, 5, 8 ]
  else if (proc_id == 2) then
    nlocal = 3
    istart = 6
    to = [ 3, 10 ]
  else if (proc_id == 3) then
    nlocal = 2
    istart = 9
    to = [ 4, 6, 9 ]
  else
    call stop_test("This test only runs on 4 processors")
  end if
  allocate(data(nlocal))
  data(:) = 0
  src_data = real(to)

  call reorder_data(par_env, src_data, to, data)
  
  passing = .true.
  do i = 1, nlocal
    if (data(i) /= (istart + (i - 1))) then
      passing = .false.
    end if
  end do

  if (.not. passing) then
    print *, proc_id, data
    call stop_test("Test failed")
  end if
  
  call fin()
  
end program test_reorder_data
