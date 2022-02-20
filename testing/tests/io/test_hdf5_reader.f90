program test_hdf5_reader

  use testing_lib
  use kinds, only : accs_int, accs_real
  use io, only: initialise_io, cleanup_io, configure_io, &
               open_file, close_file, &
               read_scalar, read_array
  
  implicit none

   character(len=:), allocatable :: test_file    !> Testfile name
   character(len=:), allocatable :: adios2_file !> ADIOS2 config file name
 
   class(io_environment), allocatable :: io_env
   class(io_process), allocatable :: test_reader
 
   real, dimension(:), allocatable :: real_var
   integer, dimension(:), allocatable :: int_var
   integer(kind=8), dimension(1) :: sel_start
   integer(kind=8), dimension(1) :: sel_count
   integer, dimension(:), allocatable :: dist
 
   integer(accs_int) :: irank !> MPI rank ID
   integer(accs_int) :: isize !> Size of MPI world
   integer(accs_int) :: i, j, k, ierror 
   integer(accs_int) :: local_start
   integer(accs_int) :: local_end
   integer(kind=8) :: long_int_attr

   integer(accs_int), parameter :: num_cells = 10

   integer(accs_int) :: sum_int
   integer(accs_int) :: loc_sum_int
 
   call init()

   irank = par_env%proc_id
   isize = par_env%num_procs
 
   test_file = "hdf5_test_file.h5"
   adios2_file = "adios2_config.xml"
 
   call initialise_io(par_env, adios2_file, io_env)
   call configure_io(io_env, "test_reader", test_reader)  
 
   call open_file(test_file, "read", test_reader)
 
   ! Test reading long integer attribute
   call read_scalar(test_reader, "adios2_schema/mesh/dimension0", long_int_attr)
   if(long_int_attr /= 10) then
    write(message,*) "FAIL: the adios2_schema/mesh/dimension0 attribute should be 10, not ", long_int_attr
    call stop_test(message)
  end if   

   ! Array to store cell range assigned to each process      
   allocate(dist(isize+1))
 
   ! Set for element to 1 and last element to world size + 1
   dist(1) = 1
   dist(isize + 1) = num_cells + 1
 
   ! Divide the total number of cells by the world size to
   ! compute the chunk sizes
   k = int(real(num_cells) / isize)
   j = 1
   do i = 1, isize
      dist(i) = j
      j = j + k
   enddo
 
   ! First and last cell index assigned to this process
   local_start = dist(irank + 1)
   local_end = dist(irank + 2) - 1

   ! Starting point for reading chunk of data
   sel_start = (/ dist(irank + 1) - 1/)
   ! How many data points will be read?
   sel_count = (/ dist(irank + 2) - dist(irank + 1)/)
 
   ! Allocate memory for real & int arrays on each MPI rank
   allocate(int_var(sel_count(1)))
   allocate(real_var(sel_count(1)))
 
   ! Read XYZ coordinates for variable "/cell/x" 
   call read_array(test_reader, "h5Floats", sel_start , sel_count, real_var)
   call read_array(test_reader, "h5Ints", sel_start , sel_count, int_var)

   loc_sum_int = 0
   do i = 1, int(sel_count(1))
      loc_sum_int = loc_sum_int + int_var(i)
   end do

   select type(par_env)
   type is (parallel_environment_mpi)
     call mpi_allreduce(loc_sum_int, sum_int, 1, MPI_INTEGER, MPI_SUM, par_env%comm, ierror)
   class default
     write(message, *) "ERROR: Unknown parallel environment!"
    call stop_test(message)
  end select

   if(sum_int /= -45) then
      write(message,*) "FAIL: Sum of integer array should be -45, not ", sum_int
      call stop_test(message)
   end if
 
   ! Close the file and ADIOS2 engine
   call close_file(test_reader)
 
   ! Finalise the ADIOS2 IO environment
   call cleanup_io(io_env)
  
   ! Deallocate memory for XYZ coordinates array
   deallocate(int_var)
   deallocate(real_var)

  call fin()

end program test_hdf5_reader
