program pi

  use ISO_FORTRAN_ENV
  use mpi

  use parallel, only: setup_parallel_environment, &
                      cleanup_parallel_environment, &
                      sync, &
                      timer

  implicit none

  integer :: comm
  integer :: rank
  integer :: numprocs
  integer :: ierr

  double precision :: step, x, s, finalsum, mypi
  double precision :: start_time, end_time
  integer(kind=int64)           :: num_steps, i, mymax, mymin
  integer                       :: argl

  num_steps = 1000000000

  call setup_parallel_environment(comm, rank, numprocs)

! Output start message

  if (rank == 0) then
    write(*,'(A)') "Calculating PI using:"
    write(*,'(A,1I16,A)') "                  ",num_steps, " slices"
    write(*,'(A,1I16,A)') "                  ",numprocs," process(es)"
  end if

! Initialise time counter and sum: set step size

  call sync(comm)

  call timer(start_time)
  s = 0d0
  step = 1.0d0 / num_steps

! Remember Fortran loops from 1
  mymin = ((rank * num_steps)/numprocs) + 1
  mymax = ((rank + 1) * num_steps)/numprocs

  do i = mymin, mymax
    x = (i - 0.5d0) * step
    s = s + 4.0d0 / (1.0d0 + x*x)
  end do

  call MPI_Reduce(s, finalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

! Evaluate PI from the final sum value, and stop the clock

  mypi = finalsum * step

  call sync(comm)
  call timer(end_time)

! output value of PI and time taken
! note cpu_time is only specified as being microsecond res

  if (rank == 0) then
    write(*,'(A,1F12.10,A)') "Obtained value of PI: ", mypi
    write(*,'(A,1F12.5,A)') "Time taken:           ", (end_time-start_time), " seconds"
  end if


  call cleanup_parallel_environment(comm)

end program pi

