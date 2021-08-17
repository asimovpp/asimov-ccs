program pi

  use ISO_FORTRAN_ENV

  use parallel, only: setup_parallel_environment, cleanup_parallel_environment, &
                      sync, timer, &
                      allreduce
  use parallel_types, only: parallel_environment_mpi, &
                            reduction_operator_mpi, &
                            set_reduction_operators

  implicit none

  ! create default parallel environment
  type(parallel_environment_mpi) :: par_env
  ! create reduction operator
  type(reduction_operator_mpi) :: rop

  double precision :: step, x, s, finalsum, mypi
  double precision :: start_time, end_time
  integer(kind=int64) :: num_steps, i, mymax, mymin

  num_steps = 1000000000

  call setup_parallel_environment(par_env)
  call set_reduction_operators(rop)

! Output start message

  if (par_env%rank == 0) then
    write(*,'(A)') "Calculating PI using:"
    write(*,'(A,1I16,A)') "                  ",num_steps, " slices"
    write(*,'(A,1I16,A)') "                  ",par_env%numprocs," process(es)"
  end if

! Initialise time counter and sum: set step size

  call sync(par_env)

  call timer(start_time)
  s = 0d0
  step = 1.0d0 / num_steps

! Remember Fortran loops from 1
  mymin = ((par_env%rank * num_steps)/par_env%numprocs) + 1
  mymax = ((par_env%rank + 1) * num_steps)/par_env%numprocs

  do i = mymin, mymax
    x = (i - 0.5d0) * step
    s = s + 4.0d0 / (1.0d0 + x*x)
  end do

  call allreduce(s, finalsum, rop%sum_op, par_env)

! Evaluate PI from the final sum value, and stop the clock

  mypi = finalsum * step

  call sync(par_env)
  call timer(end_time)

! output value of PI and time taken
! note cpu_time is only specified as being microsecond res

  if (par_env%rank == 0) then
    write(*,'(A,1F12.10,A)') "Obtained value of PI: ", mypi
    write(*,'(A,1F12.5,A)') "Time taken:           ", (end_time-start_time), " seconds"
  end if

  call cleanup_parallel_environment(par_env)

end program pi

