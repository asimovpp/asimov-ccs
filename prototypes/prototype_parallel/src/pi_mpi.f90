program pi

  use iso_fortran_env

  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, &
                      sync, timer
  use parallel_types, only: parallel_environment_mpi
  use compute, only: compute_pi

  implicit none

  ! create default parallel environment
  type(parallel_environment_mpi) :: par_env

  double precision :: mypi
  double precision :: start_time, end_time
  integer(kind=int64) :: num_steps

  num_steps = 1000000000

  call initialise_parallel_environment(par_env)

! Output start message
  if (par_env%proc_id == par_env%root) then
    write(*,'(A)') "Calculating PI using:"
    write(*,'(A,1I16,A)') "                  ",num_steps, " slices"
    write(*,'(A,1I16,A)') "                  ",par_env%num_procs," process(es)"
  end if

  call sync(par_env)
  call timer(start_time)

  call compute_pi(num_steps, par_env, mypi)

  call sync(par_env)
  call timer(end_time)

! Output value of PI and time taken
  if (par_env%proc_id == par_env%root) then
    write(*,'(A,1F12.10,A)') "Obtained value of PI: ", mypi
    write(*,'(A,1F12.5,A)') "Time taken:           ", (end_time-start_time), " seconds"
  end if

  call cleanup_parallel_environment(par_env)

 end program pi