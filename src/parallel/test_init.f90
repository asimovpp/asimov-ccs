program test_init_petsc

  use petsc
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment
  use parallel_types, only: parallel_environment

  implicit none 

  ! create default parallel environment
  class(parallel_environment), allocatable :: par_env

  call initialise_parallel_environment(par_env)

  write(*,*) "MPI_COMM_WORLD = ", MPI_COMM_WORLD
  write(*,*) "PETSC_COMM_WORLD = ", PETSC_COMM_WORLD

  call cleanup_parallel_environment(par_env)

end program test_init_petsc