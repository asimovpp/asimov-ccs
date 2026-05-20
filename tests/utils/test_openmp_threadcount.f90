!!!v Test openmp is running under multiple threads.

program test_openmp_threadcount
  use testing_lib
  use mpi
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  integer :: nthreads

  call init()

#ifdef _OPENMP
  !$omp parallel
  !$omp single
  nthreads = omp_get_num_threads()
  !$omp end single
  !$omp end parallel
#else
  nthreads = 1
#endif

  if (nthreads <= 1) then
    call stop_test("Number of omp threads not greater than 1.")
  end if

  call fin()

end program test_openmp_threadcount
