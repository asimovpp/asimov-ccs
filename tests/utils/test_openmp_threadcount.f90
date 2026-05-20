!!!v Test openmp is running under multiple threads.

program test_openmp_threadcount
  use testing_lib
  use mpi
  use omp_lib

  implicit none

  integer :: nthreads

  call init()

  !$omp parallel
  !$omp single
  nthreads = omp_get_num_threads()
  !$omp end single
  !$omp end parallel

  if (nthreads <= 1) then
    call stop_test("Number of omp threads < 1.")
  end if

  call fin()

end program test_openmp_threadcount
