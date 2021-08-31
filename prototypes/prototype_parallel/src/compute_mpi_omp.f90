!> @brief Submodule file compute_mpi_omp.smod
!>
!> @details Implementation of the compute subroutines using MPI and OpenMP
submodule (compute) compute_mpi_omp

  use parallel_types, only: parallel_environment_mpi
  use parallel, only: allreduce

  implicit none

  contains

  !> @brief Compute Pi using MPI and OpenMP
  !>
  !> @param[in] integer(kind=int64) num_steps - number of steps for the algorithm
  !> @param[in] parallel_environment par_env - the parallel environment
  !> @param[out] double precision mypi - computed value for Pi
    module subroutine compute_pi(num_steps, par_env, mypi)

    integer(kind=int64), intent(in) :: num_steps
    class(parallel_environment), intent(in) :: par_env
    double precision, intent(out) :: mypi

    double precision :: step, x, s, finalsum
    integer(kind=int64) :: i, mymax, mymin

    select type (par_env)

    type is (parallel_environment_mpi)
    
      s = 0d0
      step = 1.0d0 / num_steps
    
      mymin = ((par_env%proc_id * num_steps)/par_env%num_procs) + 1
      mymax = ((par_env%proc_id + 1) * num_steps)/par_env%num_procs

!$omp parallel 
!$omp do private(x) reduction(+:s)
      do i = mymin,mymax
        x = (i - 0.5d0) * step
        s = s + 4.0d0 / (1.0d0 + x*x)   
      end do
!$omp enddo
!$omp end parallel

      call allreduce(s, finalsum, par_env%rop%sum_op, par_env)
    
      mypi = finalsum * step  

    class default
      write(*,*) "Unsupported parallel environment"
    
    end select      

  end subroutine compute_pi

end submodule compute_mpi_omp
