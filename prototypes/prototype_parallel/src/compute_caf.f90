!> @brief Submodule file compute_caf.smod
!>
!> @details Implementation of the compute subroutines using CAF
submodule (compute) compute_caf

  use parallel_types, only: parallel_environment_caf

  implicit none

  contains

  !> @brief Compute Pi using CAF
  !>
  !> @param[in] integer(kind=int64) num_steps - number of steps for the algorithm
  !> @param[in] parallel_environment par_env - the parallel environment
  !> @param[out] double precision mypi - computed value for Pi
  module subroutine compute_pi(num_steps, par_env, mypi)

    integer(kind=int64), intent(in) :: num_steps
    class(parallel_environment), intent(in) :: par_env
    double precision, intent(out) :: mypi

    double precision :: step, x, sum
    integer(kind=int64) :: i, myid, nimg

    double precision, allocatable :: partial[:] ! co-array that holds partial sum
    allocate(partial[*])

    select type (par_env)

    type is (parallel_environment_caf)   

      write(*,*) "Using CAF compute implementation"

      myid = par_env%proc_id ! id of current image
      nimg = par_env%num_procs ! number of images

      sum = 0d0
      step = 1.0/num_steps
  
      do i=myid,num_steps,nimg
        x = (i + 0.5d0)*step
        sum = sum + 4.0 / (1.0d0 + x*x)
      end do
  
      partial = sum*step

      call co_sum(partial, 1) ! sum over all values of partial on all images
  
      mypi = partial
      
    class default
      write(*,*) "Unsupported parallel environment"
    
    end select      

  end subroutine compute_pi

end submodule compute_caf
