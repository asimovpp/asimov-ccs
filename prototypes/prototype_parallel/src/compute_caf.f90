submodule (compute) compute_caf

  use iso_fortran_env
  use parallel_types, only: parallel_environment_caf

  implicit none

  contains

  module subroutine compute_pi(num_steps, par_env, mypi)

    integer(kind=int64), intent(in) :: num_steps
    class(parallel_environment), intent(in) :: par_env
    double precision, intent(out) :: mypi

    double precision :: step, x, sum, finalsum
    double precision, save :: partial[*]

    integer(kind=int64) :: i, myid, nimg

    select type (par_env)

    type is (parallel_environment_caf)   
      myid = par_env%proc_id ! id of current image
      nimg = par_env%num_procs ! number of image

      step = 1.0/num_steps
  
      do i=myid,num_steps,nimg
        x = (i+0.5)*step
        sum = sum + 4.0/(1.0+x*x)
      end do
  
      partial = sum*step

      call co_sum(partial, 1) ! sum over all values of partial on all images
  
      mypi = partial
      
    class default
      write(*,*) "Unsupported parallel environment"
    
    end select      

  end subroutine compute_pi

end submodule compute_caf