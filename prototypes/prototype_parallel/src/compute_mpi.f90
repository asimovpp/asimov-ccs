submodule (compute) compute_mpi

  use parallel_types, only: parallel_environment_mpi
  use parallel, only: allreduce

  implicit none

  contains

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
    
      ! Remember Fortran loops from 1
      mymin = ((par_env%proc_id * num_steps)/par_env%num_procs) + 1
      mymax = ((par_env%proc_id + 1) * num_steps)/par_env%num_procs

      do i = mymin,mymax
        x = (i - 0.5d0) * step
        s = s + 4.0d0 / (1.0d0 + x*x)   
      end do

      call allreduce(s, finalsum, par_env%rop%sum_op, par_env)
    
      mypi = finalsum * step  

    class default
      write(*,*) "Unsupported parallel environment"
    
    end select      

  end subroutine compute_pi

end submodule compute_mpi