!> @brief Submodule file parallel_env_caf.smod
!>
!> @details Implementation of the parallel environment using CAF

submodule (parallel) parallel_env_caf

  implicit none

  contains

  !> @brief Create the CAF parallel environment
  !>
  !> @param[out] parallel_environment_caf par_env
  module subroutine initialise_parallel_environment(par_env)

    class(parallel_environment), intent(out) :: par_env

    select type (par_env)

    type is (parallel_environment_caf)   
      par_env%proc_id = this_image()
      par_env%num_procs = num_images()
      par_env%root=1
      
    class default
      write(*,*) "Unsupported parallel environment"
    
    end select

  end subroutine

  !> @brief Cleanup the CAF parallel environment
  !>
  !> @param[in] parallel_environment_caf par_env
  module subroutine cleanup_parallel_environment(par_env)

    class(parallel_environment), intent(in) :: par_env

    select type (par_env)

    type is (parallel_environment_caf)  
      ! Nothing to do here 
    
    class default
      write(*,*) "Unsupported parallel environment"
    
    end select

  end subroutine

end submodule parallel_env_caf