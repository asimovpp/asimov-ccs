!v Submodule for generic parallelism configuration

submodule(core) core_parallel

  use parallel, only: create_new_par_env

  implicit none

contains

  !v Subroutine to configure sub parallel environments.
  !
  ! Currently this only configures the shared environment, but it could be extended for particle
  ! environments, etc.
  module subroutine configure_parallelism(run_options, par_env, shared_env)

    type(ccs_options), intent(in) :: run_options                        !< The runtime configuration
    class(parallel_environment), allocatable, intent(in) :: par_env    !< The main parallel environment
    class(parallel_environment), allocatable, intent(out) :: shared_env !< The shared parallel environment

    call configure_parallelism_impl(run_options%parallel, par_env, shared_env)

  end subroutine configure_parallelism

  !> Internal implementation to configure parallelism based on the parallelism options set at runtime.
  subroutine configure_parallelism_impl(par_opt, par_env, shared_env)

    type(parallel_options), intent(in) :: par_opt                       !< The runtime parallelism options
    class(parallel_environment), allocatable, intent(in) :: par_env    !< The main parallel environment
    class(parallel_environment), allocatable, intent(out) :: shared_env !< The shared parallel environment

    call create_new_par_env(par_env, par_opt%split_type, par_opt%use_mpi_splitting, shared_env)

  end subroutine configure_parallelism_impl

end submodule core_parallel
