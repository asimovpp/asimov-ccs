!v Submodule for generic parallelism configuration

submodule(core) core_parallel

  use parallel, only: create_new_par_env
  
  implicit none

contains
  
  module subroutine configure_parallelism(run_options, par_env, shared_env)

    type(ccs_options), intent(in) :: run_options
    class(parallel_environment), allocatable, intent(in) ::  par_env    !< The main parallel environment
    class(parallel_environment), allocatable, intent(out) :: shared_env !< The shared parallel environment

    call create_new_par_env(par_env, run_options%split_type, run_options%use_mpi_splitting, shared_env)

  end subroutine configure_parallelism
  
end submodule core_parallel
