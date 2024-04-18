submodule(core) initialise_mesh

use ccs_base, only: mesh
use mesh_utils, only: build_mesh, read_mesh, write_mesh

implicit none

contains

  module subroutine initialise_mesh(par_env, shared_env, run_options)
    type(parallel_environment), intent(in) :: par_env
    type(parallel_environment), intent(in) :: shared_env
    type(ccs_options), intent(in) :: run_options

    call timer_register_start("Mesh read time", timer_index_build)
    if (cps /= huge(0)) then
      ! Create a cubic mesh
      if (irank == par_env%root) print *, "Building mesh"
      mesh = build_mesh(par_env, shared_env, cps, cps, cps, domain_size)
    else
      if (irank == par_env%root) print *, "Reading mesh file"
      call read_mesh(par_env, shared_env, case_name, mesh)
    end if
    call read_mesh(par_env, shared_env, case_name, mesh)
    call set_mesh_object(mesh)
    call timer_stop(timer_index_build)
  
    ! Write out mesh to file
    call timer_register_start("I/O time for mesh", timer_index_io_init)
    call write_mesh(par_env, case_path, mesh)
    call timer_stop(timer_index_io_init)
  end subroutine initialise_mesh

end submodule initialise_mesh