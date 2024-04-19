submodule(core) core_init_mesh
#include "ccs_macros.inc"

  use utils, only: exit_print, debug_print, str
  use ccs_base, only: mesh
  use parallel_types_mpi, only: parallel_environment_mpi
  use mesh_utils, only: build_mesh, build_square_mesh, read_mesh, write_mesh
  use meshing, only: set_mesh_object
  use timers, only: timer_register_start, timer_stop
  use kinds, only: ccs_int

implicit none

contains

  module subroutine initialise_mesh(par_env, shared_env, run_options)
    class(parallel_environment), intent(in), allocatable :: par_env
    class(parallel_environment), intent(in), allocatable :: shared_env
    type(ccs_options), intent(in) :: run_options
  
    integer(ccs_int):: timer_index_build
    integer(ccs_int):: timer_index_io_init

    associate (cps => run_options%cps)
    select type(par_env)
    type is (parallel_environment_mpi)
      call timer_register_start("Mesh read time", timer_index_build)
      if (run_options%init_mesh_type == build_mesh_2d) then
        ! Create a cubic mesh
        if (par_env%proc_id == par_env%root) print *, "Building mesh"
        mesh = build_square_mesh(par_env, shared_env, cps, run_options%domain_size)
      else if (run_options%init_mesh_type == build_mesh_3d) then
        ! Create a cubic mesh
        if (par_env%proc_id == par_env%root) print *, "Building mesh"
        mesh = build_mesh(par_env, shared_env, cps, cps, cps, run_options%domain_size)
      else if (run_options%init_mesh_type == read_input_mesh) then
        if (par_env%proc_id == par_env%root) print *, "Reading mesh file"
        call read_mesh(par_env, shared_env, run_options%case_name, mesh)
      else 
        call error_abort("invalid init mesh type specified")
      end if
      call set_mesh_object(mesh)
      call timer_stop(timer_index_build)
  
      ! Write out mesh to file
      call timer_register_start("I/O time for mesh", timer_index_io_init)
      call write_mesh(par_env, run_options%case_path, mesh)
      call timer_stop(timer_index_io_init)
    class default
      call error_abort("Invalid parallel environment")
    end select
    end associate
  end subroutine initialise_mesh

end submodule core_init_mesh