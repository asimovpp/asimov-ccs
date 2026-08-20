submodule(core) core_init_mesh
#include "ccs_macros.inc"

  use utils, only: exit_print
  use ccs_base, only: mesh
  use parallel, only: is_root
  use mesh_utils, only: build_mesh, build_square_mesh, read_mesh
  use meshing, only: set_mesh_object
  use profiler, only: profiler_begin_region, profiler_end_region
  use logging, only: log_unit_out

  implicit none

contains

  module subroutine initialise_mesh(par_env, shared_env, run_options)
    class(parallel_environment), intent(in), allocatable :: par_env    !< The global parallel environment
    class(parallel_environment), intent(in), allocatable :: shared_env !< The shared parallel environment
    type(ccs_options), intent(in) :: run_options                       !< Object containing relevant options for building/reading the mesh

    call profiler_begin_region('Mesh initialisation')
    if (is_root(par_env)) then
      write (log_unit_out, *) "******************************************************************************"
      write (log_unit_out, *) "* MESH INFO"
    end if

    associate (cps => run_options%mesh%cps)
      call profiler_begin_region("Mesh read time")
      select case (run_options%mesh%init_mesh_type)
      case (build_mesh_2d)
        ! Create a cubic mesh
        if (is_root(par_env)) then
          write (log_unit_out, *) "* Building 2D mesh"
        end if
        mesh = build_square_mesh(par_env, shared_env, run_options)
      case (build_mesh_3d)
        ! Create a cubic mesh
        if (is_root(par_env)) then
          write (log_unit_out, *) "* Building 3D mesh"
        end if
        mesh = build_mesh(par_env, shared_env, run_options)
      case (read_input_mesh)
        if (is_root(par_env)) then
          write (log_unit_out, *) "* Reading mesh file"
        end if
        call read_mesh(par_env, shared_env, run_options, mesh)
      case default
        call error_abort("invalid init mesh type specified")
      end select
      call set_mesh_object(mesh)
      call profiler_end_region("Mesh read time")

      if (is_root(par_env)) then
        if ((run_options%mesh%init_mesh_type == build_mesh_2d) .or. &
            (run_options%mesh%init_mesh_type == build_mesh_3d)) then
          write (log_unit_out, *) "* Cells per side: ", cps
          write (log_unit_out, '(1x, a, e10.3)') "* Domain size: ", run_options%mesh%domain_size
        end if
        write (log_unit_out, *) "* Global number of cells is ", mesh%topo%global_num_cells
        write (log_unit_out, *) "******************************************************************************"
      end if
    end associate

    call profiler_end_region('Mesh initialisation')

  end subroutine initialise_mesh

end submodule core_init_mesh
