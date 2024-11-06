submodule(core) core_init_mesh
#include "ccs_macros.inc"

  use utils, only: exit_print, str
  use ccs_base, only: mesh
  use parallel_types_mpi, only: parallel_environment_mpi
  use parallel, only: is_root
  use mesh_utils, only: build_mesh, build_square_mesh, read_mesh, write_mesh
  use meshing, only: set_mesh_object
  use timers, only: timer_register_start, timer_stop
  use kinds, only: ccs_int

implicit none

contains

  module subroutine initialise_mesh(par_env, shared_env, run_options)
    class(parallel_environment), intent(in), allocatable :: par_env    !< The global parallel environment
    class(parallel_environment), intent(in), allocatable :: shared_env !< The shared parallel environment
    type(ccs_options), intent(in) :: run_options                       !< Object containing relevant options for building/reading the mesh
  
    integer(ccs_int) :: timer_index_build
    integer(ccs_int) :: timer_index_io_init

    if (is_root(par_env)) then
      print *, "******************************************************************************"
      print *, "* MESH INFO"
    end if
    
    associate (cps => run_options%mesh%cps)
      call timer_register_start("Mesh read time", timer_index_build)
      select case (run_options%mesh%init_mesh_type)
      case (build_mesh_2d)
        ! Create a cubic mesh
        if (is_root(par_env)) then
          print *, "* Building 2D mesh"
        end if
        mesh = build_square_mesh(par_env, shared_env, run_options)
      case (build_mesh_3d) 
        ! Create a cubic mesh
        if (is_root(par_env)) then
          print *, "* Building 3D mesh"
        end if
        mesh = build_mesh(par_env, shared_env, run_options)
      case (read_input_mesh)
        if (is_root(par_env)) then
          print *, "* Reading mesh file"
        end if
        call read_mesh(par_env, shared_env, run_options, mesh)
      case default
        call error_abort("invalid init mesh type specified")
      end select
      call set_mesh_object(mesh)
      call timer_stop(timer_index_build)

      if (is_root(par_env)) then
        if ((run_options%mesh%init_mesh_type == build_mesh_2d) .or. &
            (run_options%mesh%init_mesh_type == build_mesh_3d)) then
          print *, "* Cells per side: ", cps
          write (*, '(1x, a, e10.3)') "* Domain size: ", run_options%mesh%domain_size
        end if
        print *, "* Global number of cells is ", mesh%topo%global_num_cells
        print *, "******************************************************************************"
      end if
    end associate
  
    ! Write out mesh to file
    call timer_register_start("I/O time for mesh", timer_index_io_init)
    call write_mesh(par_env, run_options%paths%case_path, mesh)
    call timer_stop(timer_index_io_init)

  end subroutine initialise_mesh

end submodule core_init_mesh
