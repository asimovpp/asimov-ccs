submodule(core) core_fields
#include "ccs_macros.inc"

  use types, only: vector_spec, field_spec
  use constants, only: cell
  use ccs_base, only: mesh
  use parallel, only: is_root
  use utils, only: set_size

implicit none

contains

  subroutine set_field_properties(par_env, run_options)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for building/reading the mesh

    integer(ccs_int) :: n_boundaries
    logical:: store_residuals, enable_cell_corrections
    type(vector_spec) :: vec_properties
    type(field_spec)  :: field_properties

  	! Read boundary conditions
    if (is_root(par_env)) then
      print *, "Read and allocate BCs"
    end if
    call get_boundary_count(run_options%config_file, n_boundaries)
    call get_store_residuals(run_options%config_file, store_residuals)
    call get_enable_cell_corrections(run_options%config_file, enable_cell_corrections)

  	! Create and initialise field vectors
    if (is_root(par_env)) then
      print *, "Initialise field vectors"
    end if
    call initialise(vec_properties)

    call set_vector_location(cell, vec_properties)
    call set_size(par_env, mesh, vec_properties)

    call set_field_config_file(run_options%config_file, field_properties)
    call set_field_n_boundaries(n_boundaries, field_properties)
    call set_field_store_residuals(store_residuals, field_properties)
    call set_field_enable_cell_corrections(enable_cell_corrections, field_properties)

    call set_field_vector_properties(vec_properties, field_properties)
  end subroutine set_field_properties

  module subroutine initialise_fields(par_env, run_options, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for building/reading the mesh
    type(fluid), intent(out) :: flow_fields                           !< The fluid fields object being initialised
    
    ! Set field properties
    call set_field_properties(par_env, run_options)
  end subroutine initialise_fields

end submodule core_fields