submodule(core) core_fields
#include "ccs_macros.inc"

  use types, only: vector_spec, field_spec
  use constants, only: face, cell, face_centred, cell_centred_central
  use ccs_base, only: mesh
  use parallel, only: is_root
  use utils, only: set_size, initialise, set_is_fluid_field_solved, add_fluid_field_to_outputlist
  use read_config, only: get_store_residuals, get_enable_cell_corrections, get_boundary_count
  use vec, only: set_vector_location
  use fields, only: set_field_config_file, set_field_n_boundaries, set_field_store_residuals, &
                    set_field_enable_cell_corrections, set_field_vector_properties, create_field, &
                    set_field_name, set_field_type

implicit none

contains

  subroutine set_field_properties(par_env, run_options, field_properties)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for building/reading the mesh
    type(field_spec), intent(out) :: field_properties                 !< The resulting field_spec object

    integer(ccs_int) :: n_boundaries
    logical:: store_residuals, enable_cell_corrections
    type(vector_spec) :: vec_properties

  	! Read boundary conditions
    if (is_root(par_env)) then
      print *, "Read and allocate BCs"
    end if
    ! XXX: these calls should probably be moved to the config reading section
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
    
    type(field_spec) :: field_properties
    
    ! Set field properties
    call set_field_properties(par_env, run_options, field_properties)

    call build_common_fields(par_env, run_options, field_properties, flow_fields)
  end subroutine initialise_fields

  subroutine build_common_fields(par_env, run_options, field_properties, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for building/reading the mesh
    type(field_spec), intent(inout) :: field_properties               !< The field spec object used to allocate the fields
    type(fluid), intent(out) :: flow_fields                           !< The fluid fields object being initialised

    integer(ccs_int) :: i
    type(vector_spec) :: vec_properties
    
    ! Expect to find u, v, w, p, p_prime, scalar
    if (is_root(par_env)) then
      print *, "Build field list"
    end if

    do i = 1, size(run_options%variable_names)
      if (is_root(par_env)) then
        print *, "Creating field ", trim(run_options%variable_names(i))
      end if
      call set_field_type(run_options%variable_types(i), field_properties)
      call set_field_name(run_options%variable_names(i), field_properties)
      call create_field(par_env, field_properties, flow_fields)
      call add_fluid_field_to_outputlist(run_options, i, flow_fields)
      call set_is_fluid_field_solved(run_options%solve(i), i, flow_fields)
    end do

    if (is_root(par_env)) then
      print *, "Built ", size(flow_fields%fields), " dynamically-defined fields"
    end if
  
    call set_vector_location(face, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call set_field_vector_properties(vec_properties, field_properties)
    call set_field_type(face_centred, field_properties)
    call set_field_name("mf", field_properties)
    call create_field(par_env, field_properties, flow_fields)
  end subroutine build_common_fields

end submodule core_fields