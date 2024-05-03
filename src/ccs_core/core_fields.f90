submodule(core) core_fields
#include "ccs_macros.inc"

  use types, only: vector_spec, field_spec
  use constants, only: face, cell, face_centred, cell_centred_central
  use bc_constants
  use boundary_conditions, only: set_bc_type, allocate_bc_arrays
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

  !> Sets field spec based on specified run options in preparation for building fields
  subroutine set_field_properties(par_env, run_options, field_properties)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for setting field properties
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
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for initialising fields
    type(fluid), intent(out) :: flow_fields                           !< The fluid fields object being initialised
    
    type(field_spec) :: field_properties
    
    ! Set field properties
    call set_field_properties(par_env, run_options, field_properties)

    ! build the fields specified in the case config file.
    call build_user_fields(par_env, run_options, field_properties, flow_fields)
    
    ! build the common fields specified.
    call build_common_fields(par_env, field_properties, flow_fields)
    
    ! Finally build any case specific fields.
    call build_case_fields(par_env, field_properties, flow_fields)
  end subroutine initialise_fields

  !> Builds the user specified fields from the case config file
  subroutine build_user_fields(par_env, run_options, field_properties, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for building fields
    type(field_spec), intent(inout) :: field_properties               !< The field spec object used to allocate the fields
    type(fluid), intent(out) :: flow_fields                           !< The fluid fields object being initialised

    integer(ccs_int) :: i
    type(vector_spec) :: vec_properties
    
    ! Expect to find u, v, w, p, p_prime, scalar
    if (is_root(par_env)) then
      print *, "Build field list"
    end if

    do i = 1, size(run_options%variable_names)
      ! Make sure we don't attempt to define mf. 
      if (trim(run_options%variable_names(i)) == 'mf') then 
        cycle
      end if

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
  end subroutine build_user_fields

  !> builds any common fields that should be inaccessible to the user.
  subroutine build_common_fields(par_env, field_properties, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(field_spec), intent(inout) :: field_properties               !< The field spec object used to allocate the fields
    type(fluid), intent(inout) :: flow_fields                           !< The fluid fields object being initialised

    type(vector_spec) :: vec_properties

    call set_vector_location(face, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call set_field_vector_properties(vec_properties, field_properties)
    call set_field_type(face_centred, field_properties)
    call set_field_name("mf", field_properties)
    call create_field(par_env, field_properties, flow_fields)
  end subroutine build_common_fields
  
  !> builds any case specific fields not specified in the case config file.
  subroutine build_case_fields(par_env, field_properties, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(field_spec), intent(inout) :: field_properties               !< The field spec object used to allocate the fields
    type(fluid), intent(inout) :: flow_fields                           !< The fluid fields object being initialised

    type(vector_spec) :: vec_properties
    character(len=:), dimension(:), allocatable :: field_names
    integer(ccs_int), dimension(:), allocatable :: field_types
    character(len=:), dimension(:), allocatable :: field_bc_types
    integer(ccs_int) :: i
    integer(ccs_int) :: bc_index
    integer(ccs_int) :: n_boundaries
    integer(ccs_int) :: new_field_index

    field_names = ["viscosity", "density"]
    field_types = [cell_centred_central, cell_centred_central]
    field_bc_types = ["neumann", "neumann"]

    call set_vector_location(cell, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call set_field_vector_properties(vec_properties, field_properties)
    new_field_index = size(flow_fields%fields)
    do i = 1, size(field_names) 
      if (.not. is_field_built(field_names(i), flow_fields)) then
        call set_field_type(field_types(i), field_properties)
        call set_field_name(field_names(i), field_properties)
        call create_field(par_env, field_properties, flow_fields)

        ! For a field that's not specified in the case config file need to set the bcs manually
        associate (bcs => flow_fields%fields(new_field_index)%ptr%bcs, &
                   existing_bcs => flow_fields%fields(1)%ptr%bcs)
          n_boundaries = size(existing_bcs%ids)
          call allocate_bc_arrays(n_boundaries, bcs)
          bcs%ids = existing_bcs%ids  ! Not strictly necessary as we could just set to 1..n_boundaries since boundary types are all going to be the same, but leave it in for safety
          do bc_index = 1, n_boundaries
            call set_bc_type(bc_index, field_bc_types(i), bcs)
          end do
          new_field_index = new_field_index + 1
        end associate
      end if
    end do
  end subroutine build_case_fields

  !> Checks whether the specified field has been built
  function is_field_built(field_name, flow_fields) result(is_built)
    character(len=*), intent(in) :: field_name    !< The name of the field being checked
    type(fluid), intent(in) :: flow_fields        !< The fluid fields object
    logical :: is_built

    integer(ccs_int) :: i

    do i = 1, size(flow_fields%fields)
      if (trim(field_name) == flow_fields%fields(i)%name) then
        is_built = .true.
        return
      end if
    end do
    is_built = .false.
  end function is_field_built

end submodule core_fields