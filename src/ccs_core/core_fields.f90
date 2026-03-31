submodule(core) core_fields
#include "ccs_macros.inc"

  use types, only: vector_spec, field_spec, field
  use constants, only: face, cell, face_centred, cell_centred_central
  use bc_constants
  use boundary_conditions, only: set_bc_type
  use ccs_base, only: mesh
  use parallel, only: is_root
  use utils, only: set_size, initialise, add_field_to_outputlist, exit_print
  use read_config, only: get_store_residuals, get_enable_cell_corrections, get_boundary_count
  use vec, only: set_vector_location
  use fields, only: set_field_config_file, set_field_n_boundaries, set_field_store_residuals, &
                    set_field_enable_cell_corrections, set_field_vector_properties, create_field, &
                    set_field_name, set_field_type, set_is_field_solved, get_field, print_field_config
  use fortran_yaml_c_interface, only: parse
  use profiler, only: profiler_begin_region, profiler_end_region

implicit none

contains

  !> Sets field spec based on specified run options in preparation for building fields
  subroutine set_field_properties(par_env, run_options, field_properties)
    class(parallel_environment), intent(in), allocatable:: par_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options                   !< Object containing relevant options for setting field properties
    type(field_spec), intent(out) :: field_properties              !< The resulting field_spec object

    integer(ccs_int) :: n_boundaries
    logical:: store_residuals, enable_cell_corrections
    type(vector_spec) :: vec_properties

    class(*), pointer :: config_file
    character(:), allocatable :: error

    ! Read boundary conditions
    if (is_root(par_env)) then
      print *, "Read and allocate BCs"
    end if
    ! XXX: these calls should probably be moved to the config reading section
    config_file => parse(run_options%paths%ccs_config_file, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if
    call get_boundary_count(config_file, n_boundaries)
    call get_store_residuals(run_options%paths%ccs_config_file, store_residuals)
    call get_enable_cell_corrections(run_options%paths%ccs_config_file, enable_cell_corrections)

    ! Create and initialise field vectors
    if (is_root(par_env)) then
      print *, "Initialise field vectors"
    end if
    call initialise(vec_properties)

    call set_vector_location(cell, vec_properties)
    call set_size(par_env, mesh, vec_properties)

    call set_field_config_file(run_options%paths%ccs_config_file, field_properties)
    call set_field_n_boundaries(n_boundaries, field_properties)
    call set_field_store_residuals(store_residuals, field_properties)
    call set_field_enable_cell_corrections(enable_cell_corrections, field_properties)

    call set_field_vector_properties(vec_properties, field_properties)
  end subroutine set_field_properties

  !> Builds required fields, including user specified fields, required common fields, and those required by the specific case.
  module subroutine initialise_fields(par_env, run_options, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env    !< The parallel environment
    type(ccs_options), intent(in) :: run_options                      !< Object containing relevant options for initialising fields
    type(fluid), intent(out) :: flow_fields                           !< The fluid fields object being initialised
    
    type(field_spec) :: field_properties
    
    call profiler_begin_region('Field initialisation')
    ! Set field properties
    call set_field_properties(par_env, run_options, field_properties)

    ! build the fields specified in the case config file.
    call build_user_fields(par_env, run_options, field_properties, flow_fields)
    
    ! build the common fields specified.
    call build_common_fields(par_env, field_properties, flow_fields)
    
    ! Finally build any case specific fields.
    call build_case_fields(par_env, run_options, field_properties, flow_fields)

    call profiler_end_region('Field initialisation')

    call print_field_config(par_env, flow_fields)

  end subroutine initialise_fields

  !> Builds the user specified fields from the case config file
  subroutine build_user_fields(par_env, run_options, field_properties, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options                   !< Object containing relevant options for building fields
    type(field_spec), intent(in) :: field_properties               !< The field spec object used to allocate the fields
    type(fluid), intent(inout) :: flow_fields                      !< The fluid fields object being initialised

    integer(ccs_int) :: i
    type(field_spec) :: my_field_properties
    
    ! Create a local copyt of field_properties
    my_field_properties = field_properties
    
    do i = 1, size(run_options%variables%variable_names)
      ! Make sure we don't attempt to define mf. 
      if (trim(run_options%variables%variable_names(i)) == 'mf') then 
        if (is_root(par_env)) then
          print *, "mf already defined in code. skipping definition in case config file"
        end if
        cycle
      end if

      call set_field_type(run_options%variables%variable_types(i), my_field_properties)
      call set_field_name(run_options%variables%variable_names(i), my_field_properties)
      call create_field(par_env, my_field_properties, flow_fields)
      call add_fluid_field_to_outputlist(run_options, i, flow_fields)
      call set_field_solver_params(run_options, i, flow_fields)
    end do

    if (is_root(par_env)) then
      print *, "Built ", size(flow_fields%fields), " dynamically-defined fields"
    end if
  end subroutine build_user_fields

  !> builds any common fields that should be inaccessible to the user.
  subroutine build_common_fields(par_env, field_properties, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env !< The parallel environment
    type(field_spec), intent(in) :: field_properties               !< The field spec object used to allocate the fields
    type(fluid), intent(inout) :: flow_fields                      !< The fluid fields object being initialised

    type(vector_spec) :: vec_properties
    type(field_spec) :: my_field_properties

    integer :: nfields_init

    ! Create local copy of field properties
    my_field_properties = field_properties

    nfields_init = size(flow_fields%fields)
    
    call set_vector_location(face, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call set_field_vector_properties(vec_properties, my_field_properties)
    call set_field_type(face_centred, my_field_properties)
    call set_field_name("mf", my_field_properties)
    call create_field(par_env, my_field_properties, flow_fields)

    if (is_root(par_env)) then
      print *, "Built ", size(flow_fields%fields) - nfields_init, " common fields"
    end if

  end subroutine build_common_fields
  
  !> builds any case specific fields not specified in the case config file.
  subroutine build_case_fields(par_env, run_options, field_properties, flow_fields)
    class(parallel_environment), intent(in), allocatable:: par_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options                   !< Object containing relevant options for building fields
    type(field_spec), intent(in) :: field_properties               !< The field spec object used to allocate the fields
    type(fluid), intent(inout) :: flow_fields                      !< The fluid fields object being initialised

    type(vector_spec) :: vec_properties
    character(len=ccs_string_len), dimension(:), allocatable :: field_names
    integer(ccs_int), dimension(:), allocatable :: field_types
    integer(ccs_int) :: i
    integer(ccs_int) :: field_index

    type(field_spec) :: my_field_properties

    ! Create local copy of field properties
    my_field_properties = field_properties

    allocate(field_names(2))
    field_names(1) = "viscosity"
    field_names(2) = "density"
    field_types = [cell_centred_central, cell_centred_central]

    call set_vector_location(cell, vec_properties)
    call set_size(par_env, mesh, vec_properties)
    call set_field_vector_properties(vec_properties, my_field_properties)
    field_index = size(flow_fields%fields) + 1
    do i = 1, size(field_names) 
      if (.not. is_field_built(field_names(i), flow_fields)) then
        call set_field_type(field_types(i), my_field_properties)
        call set_field_name(field_names(i), my_field_properties)
        call create_field(par_env, my_field_properties, flow_fields)

        call add_fluid_field_to_outputlist(run_options, field_index, flow_fields)
        field_index = field_index + 1
      end if
    end do
  end subroutine build_case_fields

  !> Checks whether the specified field has been built
  function is_field_built(field_name, flow_fields) result(is_built)
    character(len=*), intent(in) :: field_name    !< The name of the field being checked
    type(fluid), intent(in) :: flow_fields        !< The fluid fields object
    logical :: is_built

    integer(ccs_int) :: i
    class(field), pointer :: phi

    do i = 1, size(flow_fields%fields)
      call get_field(flow_fields, i, phi)
      if (trim(field_name) == phi%name) then
        is_built = .true.
        return
      end if
      nullify(phi)
    end do
    is_built = .false.
  end function is_field_built
  
  !> Adds the field specified by field index to the outputlist
  subroutine add_fluid_field_to_outputlist(run_options, field_index, flow)
    type(ccs_options), intent(in) :: run_options  !< Object containing relevant options for building/reading the mesh
    integer(ccs_int), intent(in) :: field_index   !< The index of the field being set
    type(fluid), intent(inout) :: flow            !< The fluid fields object being initialised

    class(field), pointer :: phi
    
    call get_field(flow, field_index, phi)
    if (any(phi%name == run_options%variables%output_variables)) then
      call add_field_to_outputlist(phi)
    end if
    nullify(phi)
  end subroutine add_fluid_field_to_outputlist
  

  ! Sets field solver parameters from run_options and the default values
  subroutine set_field_solver_params(run_options, field_index, flow)
    type(ccs_options), intent(in) :: run_options
    integer(ccs_int), intent(in) :: field_index   !< The index of the field being set
    type(fluid), intent(inout) :: flow            !< The fluid fields object being initialised

    class(field), pointer :: phi

    call get_field(flow, field_index, phi)
    phi%solver_parameters = run_options%solve%solver_eq_parameters(field_index)

    ! Set values to default if not specified per field
    if (phi%solver_parameters%res_target == huge(ccs_real)) then
      phi%solver_parameters%res_target = run_options%solve%default_res_target
    end if

    if (phi%solver_parameters%solver_name == "") then
      if (phi%name == "p" .OR. phi%name == "p_prime") then
          phi%solver_parameters%solver_name = run_options%solve%default_pressure_solver
      else
          phi%solver_parameters%solver_name = run_options%solve%default_solver
      end if
    end if

    if (phi%solver_parameters%precon_name == "") then
      if (phi%name == "p" .OR. phi%name == "p_prime") then
          phi%solver_parameters%precon_name = run_options%solve%default_pressure_precon
      else
          phi%solver_parameters%precon_name = run_options%solve%default_precon
      end if
    end if

    if (phi%solver_parameters%solve .and. phi%solver_parameters%relaxation_factor == huge(ccs_real)) then
      call error_abort("No values assigned to underrelaxation factor for variable "//phi%name)
    end if

    nullify(phi)

  end subroutine

end submodule core_fields
