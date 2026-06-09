submodule(core) core_fields
#include "ccs_macros.inc"

  use types, only: vector_spec, field_spec, field
  use constants, only: face, cell, face_centred, cell_centred_central
  use bc_constants
  use boundary_conditions, only: set_bc_type, translate_bcs_phi
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
  use logging, only: log_unit_out

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
      write(log_unit_out,*) "Read and allocate BCs"
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
      write(log_unit_out,*) "Initialise field vectors"
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

    call translate_bcs(par_env, run_options, flow_fields)

    call profiler_end_region('Field initialisation')

    call print_field_config(par_env, flow_fields)

  end subroutine initialise_fields

  !v Loop through fields and translate their potentially complex boundary conditions 
  ! into basic boundary conditions and also update gradient (this cannot be done before bcs are 'translated')
  subroutine translate_bcs(par_env, run_options, flow_fields)
    use fields, only: get_field, get_field_is_face_based
    use fv, only: update_gradient

    class(parallel_environment), intent(in), allocatable:: par_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options !< Object containing runtime options
    type(fluid), intent(inout) :: flow_fields    !< Object containing all the fields to process

    class(field), pointer :: phi
    real(ccs_real), dimension(:, :), allocatable :: bnd_normals
    integer(ccs_int) :: i
    logical :: is_face_based

    call get_bnd_normals(par_env, run_options, bnd_normals)

    do i=1, size(flow_fields%fields)
      call get_field(flow_fields, i, phi)

      call translate_bcs_phi(bnd_normals, phi)

      call get_field_is_face_based(phi, is_face_based)
      if (.not. is_face_based) then
        call update_gradient(phi)
      endif
    end do

  end subroutine

  !v For each boundary get normal vector
  ! Only handles flat boundaries
  subroutine get_bnd_normals(par_env, run_options, bnd_normals)
     use mpi
     use parallel_types_mpi, only: parallel_environment_mpi
     use kinds, only: CCS_MPI_PRECISION
     use types, only: cell_locator, neighbour_locator, face_locator
     use meshing, only: create_cell_locator, create_face_locator, create_neighbour_locator, &
                        get_local_num_cells, get_boundary_status, get_face_normal, count_neighbours, get_local_index

    class(parallel_environment), intent(in), allocatable:: par_env !< The parallel environment
    type(ccs_options), intent(in) :: run_options !< Object containing runtime options
    real(ccs_real), dimension(:, :), allocatable, intent(out) :: bnd_normals !< Array containing normal for each mesh boundary

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    real(ccs_real), dimension(3) :: normal
    integer(ccs_int) :: index_p, index_nb, j, local_num_cells, nnb
    integer(ccs_int) :: n_boundaries
    logical :: is_boundary
    integer :: ierr

    n_boundaries = size(run_options%mesh%bnd_names)

    allocate(bnd_normals(3, n_boundaries))
    bnd_normals(:, :) = -10.0_ccs_real

    call get_local_num_cells(local_num_cells)

    do index_p=1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)

      do j=1, nnb
        call create_face_locator(index_p, j, loc_f)
        call get_boundary_status(loc_f, is_boundary)

        if (is_boundary) then
          call create_neighbour_locator(loc_p, j, loc_nb)
          call get_local_index(loc_nb, index_nb)
          call get_face_normal(loc_f, normal)

          bnd_normals(:, -index_nb) = normal(:)
        end if

      end do
    end do

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_ALLREDUCE(MPI_IN_PLACE, bnd_normals(:, :), 3*n_boundaries, CCS_MPI_PRECISION, MPI_MAX, par_env%comm, ierr)
    class default
      call error_abort("invalid parallel environment")
    end select

  end subroutine
  

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
          write(log_unit_out,*) "mf already defined in code. skipping definition in case config file"
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
      write(log_unit_out,*) "Built ", size(flow_fields%fields), " dynamically-defined fields"
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
      write(log_unit_out,*) "Built ", size(flow_fields%fields) - nfields_init, " common fields"
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
  

  ! Sets field solver parameters from run_options
  subroutine set_field_solver_params(run_options, field_index, flow)
    type(ccs_options), intent(in) :: run_options  !< Runtime options to extract solver parameters information from
    integer(ccs_int), intent(in) :: field_index   !< The index of the field being set
    type(fluid), intent(inout) :: flow            !< The fluid fields object being initialised

    class(field), pointer :: phi

    ! Get field using the run_options equation name
    call get_field(flow, run_options%solve%solver_eq_parameters(field_index)%name, phi)
    phi%solver_parameters = run_options%solve%solver_eq_parameters(field_index)
    nullify(phi)

  end subroutine

end submodule core_fields
