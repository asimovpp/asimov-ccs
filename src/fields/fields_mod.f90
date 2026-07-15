!!!v Module file field.mod
!!!
!!!  Provides field interface

module fields
#include "ccs_macros.inc"

  use parallel_types, only: parallel_environment
  use constants, only: cell_centred_central, cell_centred_upwind, cell_centred_gamma, cell_centred_linear_upwind, face_centred
  use types, only: field, field_spec, field_ptr, fluid, &
                   vector_spec, face_field, central_field, upwind_field, gamma_field, linear_upwind_field
  use kinds, only: ccs_int

  use utils, only: update, get_scheme_name, debug_print
  use parallel, only: is_root
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use vec, only: create_vector, get_vector_data_readonly
  use timestepping, only: initialise_old_values
  use error_codes
  use logging, only: log_unit_out

  implicit none

  private

  public :: create_field
  public :: get_field_is_face_based
  public :: set_field_config_file
  public :: set_field_vector_properties
  public :: set_field_n_boundaries
  public :: set_field_store_residuals
  public :: set_field_enable_cell_corrections
  public :: set_field_name
  public :: set_field_type
  public :: get_is_field_solved
  public :: set_is_field_solved
  public :: get_field
  public :: get_field_idx
  public :: add_field
  public :: dealloc_fluid_fields
  public :: get_field_name
  public :: count_fields
  public :: print_field_config

  !> Generic interface to get a field from the flow
  interface get_field
    module procedure get_field_byname
    module procedure get_field_byidx
  end interface get_field

contains

  !> Build a field variable with data and gradient vectors + transient data and boundary arrays.
  subroutine create_field(par_env, field_properties, flow)

    use utils, only: debug_print

    implicit none

    class(parallel_environment), intent(in) :: par_env
    type(field_spec), intent(in) :: field_properties !< Field descriptor
    type(fluid), intent(inout) :: flow !< The flow field container where new field is to be constructed

    integer :: nfields
    
    call print_field(par_env, field_properties)
    
    associate (ccs_config_file => field_properties%ccs_config_file, &
               vec_properties => field_properties%vec_properties, &
               field_type => field_properties%field_type, &
               field_name => field_properties%field_name, &
               n_boundaries => field_properties%n_boundaries, &
               store_residuals => field_properties%store_residuals, &
               enable_cell_corrections => field_properties%enable_cell_corrections)
      call allocate_field(vec_properties, field_type, n_boundaries, store_residuals, flow)

      nfields = size(flow%fields)
      associate(phi => flow%fields(nfields)%ptr)
        ! XXX: ccs_config_file is host-associated from program scope.
        call read_bc_config(ccs_config_file, field_name, phi)
      
        phi%enable_cell_corrections = enable_cell_corrections
        phi%name = field_name

        !! --- Ensure data is updated/parallel-constructed ---
        ! XXX: Potential abstraction --- see update(vec), etc.
        call update(phi%values)

        if (store_residuals) then
          call update(phi%residuals)
        end if

        if (field_type /= face_centred) then
          ! Current design only computes/stores gradients at cell centres
          call update(phi%x_gradients)
          call update(phi%y_gradients)
          call update(phi%z_gradients)
        end if

        !! --- End update ---
      end associate
    end associate

  end subroutine create_field

  !> Sets the solve flag for a field
  pure subroutine set_is_field_solved(solve, phi)
    logical, intent(in) :: solve      !< flag indicating whether to solve for the given field
    type(field), intent(inout) :: phi !< Field variable

    phi%solver_parameters%solve = solve
    
  end subroutine set_is_field_solved

  !> Gets the solve flag for a field
  pure subroutine get_is_field_solved(phi, solve)
    class(field), intent(in) :: phi !< Field variable
    logical, intent(out) :: solve   !< flag indicating whether to solve for the given field

    solve = phi%solver_parameters%solve
    
  end subroutine get_is_field_solved

  !> Print a brief description (name, type) of a field as it is created
  subroutine print_field(par_env, field_properties)

    class(parallel_environment), intent(in) :: par_env
    type(field_spec), intent(in) :: field_properties

    character(len=:), allocatable :: field_name
    integer(ccs_int) :: field_type

    field_name = field_properties%field_name
    field_type = field_properties%field_type
    
    if (is_root(par_env)) then
      write(log_unit_out,*) "Create field: ", trim(field_name), " ("//get_scheme_name(field_type)//")"
    end if
    
  end subroutine print_field
  
  !> Allocate a field variable
  subroutine allocate_field(vec_properties, field_type, n_boundaries, store_residuals, flow)

    use utils, only: debug_print

    implicit none

    !! Logically vec_properties should be a field_properties variable, but this doesn't yet exist.
    type(vector_spec), intent(in) :: vec_properties !< Vector descriptor for vectors wrapped by field
    integer, intent(in) :: field_type               !< Identifier for what kind of field
    integer(ccs_int), intent(in) :: n_boundaries    !< Mesh boundary count
    logical, intent(in) :: store_residuals          !< Wether or not residual field needs to be stored (and allocated)
    type(fluid), intent(inout) :: flow              !< The flow field container where new field is to be constructed

    type(field_ptr) :: phi_ptr !< The field being constructed

    if (field_type == face_centred) then
      call dprint("Create face field")
      allocate (face_field :: phi_ptr%ptr)
    else if (field_type == cell_centred_upwind) then
      call dprint("Create upwind field")
      allocate (upwind_field :: phi_ptr%ptr)
    else if (field_type == cell_centred_gamma) then
      call dprint("Create gamma field")
      allocate (gamma_field :: phi_ptr%ptr)
    else if (field_type == cell_centred_linear_upwind) then
      call dprint("Create linear upwind field")
      allocate (linear_upwind_field :: phi_ptr%ptr)
    else if (field_type == cell_centred_central) then
      call dprint("Create central field")
      allocate (central_field :: phi_ptr%ptr)
    else
      error stop "Trying to allocate an unknown field type"
    end if

    call dprint("Create field values vector")
    call create_vector(vec_properties, phi_ptr%ptr%values)
    call get_vector_data_readonly(phi_ptr%ptr%values, phi_ptr%ptr%values_ro)

    if (store_residuals) then
      call dprint("Create residuals field vector")
      call create_vector(vec_properties, phi_ptr%ptr%residuals)
    end if

    if (field_type /= face_centred) then
      ! Current design only computes/stores gradients at cell centres
      call dprint("Create field gradients vector")
      call create_vector(vec_properties, phi_ptr%ptr%x_gradients)
      call create_vector(vec_properties, phi_ptr%ptr%y_gradients)
      call create_vector(vec_properties, phi_ptr%ptr%z_gradients)
      call get_vector_data_readonly(phi_ptr%ptr%x_gradients, phi_ptr%ptr%x_gradients_ro)
      call get_vector_data_readonly(phi_ptr%ptr%y_gradients, phi_ptr%ptr%y_gradients_ro)
      call get_vector_data_readonly(phi_ptr%ptr%z_gradients, phi_ptr%ptr%z_gradients_ro)

      ! Currently no need for old face values
      call dprint("Create field old values")
      call initialise_old_values(vec_properties, phi_ptr%ptr)
    end if

    call allocate_bc_arrays(n_boundaries, phi_ptr%ptr%bcs)
    
    allocate(phi_ptr%ptr%solver_parameters)

    call add_field(phi_ptr, flow)
    
  end subroutine allocate_field

  !> Get whether a field is face based or not
  pure subroutine get_field_is_face_based(my_field, is_face_based)
    class(field), intent(in) :: my_field
    logical, intent(out) :: is_face_based

    select type (my_field)
      type is (face_field)
        is_face_based = .true.
      class default
        is_face_based = .false.
    end select

  end subroutine

  !> Set config file used for field creation
  pure subroutine set_field_config_file(ccs_config_file, field_properties)

    character(len=*), intent(in) :: ccs_config_file
    type(field_spec), intent(inout) :: field_properties

    field_properties%ccs_config_file = ccs_config_file

  end subroutine set_field_config_file

  !> Set the vector properties used for field creation
  subroutine set_field_vector_properties(vec_properties, field_properties)

    type(vector_spec), intent(in) :: vec_properties
    type(field_spec), intent(inout) :: field_properties

    field_properties%vec_properties = vec_properties

  end subroutine set_field_vector_properties

  !> Set the number of boundaries used for field creation
  pure subroutine set_field_n_boundaries(n_boundaries, field_properties)

    integer(ccs_int), intent(in) :: n_boundaries
    type(field_spec), intent(inout) :: field_properties

    field_properties%n_boundaries = n_boundaries

  end subroutine set_field_n_boundaries

  !> Set whether or not residuals should be stored
  pure subroutine set_field_store_residuals(store_residuals, field_properties)

    logical, intent(in) :: store_residuals
    type(field_spec), intent(inout) :: field_properties

    field_properties%store_residuals = store_residuals

  end subroutine set_field_store_residuals

  !> Set whether or not cell shape corrections should be used
  pure subroutine set_field_enable_cell_corrections(enable_cell_corrections, field_properties)

    logical, intent(in) :: enable_cell_corrections 
    type(field_spec), intent(inout) :: field_properties

    field_properties%enable_cell_corrections = enable_cell_corrections

  end subroutine set_field_enable_cell_corrections

  !> Set the name of a field to be created
  pure subroutine set_field_name(name, field_properties)

    character(len=*), intent(in) :: name
    type(field_spec), intent(inout) :: field_properties

    field_properties%field_name = name

  end subroutine set_field_name

  !> Set the type of field to be created
  pure subroutine set_field_type(field_type, field_properties)

    integer(ccs_int), intent(in) :: field_type
    type(field_spec), intent(inout) :: field_properties

    field_properties%field_type = field_type

  end subroutine set_field_type

  !> Gets the field from the fluid structure specified by field_name
  subroutine get_field_byname(flow, field_name, flow_field)
    type(fluid), intent(in) :: flow                   !< the structure containing all the fluid fields
    character(len=*), intent(in) :: field_name
    class(field), pointer, intent(out) :: flow_field  !< the field of interest

    integer(ccs_int) :: i

    logical :: found = .false.

    do i = 1, size(flow%fields)
      call get_field_byidx(flow, i, flow_field)
      if (trim(flow_field%name) == trim(field_name)) then
        found = .true.
        exit
      else
        found = .false.
        nullify (flow_field)
      end if
    end do

    if (.not. found) then
      error stop field_not_found ! Field name not found
    end if

  end subroutine get_field_byname

  ! Deallocates fluid arrays
  subroutine dealloc_fluid_fields(flow)
    type(fluid), intent(inout) :: flow  !< The fluid structure to deallocate

    deallocate (flow%fields)
  end subroutine dealloc_fluid_fields

  !> Get the count of stored fields
  pure subroutine count_fields(flow, nfields)

    type(fluid), intent(in) :: flow          !< The flowfield
    integer(ccs_int), intent(out) :: nfields !< The count of fields

    nfields = size(flow%fields)

  end subroutine count_fields

  !< Sets the pointer to the field and the corresponding field name in the fluid structure
  subroutine add_field(flow_field_ptr, flow)
    type(field_ptr), target, intent(in) :: flow_field_ptr !< the field
    type(fluid), intent(inout) :: flow                    !< the fluid structure

    logical, save :: first_call = .true.

    ! Handle the case when a program body is called in a loop (e.g. as part of a convergence test)
    if (.not. allocated(flow%fields)) then
      first_call = .true.
    end if

    if (first_call) then
      allocate (flow%fields(1))
      flow%fields(1) = flow_field_ptr
      first_call = .false.
    else
      flow%fields = [flow%fields, flow_field_ptr]
    end if

  end subroutine add_field
  !> Gets the field from the fluid structure specified by field_index
  subroutine get_field_byidx(flow, field_index, flow_field)
    type(fluid), intent(in) :: flow                   !< the structure containing all the fluid fields
    integer, intent(in) :: field_index
    class(field), pointer, intent(out) :: flow_field  !< the field of interest

    if (field_index > size(flow%fields)) then
      error stop field_index_exceeded ! Field index exceeds number of flow fields
    end if

    if (field_index <= 0) then
      error stop "Field index less than or equal to zero"
    end if

    flow_field => flow%fields(field_index)%ptr

  end subroutine get_field_byidx

  !> Get the name of the i'th field
  subroutine get_field_name(flow, idx, field_name)

    type(fluid), intent(in) :: flow                          !< The flowfield
    integer(ccs_int), intent(in) :: idx                      !< The field counter
    character(len=:), allocatable, intent(out) :: field_name !< The field name

    class(field), pointer :: phi

    call get_field(flow, idx, phi)
    field_name = phi%name
    nullify (phi)

  end subroutine get_field_name

  !> Get the field idx from of flow_field in flow
  pure subroutine get_field_idx(flow, flow_field, idx)
    type(fluid), intent(in) :: flow                          !< The flowfield
    class(field), intent(in) :: flow_field  !< the field of interest
    integer(ccs_int), intent(out) :: idx                      !< The field counter
    integer(ccs_int) :: nfields, ifield

    idx = -1
    call count_fields(flow, nfields)

    do ifield = 1, nfields
      if (trim(flow_field%name) == trim(flow%fields(ifield)%ptr%name)) then
        idx = ifield
        return
      end if
    end do

  end subroutine get_field_idx


  !v Outputs each field configuration
  subroutine print_field_config(par_env, flow)

    use logging, only: log_unit_out
    use kinds, only: CCS_PRECISION_STR
    use parallel, only: is_root
    use constants, only: L2, Linfty

    class(parallel_environment), intent(in) :: par_env !< Parallel environment
    type(fluid), intent(in) :: flow !< fluid containing the list of fields to solve and display info of

    class(field), pointer :: phi

    integer :: ifield, nfields
    call count_fields(flow, nfields)

    if (is_root(par_env)) then
        write(log_unit_out,*)  " "
        write(log_unit_out,*)  "******************************************************************************"
        write(log_unit_out,*)  "* SOLVER CONFIGURATION"
        write(log_unit_out,*)  "******************************************************************************"
        write(log_unit_out,*)  "* Precision: ", CCS_PRECISION_STR
        do ifield=1, nfields
          call get_field(flow, ifield, phi)
          if (phi%solver_parameters%solve) then
            write(log_unit_out,*) "************ ", trim(phi%name), " ************"
            write(log_unit_out,'(A, E10.2)')  "   Residuals target: ", phi%solver_parameters%res_target
            select case(phi%solver_parameters%res_norm)
              case (L2)
                write(log_unit_out,*) "  Residuals norm: L2"
              case (Linfty)
                write(log_unit_out,*) "  Residuals norm: Linfty"
            end select
            write(log_unit_out,'(A, F4.2)')  "   Relaxation factor: ", phi%solver_parameters%relaxation_factor
            write(log_unit_out,*) "  Solver: ", phi%solver_parameters%solver_name
            write(log_unit_out,*) "  Preconditioner: ", phi%solver_parameters%precon_name
          end if
        end do
        write(log_unit_out,*)  " "
    end if

  end subroutine

end module fields
