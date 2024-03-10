!!!v Module file field.mod
!!!
!!!  Provides field interface

module fields
#include "ccs_macros.inc"

  use constants, only: cell_centred_central, cell_centred_upwind, cell_centred_gamma, cell_centred_linear_upwind, face_centred
  use types, only: field, field_spec, field_ptr, fluid, &
                   vector_spec, face_field, central_field, upwind_field, gamma_field, linear_upwind_field
  use kinds, only: ccs_int

  use utils, only: update, add_field
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use vec, only: create_vector, get_vector_data_readonly
  use fv, only: update_gradient
  use timestepping, only: initialise_old_values

  implicit none

  private

  public :: create_field
  public :: set_field_config_file
  public :: set_field_vector_properties
  public :: set_field_n_boundaries
  public :: set_field_store_residuals
  public :: set_field_enable_cell_corrections
  public :: set_field_name
  public :: set_field_type

contains

  !> Build a field variable with data and gradient vectors + transient data and boundary arrays.
  subroutine create_field(field_properties, flow)

    use utils, only: debug_print

    implicit none

    type(field_spec), intent(in) :: field_properties !< Field descriptor
    type(fluid), intent(inout) :: flow !< The flow field container where new field is to be constructed

    integer :: nfields
    
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

          call update_gradient(phi)
        end if

        phi%name = field_name
        !! --- End update ---
      end associate
    end associate

  end subroutine create_field

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

    call add_field(phi_ptr, flow)
    
  end subroutine allocate_field

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

end module fields
