!!!v Module file field.mod
!!!
!!!  Provides field interface

module fields
#include "ccs_macros.inc"
  
  use constants, only: cell_centred_central, cell_centred_upwind, face_centred
  use types, only: field, vector_spec, face_field, central_field, upwind_field
  use kinds, only: ccs_int

  use utils, only: update
  use boundary_conditions, only: read_bc_config, allocate_bc_arrays
  use vec, only: create_vector
  use fv, only: update_gradient
  use timestepping, only: initialise_old_values
  
  implicit none

  private

  public :: create_field
  
contains

  !> Build a field variable with data and gradient vectors + transient data and boundary arrays.
  subroutine create_field(ccs_config_file, vec_properties, field_type, field_name, n_boundaries, phi)

    use utils, only : debug_print
    
    implicit none
    
    !! Logically vec_properties should be a field_properties variable, but this doesn't yet exist.
    character(len=*), intent(in) :: ccs_config_file !< String naming the config file
    type(vector_spec), intent(in) :: vec_properties !< Vector descriptor for vectors wrapped by field
    integer, intent(in) :: field_type               !< Identifier for what kind of field
    character(len=*), intent(in) :: field_name      !< Field name -- should match against boundary conditions, etc.
    integer(ccs_int), intent(in) :: n_boundaries    !< Mesh boundary count
    class(field), allocatable, intent(out) :: phi   !< The field being constructed

    call allocate_field(vec_properties, field_type, n_boundaries, phi)

    ! XXX: ccs_config_file is host-associated from program scope.
    call read_bc_config(ccs_config_file, field_name, phi)

    !! --- Ensure data is updated/parallel-constructed ---
    ! XXX: Potential abstraction --- see update(vec), etc.
    call update(phi%values)
    if (field_type /= face_centred) then
       ! Current design only computes/stores gradients at cell centres
       call update(phi%x_gradients)
       call update(phi%y_gradients)
       call update(phi%z_gradients)

       call update_gradient(vec_properties%mesh, phi)
    end if
    !! --- End update ---
    
  end subroutine create_field

  !> Allocate a field variable
  subroutine allocate_field(vec_properties, field_type, n_boundaries, phi)

    use utils, only : debug_print
    
    implicit none
    
    !! Logically vec_properties should be a field_properties variable, but this doesn't yet exist.
    type(vector_spec), intent(in) :: vec_properties !< Vector descriptor for vectors wrapped by field
    integer, intent(in) :: field_type               !< Identifier for what kind of field
    integer(ccs_int), intent(in) :: n_boundaries    !< Mesh boundary count
    class(field), allocatable, intent(out) :: phi   !< The field being constructed
    
    if (field_type == face_centred) then
       call dprint("Create face field")
       allocate (face_field :: phi)
    else if (field_type == cell_centred_upwind) then
       call dprint("Create upwind field")
       allocate (upwind_field :: phi)
    else if (field_type == cell_centred_central) then
       call dprint("Create central field")
       allocate (central_field :: phi)
    end if

    call dprint("Create field values vector")
    call create_vector(vec_properties, phi%values)

    if (field_type /= face_centred) then
       ! Current design only computes/stores gradients at cell centres
       call dprint("Create field gradients vector")
       call create_vector(vec_properties, phi%x_gradients)
       call create_vector(vec_properties, phi%y_gradients)
       call create_vector(vec_properties, phi%z_gradients)

       ! Currently no need for old face values
       call dprint("Create field old values")
       call initialise_old_values(vec_properties, phi)
    end if

    call allocate_bc_arrays(n_boundaries, phi%bcs)

  end subroutine allocate_field
  
end module fields
