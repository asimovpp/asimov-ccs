!v boundary conditions module
!
!  Various BC related functionality. Need to expand.

module boundary_conditions
#include "ccs_macros.inc"

  use utils, only: exit_print, debug_print, str
  use types, only: bc_config, field, bc_profile
  use kinds, only: ccs_int, ccs_real
  use fortran_yaml_c_interface, only: parse
  use read_config, only: get_bc_field
  use bc_constants
  use error_codes

  implicit none

  private
  public :: read_bc_config
  public :: allocate_bc_arrays
  public :: get_bc_index
  public :: set_bc_real_value
  public :: set_bc_type
  public :: set_bc_id
  public :: set_bc_profile

contains

  !> Reads config file and assigns data to BC structure
  subroutine read_bc_config(filename, bc_field, phi)
    character(len=*), intent(in) :: filename !< name of the config file
    character(len=*), intent(in) :: bc_field !< string denoting which field we want to read in
    class(field), intent(inout) :: phi       !< the bc struct of the corresponding field

    class(*), pointer :: config_file
    character(:), allocatable :: error

    config_file => parse(filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    call dprint("reading bc config " // bc_field)
    call get_bc_field(config_file, "name", phi)
    call get_bc_field(config_file, "type", phi, required=.false.)
    call get_bc_field(config_file, "value", phi, required=.false.)
    call get_bc_field(config_file, bc_field, phi)
  end subroutine read_bc_config

  !> Sets the appropriate integer values for strings with given by the key-value pair attribute, value
  pure subroutine set_bc_type(boundary_index, bc_type, bcs)
    integer(ccs_int), intent(in) :: boundary_index !< Index of the boundary within bcs struct arrays
    character(len=*), intent(in) :: bc_type        !< string giving the bc type
    type(bc_config), intent(inout) :: bcs          !< bcs struct

    select case (bc_type)
    case ("periodic")
      bcs%bc_types(boundary_index) = bc_type_periodic
    case ("sym")
      bcs%bc_types(boundary_index) = bc_type_sym
    case ("dirichlet")
      bcs%bc_types(boundary_index) = bc_type_dirichlet
    case ("neumann")
      bcs%bc_types(boundary_index) = bc_type_neumann
    case ("extrapolate")
      bcs%bc_types(boundary_index) = bc_type_extrapolate
    case ("wall")
      bcs%bc_types(boundary_index) = bc_type_wall
    case ("profile")
      bcs%bc_types(boundary_index) = bc_type_profile
    case default
      error stop invalid_bc_name ! Invalid BC type string received
    end select

  end subroutine set_bc_type

  !> Sets the bc struct's id field to the appropriate integer value
  pure subroutine set_bc_id(boundary_index, name, bcs)
    integer(ccs_int), intent(in) :: boundary_index !< index of the boundary within the bc struct's arrays
    character(len=*), intent(in) :: name           !< string giving the bc name
    type(bc_config), intent(inout) :: bcs          !< the bcs struct

    ! XXX: in the general case this mapping should be read in from the mesh file
    select case (name)
    case ("left")
      bcs%ids(boundary_index) = 1
    case ("right")
      bcs%ids(boundary_index) = 2
    case ("bottom")
      bcs%ids(boundary_index) = 3
    case ("top")
      bcs%ids(boundary_index) = 4
    case ("back")
      bcs%ids(boundary_index) = 5
    case ("front")
      bcs%ids(boundary_index) = 6
    case ("other1")
      bcs%ids(boundary_index) = 7
    case ("other2")
      bcs%ids(boundary_index) = 8
    case ("other3")
      bcs%ids(boundary_index) = 9
    case ("other4")
      bcs%ids(boundary_index) = 10
    case ("other5")
      bcs%ids(boundary_index) = 11
    case default
      error stop invalid_bc_id ! Invalid BC name ID received
    end select

  end subroutine set_bc_id

  !> Sets the bc struct's value field to the given real value
  pure subroutine set_bc_real_value(boundary_index, val, bcs)
    integer(ccs_int), intent(in) :: boundary_index !< index of the boundary within the bc struct's arrays
    real(ccs_real), intent(in) :: val              !< the value to set
    type(bc_config), intent(inout) :: bcs          !< the bcs struct

    bcs%values(boundary_index) = val
  end subroutine set_bc_real_value

  !> Allocates arrays of the appropriate size for the name, type and value of the bcs
  pure subroutine allocate_bc_arrays(n_boundaries, bcs)
    integer(ccs_int), intent(in) :: n_boundaries !< the number of boundaries
    type(bc_config), intent(inout) :: bcs        !< the bc struct

    if (.not. allocated(bcs%ids)) then
      allocate (bcs%ids(n_boundaries))
    end if
    if (.not. allocated(bcs%bc_types)) then
      allocate (bcs%bc_types(n_boundaries))
    end if
    if (.not. allocated(bcs%values)) then
      allocate (bcs%values(n_boundaries))
    end if
    if (.not. allocated(bcs%profiles)) then
      allocate (bcs%profiles(n_boundaries))
    end if
  end subroutine allocate_bc_arrays

  !> Gets the index of the given boundary condition within the bc struct arrays
  pure subroutine get_bc_index(phi, index_nb, index_bc)
    class(field), intent(in) :: phi           !< The field whose bc we're getting
    integer(ccs_int), intent(in) :: index_nb  !< The index of the neighbouring boundary cell
    integer(ccs_int), intent(out) :: index_bc !< The index of the appropriate boundary in the bc struct

    ! Local variable
    integer(ccs_int), dimension(1) :: index_tmp ! The intrinsic returns a rank-1 array ...

    index_tmp = findloc(phi%bcs%ids, -index_nb)
    if (index_tmp(1) == 0) then
      error stop bc_index_not_found ! BC index not found
    end if
    
    index_bc = index_tmp(1)
  end subroutine get_bc_index

  !> Set boundary condition profile to the index_bc boundary
  pure subroutine set_bc_profile(phi, profile, index_bc)
    class(field), intent(inout) :: phi           !< The field whose profile we are setting
    type(bc_profile), intent(in) :: profile      !< BC profile
    integer(ccs_int), intent(in) :: index_bc     !< The index of the appropriate boundary in the bc struct

    phi%bcs%profiles(index_bc) = profile

  end subroutine

end module boundary_conditions
