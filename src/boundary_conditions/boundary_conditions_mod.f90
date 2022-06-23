!>  boundary conditions module
!
!>  Various BC related functionality. Need to expand.

module boundary_conditions
#include "ccs_macros.inc"

  use utils, only: exit_print, debug_print, str
  use types, only: bc_config
  use kinds, only: ccs_int, ccs_real
  
  implicit none

  private
  public :: read_bc_config
  public :: set_bc_attribute
  public :: allocate_bc_field

  interface set_bc_attribute
    module procedure set_bc_real_attribute
    module procedure set_bc_string_attribute
  end interface

  contains

  !>  Wrapper for reading config file and assigning data to BC structure
  !
  !> @param[in] filename - name of config file
  !> @param[out] bcs     - boundary conditions structure
  subroutine read_bc_config(filename, bc_field, bcs) 
    use yaml, only: parse, error_length
    use read_config, only: get_bc_field_data
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: bc_field
    type(bc_config), intent(inout) :: bcs

    class(*), pointer :: config_file
    character(len=error_length) :: error

    config_file => parse(filename, error=error)
    if (error/='') then
      call error_abort(trim(error))
    endif

    call get_bc_field_data(config_file, "name", bcs)
    call get_bc_field_data(config_file, "type", bcs)
    call get_bc_field_data(config_file, bc_field, bcs)
  end subroutine read_bc_config
  
  subroutine set_bc_string_attribute(boundary_index, attribute, value, bcs)
    integer(ccs_int), intent(in) :: boundary_index
    character(len=*), intent(in) :: attribute
    character(len=*), intent(in) :: value
    type(bc_config), intent(inout) :: bcs

    select case (attribute)
    case ("name")
      call dprint("name size " // str(size(bcs%name)) // " " // str(boundary_index) // " " // value)
      select case (value)
      case ("left")
        bcs%name(boundary_index) = -1
      case ("right")
        bcs%name(boundary_index) = -2
      case ("bottom")
        bcs%name(boundary_index) = -3
      case ("top")
        bcs%name(boundary_index) = -4
      case ("jet")
        bcs%name(boundary_index) = -5
      case ("coflow")
        bcs%name(boundary_index) = -6
      case ("outflow")
        bcs%name(boundary_index) = -7
      case ("atmos")
        bcs%name(boundary_index) = -8
      end select
    case ("type")
      select case (value)
      case ("inlet")
        bcs%bc_type(boundary_index) = 1
      case ("outlet")
        bcs%bc_type(boundary_index) = 2
      case ("symp")
        bcs%bc_type(boundary_index) = 3
      case ("wall")
        bcs%bc_type(boundary_index) = 4
      end select
    case default
      call dprint("invalid bc attribute")
    end select

  end subroutine set_bc_string_attribute
  
  subroutine set_bc_real_attribute(boundary_index, value, bcs)
    integer(ccs_int), intent(in) :: boundary_index
    real(ccs_real), intent(in) :: value
    type(bc_config), intent(inout) :: bcs

    bcs%value(boundary_index) = value
  end subroutine set_bc_real_attribute

  subroutine allocate_bc_field(field, n_boundaries, bcs)
    character(len=*), intent(in) :: field
    integer(ccs_int), intent(in) :: n_boundaries
    type(bc_config), intent(inout) :: bcs

    select case (field)
    case ("name")
      if (.not. allocated(bcs%name)) then
        allocate(bcs%name(n_boundaries))
        call dprint("allocated name")
      end if
    case ("type")
      if (.not. allocated(bcs%bc_type)) then
        allocate(bcs%bc_type(n_boundaries))
        call dprint("allocated type")
      end if
    case default
      if (.not. allocated(bcs%value)) then
        allocate(bcs%value(n_boundaries))
        call dprint("allocated value")
      end if
    end select
  end subroutine allocate_bc_field
end module boundary_conditions
