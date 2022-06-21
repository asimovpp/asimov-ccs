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
  subroutine read_bc_config(filename, bcs) 
    use yaml, only: parse, error_length
    use read_config, only: get_bc_field_data
    character(len=*), intent(in) :: filename
    type(bc_config), intent(inout) :: bcs

    class(*), pointer :: config_file
    character(len=error_length) :: error

    config_file => parse(filename, error=error)
    if (error/='') then
      call error_abort(trim(error))
    endif

    call get_bc_field_data(config_file, "name", bcs)
    call get_bc_field_data(config_file, "type", bcs)
    call get_bc_field_data(config_file, "u", bcs)
    call get_bc_field_data(config_file, "v", bcs)
    call get_bc_field_data(config_file, "w", bcs)
    call get_bc_field_data(config_file, "den", bcs)
    call get_bc_field_data(config_file, "T", bcs)
    call get_bc_field_data(config_file, "fmix", bcs)
    call get_bc_field_data(config_file, "pgr", bcs)
    call get_bc_field_data(config_file, "ti", bcs)
    call get_bc_field_data(config_file, "ls", bcs)
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
      case ("jet")
        bcs%name(boundary_index) = 1
      case ("coflow")
        bcs%name(boundary_index) = 2
      case ("outflow")
        bcs%name(boundary_index) = 3
      case ("atmos")
        bcs%name(boundary_index) = 4
      end select
    case ("type")
      select case (value)
      case ("inlet")
        bcs%bc_type(boundary_index) = 1
      case ("outlet")
        bcs%bc_type(boundary_index) = 2
      case ("symp")
        bcs%bc_type(boundary_index) = 3
      end select
    case default
      call dprint("invalid bc attribute")
    end select

  end subroutine set_bc_string_attribute
  
  subroutine set_bc_real_attribute(boundary_index, attribute, value, bcs)
    integer(ccs_int), intent(in) :: boundary_index
    character(len=*), intent(in) :: attribute
    real(ccs_real), intent(in) :: value
    type(bc_config), intent(inout) :: bcs

    select case (attribute)
    case ("u")
      bcs%u(boundary_index) = value
    case ("v")
      bcs%v(boundary_index) = value
    case ("w")
      bcs%w(boundary_index) = value
    case ("den")
      bcs%den(boundary_index) = value
    case ("T")
      bcs%T(boundary_index) = value
    case ("fmix")
      bcs%fmix(boundary_index) = value
    case ("pgr")
      bcs%pgr(boundary_index) = value
    case ("ti")
      bcs%ti(boundary_index) = value
    case ("ls")
      bcs%ls(boundary_index) = value
    case default
      call dprint("invalid bc attribute " // attribute)
    end select

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
    case ("u")
      if (.not. allocated(bcs%u)) then
        allocate(bcs%u(n_boundaries))
        call dprint("allocated u")
      end if
    case ("v")
      if (.not. allocated(bcs%v)) then
        allocate(bcs%v(n_boundaries))
        call dprint("allocated v")
      end if
    case ("w")
      if (.not. allocated(bcs%w)) then
        allocate(bcs%w(n_boundaries))
        call dprint("allocated w")
      end if
    case ("den")
      if (.not. allocated(bcs%den)) then
        allocate(bcs%den(n_boundaries))
        call dprint("allocated den")
      end if
    case ("T")
      if (.not. allocated(bcs%T)) then
        allocate(bcs%T(n_boundaries))
        call dprint("allocated T")
      end if
    case ("fmix")
      if (.not. allocated(bcs%fmix)) then
        allocate(bcs%fmix(n_boundaries))
        call dprint("allocated fmix")
      end if
    case ("pgr")
      if (.not. allocated(bcs%pgr)) then
        allocate(bcs%pgr(n_boundaries))
        call dprint("allocated pgr")
      end if
    case ("ti")
      if (.not. allocated(bcs%ti)) then
        allocate(bcs%ti(n_boundaries))
        call dprint("allocated ti")
      end if
    case ("ls")
      if (.not. allocated(bcs%ls)) then
        allocate(bcs%ls(n_boundaries))
        call dprint("allocated ls")
      end if
    end select
  end subroutine allocate_bc_field
end module boundary_conditions
