!>  boundary conditions module
!
!>  Various BC related functionality. Need to expand.

module boundary_conditions
#include "ccs_macros.inc"

  use utils, only: exit_print
  use types, only: bc_config, field
  use kinds, only: ccs_int, ccs_real
  use yaml, only: parse, error_length
  use read_config, only: get_bc_field_data
  use bc_constants
  
  implicit none

  private
  public :: read_bc_config
  public :: set_bc_attribute
  public :: allocate_bc_field
  public :: get_bc_index

  interface set_bc_attribute
    module procedure set_bc_real_attribute
    module procedure set_bc_string_attribute
  end interface

  contains

  !>  Reads config file and assigns data to BC structure
  subroutine read_bc_config(filename, bc_field, bcs) 
    character(len=*), intent(in) :: filename  !< name of the config file
    character(len=*), intent(in) :: bc_field  !< string denoting which field we want to read in
    type(bc_config), intent(inout) :: bcs     !< the bc struct of the corresponding field

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
  
  !> Sets the appropriate integer values for strings with given by the key-value pair attribute, value
  subroutine set_bc_string_attribute(boundary_index, attribute, value, bcs)
    integer(ccs_int), intent(in) :: boundary_index  !< Index of the boundary within bcs struct arrays 
    character(len=*), intent(in) :: attribute       !< string giving the attribute name
    character(len=*), intent(in) :: value           !< string giving the attribute value
    type(bc_config), intent(inout) :: bcs           !< bcs struct

    select case (attribute)
    case ("name")
      select case (value)
      case ("left")
        bcs%names(boundary_index) = bc_region_left
      case ("right")
        bcs%names(boundary_index) = bc_region_right
      case ("bottom")
        bcs%names(boundary_index) = bc_region_bottom
      case ("top")
        bcs%names(boundary_index) = bc_region_top
      case ("jet")
        bcs%names(boundary_index) = bc_region_jet
      case ("coflow")
        bcs%names(boundary_index) = bc_region_coflow
      case ("outflow")
        bcs%names(boundary_index) = bc_region_outflow
      case ("atmos")
        bcs%names(boundary_index) = bc_region_atmos
      end select
    case ("type")
      select case (value)
      case ("periodic")
        bcs%bc_types(boundary_index) = bc_type_periodic
      case ("sym")
        bcs%bc_types(boundary_index) = bc_type_sym
      case ("dirichlet")
        bcs%bc_types(boundary_index) = bc_type_dirichlet
      case ("const_grad")
        bcs%bc_types(boundary_index) = bc_type_const_grad
      case ("inlet")
        bcs%bc_types(boundary_index) = bc_type_inlet
      case ("outlet")
        bcs%bc_types(boundary_index) = bc_type_outlet
      case ("symp")
        bcs%bc_types(boundary_index) = bc_type_symp
      case ("wall")
        bcs%bc_types(boundary_index) = bc_type_wall
      end select
    case default
      call error_abort("invalid bc attribute")
    end select

  end subroutine set_bc_string_attribute
  
  !> Sets the bc struct's value field to the given real value
  subroutine set_bc_real_attribute(boundary_index, value, bcs)
    integer(ccs_int), intent(in) :: boundary_index  !< index of the boundary within the bc struct's arrays
    real(ccs_real), intent(in) :: value             !< the value to set 
    type(bc_config), intent(inout) :: bcs           !< the bcs struct

    bcs%values(boundary_index) = value
  end subroutine set_bc_real_attribute

  !> Allocates arrays of the appropriate size for the name, type and value of the bcs
  subroutine allocate_bc_field(field, n_boundaries, bcs)
    character(len=*), intent(in) :: field         !< string denoting which array to allocate within the bc struct
    integer(ccs_int), intent(in) :: n_boundaries  !< the number of boundaries 
    type(bc_config), intent(inout) :: bcs         !< the bc struct

    select case (field)
    case ("name")
      if (.not. allocated(bcs%names)) then
        allocate(bcs%names(n_boundaries))
      end if
    case ("type")
      if (.not. allocated(bcs%bc_types)) then
        allocate(bcs%bc_types(n_boundaries))
      end if
    case default
      if (.not. allocated(bcs%values)) then
        allocate(bcs%values(n_boundaries))
      end if
    end select
  end subroutine allocate_bc_field

  !> Gets the index of the given boundary condition within the bc struct arrays
  subroutine get_bc_index(phi, index_nb, index_bc)
    class(field), intent(in) :: phi             !< The field whose bc we're getting
    integer(ccs_int), intent(in) :: index_nb    !< The index of the neighbouring boundary cell
    integer(ccs_int), intent(out) :: index_bc   !< The index of the appropriate boundary in the bc struct

    ! XXX: There might be a better way of doing this on the assumption that the boundary condition labels are ordered.
    do index_bc = 1, size(phi%bcs%names)
      if (phi%bcs%names(index_bc) == index_nb) then
        exit
      end if
    end do
  end subroutine get_bc_index
end module boundary_conditions
