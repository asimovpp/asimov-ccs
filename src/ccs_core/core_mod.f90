!v Module file core.mod
!
! The core module implements the core functionality to define and run a CCS case.

module core

  use kinds, only: ccs_int, ccs_real
  use types, only: fluid

  implicit none

  private

  public :: core_initialise_flow
  public :: core_initialise_mass_flux

  interface

    !> Initialise every cell based field prompting values from get_init_flow
    module subroutine core_initialise_flow(flow_fields, get_init_flow)
      type(fluid), intent(inout) :: flow_fields
      interface 
        pure subroutine get_init_flow(loc_p, field_name, init_val)
          use kinds, only: ccs_real
          use types, only: cell_locator
          type(cell_locator), intent(in) :: loc_p
          character(len=*), intent(in) :: field_name
          real(ccs_real), intent(inout) :: init_val
        end subroutine
      end interface
    end subroutine

    !> Initialise mass flux field using calling get_init_mass_flux
    module subroutine core_initialise_mass_flux(flow_fields, get_init_mass_flux)
      type(fluid), intent(inout) :: flow_fields
      interface
        pure subroutine get_init_mass_flux(loc_f, init_val)
          use kinds, only: ccs_real
          use types, only: face_locator
          type(face_locator), intent(in) :: loc_f
          real(ccs_real), intent(inout) :: init_val
        end subroutine
      end interface
    end subroutine

  end interface
  
end module core
