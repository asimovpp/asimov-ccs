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
        pure subroutine get_init_flow(x_p, field_name, init_val)
          use kinds, only: ccs_real
          use constants, only: ndim
          real(ccs_real), dimension(ndim), intent(in) :: x_p
          character(len=*), intent(in) :: field_name
          real(ccs_real), intent(inout) :: init_val
        end subroutine
      end interface
    end subroutine

    !> Initialise mass flux field using calling get_init_mass_flux
    module subroutine core_initialise_mass_flux(flow_fields, get_init_mass_flux)
      type(fluid), intent(inout) :: flow_fields
      interface
        pure subroutine get_init_mass_flux(index_f, init_flux)
          use kinds, only: ccs_int, ccs_real
          integer(ccs_int), intent(in) :: index_f
          real(ccs_real), intent(inout) :: init_flux
        end subroutine
      end interface
    end subroutine

  end interface
  
end module core
