!v Module file io_visualisation_mod.f90
!
!  Provides an interface to visualisation-related IO functions.
module io_visualisation

  use types, only: io_environment, io_process, ccs_mesh, fluid
  use parallel_types, only: parallel_environment
  use kinds, only: ccs_int, ccs_real

  implicit none

  private

  public :: reset_io_visualisation
  public :: reset_io_visualisation_module
  public :: write_solution
  public :: write_fields
  public :: write_xdmf

  interface

    !> Reset module to its original stage
    module subroutine reset_io_visualisation_module()
    end subroutine

    !> Reset module to its original stage
    module subroutine reset_io_visualisation()
    end subroutine

    !> Write the flow solution for the current time-step to file
    module subroutine write_solution(par_env, case_name, mesh, flow, step, maxstep, dt)
      class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
      character(len=:), allocatable, intent(in) :: case_name                   !< The case name
      type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
      type(fluid), intent(inout) :: flow                                       !< The flow variables
      integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
      integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count
      real(ccs_real), optional, intent(in) :: dt                               !< The time-step size
    end subroutine

    !> Write the field data to file
    module subroutine write_fields(par_env, case_name, mesh, flow, step, maxstep)
      class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
      character(len=:), allocatable, intent(in) :: case_name                   !< The case name
      type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
      type(fluid), intent(inout) :: flow                                       !< The flow variables
      integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
      integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count
    end subroutine

    !> Write the XML descriptor of the field data and grid
    module subroutine write_xdmf(par_env, case_name, flow, step, maxstep, dt)
      class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
      character(len=:), allocatable, intent(in) :: case_name                   !< The case name
      type(fluid), intent(inout) :: flow                                       !< The flow variables
      integer(ccs_int), optional, intent(in) :: step                           !< The current time-step count
      integer(ccs_int), optional, intent(in) :: maxstep                        !< The maximum time-step count
      real(ccs_real), optional, intent(in) :: dt                               !< The time-step size
    end subroutine

  end interface

end module io_visualisation
