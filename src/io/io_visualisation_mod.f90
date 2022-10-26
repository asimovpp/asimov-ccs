!v Module file io_visualisation_mod.f90
!
!  Provides an interface to visualisation-related IO functions.
module io_visualisation

  use types, only: io_environment, io_process, ccs_mesh, output_list
  use parallel_types, only: parallel_environment
  use kinds, only: ccs_int, ccs_real

  implicit none

  private

  public :: write_solution
  public :: write_fields
  public :: write_xdmf

  interface

  !> Write the flow solution for the current time-step to file
  module subroutine write_solution(par_env, case_name, step, maxstep, dt, mesh, output_field_list)
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    integer(ccs_int), intent(in) :: step                                     !< The current time-step count
    integer(ccs_int), intent(in) :: maxstep                                  !< The maximum time-step count
    real(ccs_real), intent(in) :: dt                                         !< The time-step size
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list      !< List of fields to output
  end subroutine

  !> Write the field data to file
  module subroutine write_fields(par_env, case_name, step, maxstep, dt, mesh, output_field_list)
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    integer(ccs_int), intent(in) :: step                                     !< The current time-step count
    integer(ccs_int), intent(in) :: maxstep                                  !< The maximum time-step count
    real(ccs_real), intent(in) :: dt                                         !< The time-step size
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list      !< List of fields to output
  end subroutine

  !> Write the XML descriptor of the field data and grid
  module subroutine write_xdmf(par_env, case_name, step, maxstep, dt, mesh, output_field_list)
    class(parallel_environment), allocatable, target, intent(in) :: par_env  !< The parallel environment
    character(len=:), allocatable, intent(in) :: case_name                   !< The case name
    integer(ccs_int), intent(in) :: step                                     !< The current time-step count
    integer(ccs_int), intent(in) :: maxstep                                  !< The maximum time-step count
    real(ccs_real), intent(in) :: dt                                         !< The time-step size
    type(ccs_mesh), intent(in) :: mesh                                       !< The mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list      !< List of fields to output
  end subroutine

  end interface

end module io_visualisation