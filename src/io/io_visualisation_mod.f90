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

  module subroutine write_solution(par_env, case_name, step, maxstep, dt, mesh, output_field_list)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    character(len=:), allocatable, intent(in) :: case_name
    integer(ccs_int), intent(in) :: step
    integer(ccs_int), intent(in) :: maxstep
    real(ccs_real), intent(in) :: dt
    type(ccs_mesh), intent(in) :: mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list
  end subroutine

  module subroutine write_fields(par_env, case_name, step, maxstep, dt, mesh, output_field_list)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    character(len=:), allocatable, intent(in) :: case_name
    integer(ccs_int), intent(in) :: step
    integer(ccs_int), intent(in) :: maxstep
    real(ccs_real), intent(in) :: dt
    type(ccs_mesh), intent(in) :: mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list
  end subroutine

  module subroutine write_xdmf(par_env, case_name, step, maxstep, dt, mesh, output_field_list)
    class(parallel_environment), allocatable, target, intent(in) :: par_env
    character(len=:), allocatable, intent(in) :: case_name
    integer(ccs_int), intent(in) :: step
    integer(ccs_int), intent(in) :: maxstep
    real(ccs_real), intent(in) :: dt
    type(ccs_mesh), intent(in) :: mesh
    type(output_list), dimension(:), intent(inout) :: output_field_list
  end subroutine

  end interface

end module io_visualisation