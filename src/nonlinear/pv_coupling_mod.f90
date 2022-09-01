!v Module file pv_coupling.mod
!
!  An interface to pressure-velocity coupling methods (SIMPLE, etc)

module pv_coupling

  use kinds, only: ccs_int
  use types, only: field, ccs_mesh
  use parallel_types, only: parallel_environment

  implicit none

  private

  public :: solve_nonlinear

  interface

    module subroutine solve_nonlinear(par_env, mesh, cps, it_start, it_end, u, v, w, p, p_prime, mf)
      class(parallel_environment), allocatable, intent(in) :: par_env
      type(ccs_mesh), intent(in) :: mesh
      integer(ccs_int), intent(in) :: cps, it_start, it_end
      class(field), intent(inout) :: u, v, w, p, p_prime, mf
    end subroutine solve_nonlinear

  end interface

end module pv_coupling
