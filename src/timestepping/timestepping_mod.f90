!v Module file timestepping_mod.f90
!
!  Provides the interface for timestepping methods.
module timestepping

  use kinds, only: ccs_real, ccs_int
  use types, only: field, ccs_mesh, ccs_vector, ccs_matrix, vector_spec

  implicit none

  private

  public :: apply_timestep
  public :: set_timestep
  public :: get_timestep
  public :: update_old_values
  public :: activate_timestepping
  public :: initialise_old_values

  interface
    !> Apply one timestep correction
    module subroutine apply_timestep(mesh, phi, diag, M, b)
      type(ccs_mesh), intent(in) :: mesh !< mesh object
      class(field), intent(in) :: phi !< flow variable
      class(ccs_vector), intent(inout) :: diag !< preallocated vector with the same size as M diagonal
      class(ccs_matrix), intent(inout) :: M !< equation system
      class(ccs_vector), intent(inout) :: b !< rhs vector
    end subroutine

    !> Set timestep size
    module subroutine set_timestep(timestep)
      real(ccs_real), intent(in) :: timestep
    end subroutine

    !> Get timestep size
    module function get_timestep() result(timestep)
      real(ccs_real) :: timestep
    end function

    !> Place current field values into old field values
    module subroutine update_old_values(x)
      class(field), intent(inout) :: x
    end subroutine

    !v Enable timestepping
    !
    !  If this is not called at the start of a program,
    !  calls to apply_timestep will return without doing anything.
    module subroutine activate_timestepping()
    end subroutine

    !> Create vectors for storing one or more previous timesteps.
    module subroutine initialise_old_values(vec_properties, x)
      type(vector_spec), intent(in) :: vec_properties
      class(field), intent(inout) :: x
    end subroutine

  end interface

end module timestepping
