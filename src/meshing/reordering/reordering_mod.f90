!v Module file reordering_mod.f90
!
!  Provides the interface for mesh reordering methods.
module reordering

  use kinds, only: ccs_int
  use types, only: ccs_mesh
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: reorder_cells
  public :: compute_bandwidth
  public :: disable_reordering

  interface

    !> Get and apply new order to mesh cells.
    module subroutine reorder_cells(par_env, shared_env, mesh)
      class(parallel_environment), intent(in) :: par_env !< The parallel environment
      class(parallel_environment), intent(in) :: shared_env !< The parallel environment
      type(ccs_mesh), intent(inout) :: mesh
    end subroutine

    !> Given a mesh reordering, apply it.
    module subroutine apply_reordering(new_indices, par_env, shared_env, mesh)
      integer(ccs_int), dimension(:), intent(in) :: new_indices
      class(parallel_environment), intent(in) :: par_env !< The parallel environment
      class(parallel_environment), intent(in) :: shared_env !< The parallel environment
      type(ccs_mesh), intent(inout) :: mesh
    end subroutine

    !> Calculate and print the bandwidth of a mesh
    module subroutine compute_bandwidth(mesh)
      type(ccs_mesh), intent(in) :: mesh
    end subroutine

    !> Generate a mesh cell reordering mapping.
    module subroutine get_reordering(mesh, new_indices)
      type(ccs_mesh), intent(in) :: mesh
      integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices
    end subroutine

    !> Disable reordering, turning calls to "reorder_cells" no-ops.
    module subroutine disable_reordering()
    end subroutine

  end interface

end module reordering
