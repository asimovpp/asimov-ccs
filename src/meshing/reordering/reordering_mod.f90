!v Module file reordering_mod.f90
!
!  Provides the interface for mesh reordering methods.
module reordering

  use kinds, only: ccs_int, ccs_real
  use types, only: ccs_mesh
  use parallel_types, only: parallel_environment

  implicit none

  private
  public :: reorder_cells
  public :: disable_reordering
  public :: print_bandwidth

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

    !> Print statistics on the local matrix bandwidth
    module subroutine print_bandwidth(par_env, mesh)
      class(parallel_environment), allocatable, target, intent(in) :: par_env !< The parallel environment
      type(ccs_mesh), intent(in) :: mesh
    end subroutine print_bandwidth
    
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
