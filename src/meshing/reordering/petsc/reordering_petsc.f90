submodule(reordering) reordering_petsc

  use kinds, only: ccs_real, ccs_err
  use types, only: cell_locator, neighbour_locator

  implicit none

contains

  !v Determine how the mesh should be reordered using PETSc reordering
  module subroutine get_reordering(mesh, new_indices)
#include "petsc/finclude/petscmat.h"

    use petsc, only: PETSC_DETERMINE, PETSC_NULL_INTEGER, INSERT_VALUES
    use petscmat
    use petscis, only: tIS, ISGetIndicesF90, ISRestoreIndicesF90, ISDestroy

    type(ccs_mesh), intent(in) :: mesh                                      !< the mesh to be reordered
    integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices !< new indices in "to(from)" format

    type(tMat) :: M
    integer(ccs_err) :: ierr
    type(tIS) :: rperm, cperm
    integer(ccs_int), pointer :: row_indices(:)

    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i, j
    integer(ccs_int) :: ctr, nnb
    integer(ccs_int), dimension(7) :: idx ! Hardcoded to 7-point stencil
    real(ccs_real), dimension(7) :: row   ! Hardcoded to 7-point stencil
    logical :: cell_local
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: idx_new

    ! First build adjacency matrix for local cells
    call MatCreate(MPI_COMM_SELF, M, ierr)
    call MatSetFromOptions(M, ierr)
    call get_local_num_cells(mesh, local_num_cells)
    call MatSetSizes(M, local_num_cells, local_num_cells, &
                     PETSC_DETERMINE, PETSC_DETERMINE, ierr)
    call MatSeqAIJSetPreallocation(M, 7, PETSC_NULL_INTEGER, ierr)
    do i = 1, local_num_cells
      row(:) = 0.0
      idx(:) = 0
      ctr = 1

      row(ctr) = 1.0
      idx(ctr) = i
      ctr = ctr + 1

      call set_cell_location(mesh, i, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call set_neighbour_location(loc_p, j, loc_nb)
        call get_local_status(loc_nb, cell_local)
        if (cell_local) then
          call get_local_index(loc_nb, idx(ctr))
          row(ctr) = 1.0
          ctr = ctr + 1
        end if
      end do
      idx = idx - 1 ! F->C
      call MatSetValues(M, 1, i - 1, 7, idx, row, INSERT_VALUES, ierr)
    end do
    call MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr)

    ! Get index sets for reordering
    call MatGetOrdering(M, MATORDERINGRCM, rperm, cperm, ierr)
    call MatDestroy(M, ierr)
    call ISDestroy(cperm, ierr)

    ! Fill local indices in original ordering -> destination, i.e. to(i) => new index of cell i.
    allocate (new_indices(local_num_cells))

    call ISGetIndicesF90(rperm, row_indices, ierr)
    if (local_num_cells >= 1) then
      do i = 1, local_num_cells
        idx_new = row_indices(i) + 1 ! C->F
        new_indices(idx_new) = i
      end do
    end if

    call ISRestoreIndicesF90(rperm, row_indices, ierr)
    call ISDestroy(rperm, ierr)
  end subroutine get_reordering

end submodule
