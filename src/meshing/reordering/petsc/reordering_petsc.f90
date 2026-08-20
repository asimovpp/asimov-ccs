submodule(reordering) reordering_petsc
#include "ccs_macros.inc"

  use utils, only: debug_print
  use kinds, only: ccs_real, ccs_err
  use types, only: cell_locator, neighbour_locator
  use meshing, only: create_cell_locator, create_neighbour_locator, &
                     count_neighbours, &
                     get_local_status, get_local_index, &
                     get_local_num_cells

  implicit none

contains

  !v Determine how the mesh should be reordered using PETSc reordering
  module subroutine get_reordering(new_indices)
#include "petsc/finclude/petscmat.h"

    use mpi
    use petscmat
#if PETSC_VERSION_GE(3,23,0)
    use petsc, only: PETSC_NULL_INTEGER_ARRAY, PETSC_DETERMINE, INSERT_VALUES
    use petscis, only: ISGetIndices, ISRestoreIndices, tIS, ISDestroy
#else
    use petsc, only: PETSC_NULL_INTEGER, PETSC_DETERMINE, INSERT_VALUES
    use petscis, only: ISGetIndicesF90, ISRestoreIndicesF90, tIS, ISDestroy
#endif

    integer(ccs_int), dimension(:), allocatable, intent(out) :: new_indices !< new indices in "to(from)" format

    type(tMat) :: M
    integer(ccs_err) :: ierr
    type(tIS) :: rperm, cperm
    integer(ccs_int), pointer :: row_indices(:)

    integer(ccs_int) :: local_num_cells

    integer(ccs_int) :: i, j
    integer(ccs_int) :: ctr, nnb, max_nb
    integer(ccs_int), allocatable, dimension(:) :: idx
    real(ccs_real), allocatable, dimension(:) :: row
    logical :: cell_local
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb

    integer(ccs_int) :: idx_new

    call dprint("Reordering with PETSc.")

    call get_local_num_cells(local_num_cells)

    max_nb = 0
    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)
      max_nb = max(max_nb, nnb)
    end do

    allocate (idx(max_nb + 1))
    allocate (row(max_nb + 1))

    ! First build adjacency matrix for local cells
    call MatCreate(MPI_COMM_SELF, M, ierr)
    call MatSetFromOptions(M, ierr)
    call MatSetSizes(M, local_num_cells, local_num_cells, &
                     PETSC_DETERMINE, PETSC_DETERMINE, ierr)
#if PETSC_VERSION_GE(3,23,0)
    call MatSeqAIJSetPreallocation(M, max_nb, PETSC_NULL_INTEGER_ARRAY, ierr)
#else
    call MatSeqAIJSetPreallocation(M, max_nb, PETSC_NULL_INTEGER, ierr)
#endif

    do i = 1, local_num_cells
      row(:) = 0.0
      idx(:) = 0
      ctr = 1

      row(ctr) = 1.0
      idx(ctr) = i
      ctr = ctr + 1

      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_status(loc_nb, cell_local)
        if (cell_local) then
          call get_local_index(loc_nb, idx(ctr))
          row(ctr) = 1.0
          ctr = ctr + 1
        end if
      end do
      idx = idx - 1 ! F->C

#if PETSC_VERSION_GE(3,23,0)
      call MatSetValues(M, 1, [i - 1], max_nb, idx, row, INSERT_VALUES, ierr)
#else
      call MatSetValues(M, 1, i - 1, max_nb, idx, row, INSERT_VALUES, ierr)
#endif
    end do
    call MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY, ierr)

    deallocate (idx)
    deallocate (row)

    ! Get index sets for reordering
    call MatGetOrdering(M, MATORDERINGRCM, rperm, cperm, ierr)
    call MatDestroy(M, ierr)
    call ISDestroy(cperm, ierr)

    ! Fill local indices in original ordering -> destination, i.e. to(i) => new index of cell i.
    allocate (new_indices(local_num_cells))

#if PETSC_VERSION_GE(3,23,0)
    call ISGetIndices(rperm, row_indices, ierr)
#else
    call ISGetIndicesF90(rperm, row_indices, ierr)
#endif
    if (local_num_cells >= 1) then
      do i = 1, local_num_cells
        idx_new = row_indices(i) + 1 ! C->F
        new_indices(idx_new) = i
      end do
    end if

#if PETSC_VERSION_GE(3,23,0)
    call ISRestoreIndices(rperm, row_indices, ierr)
#else
    call ISRestoreIndicesF90(rperm, row_indices, ierr)
#endif
    call ISDestroy(rperm, ierr)
  end subroutine get_reordering

end submodule
