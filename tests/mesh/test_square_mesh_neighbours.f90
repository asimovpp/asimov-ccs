!v Test that cells have correct numbers of neighbours when building square mesh
!
!  for any mesh with >1 cell, every cell must have at least 1 neighbour.
program test_square_mesh_neighbours

  use testing_lib

  use meshing, only: create_cell_locator, create_neighbour_locator, count_neighbours, &
                     get_boundary_status, get_local_num_cells, get_local_index
  use meshing, only: set_mesh_object, nullify_mesh_object
  use mesh_utils, only: build_square_mesh

  implicit none

  type(cell_locator) :: loc_p

  integer(ccs_int) :: n
  real(ccs_real) :: l

  integer(ccs_int) :: local_num_cells
  integer(ccs_int) :: i

  integer(ccs_int) :: nnb
  integer(ccs_int) :: j

  integer(ccs_int) :: index_nb
  
  type(neighbour_locator) :: loc_nb
  logical :: is_boundary
  integer(ccs_int) :: boundary_ctr
  integer(ccs_int) :: global_boundary_ctr
  integer(ccs_int) :: expected_boundary_ctr

  integer(ccs_int), dimension(5) :: m = (/4, 8, 12, 16, 20/)
  integer(ccs_int) :: n_v, n_e
  integer(ccs_int) :: mctr

  call init()

  ! XXX: use smaller size than 2D test - 20^3 ~= 100^2
  do mctr = 2, size(m)
    n = m(mctr)

    l = parallel_random(par_env)
    mesh = build_square_mesh(par_env, shared_env, n, l)
    call set_mesh_object(mesh)

    boundary_ctr = 0
    call get_local_num_cells(local_num_cells)
    do i = 1, local_num_cells

      call create_cell_locator(i, loc_p)
      call count_neighbours(loc_p, nnb)

      ! In the case of a cell at the end of a chain of cells it should have 1 interior neighbour
      ! and 1 boundary/external neighbour - c.f. 1D boundary cell.
      ! Even in the limit of single 1D cell should have 2 boundary neighbours.
      call assert_gt(nnb, 1, "FAIL: cell should have 2 or more neighbours ")
      call assert_lt(nnb, 5, "FAIL: cell should have at most 4 neighbours ")

      ! Loop over neighbours
      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)

        ! Check for zero neighbour index. This indicates a cell was not linked as a neighbour. For
        ! the build mesh case we should always be able to fill our neighbours.
        call get_local_index(loc_nb, index_nb)
        call assert_neq(index_nb, 0, "All neighbours should be filled!")
        
        call get_boundary_status(loc_nb, is_boundary)
        if (is_boundary) then
          ! Boundary neighbour/face
          boundary_ctr = boundary_ctr + 1
        else
          call test_mesh_internal_neighbours(loc_nb)
        end if
      end do

    end do

    ! Check total boundary neighbours
    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(boundary_ctr, global_boundary_ctr, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
    class default
      write (message, *) "ERROR: Unknown parallel environment!"
      call stop_test(message)
    end select

    expected_boundary_ctr = 4 * n ! XXX: specific to 2D Cartesian mesh. This just counts the boundary neighbours on the perimeter of the square mesh.
    n_v = 4
    n_e = 4 * n - 8
    ! each cell on the edge excluding square vertices
    call assert_eq(global_boundary_ctr, expected_boundary_ctr, "FAIL: mesh boundary count is incorrect")

    call nullify_mesh_object()
  end do

  call fin()

contains

  subroutine test_mesh_internal_neighbours(loc_nb)

    use meshing, only: count_neighbours, get_local_index, get_boundary_status, get_local_status

    type(neighbour_locator), intent(in) :: loc_nb

    integer(ccs_int) :: index_nb
    type(cell_locator) :: cell_loc_nb
    integer(ccs_int) :: nnb
    logical :: found_parent
    integer(ccs_int) :: j
    type(neighbour_locator) :: loc_nb_nb
    logical :: is_boundary
    logical :: is_local

    associate (parent_idx => loc_nb%index_p)
      call get_local_index(loc_nb, index_nb)

      ! Neighbour index should not be its parents
      if (index_nb == parent_idx) then
        write (message, *) "FAIL: Neighbour has same index ", index_nb, " as parent cell ", parent_idx
        call stop_test(message)
      end if

      call get_local_status(loc_nb, is_local)
      if (is_local) then
        ! Parent should be in neighbour's neighbour list
        call create_cell_locator(index_nb, cell_loc_nb)
        call count_neighbours(cell_loc_nb, nnb)
        found_parent = .false.
        do j = 1, nnb
          call create_neighbour_locator(cell_loc_nb, j, loc_nb_nb)
          call get_boundary_status(loc_nb_nb, is_boundary)
          if (.not. is_boundary) then ! We are looking for parent cell - by definition not a boundary!
            call get_local_index(loc_nb_nb, index_nb)
            if (index_nb == parent_idx) then
              found_parent = .true.
              exit
            end if
          end if
        end do
        if (.not. found_parent) then
          write (message, *) "FAIL: Couldn't find cell in neighbour's neighbour list!"
          call stop_test(message)
        end if
      end if

    end associate

  end subroutine test_mesh_internal_neighbours

end program test_square_mesh_neighbours
