!v Test that cells have correct numbers of neighbours
!
!  for any mesh with >1 cell, every cell must have at least 1 neighbour.
program test_mesh_neighbours

  use testing_lib

  use ccs_base, only: bnd_names_default
  use core
  
  use meshing, only: create_cell_locator, create_neighbour_locator, count_neighbours, &
                     get_boundary_status, get_local_num_cells, get_local_index
  use meshing, only: set_mesh_object, nullify_mesh_object
  use mesh_utils, only: build_mesh, test_mesh_internal_neighbours

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
  integer(ccs_int) :: mctr

  type(ccs_options) :: run_options

  call init()

  ! XXX: use smaller size than 2D test - 20^3 ~= 100^2
  do mctr = 2, size(m)
    n = m(mctr)

    l = parallel_random(par_env)
    run_options%mesh%bnd_names = bnd_names_default
    run_options%mesh%cps = n
    run_options%mesh%domain_size = l
    mesh = build_mesh(par_env, shared_env, run_options)
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
      call assert_lt(nnb, 7, "FAIL: cell should have at most 6 neighbours ")

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
      call stop_test("ERROR: Unknown parallel environment!")
    end select

    expected_boundary_ctr = 6 * n * n ! XXX: specific to 3D Cartesian mesh. For a cube this just counts the surface area in terms of cells.
    ! each cell on the edge excluding cube vertices, and 8 for each cell on a face
    ! excluding cube vertices and edges
    call assert_eq(global_boundary_ctr, expected_boundary_ctr, "FAIL: mesh boundary count is incorrect")

    call nullify_mesh_object()
  end do

  call fin()

end program test_mesh_neighbours
