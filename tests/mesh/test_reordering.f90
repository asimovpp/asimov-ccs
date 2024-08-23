!> @brief Test the RCM reordering of the mesh.
program test_rcm

  use mpi

  use testing_lib

  use constants, only: ndim
  use meshing, only: create_cell_locator, create_face_locator, create_vert_locator, get_centre, &
                     get_local_num_cells, get_global_index, create_neighbour_locator, get_boundary_status, &
                     get_natural_index
  use meshing, only: set_mesh_object, nullify_mesh_object
  use mesh_utils, only: build_square_topology, partition_stride
  use partitioning, only: compute_partitioner_input, compute_connectivity
  use reordering, only: reorder_cells

  implicit none

  real(ccs_real) :: l
  integer(ccs_int) :: n

  integer(ccs_int), dimension(5) :: m = (/2, 4, 8, 16, 20/)
  integer(ccs_int) :: mctr

  integer(ccs_int) :: local_num_cells
  integer(ccs_int) :: nnb
  integer(ccs_int), dimension(:, :), allocatable :: global_neighbours_ref

  type(ccs_mesh), allocatable :: tmp_mesh

  call init()

  ! XXX: use smaller size than 2D test - 20^3 ~= 100^2
  do mctr = 1, size(m)
    allocate(tmp_mesh)
    mesh = tmp_mesh
    n = m(mctr)

    l = parallel_random(par_env)

    ! Build topology
    call build_square_topology(par_env, shared_env, n, mesh)
    call compute_partitioner_input(par_env, shared_env, mesh)

    call set_mesh_object(mesh)
    call partition_stride(par_env, shared_env, roots_env, mesh)
    call compute_connectivity(par_env, shared_env, mesh)

    call get_local_num_cells(local_num_cells)

    call get_reference_connectivity(global_neighbours_ref)

    ! Reorder and compare connectivity
    call mock_geo(mesh)
    call reorder_cells(par_env, shared_env, mesh)

    call test_global_connectivity()
    call nullify_mesh_object()
    deallocate(tmp_mesh)

  end do

  call fin()

contains

  subroutine test_global_connectivity()

    integer(ccs_int) :: i, j

    integer(ccs_int) :: natural_index_p
    integer(ccs_int) :: natural_index_nb
    integer(ccs_int) :: index_nb

    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    logical :: is_boundary

    character(len=:), allocatable :: msg

    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_natural_index(loc_p, natural_index_p)

      do j = 1, nnb
        index_nb = mesh%topo%nb_indices(j, i)
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        if (.not. is_boundary) then
          call get_natural_index(loc_nb, natural_index_nb)

          ! Prepare message in case of failure
          msg = "FAIL: neighbour natural index doesn't match reference at " // str(j) //  &
                " " // str(natural_index_p)
          call assert_eq(natural_index_nb, global_neighbours_ref(j, natural_index_p), msg)
        end if
      end do
    end do

  end subroutine test_global_connectivity

  subroutine get_reference_connectivity(global_neighbours_ref)

    ! Declaration ensures deallocation on entry
    integer(ccs_int), dimension(:, :), allocatable, intent(out) :: global_neighbours_ref

    integer(ccs_int) :: i, j

    integer(ccs_int) :: global_index_p
    integer(ccs_int) :: global_index_nb
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    logical :: is_boundary

    nnb = 4 ! Constant for square mesh
    allocate (global_neighbours_ref(nnb, mesh%topo%global_num_cells))

    global_neighbours_ref(:, :) = 0

    do i = 1, local_num_cells
      call create_cell_locator(i, loc_p)
      call get_global_index(loc_p, global_index_p)

      do j = 1, nnb
        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_boundary_status(loc_nb, is_boundary)
        if (.not. is_boundary) then
          call get_global_index(loc_nb, global_index_nb)

          global_neighbours_ref(j, global_index_p) = global_index_nb
        end if
      end do
    end do

    call MPI_Allreduce(MPI_IN_PLACE, global_neighbours_ref, mesh%topo%global_num_cells, &
                       MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

  end subroutine get_reference_connectivity

  subroutine mock_geo(mesh)

    type(ccs_mesh), intent(inout) :: mesh

    associate (total_num_cells => mesh%topo%total_num_cells)
      if (associated(mesh%geo%volumes)) then
        call destroy_shared_array(shared_env, mesh%geo%volumes, mesh%geo%volumes_window)
      end if
      call create_shared_array(shared_env, total_num_cells, mesh%geo%volumes, mesh%geo%volumes_window)

      if (associated(mesh%geo%x_p)) then
        call destroy_shared_array(shared_env, mesh%geo%x_p, mesh%geo%x_p_window)
      end if
      call create_shared_array(shared_env, [ndim, total_num_cells], mesh%geo%x_p, mesh%geo%x_p_window)
    end associate

  end subroutine

end program test_rcm
