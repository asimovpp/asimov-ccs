!> @brief Test the RCM reordering of the mesh.
program test_rcm

  use mpi
  
  use testing_lib

  use constants, only: ndim
  use meshing, only: set_cell_location, set_face_location, set_vert_location, get_centre, &
       get_local_num_cells, get_global_index, set_neighbour_location, get_boundary_status
  use mesh_utils, only: build_square_topology, partition_stride
  use partitioning, only: compute_partitioner_input, compute_connectivity
  use reordering, only: reorder_cells

  implicit none

  type(ccs_mesh), allocatable :: mesh

  real(ccs_real) :: l
  integer(ccs_int) :: n

  integer(ccs_int), dimension(5) :: m = (/ 2, 4, 8, 16, 20 /)
  integer(ccs_int) :: mctr
 
    integer(ccs_int) :: local_num_cells
  integer(ccs_int) :: nnb
  integer(ccs_int), dimension(:, :), allocatable :: global_neighbours_ref
  
  call init()

  ! XXX: use smaller size than 2D test - 20^3 ~= 100^2
  do mctr = 1, size(m)
    
    allocate(mesh)

    n = m(mctr)

    l = parallel_random(par_env)

    ! Build topology
    call build_square_topology(par_env, n, mesh)
    call compute_partitioner_input(par_env, mesh)
    call partition_stride(par_env, mesh)
    call compute_connectivity(par_env, mesh)

    call get_local_num_cells(mesh, local_num_cells)

    call get_reference_connectivity(global_neighbours_ref)

    ! Reorder and compare connectivity
    call mock_geo(mesh)
    call reorder_cells(par_env, mesh)

    call test_global_connectivity()
    
    deallocate(mesh)

  end do

  call fin()

contains

  subroutine test_global_connectivity()

    integer(ccs_int) :: i, j

    integer(ccs_int) :: natural_index_p
    integer(ccs_int) :: natural_index_nb
    integer(ccs_int) :: index_nb
    
    do i = 1, local_num_cells
       natural_index_p = mesh%topo%natural_indices(i)

       do j = 1, nnb
          index_nb = mesh%topo%nb_indices(j, i)
          print *, i, j, index_nb
          if (index_nb > 0) then
             natural_index_nb = mesh%topo%natural_indices(index_nb)

             if (natural_index_nb /= global_neighbours_ref(j, natural_index_p)) then
                write(message, *) "FAIL: neighbour natural index doesn't match reference at ", j, &
                     " ", natural_index_p
                call stop_test(message)
             end if
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
    allocate(global_neighbours_ref(nnb, mesh%topo%global_num_cells))

    global_neighbours_ref(:, :) = 0

    do i = 1, local_num_cells
       call set_cell_location(mesh, i, loc_p)
       call get_global_index(loc_p, global_index_p)
       
       do j = 1, nnb
          call set_neighbour_location(loc_p, j, loc_nb)
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

    associate(total_num_cells => mesh%topo%total_num_cells)
      if (allocated(mesh%geo%volumes)) then
        deallocate(mesh%geo%volumes)
      end if
      allocate (mesh%geo%volumes(total_num_cells))

      if (allocated(mesh%geo%x_p)) then
        deallocate(mesh%geo%x_p)
      end if
      allocate (mesh%geo%x_p(ndim, total_num_cells))
    end associate

  end subroutine

end program test_rcm
