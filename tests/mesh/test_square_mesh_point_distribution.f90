!v Test the parallel distribution of a mesh's points.
!
!  A mesh of size N should have its points distributed in parallel such that
!  \f$\sum_p n_p = N\f$.
program test_square_mesh_point_distribution

  use testing_lib
  use mesh_utils, only: build_square_mesh
  use meshing, only: get_local_num_cells, get_global_num_cells

  implicit none

  type(ccs_mesh) :: mesh

  integer(ccs_int) :: n
  integer(ccs_int) :: nlocal
  integer(ccs_int) :: n_expected
  integer(ccs_int) :: n_global

  integer(ccs_int), dimension(7) :: m = (/4, 8, 16, 20, 40, 80, 100/)
  integer(ccs_int) :: mctr

  integer(ccs_int) :: global_num_cells

  call init()

  do mctr = 1, size(m)
    n = m(mctr)
    mesh = build_square_mesh(par_env, shared_env, n, 1.0_ccs_real)

    call get_local_num_cells(mesh, nlocal)
    if (nlocal < 0) then
      ! XXX: Zero cells on a PE is not necessarily invalid...
      ! ? exit
      ! select type(par_env)
      ! type is(parallel_environment_mpi)
      !   call MPI_Allreduce(nlocal, n_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
      ! class default
      !   write (message,*) "ERROR: Unknown parallel environment!"
      !   call stop_test(message)
      ! end select
    end if

    n_expected = n**2

    if (nlocal > n_expected) then
      write (message, *) "FAIL: Local number of cells ", nlocal, &
        " exceeds requested total!", n
      call stop_test(message)
    end if

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(nlocal, n_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
    class default
      write (message, *) "ERROR: Unknown parallel environment!"
      call stop_test(message)
    end select

    if (n_global /= n_expected) then
      write (message, *) "FAIL: expected ", n_expected, " got ", n_global, &
        " (test_mesh:test_mesh_point_distribution/1)"
      call stop_test(message)
    end if

    call get_global_num_cells(mesh, global_num_cells)
    call assert_eq(n_expected, global_num_cells, &
                   "(test_mesh:test_mesh_point_distribution/2)")

  end do

  call fin()

end program test_square_mesh_point_distribution
