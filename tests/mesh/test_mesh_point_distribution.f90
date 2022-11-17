!v Test the parallel distribution of a mesh's points.
!
!  A mesh of size N should have its points distributed in parallel such that
!  \f$\sum_p n_p = N\f$.
program test_mesh_point_distribution

  use testing_lib
  use mesh_utils, only: build_mesh

  implicit none

  type(ccs_mesh) :: mesh

  integer(ccs_int) :: nx, ny, nz

  integer(ccs_int) :: n_expected
  integer(ccs_int) :: n_global

  call init()

  nx = 4
  ny = 4
  nz = 4

  mesh = build_mesh(par_env, nx, ny, nz, 1.0_ccs_real)
  associate (nlocal => mesh%topo%local_num_cells)
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

    n_expected = nx * ny * nz

    if (nlocal > n_expected) then
      write (message, *) "FAIL: Local number of cells ", nlocal, &
        " exceeds requested total!"
      call stop_test(message)
    end if

    select type (par_env)
    type is (parallel_environment_mpi)
      call MPI_Allreduce(nlocal, n_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
    class default
      write (message, *) "ERROR: Unknown parallel environment!"
      call stop_test(message)
    end select
  end associate

  if (n_global /= n_expected) then
    write (message, *) "FAIL: expected ", n_expected, " got ", n_global, &
      " (test_mesh:test_mesh_point_distribution/1)"
    call stop_test(message)
  end if

  call assert_equal(n_expected, mesh%topo%global_num_cells, &
                    '("FAIL: expected ", i0, " got ", i0, " (test_mesh:test_mesh_point_distribution/2)")')

  call fin()

end program test_mesh_point_distribution
