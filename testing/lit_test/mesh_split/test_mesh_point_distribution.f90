!> @brief Test the parallel distribution of a mesh's points.
!
!> @description A mesh of size N should have its points distributed in parallel such that
!!              \f$\sum_p n_p = N\f$.
program test_mesh_point_distribution

  use testing_lib
  use mesh_utils, only : build_square_mesh
  
  type(mesh) :: square_mesh

  integer(accs_int) :: n

  integer(accs_int) :: n_expected
  integer(accs_int) :: n_global
  
  call init()
  
  do n = 1, 100
    square_mesh = build_square_mesh(n, 1.0_accs_real, par_env)

    if (square_mesh%nlocal < 0) then
      ! XXX: Zero cells on a PE is not necessarily invalid...
      ! ? exit
      select type(par_env)
      type is(parallel_environment_mpi)
        call MPI_Allreduce(square_mesh%nlocal, n_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
      class default
        write (message,*) "ERROR: Unknown parallel environment!"
        call stop_test(message)
      end select
    end if
    
    n_expected = n**2

    select type(par_env)
    type is(parallel_environment_mpi)
      call MPI_Allreduce(square_mesh%nlocal, n_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
    class default
      write (message,*) "ERROR: Unknown parallel environment!"
      call stop_test(message)
    end select
    
    if (n_global /= n_expected) then
      write (message,*) "FAIL: expected ", n_expected, " got ", n_global, " (test_mesh:test_mesh_point_distribution/1)"
      call stop_test(message)
    end if

    if (square_mesh%n /= n_expected) then
      write (message,*) "FAIL: expected ", n_expected, " got ", square_mesh%n, " (test_mesh:test_mesh_point_distribution/2)"
      call stop_test(message)
    end if

  end do

  call fin()

end program test_mesh_point_distribution
