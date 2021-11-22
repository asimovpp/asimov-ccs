!> @brief Test program for mesh_mod.f90
program test_mesh

  use MPI

  use kinds
  use types
  use parallel
  use parallel_types
  use parallel_types_mpi

  implicit none

  class(parallel_environment), allocatable, target :: par_env
  integer(accs_err) :: ierr

  ! XXX: This would be a good candidate for a testing library
  integer, parameter :: line_length=72
  integer :: ctr
  
  call init()
  call test_runner()
  call fin()
  
contains

  subroutine test_runner()

    call test_mesh_point_distribution()
    call test_square_mesh_cell_centres()
    
  end subroutine test_runner
  
  !> @brief Test the cell centres of a square mesh.
  !
  !> @description The cell centres of a mesh should all fall within the meshed domain, for a square
  !!              mesh \f$x\in[0,1]^d\f$.
  subroutine test_square_mesh_cell_centres()

    use mesh_utils, only : build_square_mesh

    type (mesh) :: square_mesh

    logical :: passing
    logical :: global_passing
    
    real(accs_real) :: l
    integer(accs_int) :: n

    integer(accs_int) :: i
    
    do n = 1, 100
      call random_number(l)
      square_mesh = build_square_mesh(n, l, par_env)

      do i = 1, square_mesh%nlocal
        associate(x => square_mesh%xc(1, i), &
             y => square_mesh%xc(2, i))
          if ((x > l) .or. (x < 0_accs_real) &
               .or. (y > l) .or. (y < 0_accs_real)) then
            print *, "FAIL: expected 0 <= x,y <= ", l, " got ", x, " ", y
            passing = .false.
          end if
        end associate
      end do
      
    end do

    select type(par_env)
    type is(parallel_environment_mpi)
      call MPI_Allreduce(passing, global_passing, 1, MPI_LOGICAL, MPI_LAND, par_env%comm, ierr)
    class default
      print *, "ERROR: Unknown parallel environment!"
      stop
    end select

    ! XXX: This would be a good candidate for a testing library
    if (global_passing) then
      if (par_env%proc_id == par_env%root) then
        if (ctr > line_length) then
          print *, " "
          write(*,"(A)",advance="no") " "
        end if
        write(*,"(A)",advance="no") "."
      end if
      ctr = ctr + 1
    else
      ctr = 0
    end if
    
  end subroutine test_square_mesh_cell_centres

  !> @brief Test the parallel distribution of a mesh's points.
  !
  !> @description A mesh of size N should have its points distributed in parallel such that
  !!              \f$\sum_p n_p = N\f$.
  subroutine test_mesh_point_distribution()

    use mesh_utils, only : build_square_mesh
    
    type(mesh) :: square_mesh

    logical :: passing = .true.
    integer(accs_int) :: n

    integer(accs_int) :: n_expected
    integer(accs_int) :: n_global
    
    do n = 1, 100
      square_mesh = build_square_mesh(n, 1.0_accs_real, par_env)
      n_expected = n**2

      select type(par_env)
      type is(parallel_environment_mpi)
        call MPI_Allreduce(square_mesh%nlocal, n_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
      class default
        print *, "ERROR: Unknown parallel environment!"
        stop
      end select
      
      if (n_global /= n_expected) then
        print *, "FAIL: expected ", n_expected, " got ", n_global, " (test_mesh:test_mesh_point_distribution/1)"
        passing = .false.
      end if

      if (square_mesh%n /= n_expected) then
        print *, "FAIL: expected ", n_expected, " got ", n_global, " (test_mesh:test_mesh_point_distribution/2)"
        passing = .false.
      end if

      if (.not. passing) then
        exit
      end if
    end do

    if (passing) then
      if (par_env%proc_id == par_env%root) then
        if (ctr > line_length) then
          print *, " "
          write(*,"(A)",advance="no") " "
        end if
        write(*,"(A)",advance="no") "."
      end if
      ctr = ctr + 1
    else
      ctr = 0
    end if
    
  end subroutine test_mesh_point_distribution

  !> @brief Test initialisation
  !
  !> @description Performs initialisation for the test (setting up parallel environment, etc.)
  subroutine init()

    call initialise_parallel_environment(par_env)

    ! XXX: This would be a good candidate for a testing library
    if (par_env%proc_id == par_env%root) then
      print *, "Testing: mesh"
      write(*,"(A)",advance="no") " "
    end if

    ctr = 0
    
  end subroutine init

  !> @brief Test finalisation
  !
  !> @description Performs finalisation for the test (tearing down parallel environment, etc.)
  subroutine fin()

    call cleanup_parallel_environment(par_env)

    ! XXX: This would be a good candidate for a testing library
    if (par_env%proc_id == par_env%root) then
      print *, " "
      print *, "----------------------------------"
    end if

  end subroutine fin
  
end program test_mesh
