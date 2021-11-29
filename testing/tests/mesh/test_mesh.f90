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
  integer :: real_type

  logical :: test_suite_passing = .true.

  if (kind(0.0_accs_real) == kind(0.0d0)) then
    real_type = MPI_DOUBLE
  else
    real_type = MPI_FLOAT
  end if
  
  call init()
  call test_runner()
  call fin()
  
contains

  subroutine test_runner()

    call test_mesh_point_distribution()
    call test_square_mesh_cell_centres()
    call test_mesh_square_mesh_volume()
    
  end subroutine test_runner
  
  !> @brief Test the cell centres of a square mesh.
  !
  !> @description The cell centres of a mesh should all fall within the meshed domain, for a square
  !!              mesh \f$x\in[0,1]^d\f$.
  subroutine test_square_mesh_cell_centres()

    use mesh_utils, only : build_square_mesh

    type (mesh) :: square_mesh

    logical :: passing = .true.
    logical :: global_passing
    
    real(accs_real) :: l
    integer(accs_int) :: n

    integer(accs_int) :: i
    
    do n = 1, 100
      l = parallel_random(par_env)
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
      test_suite_passing = .false.
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

      if (square_mesh%nlocal < 0) then
        ! XXX: Zero cells on a PE is not necessarily invalid...
        passing = .false.
        select type(par_env)
        type is(parallel_environment_mpi)
          call MPI_Allreduce(square_mesh%nlocal, n_global, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
        class default
          print *, "ERROR: Unknown parallel environment!"
          stop
        end select
      end if
      
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
      test_suite_passing = .false.
    end if
    
  end subroutine test_mesh_point_distribution

  !> @brief Test the square mesh generator creates a correctly-sized mesh.
  !
  !> @description A "square" domain of side L should result in a mesh of volume L^d, this can be
  !> verified by summing the volumes of all cells.
  subroutine test_mesh_square_mesh_volume()

    use mesh_utils, only : build_square_mesh
    
    type(mesh) :: square_mesh

    logical :: passing = .true.

    integer(accs_int) :: n
    real(accs_real) :: l
    real(accs_real) :: vol
    real(accs_real) :: vol_global
    real(accs_real) :: expected_vol

    integer(accs_int) :: i
    
    do n = 1, 100
      l = parallel_random(par_env)
      square_mesh = build_square_mesh(n, l, par_env)
      expected_vol = l**2 ! XXX: Currently the square mesh is hard-coded 2D...

      vol = 0_accs_real
      do i = 1, square_mesh%nlocal
        vol = vol + square_mesh%vol
      end do
      
      select type(par_env)
      type is(parallel_environment_mpi)
        call MPI_Allreduce(vol, vol_global, 1, real_type, MPI_SUM, par_env%comm, ierr)
      class default
        print *, "ERROR: Unknown parallel environment!"
        stop
      end select

      ! XXX: This would be a good candidate for a testing library
      if (abs(expected_vol - vol_global) > 1.0e-8) then
        print *, "FAIL: expected ", expected_vol, " got ", vol_global
        print *, square_mesh%h, l/n
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
      test_suite_passing = .false.
    end if
    
  end subroutine test_mesh_square_mesh_volume

  !> @brief Test initialisation
  !
  !> @description Performs initialisation for the test (setting up parallel environment, etc.)
  subroutine init()

    integer, allocatable :: seed(:)
    integer :: n

    call initialise_parallel_environment(par_env)

    ! XXX: This would be a good candidate for a testing library
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    if (par_env%proc_id == par_env%root) then
      print *, "Using seed: ", seed
      print *, "----------------------------------"
    end if
    deallocate(seed)

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

    logical :: global_test_suite_passing
    logical :: i_am_root = .false.
    
    select type(par_env)
    type is(parallel_environment_mpi)
      call MPI_Allreduce(test_suite_passing, global_test_suite_passing, 1, MPI_LOGICAL, MPI_LAND, par_env%comm, ierr)
    class default
      print *, "ERROR: Unknown parallel environment!"
      stop
    end select

    if (par_env%proc_id == par_env%root) then
      i_am_root = .true.
    end if
    call cleanup_parallel_environment(par_env)

    ! XXX: This would be a good candidate for a testing library
    if (i_am_root) then
      print *, " "
      print *, "----------------------------------"

      if (.not. global_test_suite_passing) then
        print *, "Test suite test_mesh failed!"
        stop -1
      end if
    end if

  end subroutine fin

  !> @brief Helper function to get a random number in parallel
  !
  !> @description Generates a random number and broadcasts the value on the root of the parallel
  !! environment, ensuring a uniform value is used.
  !
  !> @note Does this belong in the parallel module?
  real(accs_real) function parallel_random(par_env)

    class(parallel_environment), intent(in) :: par_env

    call random_number(parallel_random)

    select type(par_env)
    type is(parallel_environment_mpi)
      call MPI_Bcast(parallel_random, 1, real_type, par_env%root, par_env%comm, ierr)
    class default
      print *, "ERROR: Unknown parallel environment!"
      stop
    end select
    
  end function parallel_random
  
end program test_mesh
