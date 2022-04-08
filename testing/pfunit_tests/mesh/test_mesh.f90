!> @brief Test program for mesh_mod.f90
program test_mesh

  use MPI

  use constants, only : ndim
  
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

  integer(accs_int), parameter :: nmax = 100
  real(accs_real), parameter :: eps = epsilon(0.0_accs_real)

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

    call test_square_mesh_indices()
    call test_mesh_point_distribution()
    call test_square_mesh_centres()
    call test_mesh_square_mesh_volume()
    call test_mesh_square_mesh_closed()
    call test_mesh_neighbours()
    
  end subroutine test_runner

  !> @brief Test the indexing of cells
  subroutine test_square_mesh_indices()

    use mesh_utils, only : build_square_mesh, global_index

    type (mesh) :: square_mesh

    logical :: passing = .true.
    logical :: global_passing
    
    real(accs_real) :: l
    integer(accs_int) :: n

    integer(accs_int) :: i

    type(cell_locator) :: cell_location
    integer(accs_int) :: idxg

    do n = 1, nmax
      l = parallel_random(par_env)
      square_mesh = build_square_mesh(n, l, par_env)

      associate(nlocal => square_mesh%nlocal, &
           nglobal => square_mesh%n)
        do i = 1, nlocal
          call set_cell_location(square_mesh, i, cell_location)
          call global_index(cell_location, idxg)
          if ((idxg < 1) .or. (idxg > nglobal)) then
            if (idxg /= -1) then
              print *, "FAIL: expected global index 1 <= idx <= ", nglobal, " got ", idxg
              passing = .false.
            end if
            exit
          end if
        end do
      end associate
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
    
  end subroutine test_square_mesh_indices
  
  !> @brief Test the cell/face centres of a square mesh.
  !
  !> @description The cell/face centres of a mesh should all fall within the meshed domain, for a
  !!              square mesh \f$x\in[0,1]^d\f$.
  subroutine test_square_mesh_centres()

    use mesh_utils, only : build_square_mesh, centre

    type (mesh) :: square_mesh

    logical :: passing = .true.
    logical :: global_passing
    
    real(accs_real) :: l
    integer(accs_int) :: n

    integer(accs_int) :: i
    integer(accs_int) :: j

    type(cell_locator) :: cell_location
    real(accs_real), dimension(ndim) :: cc
    type(face_locator) :: face_location
    real(accs_real), dimension(ndim) :: fc
    
    do n = 1, nmax
      l = parallel_random(par_env)
      square_mesh = build_square_mesh(n, l, par_env)

      do i = 1, square_mesh%nlocal
        call set_cell_location(square_mesh, i, cell_location)
        call centre(cell_location, cc)
        associate(x => cc(1), y => cc(2))
          if ((x > l) .or. (x < 0_accs_real) &
               .or. (y > l) .or. (y < 0_accs_real)) then
            print *, "FAIL: expected cell centre 0 <= x,y <= ", l, " got ", x, " ", y
            passing = .false.
          end if
        end associate

        associate(nnb => square_mesh%nnb(i))
          do j = 1, nnb
            call set_face_location(square_mesh, i, j, face_location)
            call centre(face_location, fc)
            associate(x => fc(1), y => fc(2))
              if ((x > (l + eps)) .or. (x < (0.0_accs_real - eps)) &
                   .or. (y > (l + eps)) .or. (y < (0.0_accs_real - eps))) then
                print *, "FAIL: expected face centre 0 <= x,y <= ", l, " got ", x, " ", y
                passing = .false.
              end if
            end associate
          end do
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
    
  end subroutine test_square_mesh_centres

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
    
    do n = 1, nmax
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
  !>              verified by summing the volumes of all cells.
  subroutine test_mesh_square_mesh_volume()

    use mesh_utils, only : build_square_mesh, volume
    
    type(mesh) :: square_mesh

    logical :: passing = .true.

    integer(accs_int) :: n
    real(accs_real) :: l
    real(accs_real) :: vol
    real(accs_real) :: vol_global
    real(accs_real) :: expected_vol

    integer(accs_int) :: i
    type(cell_locator) :: cell_location
    real(accs_real) :: V
    
    do n = 1, nmax
      l = parallel_random(par_env)
      square_mesh = build_square_mesh(n, l, par_env)
      expected_vol = l**2 ! XXX: Currently the square mesh is hard-coded 2D...

      vol = 0.0_accs_real
      do i = 1, square_mesh%nlocal
        call set_cell_location(square_mesh, i, cell_location)
        call volume(cell_location, V)
        vol = vol + V
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

  !> @brief Test the square mesh generator creates a closed mesh.
  !
  !> @description A valid mesh should be "closed" that is the surface integral should be zero. The
  !!              same is true of mesh cells.
  subroutine test_mesh_square_mesh_closed()

    use mesh_utils, only : build_square_mesh, face_normal, face_area
    
    type(mesh), target :: square_mesh
    type(face_locator) :: face_location
    
    logical :: passing = .true.
    logical :: global_passing
    
    integer(accs_int) :: n
    real(accs_real) :: l

    real(accs_real), dimension(ndim) :: S
    real(accs_real), dimension(ndim) :: norm
    real(accs_real) :: A

    integer(accs_int) :: i, j
    integer :: ndim

    do n = 1, nmax
      l = parallel_random(par_env)
      square_mesh = build_square_mesh(n, l, par_env)
      ndim = size(square_mesh%nf, 1)
      
      ! Loop over cells
      do i = 1, square_mesh%nlocal
        S(:) = 0.0_accs_real

        ! Loop over neighbours/faces
        do j = 1, square_mesh%nnb(i)
          
          call set_face_location(square_mesh, i, j, face_location)
          call face_area(face_location, A)
          call face_normal(face_location, norm)
          S(:) = S(:) + norm(:) * A
        end do

        ! Loop over axes
        do j = 1, ndim
          if (abs(S(j) - 0.0_accs_real) > epsilon(0.0_accs_real)) then
            print *, "FAIL: expected", 0.0_accs_real, " got ", S(j)
            passing = .false.
            exit
          end if
        end do
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
    
  end subroutine test_mesh_square_mesh_closed

  !> @brief Test that cells have correct numbers of neighbours
  !
  !> @description for any mesh with >1 cell, every cell must have at least 1 neighbour.
  subroutine test_mesh_neighbours()

    use mesh_utils, only : build_square_mesh, count_neighbours, boundary_status
    
    type(mesh), target :: square_mesh
    type(cell_locator) :: cell_location
    
    logical :: passing = .true.
    logical :: global_passing
    
    integer(accs_int) :: n
    real(accs_real) :: l

    integer(accs_int) :: i

    integer(accs_int) :: nnb
    integer(accs_int) :: j

    type(neighbour_locator) :: nb_location
    logical :: is_boundary
    integer(accs_int) :: boundary_ctr
    integer(accs_int) :: global_boundary_ctr
    integer(accs_int) :: expected_boundary_ctr
    
    do n = 1, nmax
      l = parallel_random(par_env)
      square_mesh = build_square_mesh(n, l, par_env)

      boundary_ctr = 0
      do i = 1, square_mesh%nlocal

        call set_cell_location(square_mesh, i, cell_location)
        call count_neighbours(cell_location, nnb)

        if (nnb < 2) then
          ! In the case of a cell at the end of a chain of cells it should have 1 interior neighbour
          ! and 1 boundary/external neighbour - c.f. 1D boundary cell.
          ! Even in the limit of single 1D cell should have 2 boundary neighbours.
          print *, "FAIL: cell should have 2 or more neighbours, got ", nnb
          passing = .false.
        else if (nnb > 4) then
          ! XXX: specific to 2D Cartesian mesh
          print *, "FAIL: cell should have at most ", 4, " neighbours, got ", nnb
          passing = .false.
        end if

        ! Loop over neighbours
        do j = 1, nnb
          call set_neighbour_location(cell_location, j, nb_location)
          call boundary_status(nb_location, is_boundary)
          if (is_boundary) then
            ! Boundary neighbour/face
            boundary_ctr = boundary_ctr + 1
          else
            if (.not. test_mesh_internal_neighbours(nb_location)) then
              passing = .false.
            end if
          end if
        end do
        
      end do

      ! Check total boundary neighbours
      select type(par_env)
      type is (parallel_environment_mpi)
        call MPI_Allreduce(boundary_ctr, global_boundary_ctr, 1, MPI_INT, MPI_SUM, par_env%comm, ierr)
      class default
        print *, "ERROR: Unknown parallel environment!"
        stop
      end select

      expected_boundary_ctr = 4 * n ! XXX: specific to 2D Cartesian mesh (square mesh has 2^d sides
                                    !      of length n)
      if (global_boundary_ctr /= expected_boundary_ctr) then
        passing = .false.
        if (par_env%proc_id == par_env%root) then
          print *, "FAIL: mesh boundary count is incorrect, expected ", expected_boundary_ctr, &
               " got ", global_boundary_ctr
        end if
      end if
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
    
  end subroutine test_mesh_neighbours

  logical function test_mesh_internal_neighbours(nb_location)

    use mesh_utils, only : count_neighbours, local_index, boundary_status
    
    type(neighbour_locator), intent(in) :: nb_location

    logical :: passing
    
    integer(accs_int) :: nbidx
    ! type(cell_locator) :: nb_cell_location
    ! integer(accs_int) :: nnb
    ! logical :: found_parent
    ! integer(accs_int) :: j
    ! type(neighbour_locator) :: nbnb_location
    ! logical :: is_boundary
    
    passing = .true.

    associate(mesh => nb_location%mesh, &
         i => nb_location%cell_idx)
      call local_index(nb_location, nbidx)
            
      ! Neighbour index should not be its parents
      if (nbidx == i) then
        print *, "FAIL: Neighbour has same index ", nbidx, " as parent cell ", i
        passing = .false.
      end if

      ! if (passing) then ! Doesn't make sense to continue if part1 fails
      !   ! Parent should be in neighbour's neighbour list
      !   call set_cell_location(mesh, nbidx, nb_cell_location)
      !   ! call count_neighbours(nb_cell_location, nnb)
      !   ! found_parent = .false.
      !   ! do j = 1, nnb
      !   !   call set_neighbour_location(nb_cell_location, j, nbnb_location)
      !   !   call boundary_status(nbnb_location, is_boundary)
      !   !   if (.not. is_boundary) then ! We are looking for parent cell - by definition not a boundary!
      !   !     call local_index(nbnb_location, nbidx)
      !   !     if (nbidx == i) then
      !   !       found_parent = .true.
      !   !       exit
      !   !     end if
      !   !   end if
      !   ! end do
      !   ! if (.not. found_parent) then
      !   !   print *, "FAIL: Couldn't find cell in neighbour's neighbour list!"
      !   !   passing = .false.
      !   ! end if
      ! end if
      
    end associate

    test_mesh_internal_neighbours = passing
    
  end function test_mesh_internal_neighbours
  
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
