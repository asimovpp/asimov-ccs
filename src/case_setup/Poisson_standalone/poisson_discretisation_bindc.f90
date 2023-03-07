
program poisson_discretisation_bindc

  use iso_c_binding
  
  implicit none

  interface
    subroutine discretise_poisson_kernel(nrows, nnz_pr, h, &
         mesh_neighbours, mesh_face_areas, csr_values) bind(c)
      use iso_c_binding

      integer(c_int) :: nrows, nnz_pr
      real(c_double) :: h
      integer(c_int), dimension(nrows * nnz_pr) :: mesh_neighbours
      real(c_double), dimension(nrows * nnz_pr) :: mesh_face_areas
      real(c_double), dimension(nrows * (nnz_pr + 1)) :: csr_values
    end subroutine discretise_poisson_kernel
  end interface

  integer, parameter :: ccs_int = kind(1)
  integer, parameter :: ccs_real = kind(0.0d0)

  integer(ccs_int) :: ncells = 100
  real(ccs_real) :: h = 0.1_ccs_real
  
  call discretise_poisson()
  
contains

  subroutine discretise_poisson()

    integer(ccs_int), dimension(:), allocatable :: flat_neighbours
    real(ccs_real), dimension(:), allocatable :: flat_areas
    real(ccs_real), dimension(:), allocatable :: csr_values

    integer :: nout, step, i, j, csr_idx
    
    ! Fake a flat mesh data structure
    call flatten_mesh(ncells, flat_neighbours, flat_areas)
    
    allocate(csr_values((2 + 1) * ncells)) ! 1 entry per neighbour + diagonal coefficient
    
    call discretise_poisson_kernel(ncells, 2, h, &
         flat_neighbours, flat_areas, csr_values)

    nout = 10
    step = ncells / nout
    do i = 1, ncells, step
      do j = 1, (2 + 1)
        csr_idx = (2 + 1) * i + (j - 1)
        print *, i, j, csr_values(csr_idx)
      end do
    end do
    
    deallocate(csr_values)
    deallocate(flat_neighbours)
    deallocate(flat_areas)
    
  end subroutine discretise_poisson
  
  ! Fake a flat mesh data structure.
  subroutine flatten_mesh(ncells, flat_neighbours, flat_areas)

    integer(ccs_int), intent(in) :: ncells
    
    integer(ccs_int), dimension(:), allocatable, intent(out) :: flat_neighbours
    real(ccs_real), dimension(:), allocatable, intent(out) :: flat_areas

    integer(ccs_int) :: index_p, index_nb
    integer(ccs_int) :: j, ctr
    real(ccs_real) :: A
    
    allocate(flat_neighbours(2 * ncells))
    allocate(flat_areas(2 * ncells))

    A = 1.0_ccs_real
    
    ctr = 1
    do index_p = 1, ncells
      
      do j = 1, 2 ! XXX: Don't leave hardcoded!

        if ((index_p == 1) .and. (j == 1)) then
          index_nb = -1
        else if ((index_p == ncells) .and. (j == 2)) then
          index_nb = -2
        else
          index_nb = index_p + (-1**j)
        end if
        flat_neighbours(ctr) = index_nb
        
        flat_areas(ctr) = A

        ctr = ctr + 1
      end do
    end do
    
  end subroutine flatten_mesh
  
end program poisson_discretisation_bindc
