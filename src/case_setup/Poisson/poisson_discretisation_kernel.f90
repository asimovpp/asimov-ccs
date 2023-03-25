! @build kernel

module foo

  use iso_c_binding

  integer, parameter :: ccs_int = c_int
  integer, parameter :: ccs_real = c_double

contains

  subroutine discretise_poisson_kernel(nrows, nnz_pr, h, mesh_neighbours, mesh_face_areas, csr_values) bind(c)

    use iso_c_binding

    integer(ccs_int), intent(in) :: nrows
    integer(ccs_int), intent(in) :: nnz_pr

    real(ccs_real), intent(in) :: h
    integer(ccs_int), dimension(nrows * nnz_pr), intent(in) :: mesh_neighbours
    real(ccs_real), dimension(nrows * nnz_pr), intent(in) :: mesh_face_areas

    real(ccs_real), dimension(nrows * (nnz_pr + 1)), intent(out) :: csr_values ! nnz_pr + 1 non-zeros per row (diagonal coefficient!)

    integer(ccs_int) :: i, j, face_idx, csr_idx
    integer(ccs_int) :: index_nb

    real(ccs_real) :: coeff_p, coeff_f, coeff_nb
    real(ccs_real) :: A
 
    print *, " +++ Entering kernel +++"
    print *, " - NROWS = ", nrows
    print *, " - NNZPR = ", nnz_pr
    print *, " - h     = ", h

    !$acc data copyin(mesh_neighbours(:), mesh_face_areas(:)) copyout(csr_values(:))
    !$acc parallel loop
    do i = 1, nrows
      coeff_p = 0.0_ccs_real

      !$acc loop seq
      do j = 1, nnz_pr
        face_idx = nnz_pr * (i - 1) + j
        csr_idx = (nnz_pr + 1) * (i - 1) + j
        index_nb = mesh_neighbours(face_idx)

        if (index_nb > 0) then
          ! Interior face
          A = mesh_face_areas(face_idx)
          coeff_f = (1.0_ccs_real / h) * A

          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f

          csr_values(csr_idx) = coeff_nb
        else
          csr_values(csr_idx) = 0.0_ccs_real
        end if
      end do

      csr_idx = ((nnz_pr + 1) * (i - 1) + nnz_pr) + 1
      csr_values(csr_idx) = coeff_p
    end do
    !$acc end data

    print *, " +++ Leaving kernel +++"   
    
  end subroutine discretise_poisson_kernel

  subroutine apply_dirichlet_bcs_kernel(nrows, nnz_pr, h, mesh_neighbours, mesh_face_areas, mesh_face_yloc, diag_values, rhs_values) bind(c)

    use iso_c_binding

    integer(ccs_int), intent(in) :: nrows
    integer(ccs_int), intent(in) :: nnz_pr

    real(ccs_real), intent(in) :: h
    integer(ccs_int), dimension(nrows * nnz_pr), intent(in) :: mesh_neighbours
    real(ccs_real), dimension(nrows * nnz_pr), intent(in) :: mesh_face_areas
    real(ccs_real), dimension(nrows * nnz_pr), intent(in) :: mesh_face_yloc ! y coordinate of each face, needed for eval_solution

    real(ccs_real), dimension(nrows), intent(out) :: diag_values
    real(ccs_real), dimension(nrows), intent(out) :: rhs_values

    integer(ccs_int) :: i, j, face_idx
    integer(ccs_int) :: index_nb

    real(ccs_real) ::  coeff, r, boundary_coeff, boundary_val
    real(ccs_real) :: A



    do i = 1, nrows
        !if (minval(mesh%topo%nb_indices(:, i)) < 0) then
      if (minval(mesh_neighbours(nnz_pr*(i-1)+1: nnz_pr*(i-1)+nnz_pr)) < 0) then  ! todo: maybe to remove?

        coeff = 0.0_ccs_real
        r = 0.0_ccs_real

        do j = 1, nnz_pr

          face_idx = nnz_pr * (i - 1) + j
          index_nb = mesh_neighbours(face_idx)

          if (index_nb .le. 0) then
            A = mesh_face_areas(face_idx)
            boundary_coeff = (2.0 / h) * A
            boundary_val = mesh_face_yloc(face_idx)  ! the solution is the y coordinate, see poisson.f90:eval_solution

            ! Coefficient
            coeff = coeff - boundary_coeff

            ! RHS vector
            r = r - boundary_val * boundary_coeff
          end if

        end do

        diag_values(i) = coeff
        rhs_values(i) = r

      end if
    end do


  end subroutine apply_dirichlet_bcs_kernel
end module foo
