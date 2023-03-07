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
    
    do i = 1, nrows
      coeff_p = 0.0_ccs_real

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
    
  end subroutine discretise_poisson_kernel

end module foo
