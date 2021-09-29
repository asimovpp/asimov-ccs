!> @brief Case-specific boundary conditions for the Poisson case.
!
!> @details

subroutine apply_dirichlet_bcs(M, b)

  use constants, only : add_mode
  use mat, only : set_eqn
  use types, only : vector_values, matrix_values, matrix, vector
  use utils, only : set_values
  use kinds, only: accs_int, accs_real

  class(matrix), intent(inout) :: M
  class(vector), intent(inout) :: b

  integer(accs_int) :: i, j
  real(accs_real) :: boundary_coeff, boundary_val

  type(vector_values) :: vec_values
  type(matrix_values) :: mat_coeffs

  allocate(mat_coeffs%rglob(1), mat_coeffs%cglob(1), mat_coeffs%val(1))
  allocate(vec_values%idx(1), vec_values%val(1))
  mat_coeffs%mode = add_mode
  vec_values%mode = add_mode

  associate(row=>mat_coeffs%rglob, &
       col=>mat_coeffs%cglob, &
       coeff=>mat_coeffs%val, &
       idx=>vec_values%idx, &
       val=>vec_values%val, &
       idx_global=>square_mesh%idx_global)

    do i = 1, square_mesh%nlocal
       if (minval(square_mesh%nbidx(:, i)) < 0) then
          coeff(1) = 0.0_accs_real
          val(1) = 0.0_accs_real
          do j = 1, square_mesh%nnb(i)
             associate(nbidx=>square_mesh%nbidx(j, i))

               if (nbidx < 0) then
                  associate(h=>square_mesh%h, &
                       Af=>square_mesh%Af)

                    !! Boundary face
                    boundary_coeff = (2.0 / h) * Af

                    if ((nbidx == -1) .or. (nbidx == -2)) then
                       !! Left or right boundary
                       boundary_val = rhs_val(idx_global(i))
                    else if (nbidx == -3) then
                       !! Bottom boundary
                       boundary_val = rhs_val(0, -0.5_accs_real)
                    else if (nbidx == -4) then
                       !! Top boundary
                       boundary_val = rhs_val(idx_global(i), 0.5_accs_real)
                    else
                       print *, "ERROR: invalid/unsupported BC ", nbidx
                       stop
                    end if

                    ! Coefficient
                    row(1) = idx_global(i) - 1
                    col(1) = idx_global(i) - 1
                    coeff(1) = coeff(1) - boundary_coeff

                    ! RHS vector
                    idx(1) = idx_global(i) - 1
                    val(1) = val(1) - boundary_val * boundary_coeff

                  end associate
               end if

             end associate
          end do
          call set_values(mat_coeffs, M)
          call set_values(vec_values, b)
       end if
    end do

  end associate

  deallocate(vec_values%idx)
  deallocate(vec_values%val)

end subroutine apply_dirichlet_bcs
