program test_distorted_kernels
  use testing_lib

  use kinds, only: ccs_real
  use fv_kernels
  use distorted_refs

  implicit none

  real(ccs_real), parameter :: rtol = 1.0e-6_ccs_real
  real(ccs_real), parameter :: atol = 1.0e-8_ccs_real

  type(diffusion_kernel) :: diffusion
  type(upwind_advection_kernel) :: upwind
  type(luds_advection_kernel) :: luds
  type(cd_advection_kernel) :: cd
  real(ccs_real), dimension(2) :: coeffs
  real(ccs_real), dimension(2) :: expected
  real(ccs_real) :: explicit
  integer :: icase
  integer :: ifield
  character(len=32) :: msg

  call init()

  ! The reference stencils isolate skewness, non-orthogonality, uneven
  ! spacing, and the combined distorted case.
  do icase = 1, n_cases
    coeffs = diffusion%eval_coeffs(diff_corr(icase))
    expected = [diff_corr(icase), -diff_corr(icase)]
    write (msg, '(a, i0)') "diff coeff c", icase
    call assert_close(coeffs, expected, rtol, atol, &
                      trim(msg))

    coeffs = cd%eval_coeffs(adv_flux(icase))
    write (msg, '(a, i0)') "cd coeff c", icase
    call assert_close(coeffs, [adv_flux(icase), 0.0_ccs_real], rtol, atol, &
                      trim(msg))

    do ifield = 1, n_fields
      explicit = diffusion%eval_explicit(phi_val(:, ifield, icase), diff_corr(icase), lf(icase), &
                                         rvec(:, :, icase), grad_phi(:, :, ifield, icase))
      write (msg, '(a, i0, a, i0)') "diff exp c", icase, " f", ifield
      call assert_close(explicit, diff_exp(ifield, icase), rtol, atol, &
                        trim(msg))

      explicit = cd%eval_explicit(phi_val(:, ifield, icase), adv_flux(icase), lf(icase), &
                                  rvec(:, :, icase), grad_phi(:, :, ifield, icase))
      write (msg, '(a, i0, a, i0)') "cd exp c", icase, " f", ifield
      call assert_close(explicit, cd_exp(ifield, icase), rtol, atol, &
                        trim(msg))

      explicit = upwind%eval_explicit(phi_val(:, ifield, icase), adv_flux(icase), lf(icase), &
                                      rvec(:, :, icase), grad_phi(:, :, ifield, icase))
      write (msg, '(a, i0, a, i0)') "upwind c", icase, " f", ifield
      call assert_close(explicit, upwind_exp(ifield, icase), rtol, atol, &
                        trim(msg))

      explicit = luds%eval_explicit(phi_val(:, ifield, icase), adv_flux(icase), lf(icase), &
                                    rvec(:, :, icase), grad_phi(:, :, ifield, icase))
      write (msg, '(a, i0, a, i0)') "luds c", icase, " f", ifield
      call assert_close(explicit, luds_exp(ifield, icase), rtol, atol, &
                        trim(msg))
    end do
  end do

  call fin()

end program test_distorted_kernels
