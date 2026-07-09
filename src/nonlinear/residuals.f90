!v Module for handling residuals
!
! Computes them, compute their norms, and prints them
! @dont_fail_linter

module residuals
#include "ccs_macros.inc"

  use constants, only: L2, Linfty
  use kinds, only: ccs_int, ccs_real
  use types, only: fluid, field, ccs_vector, ccs_matrix, ccs_residuals
  use parallel_types, only: parallel_environment
  use parallel_types_mpi, only: parallel_environment_mpi
  use utils, only: exit_print
  use fields, only: get_field, count_fields, get_field_idx, get_is_field_solved
  use mat, only: mat_vec_product
  use vec, only: vec_aypx
  use meshing, only: get_global_num_cells
  use timestepping, only: get_current_step, get_current_time
  use kinds, only: CCS_MPI_PRECISION
  use logging, only: log_unit_out

contains

  !v Allocate residuals arrays
  subroutine init_residuals(flow, residuals)
    type(fluid), intent(in) :: flow                              !< The structure containting all the fluid fields
    type(ccs_residuals), intent(inout) :: residuals
    integer(ccs_int) :: nfields

    call count_fields(flow, nfields)

    allocate (residuals%L2(nfields))
    allocate (residuals%Linfty(nfields))

    residuals%L2 = 0.0_ccs_real
    residuals%Linfty = 0.0_ccs_real

  end subroutine

  !v Computes residuals from a matric, vector and rhs.
  ! Residuals are stored in 'res' variable and their norms in the residuals data structure
  subroutine compute_residuals(flow, M, phi, rhs, res, residuals)
    type(fluid), intent(in) :: flow                              !< The structure containting all the fluid fields
    class(ccs_matrix), intent(in) :: M
    class(field), intent(in) :: phi
    class(ccs_vector), intent(in) :: rhs
    class(ccs_vector), intent(inout) :: res
    type(ccs_residuals), intent(inout) :: residuals
    integer(ccs_int) :: ifield

    call mat_vec_product(M, phi%values, res)
    call vec_aypx(rhs, -1.0_ccs_real, res)

    call get_field_idx(flow, phi, ifield)
    call normalise_residuals(res, ifield, residuals)

  end subroutine

  !v Computes normalised residuals and stores them in residuals structure.
  ! Residuals are normalised by cell volumes**(2/3) and their L2 (squared) and Linfinity norms are stored
  subroutine normalise_residuals(res, ifield, residuals)
    use types, only: ccs_vector, cell_locator
    use vec, only: get_vector_data_readonly, restore_vector_data_readonly
    use meshing, only: get_local_num_cells, create_cell_locator, get_volume
    use profiler, only: profiler_begin_region, profiler_end_region

    class(ccs_vector), intent(inout) :: res !< residuals vector
    type(ccs_residuals), intent(inout) :: residuals
    integer(ccs_int), intent(in) :: ifield
    real(ccs_real) :: L2sq     !< output L2 norm squared
    real(ccs_real) :: Linfty   !< output Linfinity norm
    real(ccs_real), dimension(:), pointer :: res_data
    integer(ccs_int) :: local_num_cells
    integer(ccs_int) :: index_p
    real(ccs_real) :: V, normalised_res
    type(cell_locator) :: loc_p

    call get_local_num_cells(local_num_cells)
    L2sq = 0.0_ccs_real
    Linfty = 0.0_ccs_real
    call get_vector_data_readonly(res, res_data)

    call profiler_begin_region("[OMP] normalise residuals loop")

    !$omp parallel do default(none) schedule(static) &
    !$omp shared(local_num_cells, res_data)          &
    !$omp private(index_p, loc_p, V, normalised_res) &
    !$omp reduction(+:L2sq) reduction(max:Linfty)
    do index_p = 1, local_num_cells
      call create_cell_locator(index_p, loc_p)
      call get_volume(loc_p, V)
      normalised_res = res_data(index_p) / (V**(2.0_ccs_real / 3.0_ccs_real))
      L2sq = L2sq + normalised_res**2
      Linfty = max(abs(normalised_res), Linfty)
    end do
    !$omp end parallel do

    call profiler_end_region("[OMP] normalise residuals loop")

    residuals%L2(ifield) = L2sq
    residuals%Linfty(ifield) = Linfty

    call restore_vector_data_readonly(res, res_data)

  end subroutine

  !v Check if a particular residual is below the res_target
  logical function is_converged(residuals, phi, ifield) result(converged)
    type(ccs_residuals), intent(in) :: residuals !< residuals object 
    class(field), intent(in) :: phi !< field to get the residual target from
    integer(ccs_int), intent(in) :: ifield !< index of the residual in 'residuals' object
    real(ccs_real) :: res_target

    res_target = phi%solver_parameters%res_target

    select case (phi%solver_parameters%res_norm)
    case (L2)
      converged = residuals%L2(ifield) <= res_target
    case (Linfty)
      converged = residuals%Linfty(ifield) <= res_target
    end select

  end function

  !v Get the maximum residuals for a specific norm
  real(ccs_real) function get_max_residuals(residuals, norm) result(max_value)
    type(ccs_residuals), intent(in) :: residuals
    integer, intent(in) :: norm

    select case (norm)
    case (L2)
      max_value = maxval(residuals%L2)
    case (Linfty)
      max_value = maxval(residuals%Linfty)
    end select

  end function

  !v Computes global residuals from local quantities
  ! Requires reduction operation
  subroutine compute_global_residuals(par_env, flow, residuals)
    use mpi
    use parallel_types_mpi, only: parallel_environment_mpi

    class(parallel_environment), allocatable, intent(in) :: par_env !< The parallel environment
    type(fluid), intent(in) :: flow                              !< The structure containting all the fluid fields
    type(ccs_residuals), intent(inout) :: residuals
    integer(ccs_int) :: nfields, ifield, global_num_cells
    integer :: ierr

    call count_fields(flow, nfields)
    call get_global_num_cells(global_num_cells)

    select type (par_env)
    type is (parallel_environment_mpi)
      do ifield = 1, nfields
        ! Computes RMS from L2 squared
        call MPI_ALLREDUCE(MPI_IN_PLACE, residuals%L2(ifield), 1, CCS_MPI_PRECISION, MPI_SUM, par_env%comm, ierr)
        residuals%L2(ifield) = sqrt(residuals%L2(ifield) / real(global_num_cells, kind=ccs_real))

        ! Linf
        call MPI_ALLREDUCE(MPI_IN_PLACE, residuals%Linfty(ifield), 1, CCS_MPI_PRECISION, MPI_MAX, par_env%comm, ierr)
      end do
    class default
      call error_abort("invalid parallel environment")
    end select

  end subroutine

  !v Prints to stdout and to residuals.log file residuals norms for current iteration
  subroutine print_residuals(par_env, flow, itr, residuals)

    use parallel, only: is_root

    class(parallel_environment), allocatable, intent(in) :: par_env !< The parallel environment
    type(fluid), intent(inout) :: flow                              !< Container for flow fields
    integer(ccs_int), intent(in) :: itr                             !< Iteration count
    type(ccs_residuals), intent(in) :: residuals

    ! Local variables
    integer :: io_unit
    integer(ccs_int) :: step                            !< The current time-step
    real(ccs_real) :: time                              !< The current time
    integer(ccs_int) :: ifield
    integer(ccs_int) :: nfields              ! Number of variables (u,v,w,p,etc)
    integer(ccs_int) :: i
    character(len=60) :: prefix           ! prefix for residual norms
    logical, save :: first_time = .true.  ! Whether first time this subroutine is called

    class(field), pointer :: phi
    logical :: phi_sol


    ! Print residuals
    if (is_root(par_env)) then
      call count_fields(flow, nfields)
      call get_current_step(step)
      call get_current_time(time)

      if (first_time) then
        ! Write header
        open (newunit=io_unit, file="residuals.log", status="replace", form="formatted")

        write (log_unit_out, *)
        if (step >= 0) then
          write (log_unit_out, '(a5, 1x, a12, 1x, a6, 2x)', advance='no') 'Step', 'time', 'Iter'

          write (io_unit, '(a6, 1x, a12, 1x, a6)', advance='no') '#step', 'time', 'iter'
        else
          write (log_unit_out, '(a5)', advance='no') 'Iter'

          write (io_unit, '(a6)', advance='no') '#iter'
        end if

        do i = 1, 2
          if (i == 1) then
            prefix = "L2_"
          else
            prefix = "Linf_"
          end if
          do ifield = 1, nfields
            call get_field(flow, ifield, phi)
            call get_is_field_solved(phi, phi_sol)

            if (phi_sol) then
              write (log_unit_out, '(1x,a12)', advance='no') phi%name
              write (io_unit, '(1x,a12)', advance='no') trim(prefix) // phi%name
            end if
          end do
        end do
        write (log_unit_out, *)
        write (io_unit, *)
        first_time = .false.
      else
        open (newunit=io_unit, file="residuals.log", status="old", form="formatted", position="append")
      end if

      ! Write step, iteration and residuals
      if (step >= 0) then
        write (log_unit_out, '(i5,1x,e12.4,1x,i6,1x)', advance='no') step, time, itr
        write (io_unit, '(i6,1x,e12.4,1x,i6,1x)', advance='no') step, time, itr
      else
        write (log_unit_out, '(i5,1x,e12.4,1x)', advance='no') step, time, itr
        write (io_unit, '(i6,1x,e12.4,1x)', advance='no') step, time, itr
      end if

      do ifield = 1, nfields
        call get_field(flow, ifield, phi)
        call get_is_field_solved(phi, phi_sol)

        if (phi_sol) then
          write (log_unit_out, '(e12.4,1x)', advance='no') residuals%L2(ifield)
          write (io_unit, '(e12.4,1x)', advance='no') residuals%L2(ifield)
        end if
      end do

      do ifield = 1, nfields
        call get_field(flow, ifield, phi)
        call get_is_field_solved(phi, phi_sol)

        if (phi_sol) then
          write (log_unit_out, '(e12.4,1x)', advance='no') residuals%Linfty(ifield)
          write (io_unit, '(e12.4,1x)', advance='no') residuals%Linfty(ifield)
        end if
      end do

      write (log_unit_out, *)
      write (io_unit, *)
      close (io_unit)
    end if

  end subroutine

end module residuals
