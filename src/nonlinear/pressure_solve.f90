!v Module with functions assiciated with the Pressure solve
!

module pressure_solve
#include "ccs_macros.inc"

  use profiler
  use kinds, only: ccs_int, ccs_real
  use fv_equations, only: poisson_equation
  use types, only: vector_spec, ccs_vector, ccs_matrix, equation_system, &
                   linear_solver, vector_values, cell_locator, ccs_residuals, &
                   matrix_values, matrix_values_spec, field
  use fv, only: update_gradient
  use vec, only: scale_vec, get_vector_data_readonly, restore_vector_data_readonly, &
                 create_vector_values, zero_vector, vec_aypx, &
                 vec_sum, vec_shift
  use mat, only: set_matrix_values_spec_nrows, set_matrix_values_spec_ncols, create_matrix_values, &
                 check_operator_symmetry
  use utils, only: update, finalise, set_values, zero, clear_entries, set_entry, set_row, set_mode, &
                   str, exit_print, debug_print
  use core, only: ccs_options
  use bc_constants
  use solver, only: create_solver, solve, set_equation_system, axpy, set_solver_method, set_solver_precon
  use constants, only: add_mode, insert_mode, ndim
  use meshing, only: create_cell_locator, get_local_num_cells, get_global_num_cells, get_max_faces
  use parallel, only: is_root
  use parallel_types, only: parallel_environment
 

  implicit none

  contains
  !v Solves the pressure correction equation
  !
  !  Solves the pressure correction equation formed by the mass-imbalance.
  subroutine calculate_pressure_correction(par_env, run_options, invA, M, vec, lin_sys, p_prime, lin_solver)


    ! Arguments
    class(parallel_environment), allocatable, intent(in) :: par_env !< the parallel environment
    type(ccs_options), intent(in) :: run_options                    !< runtime options
    class(ccs_vector), intent(inout) :: invA                        !< inverse diagonal momentum coefficients
    class(ccs_matrix), allocatable, intent(inout) :: M              !< matrix object
    class(ccs_vector), allocatable, intent(inout) :: vec            !< the RHS vector
    type(equation_system), intent(inout) :: lin_sys                 !< linear system object
    class(field), intent(inout) :: p_prime                          !< the pressure correction field
    class(linear_solver), allocatable, intent(inout) :: lin_solver  !< Linear solver that is being reused

    ! Local variables
    type(matrix_values) :: mat_coeffs
    type(vector_values) :: vec_values
    type(cell_locator) :: loc_p
    type(matrix_values_spec) :: mat_val_spec

    integer(ccs_int) :: local_num_cells, global_num_cells
    integer(ccs_int) :: index_p
    integer(ccs_int) :: max_faces
    integer(ccs_int) :: max_stencil_width

    real(ccs_real), dimension(:), pointer :: invA_data

    real(ccs_real) :: vec_mean

    type(poisson_equation) :: pois_eqn

    call profiler_begin_region("Building coefficients")
    ! First zero matrix
    call zero(M)

    ! To ensure compatibility of the Poisson equation with Neumann conditions subtract the mean of the RHS from the RHS
    if (all(p_prime%bcs%bc_types == bc_type_neumann)) then
      ! Enforce zero forcing integral
      call vec_sum(vec, vec_mean)
      call get_global_num_cells(global_num_cells)
      vec_mean = vec_mean / global_num_cells
      call vec_shift(-vec_mean, vec)
    end if

    ! The computed mass imbalance is +ve, to have a +ve diagonal coefficient we need to negate this.
    call dprint("P': negate RHS")
    call scale_vec(-1.0_ccs_real, vec)
    call update(vec)

    call zero_vector(p_prime%x_gradients)
    call update(p_prime%x_gradients)
    call zero_vector(p_prime%y_gradients)
    call update(p_prime%y_gradients)
    call zero_vector(p_prime%z_gradients)
    call update(p_prime%z_gradients)
    call zero_vector(p_prime%values)
    call update(p_prime%values)

    call update_gradient(p_prime)

    call create_vector_values(1_ccs_int, vec_values)
    call set_mode(add_mode, vec_values)

    call update(M)

    call dprint("P': get invA")
    call get_vector_data_readonly(invA, invA_data)

    call dprint("P': setup equation buffers")
    call get_max_faces(max_faces)
    max_stencil_width = max_faces + 1
    call set_matrix_values_spec_nrows(1_ccs_int, mat_val_spec)
    call set_matrix_values_spec_ncols(max_stencil_width, mat_val_spec)
    call create_matrix_values(mat_val_spec, mat_coeffs)
    call set_mode(insert_mode, mat_coeffs)

    call pois_eqn%init(max_faces)

    ! Loop over cells
    call dprint("P': cell loop")
    call get_local_num_cells(local_num_cells)
    do index_p = 1, local_num_cells
      call clear_entries(mat_coeffs)
      call clear_entries(vec_values)

      call create_cell_locator(index_p, loc_p)

      call pois_eqn%gather(p_prime, loc_p, invA_data)
      call pois_eqn%apply()
      call pois_eqn%flush_row(mat_coeffs, vec_values)

      call set_values(mat_coeffs, M)
      call set_values(vec_values, vec)
    end do

    deallocate (mat_coeffs%global_row_indices)
    deallocate (mat_coeffs%global_col_indices)
    deallocate (mat_coeffs%values)
    deallocate (vec_values%global_indices)
    deallocate (vec_values%values)

    call dprint("P': restore invA")
    call restore_vector_data_readonly(invA, invA_data)

    ! Assembly of coefficient matrix and source vector
    call dprint("P': assemble matrix, RHS")
    call update(M)
    call update(vec)
    call finalise(M)

    call profiler_end_region("Building coefficients")

    ! Create linear solver
    call dprint("P': create lin sys")
    if (allocated(p_prime%values%name)) then
      call set_equation_system(par_env, vec, p_prime%values, M, lin_sys, p_prime%values%name)
    else
      call set_equation_system(par_env, vec, p_prime%values, M, lin_sys)
    end if

    ! Verify operator is symmetric
    if (run_options%solve%debug) then
      call check_operator_symmetry(M)
    end if

    call create_solver(lin_sys, lin_solver)

    ! Customise linear solver
    call set_solver_method(p_prime%solver_parameters%solver_name, lin_solver)
    call set_solver_precon(p_prime%solver_parameters%precon_name, lin_solver)

    ! Solve the linear system
    call dprint("P': solve")
    call solve(lin_solver)

  end subroutine calculate_pressure_correction

end module