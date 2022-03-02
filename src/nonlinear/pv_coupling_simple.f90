!> @brief Submodule file pv_coupling_simple.smod
!
!> @details Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule (pv_coupling) pv_coupling_simple

  use fv, only : compute_fluxes

  implicit none

  contains

  !> @ brief Solve Navier-Stokes equations using the SIMPLE algorithm
  !
  !> @param[in,out] u, v   - arrays containing velocity fields in x, y directions
  !> @param[in,out] p      - array containing pressure values
  !> @param[in,out] pp     - array containing pressure-correction values
  !> @param[in,out] M      - coefficient matrix
  !> @param[in,out] sol    - solution vector
  !> @param[in,out] source - source vector
  !> @param[in]     cell_mesh - the mesh

  module subroutine solve_nonlinear(par_env, cell_mesh, u, v, p, pp)

    ! Arguments
    class(parallel_environment), intent(in) :: par_env
    type(mesh), intent(in) :: cell_mesh
    class(field), intent(inout) :: u, v, p, pp
    
    ! Local variables
    integer(accs_int) :: i
    class(vector), allocatable, target :: source
    class(vector), allocatable :: sol
    class(matrix), allocatable, target :: M
    

    type(vector_init_data) :: vec_sizes
    type(matrix_init_data) :: mat_sizes
    type(linear_system)    :: lin_system

    ! Initialise linear system
    call initialise(mat_sizes)
    call initialise(vec_sizes)
    call initialise(lin_system)

    ! Create coefficient matrix
    call set_global_size(mat_sizes, cell_mesh%nglobal, cell_mesh%nglobal, par_env)
    call set_nnz(mat_sizes, 5)
    call create_matrix(mat_sizes, M)

    ! Create RHS and solution vectors
    call set_global_size(vec_sizes, cell_mesh%nglobal, par_env)
    call create_vector(vec_sizes, source)
    call create_vector(vec_sizes, sol)

    outerloop: do i = it_start, it_end
      ! Solve momentum equation with guessed pressure and velocity fields
      call calculate_velocity(cell_mesh, bcs, cps, M, source, lin_system u, v)

      ! Calculate pressure correction from mass imbalance
      call calculate_pressure_correction()

      ! Update pressure field with pressure correction
      call update_pressure()

      ! Update velocity with velocity correction
      call update_velocity(u, v, pp)

      ! Todo:
      !call calculate_scalars()

      call check_convergence()
      if (converged) then
        exit outerloop
      endif

    end do outerloop

  end subroutine solve_nonlinear

  subroutine calculate_velocity(cell_mesh, bcs, cps, M, vec, lin_sys, u, v)

    ! Arguments
    type(mesh), intent(in)        :: cell_mesh
    type(bc_config), intent(in)   :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(inout)  :: M
    class(vector), intent(inout)  :: vec
    type(linear_system), intent(inout) :: lin_sys
    type(field), intent(inout)    :: u, v

    ! Local variables
    class(linear_solver), allocatable :: lin_solver
    
    ! Calculate fluxes and populate coefficient matrix
    call compute_fluxes(u, v, cell_mesh, bcs, cps, M, vec)

    ! Assembly of coefficient matrix and source vector
    call update(M)
    call update(vec)

    ! Create linear solver
    call set_linear_system(lin_sys, vec, u%vec, M, par_env)
    call create_solver(lin_sys, lin_solver)

    ! Solve the linear system
    call solve(lin_solver)

    ! Clean up
    deallocate(lin_solver)

  end subroutine calculate_velocity


  subroutine calculate_pressure_correction(cell_mesh, u, v, pp)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(field), intent(in) :: u, v
    class(field), intent(out) :: pp

    ! Local variables


    mat_coeffs%mode = insert_mode
    
    ! Loop over cells
    do icell = 1, cell_mesh%nlocal
      call set_cell_location(cell_location, cell_mesh, icell)
      call get_global_index(cell_location, idxg)
      call count_neighbours(cell_location, nnb)

      allocate(mat_coeffs%rglob(1))
      allocate(mat_coeffs%cglob(1 + nnb))
      allocate(mat_coeffs%val(1 + nnb))

      row = idxg
      coeff_p = 0.0_accs_real
      r = 0.0_accs_real

      ! Loop over faces
      do jface = 1, nnb
        call set_neighbour_location(nb_location, cell_location, jface)
        call get_boundary_status(nb_location, is_boundary)

        if (.not. is_boundary) then
          ! Interior face
          call set_face_location(face_location, cell_mesh, icell, jface)
          call get_face_area(face_location, A)
          coeff_f = (1.0 / cell_mesh%h) * A

          call get_global_index(nb_location, nbidxg)

          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = nbidxg

          ! RHS vector
          r = r + calc_mass_flux(u, v, ngb_idx, self_idx, A, face_normal, bc_flag)

        else
          col = -1
          coeff_nb = 0.0_accs_real
        endif
        call pack_entries(mat_coeffs, 1, jface+1, row, col, coeff_nb)
        call pack_entries(vec_values, 1, idx, r)

      end do

      ! Add the diagonal entry
      col = row
      call pack_entries(mat_coeffs, 1, 1, row, col, coeff_p)

      ! Set the values
      call set_values(mat_coeffs, M)
      call set_values(vec_values, b)

      deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)

    end do

    call update(M)
    call update(b)

    ! Create linear solver
    call set_linear_system(lin_system, b, p, M, par_env)
    call create_solver(lin_system, lin_solver)
    call solve(lin_solver)


  end subroutine calculate_pressure_correction


  subroutine update_pressure(cell_mesh, pp, p)

    ! Arguments
    class(mesh), intent(in) :: cell_mesh
    class(field), intent(in) :: pp
    class(field), intent(inout) :: p


    ! Local variables
    integer(accs_int) :: icell
    real(accs_real) :: alpha   !< Under-relaxation factor


    ! Loop over cells
    do icell = 1, cell_mesh%nlocal
      p(icell) = p(icell) + alpha*pp(icell)
    end do
    
     
  end subroutine update_pressure

  subroutine update_velocity()

    ! Arguments


    ! Local variables


    ! Loop over cells
    do icell = 1, cell_mesh%nlocal
      
    end do




  end subroutine update_velocity


  subroutine calculate_scalars()


  end subroutine calculate_scalars


  subroutine check_convergence(converged)

    ! Arguments
    logical, intent(inout) :: converged

  end subroutine check_convergence

end submodule pv_coupling_simple