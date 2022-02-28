!> @brief Submodule file pv_coupling_simple.smod
!
!> @details Implementation of the SIMPLE algorithm for pressure-velocity coupling.

submodule (pv_coupling) pv_coupling_simple

  use fv, only : compute_fluxes

  implicit none

  contains

  module subroutine solve_nonlinear(u, v, p, pp, M, solution, source, cell_mesh)

    ! Arguments
    class(field), intent(inout) :: u, v, p
    class(matrix), intent(inout) :: M
    class(vector), intent(inout) :: solution, source
    type(mesh), intent(in) :: cell_mesh

    ! Local variables
    integer(accs_int) :: i

    outerloop: do i = it_start, it_end
      call calculate_momentum(u, v, cell_mesh, bcs, cps, M, source)

      call calculate_pressure_correction()

      call update_pressure()

      ! Todo:
      !call calculate_scalars()

      call check_convergence()
      if (converged) then
        exit outerloop
      endif

    end do outerloop

  end subroutine

  subroutine calculate_momentum(u, v, cell_mesh, bcs, cps, M, vec)

    ! Arguments
    type(field), intent(inout)    :: u, v
    type(mesh), intent(in)        :: cell_mesh
    type(bc_config), intent(in)   :: bcs
    integer(accs_int), intent(in) :: cps
    class(matrix), intent(inout)  :: M
    class(vector), intent(inout)  :: vec

    ! Local variables
    
    ! Calculate fluxes
    call compute_fluxes(u, v, cell_mesh, bcs, cps, M, vec)

    ! Assembly of coefficient matrix and source vector
    call update(M)
    call update(vec)

    ! Create linear solver
    call set_linear_system(scalar_linear_system, vec, scalar, M, par_env)
    call create_solver(scalar_linear_system, scalar_solver)

    ! Solve the linear system
    call solve(scalar_solver)

  end subroutine calculate_momentum


  subroutine calculate_pressure_correction()

    mat_coeffs%mode = insert_mode
    
    ! Loop over cells
    do i = 1, cell_mesh%nlocal
      call set_cell_location(cell_location, cell_mesh, i)
      call get_global_index(cell_location, idxg)
      call count_neighbours(cell_location, nnb)

      allocate(mat_coeffs%rglob(1))
      allocate(mat_coeffs%cglob(1 + nnb))
      allocate(mat_coeffs%val(1 + nnb))

      row = idxg
      coeff_p = 0.0_accs_real

      ! Loop over faces
      do j = 1, nnb
        call set_neighbour_location(nb_location, cell_location, j)
        call get_boundary_status(nb_location, is_boundary)

        if (.not. is_boundary) then
          ! Interior face
          call set_face_location(face_location, cell_mesh, i, j)
          call get_face_area(face_location, A)
          coeff_f = (1.0 / cell_mesh%h) * A

          call get_global_index(nb_location, nbidxg)

          coeff_p = coeff_p - coeff_f
          coeff_nb = coeff_f
          col = nbidxg
        else
          col = -1
          coeff_nb = 0.0_accs_real
        endif
        call pack_entries(mat_coeffs, 1, j+1, row, col, coeff_nb)

      end do

      ! Add the diagonal entry
      col = row
      call pack_entries(mat_coeffs, 1, 1, row, col, coeff_p)

      ! Set the values
      call set_values(mat_coeffs, M)

      deallocate(mat_coeffs%rglob, mat_coeffs%cglob, mat_coeffs%val)

    end do

  end subroutine calculate_pressure_correction


  subroutine update_pressure()


  end subroutine update_pressure


  subroutine calculate_scalars()


  end subroutine calculate_scalars


  subroutine check_convergence(converged)

    ! Arguments
    logical, intent(inout) :: converged

  end subroutine check_convergence

end submodule pv_coupling_simple