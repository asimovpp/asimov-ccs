submodule (timestepping) timestepping_first_order

  !use ...
  
  implicit none
  
  real(ccs_real) :: dt

contains

  module subroutine apply_timestep(mesh, phi, diag, M, b)
    use kinds, only: ccs_int
    use mat, only : set_matrix_diagonal, get_matrix_diagonal
    use vec, only: get_vector_data, restore_vector_data
    use utils, only: update, finalise
    
    type(ccs_mesh), intent(in) :: mesh
    class(field), intent(in) :: phi
    class(ccs_vector), intent(inout) :: diag
    class(ccs_matrix), intent(inout) :: M
    class(ccs_vector), intent(inout) :: b

    real(ccs_real), dimension(:), pointer :: diag_data
    real(ccs_real), dimension(:), pointer :: b_data
    real(ccs_real), dimension(:), pointer :: phi_data
    
    integer(ccs_int) :: i

    ! V = mesh%volumes
    dt = 0.9 / 1.0 * mesh%h
    call finalise(M)
    call get_matrix_diagonal(M, diag)
   
    call get_vector_data(phi%old_values, phi_data)
    call get_vector_data(diag, diag_data)
    call update(b)
    call get_vector_data(b, b_data)

    
    do i = 1, mesh%nlocal
    ! A = A + V/dt
      diag_data(i) = diag_data(i) + mesh%volumes(i) / dt

    ! b = b + V/dt * phi_old
      b_data(i) = b_data(i) + mesh%volumes(i) / dt * phi_data(i)
    end do
    call restore_vector_data(phi%old_values, phi_data)
    call restore_vector_data(diag, diag_data)
    call restore_vector_data(b, b_data)
    call set_matrix_diagonal(diag, M)
  end subroutine apply_timestep

end submodule 
