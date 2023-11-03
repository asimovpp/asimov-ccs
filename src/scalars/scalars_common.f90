!v Submodule file scalars_common.smod
!
!  Implementation of the scalar transport subroutines

submodule(scalars) scalars_common
#include "ccs_macros.inc"
  use constants, only: field_u, field_v, field_w, field_p, field_p_prime, field_mf, field_viscosity, field_density

  use kinds, only: ccs_int, ccs_real !< added here
  use types, only: ccs_matrix, ccs_vector, &
       vector_spec, matrix_spec, &
       linear_solver, equation_system

  use fv, only: compute_fluxes, update_gradient
  use timestepping, only: update_old_values, get_current_step, apply_timestep

  use vec, only: create_vector, get_vector_data, restore_vector_data !< added here
  use mat, only: create_matrix, set_nnz
  use solver, only: create_solver, solve, set_equation_system

  use meshing, only: get_max_faces
  use utils, only: get_field, update, initialise, finalise, set_size, debug_print, &
       zero
  
  logical, save :: first_call = .true.
  integer(ccs_int), save :: previous_step = -1

  !> List of fields not to be updated as transported scalars
  integer(ccs_int), dimension(8), parameter :: skip_fields = &
       (/ field_u, field_v, field_w, &
          field_p, field_p_prime, &
          field_mf, field_viscosity, field_density /)

contains
  
  !> Subroutine to perform scalar transport for all scalar fields.
  module subroutine update_scalars(par_env, mesh, flow)
    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    type(ccs_mesh), intent(in) :: mesh                              !< the mesh
    type(fluid), intent(inout) :: flow                              !< The structure containting all the fluid fields

    class(ccs_matrix), allocatable :: M
    class(ccs_vector), allocatable :: rhs
    class(ccs_vector), allocatable :: D
    
    integer(ccs_int) :: nfields  ! Number of variables in the flowfield
    integer(ccs_int) :: s        ! Scalar field counter
    integer(ccs_int) :: field_id ! The field's numeric identifier
    class(field), pointer :: phi ! The scalar field

    logical :: do_update
    integer(ccs_int) :: current_step

    type(vector_spec) :: vec_properties
    type(matrix_spec) :: mat_properties

    integer(ccs_int) :: max_faces

    ! Initialise equation system (reused across scalars)
    call dprint("SCALAR: init")
    call initialise(vec_properties)
    call initialise(mat_properties)

    call dprint("SCALAR: setup matrix")
    call get_max_faces(mesh, max_faces)
    call set_size(par_env, mesh, mat_properties)
    call set_nnz(max_faces + 1, mat_properties)
    call create_matrix(mat_properties, M)

    call dprint("SCALAR: setup RHS")
    call set_size(par_env, mesh, vec_properties)
    call create_vector(vec_properties, rhs)
    call create_vector(vec_properties, D)

    ! Check whether we need to update the old values
    call dprint("SCALAR: check new timestep")
    do_update = .false.
    
    call get_current_step(current_step)
    
    if (first_call) then
       first_call = .false.
       do_update = .true.
    else if (previous_step /= current_step) then
       do_update = .true.
    end if

    previous_step = current_step

    ! Transport the scalars
    call count_fields(flow, nfields)
    do s = 1, nfields
       call get_field_id(flow, s, field_id)
       if (any(skip_fields == field_id)) then
          ! Not a scalar to solve
          cycle
       end if

       call get_field(flow, field_id, phi)
       if (do_update) then
          call update_old_values(phi)
       end if

       call transport_scalar(par_env, mesh, flow, M, rhs, D, phi)
         
       ! density values are turned to exponential here
       
    end do
    
  end subroutine update_scalars

  !> Subroutine to transport a scalar field.
  subroutine transport_scalar(par_env, mesh, flow, M, rhs, D, phi)

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    type(ccs_mesh), intent(in) :: mesh                              !< the mesh
    type(fluid), intent(inout) :: flow                              !< The structure containting all the fluid fields
    class(ccs_matrix), allocatable, intent(inout) :: M
    class(ccs_vector), allocatable, intent(inout) :: rhs
    class(ccs_vector), intent(inout) :: D
    class(field), intent(inout) :: phi ! The scalar field

    class(field), pointer :: mf  ! The advecting velocity field
    class(field), pointer :: viscosity  ! viscosity
    class(field), pointer :: density ! density
    class(linear_solver), allocatable :: lin_solver
    type(equation_system) :: lin_system
    
    !print*,"inside transport_scalar"
    call initialise(lin_system)
    call zero(rhs)
    call zero(M)

    call dprint("SCALAR: compute coefficients")
    call get_field(flow, field_mf, mf)
    call get_field(flow, field_viscosity, viscosity) 
    call get_field(flow, field_density, density)
    call compute_fluxes(phi, mf, viscosity, density, mesh, 0, M, rhs)
    call apply_timestep(mesh, phi, D, M, rhs)

    call dprint("SCALAR: assemble linear system")
    call update(M)
    call update(rhs)
    call finalise(M)

    ! density values are single digit (same as i/p)  

    if (allocated(phi%values%name)) then
       call set_equation_system(par_env, rhs, phi%values, M, lin_system, phi%values%name)
    else
       call set_equation_system(par_env, rhs, phi%values, M, lin_system)
    end if

    ! density values are single digit (same as i/p)  
       
    call dprint("SCALAR: solve linear system")
    call create_solver(lin_system, lin_solver)

    ! density values are single digit (same as i/p)  
    call solve(lin_solver)
    !call get_vector_data(density%values, density_data)
    !do index_p = 1, 5
      !print*,"IUS cell=",index_p,"density=",density_data(index_p)
    !end do 
    !call restore_vector_data(density%values, density_data)
    ! density values are exponential

    call update_gradient(mesh, phi)

    deallocate(lin_solver)
    
  end subroutine transport_scalar

  !> Get the count of stored fields - probably belongs somewhere else
  subroutine count_fields(flow, nfields)

    type(fluid), intent(in) :: flow          !< The flowfield
    integer(ccs_int), intent(out) :: nfields !< The count of fields

    nfields = size(flow%field_names)
    
  end subroutine count_fields

  !> Get the numeric ID of the i'th field
  subroutine get_field_id(flow, s, field_id)

    type(fluid), intent(in) :: flow           !< The flowfield
    integer(ccs_int), intent(in) :: s         !< The field counter
    integer(ccs_int), intent(out) :: field_id !< The field ID

    field_id = flow%field_names(s)

  end subroutine get_field_id
  
end submodule scalars_common
