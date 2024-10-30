!v Submodule file scalars_common.smod
!
!  Implementation of the scalar transport subroutines

submodule(scalars) scalars_common
#include "ccs_macros.inc"

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

  implicit none
  
  logical, save :: first_call = .true.
  integer(ccs_int), save :: previous_step = -1

  !> List of fields not to be updated as transported scalars
  character (len=20), dimension(*), parameter :: skip_fields = &
  [ character(len=20) :: "u", "v", "w", &
                         "p", "p_prime", &
                         "mf", "viscosity", "density" ]

contains
  
  !> Subroutine to perform scalar transport for all scalar fields.
  module subroutine update_scalars(par_env, mesh, eval_sources, flow)
    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    type(ccs_mesh), intent(in) :: mesh                              !< the mesh
    interface
      !v Subroutine to evaluate source terms, case-specific.
      !
      !  Note this should return the integrated source.
      subroutine eval_sources(flow, phi, R, S)
        use types, only: fluid, field, ccs_vector
        type(fluid), intent(in) :: flow !< Provides access to full flow field
        class(field), intent(in) :: phi !< Field being transported
        class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
        class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
      end subroutine eval_sources
    end interface
    type(fluid), intent(inout) :: flow                              !< The structure containting all the fluid fields

    class(ccs_matrix), allocatable :: M
    class(ccs_vector), allocatable :: rhs
    class(ccs_vector), allocatable :: D
    class(ccs_vector), allocatable :: source
    
    integer(ccs_int) :: nfields  ! Number of variables in the flowfield
    integer(ccs_int) :: s        ! Scalar field counter
    character(len=:), allocatable :: field_name ! The field's name
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
    call get_max_faces(max_faces)
    call set_size(par_env, mesh, mat_properties)
    call set_nnz(max_faces + 1, mat_properties)
    call create_matrix(mat_properties, M)

    call dprint("SCALAR: setup RHS")
    call set_size(par_env, mesh, vec_properties)
    call create_vector(vec_properties, rhs)
    call create_vector(vec_properties, D)
    call create_vector(vec_properties, source)

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
       call get_field_name(flow, s, field_name)
      
       if (any(skip_fields == field_name)) then
          ! Not a scalar to solve
          cycle
       end if

       call get_field(flow, field_name, phi)
       if (do_update) then
          call update_old_values(phi)
       end if

       call transport_scalar(par_env, flow, eval_sources, M, rhs, D, source, phi)
    end do
    
  end subroutine update_scalars

  !> Subroutine to transport a scalar field.
  subroutine transport_scalar(par_env, flow, eval_sources, M, rhs, D, S, phi)

    class(parallel_environment), allocatable, intent(in) :: par_env !< parallel environment
    type(fluid), intent(inout) :: flow                              !< The structure containting all the fluid fields
    interface
      !v Subroutine to evaluate source terms, case-specific.
      !
      !  Note this should return the integrated source.
      subroutine eval_sources(flow, phi, R, S)
        use types, only: fluid, field, ccs_vector
        type(fluid), intent(in) :: flow !< Provides access to full flow field
        class(field), intent(in) :: phi !< Field being transported
        class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
        class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
      end subroutine eval_sources
    end interface
    class(ccs_matrix), allocatable, intent(inout) :: M
    class(ccs_vector), allocatable, intent(inout) :: rhs
    class(ccs_vector), intent(inout) :: D !< Working vector for equation diagonal
    class(ccs_vector), intent(inout) :: S !< Working vector for equation RHS
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
    call get_field(flow, "mf", mf)
    call get_field(flow, "viscosity", viscosity) 
    call get_field(flow, "density", density)
    call compute_fluxes(phi, mf, viscosity, density, 0, M, rhs)
    call apply_timestep(phi, D, M, rhs)

    call apply_sources(flow, phi, eval_sources, D, S, M, rhs)

    call dprint("SCALAR: assemble linear system")
    call update(M)
    call update(rhs)
    call finalise(M)

    if (allocated(phi%values%name)) then
       call set_equation_system(par_env, rhs, phi%values, M, lin_system, phi%values%name)
    else
       call set_equation_system(par_env, rhs, phi%values, M, lin_system)
    end if
     
    call dprint("SCALAR: solve linear system")
    call create_solver(lin_system, lin_solver)

    call solve(lin_solver)
    call update_gradient(phi)

    deallocate(lin_solver)
    
  end subroutine transport_scalar

  !> Get the count of stored fields - probably belongs somewhere else
  pure subroutine count_fields(flow, nfields)

    type(fluid), intent(in) :: flow          !< The flowfield
    integer(ccs_int), intent(out) :: nfields !< The count of fields

    nfields = size(flow%fields)
    
  end subroutine count_fields

  !> Get the name of the i'th field
  subroutine get_field_name(flow, s, field_name)

    type(fluid), intent(in) :: flow                          !< The flowfield
    integer(ccs_int), intent(in) :: s                        !< The field counter
    character(len=:), allocatable, intent(out) :: field_name !< The field name

    class(field), pointer :: phi
   
    call get_field(flow, s, phi)
    field_name = phi%name
    nullify(phi)

  end subroutine get_field_name
 
  !> Compute source terms and add to the equation system
  subroutine apply_sources(flow, phi, eval_sources, R, S, M, rhs)

    use fv, only: add_fixed_source, add_linear_source
    
    type(fluid), intent(in) :: flow !< Provides access to full flow field
    class(field), intent(in) :: phi !< The transported field
    interface
      !v Subroutine to evaluate source terms, case-specific.
      !
      !  Note this should return the integrated source.
      subroutine eval_sources(flow, phi, R, S)
        use types, only: fluid, field, ccs_vector
        type(fluid), intent(in) :: flow !< Provides access to full flow field
        class(field), intent(in) :: phi !< Field being transported
        class(ccs_vector), intent(inout) :: R !< Work vector (for evaluating linear/implicit sources)
        class(ccs_vector), intent(inout) :: S !< Work vector (for evaluating fixed/explicit sources)
      end subroutine eval_sources
    end interface
    class(ccs_vector), intent(inout) :: R   !< Work vector (for evaluating linear/implicit sources)
    class(ccs_vector), intent(inout) :: S   !< Work vector (for evaluating fixed/explicit sources)
    class(ccs_matrix), intent(inout) :: M   !< The equation system/matrix
    class(ccs_vector), intent(inout) :: rhs !< The equation system/RHS vector
    
    call eval_sources(flow, phi, R, S)
    ! TODO: Insert model-specific sources here
    ! @note Model-specific source subroutines should be additive to avoid overwriting user-defined
    !       sources.
    call add_fixed_source(S, rhs)
    call add_linear_source(R, M)
    
  end subroutine apply_sources
  
end submodule scalars_common
