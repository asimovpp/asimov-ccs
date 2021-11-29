!> @brief Module file read_yaml.mod
!>
!> @details Module defining IO interface for ASiMoV-CCS

module read_yaml

  use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
    
  implicit none
  
  private  

  public :: get_case_name
  public :: get_steps
  public :: get_init
  public :: get_reference_numbers
  public :: get_solve
  public :: get_solvers
  public :: get_transient
  public :: get_target_residual
  public :: get_monitor_cell
  public :: get_convection_scheme
  public :: get_blending_factor
  public :: get_output_frequency
  public :: get_plot_format
  public :: get_output_type
  public :: get_boundaries

  interface

    module subroutine error_handler(io_err)
      type (type_error), pointer, intent(inout) :: io_err
    end subroutine
      
    module subroutine get_case_name(root, title)
      class(type_dictionary), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: title
    end subroutine
      
    module  subroutine get_steps(root, steps)
      class(type_dictionary), pointer, intent(in) :: root
      integer, intent(inout) :: steps
    end subroutine
      
    module  subroutine get_init(root, init)
      class(type_dictionary), pointer, intent(in) :: root
      character(len=5), intent(inout) :: init
    end subroutine

    module subroutine get_reference_numbers(root, pressure, temperature, density, viscosity, pref_at_cell)
      class(type_dictionary), pointer, intent(in) :: root
      integer, intent(inout) :: pressure
      integer, intent(inout) :: temperature
      integer, intent(inout) :: density
      integer, intent(inout) :: pref_at_cell
      real(real_kind), intent(inout) :: viscosity
    end subroutine

    module subroutine get_solve(root, u_sol, v_sol, w_sol, p_sol)
      class(type_dictionary), pointer, intent(in) :: root
      character(len=3), intent(inout) :: u_sol
      character(len=3), intent(inout) :: v_sol
      character(len=3), intent(inout) :: w_sol
      character(len=3), intent(inout) :: p_sol
    end subroutine

    module subroutine get_solvers(root, u_solver, v_solver, w_solver, p_solver)
      class(type_dictionary), pointer, intent(in) :: root
      integer, intent(inout) :: u_solver
      integer, intent(inout) :: v_solver
      integer, intent(inout) :: w_solver
      integer, intent(inout) :: p_solver
    end subroutine

    module subroutine get_transient(root, transient_type, dt, gamma, max_sub_steps)
      class(type_dictionary), pointer, intent(in) :: root
      character(len=5), intent(inout) :: transient_type
      real(real_kind), intent(inout) :: dt
      real(real_kind), intent(inout) :: gamma
      integer, intent(inout) :: max_sub_steps
    end subroutine

    module  subroutine get_target_residual(root, residual)
      class(type_dictionary), pointer, intent(in) :: root
      real(real_kind), intent(inout) :: residual
    end subroutine

    module  subroutine get_monitor_cell(root, monitor_cell)
      class(type_dictionary), pointer, intent(in) :: root
      integer, intent(inout) :: monitor_cell
    end subroutine

    module  subroutine get_convection_scheme(root, u_conv, v_conv, w_conv)
      class(type_dictionary), pointer, intent(in) :: root
      integer, intent(inout) :: u_conv
      integer, intent(inout) :: v_conv
      integer, intent(inout) :: w_conv
    end subroutine

    module subroutine get_blending_factor(root, u_blend, v_blend, p_blend)
      class(type_dictionary), pointer, intent(in) :: root
      real(real_kind), intent(inout) :: u_blend
      real(real_kind), intent(inout) :: v_blend
      real(real_kind), intent(inout) :: p_blend
    end subroutine

    module subroutine get_output_frequency(root, output_freq, output_iter)
      class(type_dictionary), pointer, intent(in) :: root
      integer, intent(inout) :: output_freq
      integer, intent(inout) :: output_iter
    end subroutine

    module subroutine get_plot_format(root, plot_format)
      class(type_dictionary), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: plot_format
    end subroutine

    module subroutine get_output_type(root, post_type, post_vars)
      class(type_dictionary), pointer, intent(in) :: root
      character(len=:), allocatable, intent(inout) :: post_type
      character(len=2), dimension(10), intent(inout) :: post_vars    
    end subroutine

    module subroutine get_boundaries(root, bnd_region, bnd_type, bnd_vector)
      class(type_dictionary), pointer, intent(in) :: root
      character(len=16), dimension(:), allocatable, intent(inout) :: bnd_region
      character(len=16), dimension(:), allocatable, intent(inout) :: bnd_type
      real(real_kind), dimension(:,:), allocatable, intent(inout) :: bnd_vector
    end subroutine

  end interface
end module read_yaml
  