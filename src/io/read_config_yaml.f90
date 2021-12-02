!> @brief Submodule implementation file read_config_utils.mod
!> @build yaml
!>
!> @details Module implementing the interface to read YAML config file
submodule (read_config) read_config_utils

  use yaml_types, only: type_dictionary, &
                        type_error, &
                        type_node, &
                        type_list, &
                        type_list_item, &
                        type_scalar


  implicit none

  interface get_value
    module procedure get_integer_value
    module procedure get_real_value
    module procedure get_string_value
  end interface

  contains 

  subroutine get_integer_value(dict, keyword, int_val)
    class (type_dictionary), pointer, intent(in) :: dict
    character (len=*), intent(in) :: keyword
    integer, intent(out)  :: int_val

    type(type_error), pointer :: io_err

    int_val = dict%get_integer(keyword,error=io_err)
    call error_handler(io_err)  

    if(associated(io_err) .eqv. .false.) then 
      print*,keyword," = ",int_val
    end if

  end subroutine
    
  subroutine get_real_value(dict, keyword, real_val)
    class (type_dictionary), pointer, intent(in) :: dict
    character (len=*), intent(in) :: keyword
    real(accs_real), intent(out)  :: real_val

    type(type_error), pointer :: io_err

    real_val = dict%get_real(keyword,error=io_err)
    call error_handler(io_err)  

    if(associated(io_err) .eqv. .false.) then 
      print*,keyword," = ",real_val
    end if

  end subroutine

  subroutine get_string_value(dict, keyword, string_val)
    class (type_dictionary), pointer, intent(in) :: dict
    character (len=*), intent(in) :: keyword
    character (len=:), allocatable, intent(inout) :: string_val

    type(type_error), pointer :: io_err

    string_val = trim(dict%get_string(keyword,error=io_err))
    call error_handler(io_err)  

    if(associated(io_err) .eqv. .false.) then 
      print*,keyword," = ",string_val
    end if

  end subroutine

  subroutine error_handler(io_err)
    type(type_error), pointer, intent(inout) :: io_err

    if (associated(io_err)) then
      print*,trim(io_err%message)
    endif
  
  end subroutine

  !> @brief Get the name of the test case
  !
  !> @details Get the case name for the configuration file and 
  !! store it in a string.
  !
  !> @param[in] root - the entry point to the config file    
  !> @param[in,out] title - the case name string    
  module subroutine get_case_name(root, title)
    class(*), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: title

    select type(root)
    type is(type_dictionary)
    
      call get_value(root, "title", title)

    class default
      print*,"Unknown type"
    end select
  
  end subroutine
    
  !> @brief Get the number of steps
  !
  !> @details Get the maximum number of iterations 
  !! to be preformed in the current run 
  !
  !> @param[in] root - the entry point to the config file    
  !> @param[in,out] steps - the maximum number of iterations    
  module  subroutine get_steps(root, steps)
    class(*), pointer, intent(in) :: root
    integer, intent(inout) :: steps

    select type(root)
    type is(type_dictionary)

      print*,"*Number of steps: "

      call get_value(root, 'steps', steps)

    class default
      print*,"Unknown type"
    end select

  end subroutine
    
  !> @brief Get source of initial values
  !
  !> @details Get the source of the initial values - accepted
  !! values are "user", "field" or "step" 
  !
  !> @param[in] root - the entry point to the config file    
  !> @param[in,out] init - the source of the initial values    
  module  subroutine get_init(root, init)
    class(*), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: init

    select type(root)
    type is(type_dictionary)

      print*,"* Initialisation: "

      call get_value(root, "init", init)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get reference numbers
  !
  !> @details Get the reference numbers, the fluid properties 
  !! and the operating condition 
  !
  !> @param[in] root - the entry point to the config file    
  !> @param[in,out] pressure - reference pressure      
  !> @param[in,out] temperature - reference temperature      
  !> @param[in,out] density - reference density      
  !> @param[in,out] viscosity - laminar viscosity      
  !> @param[in,out] pref_at_cell - cell at which the reference pressure is set      
  module subroutine get_reference_numbers(root, pressure, temperature, density, viscosity, pref_at_cell)
    class(*), pointer, intent(in) :: root
    real(accs_real), intent(inout) :: pressure
    real(accs_real), intent(inout) :: temperature
    real(accs_real), intent(inout) :: density
    integer, intent(inout) :: pref_at_cell
    real(accs_real), intent(inout) :: viscosity

    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      dict => root%get_dictionary('reference_numbers',required=.true.,error=io_err)

      print*,"* Reference numbers: "
      
      ! Pressure
      call get_value(dict, "pressure", pressure)
      
      ! Temperature
      call get_value(dict, "temperature", temperature)

      ! Density
      call get_value(dict, "density", density)
      
      ! Pref_at_cell
      call get_value(dict, "pref_at_cell", pref_at_cell)

      ! Viscosity
      call get_value(dict, "viscosity", viscosity)
      
    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get variables to be solved
  !
  !> @details By default, all variables will be solved. Using this 
  !! "solve" keyword, the user can specifically request that 
  !! certain variables will not be solved by setting in to "off"
  !
  !> @param[in] root - the entry point to the config file    
  !> @param[in,out] u_sol - solve u on/off
  !> @param[in,out] v_sol - solve v on/off
  !> @param[in,out] w_sol - solve w on/off
  !> @param[in,out] p_sol - solve p on/off
  !
  !> @todo extend list of variables 
  module subroutine get_solve(root, u_sol, v_sol, w_sol, p_sol)

    class(*), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: u_sol
    character(len=:), allocatable, intent(inout) :: v_sol
    character(len=:), allocatable, intent(inout) :: w_sol
    character(len=:), allocatable, intent(inout) :: p_sol

    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      dict => root%get_dictionary('solve',required=.true.,error=io_err)

      print*,"* Solve (on/off): "

      ! Solve u?
      call get_value(dict, "u", u_sol)
      
      ! Solve v?
      call get_value(dict, "v", v_sol)

      ! Solve w?
      call get_value(dict, "w", w_sol)
      
      ! Solve p?
      call get_value(dict, "p", p_sol)
      
    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get solvers to be used
  !
  !> @details Get the solvers that are to be used for each of
  !! the variables. Solver types are defined by integer values
  !
  !> @param[in] root - the entry point to the config file    
  !> @param[in,out] u_solver - solver to be used for u
  !> @param[in,out] v_solver - solver to be used for v
  !> @param[in,out] w_solver - solver to be used for w
  !> @param[in,out] p_solver - solver to be used for p
  !
  !> @todo extend list of variables   
  module subroutine get_solvers(root, u_solver, v_solver, w_solver, p_solver)

    class(*), pointer, intent(in) :: root
    integer, intent(inout) :: u_solver
    integer, intent(inout) :: v_solver
    integer, intent(inout) :: w_solver
    integer, intent(inout) :: p_solver

    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      dict => root%get_dictionary('solver',required=.true.,error=io_err)

      print*,"* Solver: "

      ! Get u_solver
      call get_value(dict, "u", u_solver)
      
      ! Get v_solver
      call get_value(dict, "v", v_solver)

      ! Get w_solver
      call get_value(dict, "w", w_solver)
      
      ! Get p_solver
      call get_value(dict, "p", p_solver)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get transient status
  !
  !> @details Enables/disables unsteady solution algorithm
  !
  !> @param[in] root - the entry point to the config file   
  !> @param[in,out] transient_type - "euler" (first order) or "quad" (second order)
  !> @param[in,out] dt - time interval (seconds) between two consecutive time steps
  !> @param[in,out] gamma - euler blending factor which blends quad
  !> @param[in,out] max_sub_step - maximum number of sub-iterations at each time step
  module subroutine get_transient(root, transient_type, dt, gamma, max_sub_steps)
    class(*), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: transient_type
    real(accs_real), intent(inout) :: dt
    real(accs_real), intent(inout) :: gamma
    integer, intent(inout) :: max_sub_steps

    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      dict => root%get_dictionary('transient',required=.false.,error=io_err)

      print*,"* Transient: "
      
      ! Transient type (euler/quad)
      call get_value(dict, "type", transient_type)
      
      ! Dt
      call get_value(dict, "dt", dt)

      ! Gamma
      call get_value(dict, "gamma", gamma)
      
      ! Maximum number of sub steps
      call get_value(dict, "max_sub_steps", max_sub_steps)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get target residual
  !
  !> @details Get the convergence criterion. 
  !! The calculation will stop when the residuals (L2-norm) of the 
  !! governing equations are less than the target residual.
  !
  !> @param[in] root - the entry point to the config file   
  !> @param[in,out] residual - convergence criterion
  module  subroutine get_target_residual(root, residual)
    class(*), pointer, intent(in) :: root
    real(accs_real), intent(inout) :: residual

    select type(root)
    type is(type_dictionary)

      call get_value(root, "target_residual", residual)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get grid cell to monitor
  !
  !> @details Get the grid cell at which to monitor the values
  !! of the flow variables (U,V,W,P,TE,ED and T)
  !
  !> @param[in] root - the entry point to the config file   
  !> @param[in,out] monitor_cell - grid cell ID
  module  subroutine get_monitor_cell(root, monitor_cell)
    class(*), pointer, intent(in) :: root
    integer, intent(inout) :: monitor_cell

    select type(root)
    type is(type_dictionary)

      call get_value(root, "monitor_cell", monitor_cell)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get convection schemes 
  !
  !> @details Get convection schemes to be used for the 
  !! different variables. The convection schemes are defined
  !! by integer values.
  !
  !> @param[in] root - the entry point to the config file   
  !> @param[in,out] u_conv - convection scheme for u
  !> @param[in,out] v_conv - convection scheme for v
  !> @param[in,out] w_conv - convection scheme for w
  module  subroutine get_convection_scheme(root, u_conv, v_conv, w_conv)
    class(*), pointer, intent(in) :: root
    integer, intent(inout) :: u_conv
    integer, intent(inout) :: v_conv
    integer, intent(inout) :: w_conv

    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      print*,"* Convection scheme: "
      
      dict => root%get_dictionary('convection_scheme',required=.false.,error=io_err)

      call get_value(dict, "u", u_conv)
      call get_value(dict, "v", v_conv)
      call get_value(dict, "w", w_conv)
      
    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get blending factor values 
  !
  !> @details Get blending factors
  !
  !> @param[in] root - the entry point to the config file
  !> @param[in,out] u_blend - blending factor for u
  !> @param[in,out] v_blend - blending factor for v
  !> @param[in,out] w_blend - blending factor for w
  module subroutine get_blending_factor(root, u_blend, v_blend, w_blend)
    class(*), pointer, intent(in) :: root
    real(accs_real), intent(inout) :: u_blend
    real(accs_real), intent(inout) :: v_blend
    real(accs_real), intent(inout) :: w_blend

    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      print*,"* Blending factor: "

      dict => root%get_dictionary('blending_factor',required=.false.,error=io_err)

      call get_value(dict, "u", u_blend)
      call get_value(dict, "v", v_blend)
      call get_value(dict, "w", w_blend)

    class default
      print*,"Unknown type"
    end select

  end subroutine

    !> @brief Get relaxation factor values 
    !
    !> @details Get relaxation factors
    !
    !> @param[in] root - the entry point to the config file
    !> @param[in,out] u_relax - relaxation factor for u
    !> @param[in,out] v_relax - relaxation factor for v
    !> @param[in,out] p_relax - relaxation factor for p
  module subroutine get_relaxation_factor(root, u_relax, v_relax, p_relax)
    class(*), pointer, intent(in) :: root
    real(accs_real), intent(inout) :: u_relax
    real(accs_real), intent(inout) :: v_relax
    real(accs_real), intent(inout) :: p_relax


    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      print*,"* Relaxation factor: "

      dict => root%get_dictionary('relaxation_factor',required=.false.,error=io_err)

      call get_value(dict, "u", u_relax)
      call get_value(dict, "v", v_relax)
      call get_value(dict, "p", p_relax)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get output frequency 
  !
  !> @details Get output frequency, set with keywords "every"
  !! "iter" or both.
  !
  !> @param[in] root - the entry point to the config file
  !> @param[inout] output_freq - the frequency of writing output files
  !> @param[inout] output iter - output files are written every output_iter/steps
  module subroutine get_output_frequency(root, output_freq, output_iter)
    class(*), pointer, intent(in) :: root
    integer, intent(inout) :: output_freq
    integer, intent(inout) :: output_iter

    class(type_dictionary), pointer :: dict
    type(type_error), pointer :: io_err

    select type(root)
    type is(type_dictionary)

      print*,"Output frequency: "

      dict => root%get_dictionary('output',required=.false.,error=io_err)

      call get_value(dict, "every", output_freq)
      call get_value(dict, "iter", output_iter)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get output file format 
  !
  !> @param[in] root - the entry point to the config file
  !> @param[inout] plot_format - output format (e.g. vtk)
  module subroutine get_plot_format(root, plot_format)
    class(*), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: plot_format

    select type(root)
    type is(type_dictionary)

      call get_value(root, "plot_format", plot_format)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get output type and variables 
  !
  !> @param[in] root - the entry point to the config file
  !> @param[inout] post_type - values at cell centres or cell vertices?
  !> @param[inout] post_vars - variables to be written out
  module subroutine get_output_type(root, post_type, post_vars)
    class(*), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: post_type
    character(len=2), dimension(10), intent(inout) :: post_vars    

    class(type_dictionary), pointer :: dict
    class(type_list), pointer :: list
    class(type_list_item), pointer :: item
    type(type_error), pointer :: io_err
    integer :: index

    select type(root)
    type is(type_dictionary)

      print*,"* Output type and variables: "

      dict => root%get_dictionary('post',required=.false.,error=io_err)

      call get_value(dict, "type", post_type)

      list => dict%get_list('variables',required=.false.,error=io_err)
      call error_handler(io_err)  

      item => list%first
      index = 1
      do while(associated(item))
        select type (element => item%node)
        class is (type_scalar)
          post_vars(index) = trim(element%string)
          print*,post_vars(index)
          item => item%next
          index = index + 1
        end select
      enddo

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get boundary condntions 
  !
  !> @param[in] root - the entry point to the config file
  !> @param[inout] bnd_region - array of boundary region names
  !> @param[inout] bnd_type - array of boundary types (e.g. periodic, symmetric, ...)
  !> @param[inout] bnd_vector - array of boundary vectors
  module subroutine get_boundaries(root, bnd_region, bnd_type, bnd_vector)

    use yaml_types, only: real_kind

    class(*), pointer, intent(in) :: root
    character(len=16), dimension(:), allocatable, intent(inout) :: bnd_region
    character(len=16), dimension(:), allocatable, intent(inout) :: bnd_type
    real(accs_real), dimension(:,:), allocatable, intent(inout) :: bnd_vector
  
    type(type_error), pointer :: io_err
    integer :: num_boundaries = 0
    integer :: index = 1
    integer :: inner_index = 1
    logical :: success

    class(type_dictionary), pointer :: dict
    class(type_list), pointer :: list
    class(type_list_item), pointer :: item, inner_item

    select type(root)
    type is(type_dictionary)

      dict => root%get_dictionary('boundary', required=.false., error=io_err)
      call error_handler(io_err)  

      list => dict%get_list('region', required=.false.,error=io_err)
      call error_handler(io_err)  

      item => list%first
      do while(associated(item))
        num_boundaries = num_boundaries + 1
        item => item%next
      end do

      allocate(bnd_region(num_boundaries))
      allocate(bnd_type(num_boundaries))
      allocate(bnd_vector(3,num_boundaries))

      print*,"* Boundaries:"

      list => dict%get_list('region', required=.false.,error=io_err)
      call error_handler(io_err)  

      print*,"** Regions:"
      item => list%first
      do while(associated(item))
        select type(element => item%node)
        class is (type_scalar)
          bnd_region(index) = trim(element%string)
          print*,bnd_region(index)
          item => item%next
          index = index + 1
          end select  
      end do

      list => dict%get_list('type', required=.false.,error=io_err)
      call error_handler(io_err)  

      index = 1 
      
      print*,"** Types:"
      item => list%first
      do while(associated(item))
        select type(element => item%node)
        class is (type_scalar)
          bnd_type(index) = trim(element%string)
          print*,bnd_type(index)
          item => item%next
          index = index + 1
          end select  
      end do
    
      list => dict%get_list('vector', required=.false.,error=io_err)
      call error_handler(io_err)  

      index = 1

      print*,"** Vectors:"
      item => list%first

      do while(associated(item))

        select type(inner_list => item%node)
        type is(type_list)

          inner_item => inner_list%first
          inner_index = 1

          do while(associated(inner_item))
            select type(inner_element => inner_item%node)
            class is(type_scalar)
              inner_item => inner_item%next
              inner_index = inner_index + 1
              bnd_vector(inner_index,index) = inner_element%to_real(real(bnd_vector(inner_index,index),real_kind),success)
              print*,bnd_vector(inner_index,index)
            end select
          end do

        end select

        item => item%next
        index = index + 1

      end do

    class default
      print*,"Unknown type"
    end select

  end subroutine
    
end submodule read_config_utils
  