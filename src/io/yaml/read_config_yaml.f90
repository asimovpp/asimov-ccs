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
    class (*), pointer, intent(in) :: dict
    character (len=*), intent(in) :: keyword
    integer, intent(out)  :: int_val

    type(type_error), pointer :: io_err

    select type(dict)
    type is(type_dictionary)

      int_val = dict%get_integer(keyword,error=io_err)
      call error_handler(io_err)
      
    class default
      print*,"Unknown type"
    end select
   
    if(associated(io_err) .eqv. .true.) then 
      print*,"Error reading ",keyword
    end if

  end subroutine
    
  subroutine get_real_value(dict, keyword, real_val)
    class (*), pointer, intent(in) :: dict
    character (len=*), intent(in) :: keyword
    real(accs_real), intent(out)  :: real_val

    type(type_error), pointer :: io_err

    select type(dict)
    type is(type_dictionary)

      real_val = dict%get_real(keyword,error=io_err)
      call error_handler(io_err)  
      
    class default
      print*,"Unknown type"
    end select

    if(associated(io_err) .eqv. .true.) then 
      print*,"Error reading ",keyword
    end if

  end subroutine

  subroutine get_string_value(dict, keyword, string_val)
    class (*), pointer, intent(in) :: dict
    character (len=*), intent(in) :: keyword
    character (len=:), allocatable, intent(inout) :: string_val

    type(type_error), pointer :: io_err

    select type(dict)
    type is(type_dictionary)

      string_val = trim(dict%get_string(keyword,error=io_err))
      call error_handler(io_err)  

    class default
      print*,"Unknown type"
    end select

    if(associated(io_err) .eqv. .true.) then 
      print*,"Error reading ",keyword
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
  !> @param[in] config_file - the entry point to the config file    
  !> @param[in,out] title - the case name string    
  module subroutine get_case_name(config_file, title)
    class(*), pointer, intent(in) :: config_file
    character(len=:), allocatable, intent(inout) :: title

    call get_value(config_file, "title", title)

  end subroutine
    
  !> @brief Get the number of steps
  !
  !> @details Get the maximum number of iterations 
  !! to be preformed in the current run 
  !
  !> @param[in] config_file - the entry point to the config file    
  !> @param[in,out] steps - the maximum number of iterations    
  module  subroutine get_steps(config_file, steps)
    class(*), pointer, intent(in) :: config_file
    integer, intent(inout) :: steps

    call get_value(config_file, 'steps', steps)

  end subroutine
    
  !> @brief Get source of initial values
  !
  !> @details Get the source of the initial values - accepted
  !! values are "user", "field" or "step" 
  !
  !> @param[in] config_file - the entry point to the config file    
  !> @param[in,out] init - the source of the initial values (user or field)
  !> @param[in,out] u_init - initial value for u
  !> @param[in,out] v_init - initial value for v
  !> @param[in,out] w_init - initial value for w
  !> @param[in,out] te_init - initial value for te
  !> @param[in,out] ed_init - initial value for ed

  module  subroutine get_init(config_file, init, u_init, v_init, w_init, te_init, ed_init)
    class(*), pointer, intent(in) :: config_file
    character(len=:), allocatable, intent(inout) :: init
    integer, optional, intent(inout) :: u_init
    integer, optional, intent(inout) :: v_init
    integer, optional, intent(inout) :: w_init
    integer, optional, intent(inout) :: te_init
    integer, optional, intent(inout) :: ed_init

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('init',required=.true.,error=io_err)

      call get_value(dict, "type", init)

      if(present(u_init)) then
        call get_value(dict, "u", u_init)
      end if

      if(present(v_init)) then
        call get_value(dict, "v", v_init)
      end if

      if(present(w_init)) then
        call get_value(dict, "w", w_init)
      end if

      if(present(te_init)) then
        call get_value(dict, "te", te_init)
      end if

      if(present(ed_init)) then
        call get_value(dict, "ed", ed_init)
      end if

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get reference numbers
  !
  !> @details Get the reference numbers, the fluid properties 
  !! and the operating condition 
  !
  !> @param[in] config_file - the entry point to the config file    
  !> @param[in,out] p_ref - reference pressure 
  !> @param[in,out] p_total - total pressure 
  !> @param[in,out] temp_ref - reference temperature      
  !> @param[in,out] dens_ref - reference density      
  !> @param[in,out] visc_ref - laminar viscosity      
  !> @param[in,out] velo_ref - reference velocity      
  !> @param[in,out] leng_ref - reference length, used to define the Reynolds number of the flow      
  !> @param[in,out] pref_at_cell - cell at which the reference pressure is set      
  module subroutine get_reference_number(config_file, p_ref, p_total, temp_ref, &
                                          dens_ref, visc_ref, velo_ref, len_ref, pref_at_cell)
    class(*), pointer, intent(in) :: config_file
    real(accs_real), optional, intent(inout) :: p_ref
    real(accs_real), optional, intent(inout) :: p_total
    real(accs_real), optional, intent(inout) :: temp_ref
    real(accs_real), optional, intent(inout) :: dens_ref
    real(accs_real), optional, intent(inout) :: visc_ref
    real(accs_real), optional, intent(inout) :: velo_ref
    real(accs_real), optional, intent(inout) :: len_ref      
    integer, optional, intent(inout) :: pref_at_cell

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('reference_numbers',required=.true.,error=io_err)

      ! Pressure
      if(present(p_ref)) then
        call get_value(dict, "pressure", p_ref)
      end if

      ! Pressure_total
      if(present(p_total)) then
        call get_value(dict, "pressure_total", p_total)
      end if

      ! Temperature
      if(present(temp_ref)) then
        call get_value(dict, "temperature", temp_ref)
      end if

      ! Density
      if(present(dens_ref)) then
        call get_value(dict, "density", dens_ref)
      end if

      ! Viscosity
      if(present(visc_ref)) then
        call get_value(dict, "viscosity", visc_ref)
      end if
      
      ! Velocity
      if(present(velo_ref)) then
        call get_value(dict, "velocity", velo_ref)
      end if

      ! Length
      if(present(len_ref)) then
        call get_value(dict, "length", len_ref)
      end if

      ! Pref_at_cell
      if(present(pref_at_cell)) then
        call get_value(dict, "pref_at_cell", pref_at_cell)
      end if

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
  !> @param[in] config_file - the entry point to the config file    
  !> @param[in,out] u_sol - solve u on/off
  !> @param[in,out] v_sol - solve v on/off
  !> @param[in,out] w_sol - solve w on/off
  !> @param[in,out] p_sol - solve p on/off
  !
  !> @todo extend list of variables 
  module subroutine get_solve(config_file, u_sol, v_sol, w_sol, p_sol)

    class(*), pointer, intent(in) :: config_file
    character(len=:), allocatable, optional, intent(inout) :: u_sol
    character(len=:), allocatable, optional, intent(inout) :: v_sol
    character(len=:), allocatable, optional, intent(inout) :: w_sol
    character(len=:), allocatable, optional, intent(inout) :: p_sol

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('solve',required=.true.,error=io_err)

      ! Solve u?
      if(present(u_sol)) then
        call get_value(dict, "u", u_sol)
      end if   
      
      ! Solve v?
      if(present(v_sol)) then
        call get_value(dict, "v", v_sol)
      end if

      ! Solve w?
      if(present(w_sol)) then
        call get_value(dict, "w", w_sol)
      end if 
      
      ! Solve p?
      if(present(p_sol)) then
        call get_value(dict, "p", p_sol)
      end if

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get solvers to be used
  !
  !> @details Get the solvers that are to be used for each of
  !! the variables. Solver types are defined by integer values
  !
  !> @param[in] config_file - the entry point to the config file    
  !> @param[in,out] u_solver - solver to be used for u
  !> @param[in,out] v_solver - solver to be used for v
  !> @param[in,out] w_solver - solver to be used for w
  !> @param[in,out] p_solver - solver to be used for p
  !> @param[in,out] te_solver - solver to be used for te
  !> @param[in,out] ed_solver - solver to be used for ed
  !
  !> @todo extend list of variables   
  module subroutine get_solver(config_file, u_solver, v_solver, w_solver, p_solver, te_solver, ed_solver)

    class(*), pointer, intent(in) :: config_file
    integer, optional, intent(inout) :: u_solver
    integer, optional, intent(inout) :: v_solver
    integer, optional, intent(inout) :: w_solver
    integer, optional, intent(inout) :: p_solver
    integer, optional, intent(inout) :: te_solver
    integer, optional, intent(inout) :: ed_solver

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('solver',required=.true.,error=io_err)

      ! Get u_solver
      if(present(u_solver)) then
        call get_value(dict, "u", u_solver)
      end if
      
      ! Get v_solver
      if(present(v_solver)) then
        call get_value(dict, "v", v_solver)
      end if

      ! Get w_solver
      if(present(w_solver)) then
        call get_value(dict, "w", w_solver)
      end if
      
      ! Get p_solver
      if(present(p_solver)) then
        call get_value(dict, "p", p_solver)
      end if

      ! Get te_solver
      if(present(te_solver)) then
        call get_value(dict, "te", te_solver)
      end if

      ! Get ed_solver
      if(present(ed_solver)) then
        call get_value(dict, "ed", ed_solver)
      end if

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get transient status
  !
  !> @details Enables/disables unsteady solution algorithm
  !
  !> @param[in] config_file - the entry point to the config file   
  !> @param[in,out] transient_type - "euler" (first order) or "quad" (second order)
  !> @param[in,out] dt - time interval (seconds) between two consecutive time steps
  !> @param[in,out] euler_blend - gamma, euler blending factor which blends quad
  !> @param[in,out] max_sub_step - maximum number of sub-iterations at each time step
  module subroutine get_transient(config_file, transient_type, dt, euler_blend, max_sub_steps)
    class(*), pointer, intent(in) :: config_file
    character(len=:), allocatable, intent(inout) :: transient_type
    real(accs_real), intent(inout) :: dt
    real(accs_real), intent(inout) :: euler_blend
    integer, intent(inout) :: max_sub_steps

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('transient',required=.false.,error=io_err)

      ! Transient type (euler/quad)
      call get_value(dict, "type", transient_type)
      
      ! Dt
      call get_value(dict, "dt", dt)

      ! Gamma
      call get_value(dict, "gamma", euler_blend)
      
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
  !> @param[in] config_file - the entry point to the config file   
  !> @param[in,out] residual - convergence criterion
  module  subroutine get_target_residual(config_file, residual)
    class(*), pointer, intent(in) :: config_file
    real(accs_real), intent(inout) :: residual

    call get_value(config_file, "target_residual", residual)

  end subroutine

  !> @brief Get grid cell to monitor
  !
  !> @details Get the grid cell at which to monitor the values
  !! of the flow variables (U,V,W,P,TE,ED and T)
  !
  !> @param[in] config_file - the entry point to the config file   
  !> @param[in,out] monitor_cell - grid cell ID
  module  subroutine get_monitor_cell(config_file, monitor_cell)
    class(*), pointer, intent(in) :: config_file
    integer, intent(inout) :: monitor_cell

    call get_value(config_file, "monitor_cell", monitor_cell)

  end subroutine

  !> @brief Get convection schemes 
  !
  !> @details Get convection schemes to be used for the 
  !! different variables. The convection schemes are defined
  !! by integer values.
  !
  !> @param[in] config_file - the entry point to the config file   
  !> @param[in,out] u_conv - convection scheme for u
  !> @param[in,out] v_conv - convection scheme for v
  !> @param[in,out] w_conv - convection scheme for w
  !> @param[in,out] te_conv - convection scheme for te
  !> @param[in,out] ed_conv - convection scheme for ed
  module  subroutine get_convection_scheme(config_file, u_conv, v_conv, w_conv, te_conv, ed_conv)
    class(*), pointer, intent(in) :: config_file
    integer, optional, intent(inout) :: u_conv
    integer, optional, intent(inout) :: v_conv
    integer, optional, intent(inout) :: w_conv
    integer, optional, intent(inout) :: te_conv
    integer, optional, intent(inout) :: ed_conv

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('convection_scheme',required=.false.,error=io_err)

      if(present(u_conv)) then
        call get_value(dict, "u", u_conv)
      end if
      
      if(present(v_conv)) then
        call get_value(dict, "v", v_conv)
      end if
      
      if(present(w_conv)) then
        call get_value(dict, "w", w_conv)
      end if
      
      if(present(te_conv)) then
        call get_value(dict, "te", te_conv)
      end if
      
      if(present(ed_conv)) then
        call get_value(dict, "ed", ed_conv)
      end if

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get blending factor values 
  !
  !> @details Get blending factors
  !
  !> @param[in] config_file - the entry point to the config file
  !> @param[in,out] u_blend - blending factor for u
  !> @param[in,out] v_blend - blending factor for v
  !> @param[in,out] w_blend - blending factor for w
  !> @param[in,out] te_blend - blending factor for te
  !> @param[in,out] ed_blend - blending factor for ed
  module subroutine get_blending_factor(config_file, u_blend, v_blend, w_blend, te_blend, ed_blend)
    class(*), pointer, intent(in) :: config_file
    real(accs_real), optional, intent(inout) :: u_blend
    real(accs_real), optional, intent(inout) :: v_blend
    real(accs_real), optional, intent(inout) :: w_blend
    real(accs_real), optional, intent(inout) :: te_blend
    real(accs_real), optional, intent(inout) :: ed_blend

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('blending_factor',required=.false.,error=io_err)

      if(present(u_blend)) then
        call get_value(dict, "u", u_blend)
      end if

      if(present(v_blend)) then
        call get_value(dict, "v", v_blend)
      end if 

      if(present(w_blend)) then
        call get_value(dict, "w", w_blend)
      end if

      if(present(te_blend)) then
        call get_value(dict, "te", te_blend)
      end if

      if(present(ed_blend)) then
        call get_value(dict, "ed", ed_blend)
      end if

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get relaxation factor values 
  !
  !> @details Get relaxation factors
  !
  !> @param[in] config_file - the entry point to the config file
  !> @param[in,out] u_relax - relaxation factor for u
  !> @param[in,out] v_relax - relaxation factor for v
  !> @param[in,out] p_relax - relaxation factor for p
  !> @param[in,out] te_relax - relaxation factor for te
  !> @param[in,out] ed_relax - relaxation factor for ed
  module subroutine get_relaxation_factor(config_file, u_relax, v_relax, p_relax, te_relax, ed_relax)
    class(*), pointer, intent(in) :: config_file
    real(accs_real), optional, intent(inout) :: u_relax
    real(accs_real), optional, intent(inout) :: v_relax
    real(accs_real), optional, intent(inout) :: p_relax
    real(accs_real), optional, intent(inout) :: te_relax
    real(accs_real), optional, intent(inout) :: ed_relax

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('relaxation_factor',required=.false.,error=io_err)

      if(present(u_relax)) then
        call get_value(dict, "u", u_relax)
      end if

      if(present(v_relax)) then
        call get_value(dict, "v", v_relax)
      end if

      if(present(p_relax)) then
        call get_value(dict, "p", p_relax)
      end if

      if(present(te_relax)) then
        call get_value(dict, "te", te_relax)
      end if

      if(present(ed_relax)) then
        call get_value(dict, "ed", ed_relax)
      end if

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get output frequency 
  !
  !> @details Get output frequency, set with keywords "every"
  !! "iter" or both.
  !
  !> @param[in] config_file - the entry point to the config file
  !> @param[inout] output_freq - the frequency of writing output files
  !> @param[inout] output iter - output files are written every output_iter/steps
  module subroutine get_output_frequency(config_file, output_freq, output_iter)
    class(*), pointer, intent(in) :: config_file
    integer, intent(inout) :: output_freq
    integer, intent(inout) :: output_iter

    class(*), pointer :: dict
    type(type_error), pointer :: io_err

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('output',required=.false.,error=io_err)

      call get_value(dict, "every", output_freq)
      call get_value(dict, "iter", output_iter)

    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get output file format 
  !
  !> @param[in] config_file - the entry point to the config file
  !> @param[inout] plot_format - output format (e.g. vtk)
  module subroutine get_plot_format(config_file, plot_format)
    class(*), pointer, intent(in) :: config_file
    character(len=:), allocatable, intent(inout) :: plot_format

    call get_value(config_file, "plot_format", plot_format)

  end subroutine

  !> @brief Get output type and variables 
  !
  !> @param[in] config_file - the entry point to the config file
  !> @param[inout] post_type - values at cell centres or cell vertices?
  !> @param[inout] post_vars - variables to be written out
  module subroutine get_output_type(config_file, post_type, post_vars)
    class(*), pointer, intent(in) :: config_file
    character(len=:), allocatable, intent(inout) :: post_type
    character(len=2), dimension(10), intent(inout) :: post_vars    

    class(*), pointer :: dict
    class(type_list), pointer :: list
    class(type_list_item), pointer :: item
    type(type_error), pointer :: io_err
    integer :: idx

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('post',required=.false.,error=io_err)

      call get_value(dict, "type", post_type)

      select type(dict)
      type is(type_dictionary)
         
        list => dict%get_list('variables',required=.false.,error=io_err)
        call error_handler(io_err)  

        item => list%first
        idx = 1
        do while(associated(item))
          select type (element => item%node)
          class is (type_scalar)
            post_vars(idx) = trim(element%string)
            print*,post_vars(idx)
            item => item%next
            idx = idx + 1
           end select
        end do

      class default
        print*,"Unknown type"
      end select
     
    class default
      print*,"Unknown type"
    end select

  end subroutine

  !> @brief Get boundary conditions 
  !
  !> @param[in] config_file - the entry point to the config file
  !> @param[inout] bnd_region - array of boundary region names
  !> @param[inout] bnd_type - array of boundary types (e.g. periodic, symmetric, ...)
  !> @param[inout] bnd_vector - array of boundary vectors
  module subroutine get_boundaries(config_file, bnd_region, bnd_type, bnd_vector)

    use yaml_types, only: real_kind

    class(*), pointer, intent(in) :: config_file
    character(len=16), dimension(:), allocatable, intent(inout) :: bnd_region
    character(len=16), dimension(:), allocatable, intent(inout) :: bnd_type
    real(accs_real), optional, dimension(:,:), allocatable, intent(inout) :: bnd_vector
  
    type(type_error), pointer :: io_err
    integer :: num_boundaries = 0
    integer :: idx = 1
    integer :: inner_idx = 1
    logical :: success

    class(type_dictionary), pointer :: dict
    class(type_list), pointer :: list
    class(type_list_item), pointer :: item, inner_item

    select type(config_file)
    type is(type_dictionary)

      dict => config_file%get_dictionary('boundary', required=.false., error=io_err)
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

      list => dict%get_list('region', required=.false.,error=io_err)
      call error_handler(io_err)  

      item => list%first
      do while(associated(item))
        select type(element => item%node)
        class is (type_scalar)
          bnd_region(idx) = trim(element%string)
          print*,bnd_region(idx)
          item => item%next
          idx = idx + 1
          end select  
      end do

      list => dict%get_list('type', required=.false.,error=io_err)
      call error_handler(io_err)  

      idx = 1 
      
      item => list%first
      do while(associated(item))
        select type(element => item%node)
        class is (type_scalar)
          bnd_type(idx) = trim(element%string)
          print*,bnd_type(idx)
          item => item%next
          idx = idx + 1
          end select  
      end do

      if(present(bnd_vector)) then
    
        list => dict%get_list('vector', required=.false.,error=io_err)
        call error_handler(io_err)  

        idx = 1

        item => list%first

        do while(associated(item))

          select type(inner_list => item%node)
          type is(type_list)

            inner_item => inner_list%first
            inner_idx = 1

            do while(associated(inner_item))
              select type(inner_element => inner_item%node)
              class is(type_scalar)
                inner_item => inner_item%next
                inner_idx = inner_idx + 1
                bnd_vector(inner_idx,idx) = inner_element%to_real(real(bnd_vector(inner_idx,idx),real_kind),success)
                print*,bnd_vector(inner_idx,idx)
              end select
            end do

          end select

          item => item%next
          idx = idx + 1

        end do

      end if

    class default
      print*,"Unknown type"
    end select

  end subroutine
    
end submodule read_config_utils
