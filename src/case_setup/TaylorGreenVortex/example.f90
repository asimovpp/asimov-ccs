
program example


  use, intrinsic :: iso_fortran_env, only:  output_unit
  use yaml, only: parse, error_length
  use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
  

  implicit none
  
  class(type_node), pointer :: root
  character(len=error_length) :: error
  
  class(type_dictionary), pointer :: dict
  class(type_list), pointer :: list, list2
  class(type_list_item), pointer :: item, inner_item
  type(type_error), pointer :: io_err
  
  character(len=:), allocatable :: string

  character(len=:), allocatable :: title

  ! Number of steps
  integer :: steps

  ! Initialisation
  character(len=5) :: init = "     "

  ! Reference numbers
  integer :: pressure
  integer :: temperature
  integer :: density
  integer :: pref_at_cell
  real(real_kind)    :: viscosity

  ! Solve
  character(len=3) :: u_sol = "   "
  character(len=3) :: v_sol = "   "
  character(len=3) :: w_sol = "   "
  character(len=3) :: p_sol = "   "

  ! Solvers
  integer :: u_solver
  integer :: v_solver
  integer :: w_solver
  integer :: p_solver

  ! Unsteady solution
  character(len=5) :: transient_type = "     "  
  real(real_kind) :: dt
  real(real_kind) :: gamma
  integer :: max_sub_steps

  ! Target residual
  real(real_kind) :: residual

  ! Monitor cell
  integer :: monitor_cell

  ! Convection schemes
  integer :: u_conv
  integer :: v_conv
  integer :: w_conv

  ! Blending factors
  real(real_kind) :: u_blend
  real(real_kind) :: v_blend
  real(real_kind) :: p_blend

  ! Output frequency & iteration
  integer :: output_freq
  integer :: output_iter

  ! Plot format
  character(len=:), allocatable :: plot_format

  ! Post type & vars
  character(len=:), allocatable :: post_type
  character(len=2), dimension(10) :: post_vars = "  "

  ! Boundardies
  character(len=16), dimension(:), allocatable :: bnd_region
  character(len=16), dimension(:), allocatable :: bnd_type
  real(real_kind), dimension(:,:), allocatable :: bnd_vector

  root => parse("./tgv_config.yaml", error = error)
  if (error/='') then
    print*,trim(error)
    stop 1
  endif
  
  select type (root)
  class is (type_dictionary)

    ! Get title
    call get_case_name(root, title)
    
    ! Get reference numbers
    call get_reference_numbers(root, pressure, temperature, density, viscosity, pref_at_cell)

    ! Get steps
    call get_steps(root, steps)

    ! Variables to solve
    call get_solve(root, u_sol, v_sol, w_sol, p_sol)

    ! Solvers
    call get_solvers(root, u_solver, v_solver, w_solver, p_solver)  
    
    ! Get initilisation
    call get_init(root, init)

    ! Get unsteady solution parameters
    call get_transient(root, transient_type, dt, gamma, max_sub_steps)

    ! Get target resdiual
    call get_target_residual(root, residual)

    ! Get monitor cell ID
    call get_monitor_cell(root, monitor_cell)

    ! Get convection schemes
    call get_convection_scheme(root, u_conv, v_conv, w_conv)

    ! Get blending factors
    call get_blending_factor(root, u_blend, v_blend, p_blend)

    ! Get output frequency and iteration
    call get_output_frequency(root, output_freq, output_freq)

    ! Get plot format
    call get_plot_format(root, plot_format)

    ! Get output type and variables
    call get_output_type(root, post_type, post_vars)

    ! Get boundaries
    call get_boundaries(root, bnd_region, bnd_type, bnd_vector)

  end select

  call root%finalize()
  
  deallocate(root)


contains

  subroutine error_handler(io_err)
    type (type_error), pointer :: io_err
    if (associated(io_err)) then
      print*,trim(io_err%message)
      !stop 1
    endif
  end subroutine

  subroutine get_case_name(root, title)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: title

    title = trim(root%get_string('title',error=io_err))
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'title=',title

  end subroutine

  subroutine get_steps(root, steps)
    class(type_dictionary), pointer, intent(in) :: root
    integer :: steps

    steps = root%get_integer('steps',error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'steps=',steps

  end subroutine

  subroutine get_init(root, init)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=5), intent(inout) :: init

    init = trim(root%get_string('init',error=io_err))
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'init=',init

  end subroutine

  subroutine get_reference_numbers(root, p, temp, density, viscosity, pref)
    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: p
    integer, intent(inout) :: temp
    integer, intent(inout) :: density
    integer, intent(inout) :: pref
    real(real_kind), intent(inout) :: viscosity

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('reference_numbers',required=.true.,error=io_err)

    print*,"Reference numbers: "
    
    ! Pressure
    pressure = dict%get_integer("pressure",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"pressure= ",pressure
    
    ! Temperature
    temperature = dict%get_integer("temperature",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"temperature= ",temperature

    ! Density
    density = dict%get_integer("density",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"density= ",density
    
    ! Pref_at_cell
    pref_at_cell = dict%get_integer("pref_at_cell",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"pref_at_cell= ",pref_at_cell

    ! Viscosity
    viscosity = dict%get_real("viscosity",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"viscosity= ",viscosity
  
  end subroutine

  subroutine get_solve(root, u_sol, v_sol, w_sol, p_sol)

    class(type_dictionary), pointer, intent(in) :: root
    character(len=3), intent(inout) :: u_sol
    character(len=3), intent(inout) :: v_sol
    character(len=3), intent(inout) :: w_sol
    character(len=3), intent(inout) :: p_sol

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('solve',required=.true.,error=io_err)

    print*,"Solve (on/off): "

    ! Solve u?
    u_sol = trim(dict%get_string('u',error=io_err))
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,'u=',u_sol
    
    ! Solve v?
    v_sol = trim(dict%get_string('v',error=io_err))
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'v=',v_sol

    ! Solve w?
    w_sol = trim(dict%get_string('w',error=io_err))
    call error_handler(io_err)   
    if(associated(io_err) == .false.) print*,'w=',w_sol
    
    ! Solve p?
    p_sol = trim(dict%get_string('p',error=io_err))
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'p=',p_sol

  end subroutine

  subroutine get_solvers(root, u_solver, v_solver, w_solver, p_solver)

    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: u_solver
    integer, intent(inout) :: v_solver
    integer, intent(inout) :: w_solver
    integer, intent(inout) :: p_solver

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('solver',required=.true.,error=io_err)

    print*,"Solver: "

    ! Get u_solver
    u_solver = dict%get_integer('u',error=io_err)
    call error_handler(io_err)   
    if(associated(io_err) == .false.) print*,'u=',u_solver
    
    ! Get v_solver
    v_solver = dict%get_integer('u',error=io_err)
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'v=',v_solver

    ! Get w_solver
    w_solver = dict%get_integer('w',error=io_err)
    call error_handler(io_err)   
    if(associated(io_err) == .false.) print*,'w=',w_solver
    
    ! Get p_solver
    p_solver = dict%get_integer('p',error=io_err)
    call error_handler(io_err)    
    if(associated(io_err) == .false.) print*,'p=',p_solver

  end subroutine

  subroutine get_transient(root, transient_type, dt, gamma, max_sub_steps)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=5), intent(inout) :: transient_type
    real(real_kind), intent(inout) :: dt
    real(real_kind), intent(inout) :: gamma
    integer, intent(inout) :: max_sub_steps

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('transient',required=.false.,error=io_err)

    print*,"Transient: "
    
    ! Transient type (euler/quad)
    transient_type = dict%get_string("type",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"type= ",transient_type
    
    ! Dt
    dt = dict%get_real("dt",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"dt= ",dt

    ! Gamma
    gamma = dict%get_integer("gamma",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"gamma= ",gamma
    
    ! Maximum number of sub steps
    max_sub_steps = dict%get_integer("max_sub_steps",error = io_err)  
    call error_handler(io_err)
    if(associated(io_err) == .false.) print*,"max_sub_steps= ",max_sub_steps

  end subroutine

  subroutine get_target_residual(root, residual)
    class(type_dictionary), pointer, intent(in) :: root
    real(real_kind), intent(inout) :: residual

    residual = root%get_real('target_residual',error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'target_residual=',residual

  end subroutine

  subroutine get_monitor_cell(root, monitor_cell)
    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: monitor_cell

    monitor_cell = root%get_integer('monitor_cell',error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'monitor_cell=',monitor_cell

  end subroutine

  subroutine get_convection_scheme(root, u_conv, v_conv, w_conv)
    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: u_conv
    integer, intent(inout) :: v_conv
    integer, intent(inout) :: w_conv

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('convection_scheme',required=.false.,error=io_err)

    u_conv = dict%get_integer('u', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'u=',u_conv

    v_conv = dict%get_integer('v', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'v=',v_conv

    w_conv = dict%get_integer('w', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'w=',w_conv

  end subroutine

  subroutine get_blending_factor(root, u_blend, v_blend, p_blend)
    class(type_dictionary), pointer, intent(in) :: root
    real(real_kind), intent(inout) :: u_blend
    real(real_kind), intent(inout) :: v_blend
    real(real_kind), intent(inout) :: p_blend

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('blending_factor',required=.false.,error=io_err)

    u_blend = dict%get_real('u', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'u=',u_blend

    v_blend = dict%get_real('v', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'v=',v_blend

    p_blend = dict%get_real('p', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'p=',p_blend

  end subroutine

  subroutine get_output_frequency(root, output_freq, output_iter)
    class(type_dictionary), pointer, intent(in) :: root
    integer, intent(inout) :: output_freq
    integer, intent(inout) :: output_iter

    class(type_dictionary), pointer :: dict

    dict => root%get_dictionary('output',required=.false.,error=io_err)

    output_freq = dict%get_integer('every', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'every=',output_freq

    output_iter = dict%get_integer('iter', error=io_err)
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'iter=',output_iter

  end subroutine

  subroutine get_plot_format(root, plot_format)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: plot_format

    plot_format = trim(root%get_string('plot_format',error=io_err))
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'plot_formt=',plot_format

  end subroutine

  subroutine get_output_type(root, post_type, post_vars)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=:), allocatable, intent(inout) :: post_type
    character(len=2), dimension(10), intent(inout) :: post_vars    

    class(type_dictionary), pointer :: dict
    integer :: index

    dict => root%get_dictionary('post',required=.false.,error=io_err)

    post_type = trim(dict%get_string('type', error=io_err))
    call error_handler(io_err)  
    if(associated(io_err) == .false.) print*,'type=',post_type

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
  end subroutine

  subroutine get_boundaries(root, bnd_region, bnd_type, bnd_vector)
    class(type_dictionary), pointer, intent(in) :: root
    character(len=16), dimension(:), allocatable, intent(inout) :: bnd_region
    character(len=16), dimension(:), allocatable, intent(inout) :: bnd_type
    real(real_kind), dimension(:,:), allocatable, intent(inout) :: bnd_vector
  
    integer :: num_boundaries = 0
    integer :: index = 1
    integer :: inner_index = 1
    logical :: success

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
            bnd_vector(inner_index,index) = inner_element%to_real(bnd_vector(inner_index,index),success)
            print*,bnd_vector(inner_index,index)
          end select
        end do

      end select

      item => item%next
      index = index + 1

    end do

  end subroutine

end program
