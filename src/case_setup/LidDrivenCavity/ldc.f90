!> @brief Program file for LidDrivenCavity case
!
!> @build mpi+petsc

program ldc

  use yaml, only: parse, error_length
  use kinds, only : ccs_int, ccs_real
  use constants, only: ccsconfig
  use bc_constants, only: bc_region_left, bc_region_right, &
                          bc_region_top, bc_region_bottom, bc_region_live, &
                          bc_type_sym, bc_type_wall
  use parallel, only: initialise_parallel_environment, &
                      cleanup_parallel_environment, &
                      read_command_line_arguments, &
                      timer
  use parallel_types, only: parallel_environment

  implicit none

  character(len=:), allocatable :: case_name   !> Case name
  character(len=:), allocatable :: ccs_config_file  !> Config file for CCS

  class(parallel_environment), allocatable :: par_env

  ! Reference numbers
  real(ccs_real) :: p_ref
  real(ccs_real) :: temp_ref
  real(ccs_real) :: density
  real(ccs_real) :: viscosity
  integer(ccs_int) :: pref_at_cell

  ! Number of steps
  integer(ccs_int) :: num_steps

  ! Convection/discretisation scheme
  integer(ccs_int) :: u_conv
  integer(ccs_int) :: v_conv

  ! Relation factors
  real(ccs_real) :: u_relax
  real(ccs_real) :: v_relax
  real(ccs_real) :: p_relax

  ! Boundary conditions - hardcoded for now
  integer(ccs_int), dimension(5) :: bnd_region
  integer(ccs_int), dimension(5) :: bnd_type
  
  integer(ccs_int) :: irank !> MPI rank ID
  ! integer(ccs_int) :: isize !> Size of MPI world

  double precision :: start_time, end_time

  ! Launch MPI
  call initialise_parallel_environment(par_env) 

  irank = par_env%proc_id
!  isize = par_env%num_procs

  call read_command_line_arguments(par_env, case_name=case_name)

  ccs_config_file = case_name//ccsconfig

  call timer(start_time)

  ! Read case name from configuration file
  call read_configuration(ccs_config_file)

  ! BC regions and types - for region "top", u = 1.0
  bnd_region = (/ bc_region_left, bc_region_right, bc_region_top, bc_region_bottom, bc_region_live /)
  bnd_type = (/ bc_type_wall, bc_type_wall, bc_type_wall, bc_type_wall, bc_type_sym /)

  if(irank == 0) then
    call print_config()
  end if

  call timer(end_time)

  if(irank == 0) then
     print*, "Elapsed time: ", end_time - start_time
  end if
  
  ! Finalise MPI
  call cleanup_parallel_environment(par_env)

  contains

  ! Read YAML configuration file
  subroutine read_configuration(config_filename)

    use read_config, only: get_reference_number, get_steps, &
                            get_convection_scheme, get_relaxation_factor

    character(len=*), intent(in) :: config_filename
    
    class(*), pointer :: config_file_pointer  !> Pointer to CCS config file
    character(len=error_length) :: error

    config_file_pointer => parse(config_filename, error=error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif
    
    call get_steps(config_file_pointer, num_steps)
    call get_reference_number(config_file_pointer, p_ref=p_ref, temp_ref=temp_ref, &
                              dens_ref=density, visc_ref=viscosity, pref_at_cell=pref_at_cell)
    call get_convection_scheme(config_file_pointer, u_conv=u_conv, v_conv=v_conv)
    call get_relaxation_factor(config_file_pointer, u_relax=u_relax, v_relax=v_relax, p_relax=p_relax)

  end subroutine

  ! Print test case configuration
  subroutine print_config()

    print*,"Solving ", case_name, " case"

    print*,"++++" 
    print*,"SIMULATION LENGTH"
    print*,"Running for ",num_steps, "time steps"

    print*,"++++"   
    print*,"REFERENCE NUMBERS"
    print*,"Reference pressure ", p_ref, " set a cell ", pref_at_cell
    print*,"Reference temperature ", temp_ref
    print*,"Density ", density
    print*,"Viscocity ", viscosity

    print*,"++++" 
    print*,"DISCRETISATION SCHEMES"  
    if(u_conv == 0) then
      print*,"u is upwind" 
    else
      print*,"u is central" 
    end if
  
    if(v_conv == 0) then
      print*,"v is upwind" 
    else
      print*,"v is central" 
    end if

    print*,"++++" 
    print*,"RELAXATION FACTORS"
    print*,"u is ", u_relax 
    print*,"v is ", v_relax 
    print*,"p is ", p_relax 

    print*,"++++" 
    print*,"BOUNDARY CONDITIONS"
    print*,"Regions: ", bnd_region
    print*,"Types: ", bnd_type

  end subroutine

end program ldc
