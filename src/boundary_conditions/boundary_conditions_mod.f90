!v boundary conditions module
!
!  Various BC related functionality. Need to expand.

module boundary_conditions
#include "ccs_macros.inc"

  use utils, only: exit_print, debug_print, str
  use types, only: bc_config, field, bc_profile, fluid
  use kinds, only: ccs_int, ccs_real
  use fortran_yaml_c_interface, only: parse
  use read_config, only: get_bc_field
  use bc_constants
  use core, only: ccs_options
  use constants, only: ccs_string_len
  use error_codes

  implicit none

  private
  public :: read_bc_config
  public :: allocate_bc_arrays
  public :: get_bc_index
  public :: set_bc_real_value
  public :: set_bc_type
  public :: set_bc_id
  public :: set_bc_profile
  public :: translate_bcs

contains

  !> Reads config file and assigns data to BC structure
  subroutine read_bc_config(filename, bc_field, phi)
    character(len=*), intent(in) :: filename !< name of the config file
    character(len=*), intent(in) :: bc_field !< string denoting which field we want to read in
    class(field), intent(inout) :: phi       !< the bc struct of the corresponding field

    class(*), pointer :: config_file
    character(:), allocatable :: error

    config_file => parse(filename, error)
    if (allocated(error)) then
      call error_abort(trim(error))
    end if

    call dprint("reading bc config " // bc_field)
    call get_bc_field(config_file, "name", phi, required=.false.)
    call get_bc_field(config_file, "type", phi, required=.false.)
    call get_bc_field(config_file, "value", phi, required=.false.)
    call get_bc_field(config_file, bc_field, phi, required=.false.)

  end subroutine read_bc_config

  !v Loop through fields and translate their potentially complex boundary conditions 
  ! into basic boundary conditions
  subroutine translate_bcs(run_options, flow_fields)
    use ccs_base, only: left, right, bottom, top, back, front
    use fields, only: get_field
    use meshing, only: get_mesh_generated

    type(ccs_options), intent(in) :: run_options !< Object containing runtime options
    type(fluid), intent(inout) :: flow_fields    !< Object containing all the fields to process

    class(field), pointer :: phi
    character(len=ccs_string_len) :: bnd_name
    integer(ccs_int) :: bc_id, i, n_boundaries
    logical :: is_mom_normal, is_generated, x_mom, y_mom, z_mom

    call get_mesh_generated(is_generated)
    is_mom_normal = .false.

    do i=1, size(flow_fields%fields)
      call get_field(flow_fields, i, phi)
      n_boundaries = size(run_options%mesh%bnd_names)

      do bc_id=1, n_boundaries
        bnd_name = trim(run_options%mesh%bnd_names(bc_id))

        if (is_generated) then
          x_mom = (phi%name == 'u')
          y_mom = (phi%name == 'v')
          z_mom = (phi%name == 'w')

          is_mom_normal = (x_mom .and. ((bnd_name == left) .or. (bnd_name == right))) .or. \
                          (y_mom .and. ((bnd_name == bottom) .or. (bnd_name == top))) .or. \
                          (z_mom .and. ((bnd_name == front) .or. (bnd_name == back)))
        else 
          is_mom_normal = .false.
        end if

        call translate_bc(is_generated, is_mom_normal, bc_id, phi)

      end do
    end do

  end subroutine


  !v Translate a boundary contition set by the user into a base bc.
  subroutine translate_bc(is_generated, is_mom_normal, bc_id, phi)

    logical, intent(in) :: is_generated      !< Flag telling if the mesh has been generated (or read from file)
    logical, intent(in) :: is_mom_normal     !< Flag telling if the bc to be processed is normal to the momentum field. Unused, if phi isn't a velocity field
    integer(ccs_int), intent(in) :: bc_id    !< Id of the boundary condition to process
    class(field), intent(inout) :: phi       !< Field to process

    character(len=ccs_string_len) :: field_name
    logical :: is_momentum, is_pressure, is_extra, is_mf
     
    is_momentum = .false.
    is_pressure = .false.
    is_mf = .false.
    is_extra = .false.

    field_name = trim(phi%name)
    if (field_name == "u" .or. (field_name == "v" .or. field_name == "w")) then
      is_momentum = .true.
    elseif (field_name == "p" .or. field_name == "p_prime") then
      is_pressure = .true.
    elseif (field_name == "mf") then
      is_mf = .true.
    else
      is_extra = .true.
    end if

    associate(i => bc_id)
      select case(phi%bcs%bc_types(i))
      case(bc_type_wall)
        if (is_momentum) then
          if (is_mom_normal .and. is_generated) then
            phi%bcs%bc_types(i) = bc_type_neumann
            phi%bcs%values(i) = 0.0_ccs_real
          else
            phi%bcs%bc_types(i) = bc_type_dirichlet
            phi%bcs%values(i) = 0.0_ccs_real
          end if
        else if (is_pressure) then
          phi%bcs%bc_types(i) = bc_type_neumann
          phi%bcs%values(i) = 0.0_ccs_real
        else if (is_mf) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
          phi%bcs%values(i) = 0.0_ccs_real
        else if (is_extra) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
          phi%bcs%values(i) = 0.0_ccs_real
        end if

      case(bc_type_slip_wall)
        if (is_momentum) then
          phi%bcs%bc_types(i) = bc_type_neumann
          phi%bcs%values(i) = 0.0_ccs_real
        else if (is_pressure) then
          phi%bcs%bc_types(i) = bc_type_neumann
          phi%bcs%values(i) = 0.0_ccs_real
        else if (is_mf) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
          phi%bcs%values(i) = 0.0_ccs_real
        else if (is_extra) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
          phi%bcs%values(i) = 0.0_ccs_real
        end if

      case(bc_type_inflow)
        if (is_momentum) then
          if (is_mom_normal) then
            phi%bcs%bc_types(i) = bc_type_dirichlet
          else
            phi%bcs%bc_types(i) = bc_type_dirichlet
            phi%bcs%values(i) = 0.0_ccs_real
          end if
        else if (is_pressure) then
          phi%bcs%bc_types(i) = bc_type_neumann
          phi%bcs%values(i) = 0.0_ccs_real
        else if (is_mf) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
        else if (is_extra) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
          phi%bcs%values(i) = 0.0_ccs_real
        end if

      case(bc_type_outflow)
        if (is_momentum) then
          if (is_mom_normal) then
            phi%bcs%bc_types(i) = bc_type_extrapolate
          else
            phi%bcs%bc_types(i) = bc_type_dirichlet
            phi%bcs%values(i) = 0.0_ccs_real
          end if
        else if (is_pressure) then
          phi%bcs%bc_types(i) = bc_type_neumann
          phi%bcs%values(i) = 0.0_ccs_real
        else if (is_mf) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
        else if (is_extra) then
          phi%bcs%bc_types(i) = bc_type_dirichlet
          phi%bcs%values(i) = 0.0_ccs_real
        end if

      end select
    end associate

  end subroutine

  !> Sets the appropriate integer values for strings with given by the key-value pair attribute, value
  pure subroutine set_bc_type(boundary_index, bc_type, bcs)
    integer(ccs_int), intent(in) :: boundary_index !< Index of the boundary within bcs struct arrays
    character(len=*), intent(in) :: bc_type        !< string giving the bc type
    type(bc_config), intent(inout) :: bcs          !< bcs struct

    select case (bc_type)
    case ("dirichlet")
      bcs%bc_types(boundary_index) = bc_type_dirichlet
    case ("neumann")
      bcs%bc_types(boundary_index) = bc_type_neumann
    case ("extrapolate")
      bcs%bc_types(boundary_index) = bc_type_extrapolate
    case ("profile")
      bcs%bc_types(boundary_index) = bc_type_profile
    case ("inflow")
      bcs%bc_types(boundary_index) = bc_type_inflow
    case ("outflow")
      bcs%bc_types(boundary_index) = bc_type_outflow
    case ("wall")
      bcs%bc_types(boundary_index) = bc_type_wall
    case ("slipwall")
      bcs%bc_types(boundary_index) = bc_type_slip_wall
    !case ("periodic")
    !  bcs%bc_types(boundary_index) = bc_type_periodic
    !case ("sym")
    !  bcs%bc_types(boundary_index) = bc_type_sym
    case default
      error stop invalid_bc_name ! Invalid BC type string received
    end select

  end subroutine set_bc_type
  

  !> Sets the bc struct's id field to the appropriate integer value
  pure subroutine set_bc_id(boundary_index, name, bcs)

    use ccs_base, only: mesh
    use meshing, only: get_bc_id

    integer(ccs_int), intent(in) :: boundary_index !< index of the boundary within the bc struct's arrays
    character(len=*), intent(in) :: name           !< string giving the bc name
    type(bc_config), intent(inout) :: bcs          !< the bcs struct

    integer :: bc_id
    
    call get_bc_id(mesh, name, bc_id)
    bcs%ids(boundary_index) = bc_id
    
  end subroutine set_bc_id

  !> Sets the bc struct's value field to the given real value
  pure subroutine set_bc_real_value(boundary_index, val, bcs)
    integer(ccs_int), intent(in) :: boundary_index !< index of the boundary within the bc struct's arrays
    real(ccs_real), intent(in) :: val              !< the value to set
    type(bc_config), intent(inout) :: bcs          !< the bcs struct

    bcs%values(boundary_index) = val
  end subroutine set_bc_real_value

  !> Allocates arrays of the appropriate size for the name, type and value of the bcs
  pure subroutine allocate_bc_arrays(n_boundaries, bcs)
    integer(ccs_int), intent(in) :: n_boundaries !< the number of boundaries
    type(bc_config), intent(inout) :: bcs        !< the bc struct

    if (.not. allocated(bcs%ids)) then
      allocate (bcs%ids(n_boundaries))
    end if
    if (.not. allocated(bcs%bc_types)) then
      allocate (bcs%bc_types(n_boundaries))
    end if
    if (.not. allocated(bcs%values)) then
      allocate (bcs%values(n_boundaries))
    end if
    if (.not. allocated(bcs%profiles)) then
      allocate (bcs%profiles(n_boundaries))
    end if
  end subroutine allocate_bc_arrays

  !> Gets the index of the given boundary condition within the bc struct arrays
  pure subroutine get_bc_index(phi, index_nb, index_bc)
    class(field), intent(in) :: phi           !< The field whose bc we're getting
    integer(ccs_int), intent(in) :: index_nb  !< The index of the neighbouring boundary cell
    integer(ccs_int), intent(out) :: index_bc !< The index of the appropriate boundary in the bc struct

    ! Local variable
    integer(ccs_int), dimension(1) :: index_tmp ! The intrinsic returns a rank-1 array ...

    index_tmp = findloc(phi%bcs%ids, -index_nb)
    if (index_tmp(1) == 0) then
      error stop bc_index_not_found ! BC index not found
    end if
    
    index_bc = index_tmp(1)
  end subroutine get_bc_index

  !> Set boundary condition profile to the index_bc boundary
  pure subroutine set_bc_profile(phi, profile, index_bc)
    class(field), intent(inout) :: phi           !< The field whose profile we are setting
    type(bc_profile), intent(in) :: profile      !< BC profile
    integer(ccs_int), intent(in) :: index_bc     !< The index of the appropriate boundary in the bc struct

    phi%bcs%profiles(index_bc) = profile

  end subroutine

end module boundary_conditions
