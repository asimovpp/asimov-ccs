submodule(core) core_init_flow

  use ccs_base, only: mesh
  
  use kinds, only: ccs_int, ccs_real

  use kinds, only: ccs_real, ccs_int
  use types, only: fluid, field, cell_locator, face_locator, neighbour_locator, vector_values

  use utils, only: update, set_mode, set_row, set_entry, set_values
  use fields, only: get_field, get_field_is_face_based

  use meshing, only: create_cell_locator, get_global_index, count_neighbours, create_neighbour_locator, &
                       get_local_index, create_face_locator, get_face_normal, get_centre, get_face_area, &
                       get_face_interpolation, get_local_num_cells, get_boundary_status

  use vec, only: get_vector_data, restore_vector_data, create_vector_values
  use fv, only: compute_boundary_values

  use io_visualisation, only: read_solution

  use constants, only: insert_mode, ndim

  implicit none

  contains

  !v Initialise both cell centre values and mass fluxes by calling get_init_flow and get_init_mass_flux
  !  on every cell or face
  module subroutine initialise_flow(par_env, run_options, flow_fields, get_init_flow, get_init_mass_flux)
    class(parallel_environment), intent(in) :: par_env !< Parallel environment
    type(ccs_options), intent(in) :: run_options       !< Runtime configuration
    type(fluid), intent(inout) :: flow_fields          !< The flow
    interface
      !> User-supplied subroutine to set field values at cell centres
      pure subroutine get_init_flow(loc_p, field_name, init_val)
        use kinds, only: ccs_real
        use types, only: cell_locator
        type(cell_locator), intent(in) :: loc_p    !< Cell locator
        character(len=*), intent(in) :: field_name !< Name of field being initialised
        real(ccs_real), intent(inout) :: init_val  !< The initial value to set for the field
      end subroutine

      !> User-supplied subroutine to set the mass flux at face centres
      pure subroutine get_init_mass_flux(loc_f, init_val)
        use kinds, only: ccs_real
        use types, only: face_locator
        type(face_locator), intent(in) :: loc_f   !< Face locator
        real(ccs_real), intent(inout) :: init_val !< The initial value to set for the mass flux
      end subroutine
    end interface

    if (.not. run_options%variables%restart) then
      call initialise_cell_values(flow_fields, get_init_flow)
      call initialise_mass_flux(flow_fields, get_init_mass_flux)
    else
      call read_solution(par_env, run_options%paths%case_path, mesh, flow_fields)
    end if

  end subroutine initialise_flow

  !v Initialise the cell centred values by calling the user-supplied initialisation per cell-centred
  !  field on each cell.
  subroutine initialise_cell_values(flow_fields, get_init_flow)

    ! Arguments
    type(fluid), intent(inout) :: flow_fields !< The flow

    interface 
      !> User-supplied subroutine to set field values at cell centres
      pure subroutine get_init_flow(loc_p, field_name, init_val)
        use kinds, only: ccs_real
        use types, only: cell_locator
        type(cell_locator), intent(in) :: loc_p    !< Cell locator
        character(len=*), intent(in) :: field_name !< Name of field being initialised
        real(ccs_real), intent(inout) :: init_val  !< The initial value to set for the field
      end subroutine get_init_flow
    end interface

    ! Local variables
    class(field), pointer :: current_field
    integer(ccs_int) :: n_local, i_field
    integer(ccs_int) :: index_p
    real(ccs_real) :: init_val
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(:), pointer :: field_data
    logical :: is_face_field

    call get_local_num_cells(n_local)

    do i_field=1, size(flow_fields%fields)

      call get_field(flow_fields, i_field, current_field)

      call get_field_is_face_based(current_field, is_face_field)
      ! Skip face centered fields
      if (is_face_field) then
        continue
      end if

      call get_vector_data(current_field%values, field_data)

      do index_p = 1, n_local
        call create_cell_locator(index_p, loc_p)

        call get_init_default(current_field%name, init_val)
        call get_init_flow(loc_p, current_field%name, init_val)
        field_data(index_p) = init_val
      end do

      call restore_vector_data(current_field%values, field_data)

      call update(current_field%values)

      nullify(current_field)

    end do

  end subroutine initialise_cell_values

  !> Utility subroutine to set "sensible" default values
  pure subroutine get_init_default(field_name, init_val)

    character(len=*), intent(in) :: field_name
    real(ccs_real), intent(out) :: init_val

    select case(field_name)
    case ("density")
      init_val = 1.0_ccs_real ! Zero density will cause solver failure!
    case ("viscosity")
      init_val = 1.0e-2_ccs_real
    case default
      init_val = 0.0_ccs_real
    end select
    
  end subroutine get_init_default
  
  !v Initialise the face centred mass flux values by calling the user-supplied initialisation on
  !  mesh faces.
  subroutine initialise_mass_flux(flow_fields, get_init_mass_flux)

    type(fluid), intent(inout) :: flow_fields !< The flow field structure

    interface
      !> User-supplied subroutine to set the mass flux at face centres
      pure subroutine get_init_mass_flux(loc_f, init_val)
        use kinds, only: ccs_real
        use types, only: face_locator
        type(face_locator), intent(in) :: loc_f   !< Face locator
        real(ccs_real), intent(inout) :: init_val !< The initial value to set for the mass flux
      end subroutine
    end interface

    real(ccs_real), dimension(ndim) :: velocity

    class(field), pointer :: u, v, w
    class(field), pointer :: mf
    
    real(ccs_real), dimension(:), pointer:: mf_data

    type(cell_locator):: loc_p
    type(face_locator):: loc_f
    type(neighbour_locator):: loc_nb
    integer(ccs_int) :: n_local, index_p, index_nb, index_f
    integer(ccs_int) :: j, nnb
    real(ccs_real) :: interpol_factor
    real(ccs_real), dimension(ndim) :: face_normal, u_nb, u_p, x_f
    logical :: is_boundary

    call get_field(flow_fields, "u", u)
    call get_field(flow_fields, "v", v)
    call get_field(flow_fields, "w", w)
    call get_field(flow_fields, "mf", mf)
    call get_vector_data(mf%values, mf_data)

    ! Loop over local cells and faces
    call get_local_num_cells(n_local)
    do index_p = 1, n_local

      call create_cell_locator(index_p, loc_p)
      call count_neighbours(loc_p, nnb)
      do j = 1, nnb

        call create_neighbour_locator(loc_p, j, loc_nb)
        call get_local_index(loc_nb, index_nb)

        call get_boundary_status(loc_nb, is_boundary)
        call create_face_locator(index_p, j, loc_f)
        call get_local_index(loc_f, index_f)
        call get_face_normal(loc_f, face_normal)

        if (is_boundary) then
          ! Get boundary values
          call compute_boundary_values(mf, 1, loc_p, loc_f, face_normal, velocity(1))
          call compute_boundary_values(mf, 2, loc_p, loc_f, face_normal, velocity(2))
          call compute_boundary_values(mf, 3, loc_p, loc_f, face_normal, velocity(3))

          mf_data(index_f) = dot_product(velocity, face_normal)
          call get_init_mass_flux(loc_f, mf_data(index_f))

        else if (index_p < index_nb) then  
          call get_centre(loc_f, x_f)
          call get_face_interpolation(loc_f, interpol_factor)

          u_p = [u%values_ro(index_p), v%values_ro(index_p), w%values_ro(index_p)]
          u_nb = [u%values_ro(index_nb), v%values_ro(index_nb), w%values_ro(index_nb)]
          velocity = interpol_factor * u_p + (1-interpol_factor) * u_nb

          ! compute initial value based on current face coordinates
          mf_data(index_f) = dot_product(velocity, face_normal)

          ! Allow case to overwrite interpolated value
          call get_init_mass_flux(loc_f, mf_data(index_f))
        end if

      end do
    end do

    call restore_vector_data(mf%values, mf_data)
    call update(mf%values)

  end subroutine

end submodule
