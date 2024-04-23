submodule(core) core_init_flow

  use kinds, only: ccs_int, ccs_real

  use fields, only: get_field_is_face_based
  use kinds, only: ccs_real, ccs_int
  use types, only: fluid, field, cell_locator, face_locator, neighbour_locator, vector_values

  use utils, only: update, get_field, set_mode, set_row, set_entry, set_values

  use meshing, only: create_cell_locator, get_global_index, count_neighbours, create_neighbour_locator, &
                       get_local_index, create_face_locator, get_face_normal, get_centre, get_face_area, &
                       get_face_interpolation

  use vec, only: get_vector_data, restore_vector_data, create_vector_values
  use fv, only: compute_boundary_values

  use constants, only: insert_mode, ndim



  implicit none


  contains

  !> Initialise both cell centre values and mass fluxes
  module subroutine initialise_flow(flow_fields, get_init_flow, get_init_mass_flux)
    type(fluid), intent(inout) :: flow_fields
    interface 
      pure subroutine get_init_flow(loc_p, field_name, init_val)
        use kinds, only: ccs_real
        use types, only: cell_locator
        type(cell_locator), intent(in) :: loc_p
        character(len=*), intent(in) :: field_name
        real(ccs_real), intent(inout) :: init_val
      end subroutine

      pure subroutine get_init_mass_flux(loc_f, init_val)
        use kinds, only: ccs_real
        use types, only: face_locator
        type(face_locator), intent(in) :: loc_f
        real(ccs_real), intent(inout) :: init_val
      end subroutine
    end interface


    call initialise_cells(flow_fields, get_init_flow)
    call initialise_mass_flux(flow_fields, get_init_mass_flux)

  end subroutine


  module subroutine initialise_cells(flow_fields, get_init_flow)

    ! Arguments
    type(fluid), intent(inout) :: flow_fields

    interface 
      pure subroutine get_init_flow(loc_p, field_name, init_val)
        use kinds, only: ccs_real
        use types, only: cell_locator
        type(cell_locator), intent(in) :: loc_p
        character(len=*), intent(in) :: field_name
        real(ccs_real), intent(inout) :: init_val
      end subroutine
    end interface

    ! Local variables
    class(field), pointer :: current_field
    integer(ccs_int) :: n_local, i_field
    integer(ccs_int) :: index_p, global_index_p
    real(ccs_real) :: init_val
    type(cell_locator) :: loc_p
    type(vector_values):: values
    logical :: is_face_field

    call get_local_num_cells(n_local)

    do i_field=1, size(flow_fields%fields)

      call get_field(flow_fields, i_field, current_field)

      call get_field_is_face_based(current_field, is_face_field)
      ! Skip face centered fields
      if (is_face_field) then
        continue
      end if

      call create_vector_values(n_local, values)
      call set_mode(insert_mode, values)

      do index_p = 1, n_local
        call create_cell_locator(index_p, loc_p)
        call get_global_index(loc_p, global_index_p)

        init_val = 0.0_ccs_real
        ! Call case specific function
        call get_init_flow(loc_p, current_field%name, init_val)

        call set_row(global_index_p, values)
        call set_entry(init_val, values)
      end do

      call set_values(values, current_field%values)

      call update(current_field%values)

      nullify(current_field)

    end do

  end subroutine

  module subroutine initialise_mass_flux(flow_fields, get_init_mass_flux)

    type(fluid), intent(inout) :: flow_fields

    interface
      pure subroutine get_init_mass_flux(loc_f, init_val)
        use kinds, only: ccs_real
        use types, only: face_locator
        type(face_locator), intent(in) :: loc_f
        real(ccs_real), intent(inout) :: init_val
      end subroutine
    end interface

    real(ccs_real) :: area
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
        call get_face_area(loc_f, area)

        if (is_boundary) then
          ! Get boundary values
          call compute_boundary_values(mf, 1, loc_p, loc_f, face_normal, velocity(1))
          call compute_boundary_values(mf, 2, loc_p, loc_f, face_normal, velocity(2))
          call compute_boundary_values(mf, 3, loc_p, loc_f, face_normal, velocity(3))

          mf_data(index_f) = area * dot_product(velocity, face_normal)
          call get_init_mass_flux(loc_f, mf_data(index_f))

        else if (index_p < index_nb) then  
          call get_centre(loc_f, x_f)
          call get_face_interpolation(loc_f, interpol_factor)

          u_p = [u%values_ro(index_p), v%values_ro(index_p), w%values_ro(index_p)]
          u_nb = [u%values_ro(index_nb), v%values_ro(index_nb), w%values_ro(index_nb)]
          velocity = interpol_factor * u_p + (1-interpol_factor) * u_nb

          ! compute initial value based on current face coordinates
          mf_data(index_f) = area * dot_product(velocity, face_normal)

          ! Allow case to overwrite interpolated value
          call get_init_mass_flux(loc_f, mf_data(index_f))
        end if

      end do
    end do

    call restore_vector_data(mf%values, mf_data)
    call update(mf%values)

  end subroutine



end submodule