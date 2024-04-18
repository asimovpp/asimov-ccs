module flow_init

  use kinds, only: ccs_int
  use case_setup, only: get_init_flow, get_init_mass_flux


  use kinds, only: ccs_real, ccs_int
  use types, only: fluid, field, cell_locator, face_locator, neighbour_locator, vector_values

  use utils, only: update, get_field, set_mode, set_row, set_entry, set_values

  use meshing, only: create_cell_locator, get_global_index, count_neighbours, create_neighbour_locator, &
                       get_local_index, create_face_locator, get_face_normal, get_centre, get_face_area, &
                       get_face_interpolation

  use vec, only: get_vector_data, restore_vector_data, create_vector_values

  use constants, only: insert_mode, ndim



  implicit none

  private

  contains

  subroutine initialise_flow(flow_fields)

    ! Arguments
    type(fluid), intent(inout) :: flow_fields

    ! Local variables
    class(field), pointer :: current_field
    integer(ccs_int) :: n_local, i_field
    integer(ccs_int) :: index_p, global_index_p
    real(ccs_real) :: init_val
    type(cell_locator) :: loc_p
    real(ccs_real), dimension(ndim) :: x_p
    type(vector_values):: values

    call get_local_num_cells(n_local)

    do i_field=1, size(flow_fields%fields)

      call get_field(flow_fields, i_field, current_field)

      call create_vector_values(n_local, values)
      call set_mode(insert_mode, values)

      do index_p = 1, n_local
        call create_cell_locator(index_p, loc_p)
        call get_global_index(loc_p, global_index_p)

        call get_centre(loc_p, x_p)

        init_val = 0.0_ccs_real
        ! Call case specific function
        call get_init_flow(x_p, current_field%name, init_val)

        call set_row(global_index_p, values)
        call set_entry(init_val, values)
      end do

      call set_values(values, current_field%values)

      call update(current_field%values)

      nullify(current_field)

    end do

  end subroutine initialise_flow

  subroutine initialise_mass_flux(flow_fields)

    type(fluid), intent(inout) :: flow_fields
    real(ccs_real) :: area, density
    real(ccs_real), dimension(ndim) :: velocity

    class(field), pointer :: u, v, w
    class(field), pointer :: mf, rho
    
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
    call get_field(flow_fields, "density", rho)
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

        if (index_p > index_nb) then  

          call create_face_locator(index_p, j, loc_f)
          call get_local_index(loc_f, index_f)
          call get_face_normal(loc_f, face_normal)
          call get_centre(loc_f, x_f)
          call get_face_interpolation(loc_f, interpol_factor)

          u_p = [u%values_ro(index_p), v%values_ro(index_p), w%values_ro(index_p)]
          u_nb = [u%values_ro(index_nb), v%values_ro(index_nb), w%values_ro(index_nb)]
          velocity = interpol_factor * u_p + (1-interpol_factor) * u_nb
          density = interpol_factor * rho%values_ro(index_p) + (1-interpol_factor) * rho%values_ro(index_nb)
          call get_face_area(loc_f, area)

          ! compute initial value based on current face coordinates
          mf_data(index_f) = density * area * dot_product(velocity, face_normal)

          ! Allow case to overwrite interpolated value
          call get_init_mass_flux(index_f, mf_data(index_f))

        else if (is_boundary) then
          continue


        end if

      end do
    end do

    call restore_vector_data(mf%values, mf_data)
    call update(mf%values)

  end subroutine initialise_mass_flux



end module