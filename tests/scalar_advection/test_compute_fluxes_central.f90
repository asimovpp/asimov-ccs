!v Test that the advection flux coefficients behave as expected for central schemes.
program test_compute_fluxes
#include "ccs_macros.inc"

  use testing_lib
  use types, only: field, central_field, face_field, matrix_values_spec
  use mesh_utils, only: build_square_mesh
  use fv, only: calc_advection_coeff
  use utils, only: update, initialise, finalise, &
                   set_size, set_values, zero, &
                   set_mode, set_entry, set_col, set_row, clear_entries
  use vec, only: create_vector, get_vector_data, restore_vector_data, set_vector_location, create_vector_values
  use mat, only: create_matrix, set_nnz, create_matrix_values, set_matrix_values_spec_nrows, set_matrix_values_spec_ncols
  use solver, only: axpy, norm
  use constants, only: add_mode, face
  use bc_constants
  use meshing, only: get_global_index, get_local_index, get_boundary_status, &
                     create_cell_locator, create_neighbour_locator, create_face_locator, &
                     count_neighbours, get_local_num_cells
  use boundary_conditions, only: allocate_bc_arrays

  implicit none

  real(ccs_real), parameter :: diffusion_factor = 1.e-2_ccs_real ! XXX: temporarily hard-coded
  type(ccs_mesh) :: mesh
  type(vector_spec) :: vec_properties
  class(field), allocatable :: scalar
  class(field), allocatable :: mf
  integer(ccs_int), parameter :: cps = 5
  integer(ccs_int) :: i
  real(ccs_real), dimension(3) :: mf_values
  integer(ccs_int), parameter :: n_boundaries = 4

  call init()

  mesh = build_square_mesh(par_env, cps, 1.0_ccs_real)

  allocate (central_field :: scalar)
  allocate (face_field :: mf)

  call allocate_bc_arrays(n_boundaries, scalar%bcs)
  call allocate_bc_arrays(n_boundaries, scalar%bcs)
  call allocate_bc_arrays(n_boundaries, scalar%bcs)
  call allocate_bc_arrays(n_boundaries, scalar%bcs)
  do i = 1, n_boundaries
    scalar%bcs%ids(i) = i
  end do
  scalar%bcs%bc_types = bc_type_dirichlet
  scalar%bcs%values = 0.0_ccs_real
  scalar%bcs%values(4) = 1.0_ccs_real

  call initialise(vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, scalar%values)

  call set_vector_location(face, vec_properties)
  call set_size(par_env, mesh, vec_properties)
  call create_vector(vec_properties, mf%values)
  call update(mf%values)

  mf_values = (/-1.0_ccs_real, 0.0_ccs_real, 1.0_ccs_real/)

  do i = 1, size(mf_values)
    call set_mass_flux(mf, mf_values(i))
    call run_compute_fluxes_test(scalar, mf, mesh)
  end do

  deallocate (scalar)
  deallocate (mf)

  call fin()

contains

  !> Sets the mass flux array
  subroutine set_mass_flux(mf, mf_value)
    class(field), intent(inout) :: mf     !< The mass flux
    real(ccs_real) :: mf_value            !< The value to set the mass flux field to
    real(ccs_real), dimension(:), pointer :: mf_data

    call get_vector_data(mf%values, mf_data)
    mf_data(:) = mf_value
    call restore_vector_data(mf%values, mf_data)
    call update(mf%values)
  end subroutine set_mass_flux

  !> Tests that the flux coefficients are in a sensible range
  subroutine run_compute_fluxes_test(scalar, mf, mesh)
    class(field), intent(inout) :: scalar   !< The scalar field structure
    class(field), intent(inout) :: mf       !< The mass flux field
    type(ccs_mesh), intent(in) :: mesh      !< The mesh structure

    real(ccs_real), dimension(:), pointer :: mf_data
    
    integer(ccs_int) :: index_p, index_nb, index_f
    integer(ccs_int) :: j, nnb
    type(cell_locator) :: loc_p
    type(neighbour_locator) :: loc_nb
    type(face_locator) :: loc_f
    logical :: is_boundary

    integer(ccs_int) :: local_num_cells

    real(ccs_real) :: sgn
    real(ccs_real) :: adv_coeff, adv_coeffaP, adv_coeffaF

    select type(scalar)
    type is (central_field)
      call get_local_num_cells(mesh, local_num_cells)

      call get_vector_data(mf%values, mf_data)

      do index_p = 1, local_num_cells
        call create_cell_locator(mesh, index_p, loc_p)
        call count_neighbours(loc_p, nnb)

        do j = 1, nnb
          call create_neighbour_locator(loc_p, j, loc_nb)
          call get_boundary_status(loc_nb, is_boundary)
          call get_local_index(loc_nb, index_nb)
          
          call create_face_locator(mesh, index_p, j, loc_f)
          call get_local_index(loc_f, index_f)

          if (.not. is_boundary) then
            if (index_nb < index_p) then
              sgn = -1.0_ccs_real
            else
              sgn = 1.0_ccs_real
            end if
            call calc_advection_coeff(scalar, loc_f, sgn * mf_data(index_f), 0, adv_coeffaP, adv_coeffaF)

            call assert_ge(adv_coeffaF, 0.0_ccs_real, "Central advection coefficient should be >= 0")
            call assert_le(adv_coeffaF, 1.0_ccs_real, "Central advection coefficient should be <= 1")
          else
            sgn = 1.0_ccs_real
            call calc_advection_coeff(scalar, loc_f, sgn * mf_data(index_f), index_nb, adv_coeffaP, adv_coeffaF)
          end if
        end do
      end do

      call restore_vector_data(mf%values, mf_data)
    class default
      call stop_test("This test is only for centrally-differenced fields!")
    end select

  end subroutine run_compute_fluxes_test

end program test_compute_fluxes
