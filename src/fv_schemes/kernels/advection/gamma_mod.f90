module advection_gamma_mod
  use advection_mod
  use types
  use vec, only: get_vector_data, restore_vector_data
  use meshing, only: get_face_interpolation, get_local_index, get_distance, get_centre
  use types, only: neighbour_locator

  implicit none

  type, extends(advection_kernel) :: gamma_kernel
  contains
    procedure eval_coeffs => advect_gamma_coeffs
    procedure eval_explicit => advect_gamma_eval
    procedure get_width => get_gamma_width
    procedure get_order => get_gamma_order
  end type gamma_kernel

  interface
    !> Calculates advection coefficient for neighbouring cell using gamma discretisation
    module pure function advect_gamma_coeffs(self, phi, loc_f, mf, bc, loc_p, loc_nb) result(coeffs)
      class(advection_kernel), intent(in) :: self
      type(gamma_field), intent(inout) :: phi       !< scalar field
      type(face_locator), intent(in) :: loc_f       !< face locator
      real(ccs_real), intent(in) :: mf              !< mass flux at the face
      integer(ccs_int), intent(in) :: bc            !< flag indicating whether cell is on boundary
      type(cell_locator), intent(in) :: loc_p       !< current cell locator
      type(neighbour_locator), intent(in) :: loc_nb !< neighbour cell locator
      real(ccs_real), dimension(2) :: coeffs        !< advection coefficient for current and neighbour cells
    end function advect_gamma_coeffs

    module subroutine advect_gamma_eval(self, result)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(out) :: result
    end subroutine advect_gamma_eval

    module pure function get_gamma_width(self) result(width)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_gamma_width

    module pure function get_gamma_order(self) result(order)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_gamma_order
  end interface
end module advection_gamma_mod
