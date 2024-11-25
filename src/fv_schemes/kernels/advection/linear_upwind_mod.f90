module advection_luw_mod
  use advection_mod
  use types
  use vec, only: get_vector_data, restore_vector_data
  use meshing, only: get_face_interpolation, get_local_index, get_distance, get_centre

  implicit none

  type, extends(advection_kernel) :: luw_kernel
  contains
    procedure eval_coeffs => advect_luw_coeffs
    procedure eval_explicit => advect_luw_eval
    procedure get_width => get_luw_width
    procedure get_order => get_luw_order
  end type luw_kernel

  interface
    !> Calculates advection coefficient for neighbouring cell using Linear Upwind discretisation
    module pure function advect_luw_coeffs(self, phi, loc_f, mf, bc, loc_p, loc_nb) result(coeffs)
      class(advection_kernel), intent(in) :: self
      type(linear_upwind_field), intent(inout) :: phi     !< scalar field
      type(face_locator), intent(in) :: loc_f             !< face locator
      real(ccs_real), intent(in) :: mf                    !< mass flux at the face
      integer(ccs_int), intent(in) :: bc                  !< flag indicating whether cell is on boundary
      type(cell_locator), intent(in) :: loc_p             !< current cell locator
      type(neighbour_locator), intent(in) :: loc_nb       !< neighbour cell locator
      real(ccs_real), dimension(2) :: coeffs              !< advection coefficient for current and neighbour cells
    end function advect_luw_coeffs

    module subroutine advect_luw_eval(self, result)
      class(advection_kernel), intent(in) :: self
      real(ccs_real), intent(out) :: result
    end subroutine advect_luw_eval

    module pure function get_luw_width(self) result(width)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: width
    end function get_luw_width

    module pure function get_luw_order(self) result(order)
      class(advection_kernel), intent(in) :: self
      integer(ccs_int) :: order
    end function get_luw_order
  end interface
end module advection_luw_mod
