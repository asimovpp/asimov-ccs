submodule(meshing) meshing_operators

  implicit none

contains

  !> Returns the distance between two cell centres
  pure module subroutine get_neighbour_distance(loc_p, loc_nb, dx)
    type(cell_locator), intent(in) :: loc_p            !< The cell distance is measured from.
    type(neighbour_locator), intent(in) :: loc_nb      !< The cell distance is measured to.
    real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance

    real(ccs_real), dimension(ndim) :: xp
    real(ccs_real), dimension(ndim) :: xnb

    call get_centre(loc_p, xp)
    call get_centre(loc_nb, xnb)

    dx(:) = xnb(:) - xp(:)

  end subroutine get_neighbour_distance

  !> Returns the distance from cell to face centres
  pure module subroutine get_face_distance(loc_p, loc_f, dx)
    type(cell_locator), intent(in) :: loc_p            !< The cell distance is measured from.
    type(face_locator), intent(in) :: loc_f            !< The face distance is measured to.
    real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance

    real(ccs_real), dimension(ndim) :: xp
    real(ccs_real), dimension(ndim) :: xf

    call get_centre(loc_p, xp)
    call get_centre(loc_f, xf)

    dx(:) = xf(:) - xp(:)

  end subroutine get_face_distance

  !> Returns the distance from cell to face centres
  pure module subroutine get_face_neighbour_distance(loc_nb, loc_f, dx)
    type(neighbour_locator), intent(in) :: loc_nb      !< The cell distance is measured from.
    type(face_locator), intent(in) :: loc_f            !< The face distance is measured to.
    real(ccs_real), dimension(ndim), intent(out) :: dx !< ndim-array of the distance

    real(ccs_real), dimension(ndim) :: xnb
    real(ccs_real), dimension(ndim) :: xf

    call get_centre(loc_nb, xnb)
    call get_centre(loc_f, xf)

    dx(:) = xf(:) - xnb(:)

  end subroutine get_face_neighbour_distance

end submodule meshing_operators
