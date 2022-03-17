submodule (meshing) meshing_operators

  implicit none

contains

  module procedure get_neighbour_distance
    real(accs_real), dimension(ndim) :: xp
    real(accs_real), dimension(ndim) :: xnb

    call get_centre(loc_p, xp)
    call get_centre(loc_nb, xnb)

    dx(:) = xnb(:) - xp(:)
    
  end procedure get_neighbour_distance

  module procedure get_face_distance
    real(accs_real), dimension(ndim) :: xp
    real(accs_real), dimension(ndim) :: xf

    call get_centre(loc_p, xp)
    call get_centre(loc_f, xf)

    dx(:) = xf(:) - xp(:)
    
  end procedure get_face_distance
    
end submodule meshing_operators
