!v Module file fv.mod
!
!  An interface to finite volume implementations (CDS, UDS, etc.)

module fv

  use kinds, only: ccs_real, ccs_int
  use types, only: ccs_matrix, ccs_vector, ccs_mesh, field, upwind_field, central_field, bc_config, face_locator, cell_locator
  use constants, only: ndim

  implicit none

  private

  public :: compute_fluxes
  public :: calc_advection_coeff
  public :: calc_diffusion_coeff
  public :: calc_mass_flux
  public :: calc_cell_coords
  public :: update_gradient
  public :: compute_boundary_values

  interface calc_advection_coeff
    module procedure calc_advection_coeff_cds
    module procedure calc_advection_coeff_uds
  end interface calc_advection_coeff

  interface calc_mass_flux
    module procedure calc_mass_flux_uv
    module procedure calc_mass_flux_no_uv
  end interface calc_mass_flux

  interface

    !> Calculates advection coefficient for neighbouring cell using CDS discretisation
    module subroutine calc_advection_coeff_cds(phi, mf, bc, coeff)
      type(central_field), intent(in) :: phi !< scalar (central) field
      real(ccs_real), intent(in) :: mf       !< mass flux at the face
      integer(ccs_int), intent(in) :: bc     !< flag indicating whether cell is on boundary
      real(ccs_real), intent(out) :: coeff   !< advection coefficient to be calculated
    end subroutine calc_advection_coeff_cds

    !> Calculates advection coefficient for neighbouring cell using UDS discretisation
    module subroutine calc_advection_coeff_uds(phi, mf, bc, coeff)
      type(upwind_field), intent(in) :: phi !< scalar (upwind) field
      real(ccs_real), intent(in) :: mf      !< mass flux at the face
      integer(ccs_int), intent(in) :: bc    !< flag indicating whether cell is on boundary
      real(ccs_real), intent(out) :: coeff  !< advection coefficient to be calculated
    end subroutine calc_advection_coeff_uds

    !> Sets the diffusion coefficient
    ! XXX: why is this a function when the equivalent advection ones are subroutines?
    module function calc_diffusion_coeff(index_p, index_nb, mesh) result(coeff)
      integer(ccs_int), intent(in) :: index_p  !< the local cell index
      integer(ccs_int), intent(in) :: index_nb !< the local neigbouring cell index
      type(ccs_mesh), intent(in) :: mesh       !< the mesh structure
      real(ccs_real) :: coeff                  !< the diffusion coefficient
    end function calc_diffusion_coeff

    !> Computes fluxes and assign to matrix and RHS
    module subroutine compute_fluxes(phi, mf, mesh, component, M, vec)
      class(field), intent(in) :: phi             !< scalar field structure
      class(field), intent(in) :: mf              !< mass flux field structure (defined at faces)
      type(ccs_mesh), intent(in) :: mesh          !< the mesh being used
      integer(ccs_int), intent(in) :: component   !< integer indicating direction of velocity field component
      class(ccs_matrix), intent(inout) :: M       !< Data structure containing matrix to be filled
      class(ccs_vector), intent(inout) :: vec     !< Data structure containing RHS vector to be filled
    end subroutine

    !> Calculates mass flux across given face. Note: assumes rho = 1 and uniform grid
    module function calc_mass_flux_uv(u_field, v_field, p, dpdx, dpdy, invAu, invAv, loc_f) result(flux)
      class(field), intent(in) :: u_field               !< x velocities field
      class(field), intent(in) :: v_field               !< y velocities field
      real(ccs_real), dimension(:), intent(in) :: p     !< array containing pressure
      real(ccs_real), dimension(:), intent(in) :: dpdx  !< pressure gradients in x
      real(ccs_real), dimension(:), intent(in) :: dpdy  !< pressure gradients in y
      real(ccs_real), dimension(:), intent(in) :: invAu !< inverse momentum diagonal in x
      real(ccs_real), dimension(:), intent(in) :: invAv !< inverse momentum diagonal in y
      type(face_locator), intent(in) :: loc_f           !< face locator
      real(ccs_real) :: flux                            !< the flux across the boundary
    end function calc_mass_flux_uv

    !> Computes Rhie-Chow correction
    module function calc_mass_flux_no_uv(p, dpdx, dpdy, invAu, invAv, loc_f) result(flux)
      real(ccs_real), dimension(:), intent(in) :: p     !< array containing pressure
      real(ccs_real), dimension(:), intent(in) :: dpdx  !< pressure gradients in x
      real(ccs_real), dimension(:), intent(in) :: dpdy  !< pressure gradients in y
      real(ccs_real), dimension(:), intent(in) :: invAu !< inverse momentum diagonal in x
      real(ccs_real), dimension(:), intent(in) :: invAv !< inverse momentum diagonal in y
      type(face_locator), intent(in) :: loc_f           !< face locator
      real(ccs_real) :: flux                            !< the flux across the boundary
    end function calc_mass_flux_no_uv

    !> Calculates the row and column indices from flattened vector index. Assumes square mesh
    module subroutine calc_cell_coords(index, cps, row, col)
      integer(ccs_int), intent(in) :: index !< cell index
      integer(ccs_int), intent(in) :: cps   !< number of cells per side
      integer(ccs_int), intent(out) :: row  !< cell row within mesh
      integer(ccs_int), intent(out) :: col  !< cell column within mesh
    end subroutine calc_cell_coords

    !v Performs an update of the gradients of a field.
    !  @note This will perform a parallel update of the gradient fields to ensure halo cells are
    !  correctly updated on other PEs. @endnote
    module subroutine update_gradient(mesh, phi)
      type(ccs_mesh), intent(in) :: mesh !< the mesh
      class(field), intent(inout) :: phi !< the field whose gradients we want to update
    end subroutine update_gradient

    !> Computes the value of the scalar field on the boundary 
    module subroutine compute_boundary_values(phi, component, loc_p, loc_f, normal, bc_value, &
                                              x_gradients, y_gradients, z_gradients)
      class(field), intent(in) :: phi                         !< the field for which boundary values are being computed
      integer(ccs_int), intent(in) :: component               !< integer indicating direction of velocity field component
      type(cell_locator), intent(in) :: loc_p                 !< location of cell
      type(face_locator), intent(in) :: loc_f                 !< location of face
      real(ccs_real), dimension(ndim), intent(in) :: normal   !< boundary face normal direction
      real(ccs_real), intent(out) :: bc_value                 !< the value of the scalar field at the specified boundary
      real(ccs_real), dimension(:), optional, intent(in) :: x_gradients, y_gradients, z_gradients
    end subroutine

  end interface

end module fv
