!v Submodule file fv_discretisation.smod
!
!  Implementations of the finite volume method using the various discretisation schemes scheme
!
!  @build discretisation

submodule(fv) fv_discretisation

use vec, only: get_vector_data, restore_vector_data
use meshing, only: get_local_index
use types, only: neighbour_locator
use meshing, only: get_distance, get_centre

  implicit none

contains
  !> Calculates advection coefficient for neighbouring cell using CDS discretisation
  module subroutine calc_advection_coeff_cds(phi, mf, bc, coeff)
    type(central_field), intent(in) :: phi !< scalar field
    real(ccs_real), intent(in) :: mf       !< mass flux at the face
    integer(ccs_int), intent(in) :: bc     !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeff   !< advection coefficient to be calculated

    real(ccs_real) :: interpolation_factor

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi)
    end associate
    associate (mflux => mf)
    end associate

    if (bc == 0) then
      interpolation_factor = 0.5_ccs_real
    else
      interpolation_factor = 0.5_ccs_real !1.0_ccs_real
    end if
    coeff = interpolation_factor
  end subroutine calc_advection_coeff_cds

  !> Calculates advection coefficient for neighbouring cell using UDS discretisation
  module subroutine calc_advection_coeff_uds(phi, mf, bc, coeff)
    type(upwind_field), intent(in) :: phi !< scalar field
    real(ccs_real), intent(in) :: mf      !< mass flux at the face
    integer(ccs_int), intent(in) :: bc    !< flag indicating whether cell is on boundary
    real(ccs_real), intent(out) :: coeff  !< advection coefficient to be calculated

    ! Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc)
    end associate

    if (mf < 0.0) then
      coeff = 1.0_ccs_real
    else
      coeff = 0.0_ccs_real
    end if

  end subroutine calc_advection_coeff_uds

  !> Calculates advection coefficient for neighbouring cell using gamma discretisation
  module subroutine calc_advection_coeff_gamma(phi, mf, bc, loc_p, loc_nb, face_area, coeff)
    type(gamma_field), intent(inout) :: phi       !< scalar field
    real(ccs_real), intent(in) :: mf              !< mass flux at the face
    integer(ccs_int), intent(in) :: bc            !< flag indicating whether cell is on boundary
    type(cell_locator), intent(in) :: loc_p       !< current cell locator
    type(neighbour_locator), intent(in) :: loc_nb !< neighbour cell locator
    real(ccs_real), intent(out) :: coeff          !< advection coefficient to be calculated
    real(ccs_real) :: face_area                   !< area of the face

    real(ccs_real),dimension(:),pointer:: phi_data 
    real(ccs_real),dimension(:),pointer:: dphidx,dphidy,dphidz
    real(ccs_real),dimension(3) :: dphiF, dphiP, d
    real(ccs_real)::phiF, phiP, dphi, ddphi, phiPt, phiCDS, gamma_m, beta_m
    
    integer(ccs_int) :: index_p, index_nb

    !store values of phi filed in phi_data array
    call get_vector_data(phi%values,phi_data)
    
    !store x-gradients of phi in dphidx array
    call get_vector_data(phi%x_gradients,dphidx)

    !store y-gradients of phi in dphidx array
    call get_vector_data(phi%y_gradients,dphidy)

    !store z-gradients of phi in dphidx array
    call get_vector_data(phi%z_gradients,dphidz)

    !get the local index of current cell and neighbouring cell
    call get_local_index(loc_p, index_p)
    call get_local_index(loc_nb, index_nb)

    !Dummy usage to prevent unused argument.
    associate (scalar => phi, foo => bc)
    end associate

    beta_m=0.5 !value can be varied between 0.1 and 0.5

    if (mf < 0.0) then
      !print*,"gamma mf<0"
      phiP=phi_data(index_nb)
      phiF=phi_data(index_p)

      !Gradient of phi at cell center (current cell)
      dphiP(1)=dphidx(index_nb)
      dphiP(2)=dphidy(index_nb)
      dphiP(3)=dphidz(index_nb)

      !Gradient of phi at cell center (neighbouring cell)
      dphiF(1)=dphidx(index_p)
      dphiF(2)=dphidy(index_p)
      dphiF(3)=dphidy(index_p)

      !Gradient phi at cell face
      dphi=phiF-phiP

      !Get the distance between present and neighbouring cell centers and store it in d
      call get_distance(loc_p, loc_nb, d)
      
      !calculate the normalized phi
      ddphi=2.0*dot_product(dphiP,d) 
      phiPt=1.0-(dphi/ddphi)

      if (phiPt<=0.0 .or. phiPt>=1.0) then !UD
        print*,"gamma<0:UD"
        coeff=1.0_ccs_real
      else if (phiPt>=beta_m.and.phiPt<1.0) then !CDS
        print*,"gamma<0:CDS"
        coeff=0.5_ccs_real 
      else if (phiPt>0.0.and.phiPt<beta_m) then !Gamma
        print*,"gamma<0:CDS+UD"
        gamma_m=phiPt/beta_m
        phiCDS=0.5*(phiP+phiF)
        coeff=((gamma_m*phiCDS)+((1-gamma_m)*phiP))/(mf*face_area)
        !print*,"gamma <0, coeff=",coeff
      end if 

    else
      !print*,"gamma mf>=0"
      phiP=phi_data(index_p)
      phiF=phi_data(index_nb)

      !Gradient of phi at cell center (current cell)
      dphiP(1)=dphidx(index_p)
      dphiP(2)=dphidy(index_p)
      dphiP(3)=dphidz(index_p)

      !Gradient of phi at cell center (neighbouring cell)
      dphiF(1)=dphidx(index_nb)
      dphiF(2)=dphidy(index_nb)
      dphiF(3)=dphidy(index_nb)

      !Gradient phi at cell face
      dphi=phiF-phiP

      !Get the distance between present and neighbouring cell centers and store it in d
      call get_distance(loc_p, loc_nb, d)

      !calculate the normalized phi
      ddphi=2.0*dot_product(dphiP,d)
      phiPt=1.0-(dphi/ddphi)

      if (phiPt<=0.0 .or. phiPt>=1.0) then !UD
        print*,"gamma>=0:UD"
        coeff=0.0_ccs_real
      else if (phiPt>=beta_m.and.phiPt<1.0) then !CDS
        print*,"gamma>=0:CDS"
        coeff=0.5_ccs_real 
      else if (phiPt>0.0.and.phiPt<beta_m) then !Gamma
        print*,"gamma>=0:CDS+UD"
        gamma_m=phiPt/beta_m
        phiCDS=0.5*(phiP+phiF)
        coeff=((gamma_m*phiCDS)+((1-gamma_m)*phiP))/(mf*face_area)
        !print*,"gamma>=0, coeff=",coeff
      end if 
    end if

    ! Restore vectors
    !print*,"Restoring vector data"
    call restore_vector_data(phi%values,phi_data)
    call restore_vector_data(phi%x_gradients,dphidx)
    call restore_vector_data(phi%y_gradients,dphidy)
    call restore_vector_data(phi%z_gradients,dphidz)

    !print*,"stopping code1"
    !stop
    !print*,"stopping code2"

  end subroutine calc_advection_coeff_gamma

end submodule fv_discretisation
