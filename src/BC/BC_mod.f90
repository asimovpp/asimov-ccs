!> @brief BC module
!
!> @details Various BC related functionality. Need to expand.

module BC
  use types, only: BC_config
  use kinds, only: accs_int, accs_real
  
  implicit none

  private
  public :: read_BC_config

  contains

  subroutine read_BC_config(filename, BCs) 
    use yaml, only: parse, error_length
    character(len=*), intent(in) :: filename
    type(BC_config), intent(out) :: BCs

    class(*), pointer :: config_file
    character(len=error_length) :: error

    config_file => parse(filename, error=error)
    if (error/='') then
      print*,trim(error)
      stop 1
    endif

    call get_BCs(config_file, BCs)
  end subroutine read_BC_config
  
  subroutine get_BCs(config_file, BCs)
    use BC_constants
    use read_config, only: get_boundaries

    class(*), pointer, intent(in) :: config_file
    type(BC_config), intent(out) :: BCs
    character(len=16), dimension(:), allocatable :: region
    character(len=16), dimension(:), allocatable :: BC_type
    real(accs_real), dimension(:,:), allocatable :: BC_data
    integer(accs_int) :: i

    call get_boundaries(config_file, region, BC_type, BC_data)

    do i = 1, size(region)
      select case(region(i))
        case("left")
          BCs%region(i) = BC_region_left
        case("right")
          BCs%region(i) = BC_region_right
        case("top")
          BCs%region(i) = BC_region_top
        case("bottom")
          BCs%region(i) = BC_region_bottom
        case default
          print *, 'invalid BC region selected'
          stop 
      end select

      select case(BC_type(i))
        case("periodic")
          BCs%BC_type(i) = BC_type_periodic
        case("sym")
          BCs%BC_type(i) = BC_type_sym
        case("dirichlet")
          BCs%BC_type(i) = BC_type_dirichlet
        case("const_grad")
          BCs%BC_type(i) = BC_type_const_grad
        case default
          print *, 'invalid BC type selected'
          stop 
      end select
    end do
    BCs%endpoints(:,:) = BC_data(2:,:2)
  end subroutine get_BCs
end module BC
