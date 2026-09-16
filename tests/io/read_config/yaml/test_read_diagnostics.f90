!> Verify optional nested diagnostics configuration behavior.
!> The first argument selects a fixture under Input and the second gives its expected boolean:
!> diagnostics_enabled.yaml enables reporting; diagnostics_disabled.yaml disables it;
!> diagnostics_absent.yaml and diagnostics_missing_key.yaml check that an absent dictionary
!> or absent setting defaults to false. diagnostics_malformed.yaml supplies an invalid boolean
!> and must terminate with a configuration error, checked by the LIT driver.
program test_read_diagnostics

  use fortran_yaml_c_interface, only: parse
  use read_config, only: get_diagnostics
  use testing_lib

  implicit none

  character(len=256) :: config_path
  character(len=16) :: expected_text
  class(*), pointer :: config_file
  character(len=:), allocatable :: error
  logical :: boundary_mass_fluxes
  logical :: expected

  call init()

  call get_command_argument(1, config_path)
  call get_command_argument(2, expected_text)
  expected = trim(expected_text) == "true"

  config_file => parse(trim(config_path), error)
  if (allocated(error)) call stop_test(trim(error))

  call get_diagnostics(config_file, boundary_mass_fluxes)
  if (boundary_mass_fluxes .neqv. expected) then
    call stop_test("Nested boundary mass-flux diagnostic did not match expectation")
  end if

  call fin()

end program test_read_diagnostics
