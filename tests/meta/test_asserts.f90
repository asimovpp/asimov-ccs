!> @brief Test the asserts from the testing lib. 
program test_asserts

  use testing_lib

  implicit none

  integer(ccs_int), dimension(3) :: vec_a
  integer(ccs_int), dimension(3) :: vec_b
  integer(ccs_int), dimension(3) :: vec_c
  
  real(ccs_real), dimension(3) :: rvec_a
  real(ccs_real), dimension(3) :: rvec_b
  real(ccs_real), dimension(3) :: rvec_c

  character(len=1024) :: outmsg
  integer :: i

  integer, parameter :: n_tests = 14

  logical, dimension(n_tests) :: res
  
  call init()

  vec_a = (/ 11, 22, 33 /)
  vec_b = (/ 11, 22, 33 /)
  vec_c = (/ 11, 22, 44 /)
  
  rvec_a = (/ 11.1_ccs_real, 22.2_ccs_real, 33.3_ccs_real /)
  rvec_b = (/ 11.1_ccs_real, 22.2_ccs_real, 33.3_ccs_real /)
  rvec_c = (/ 11.1_ccs_real, 22.2_ccs_real, 44.4_ccs_real /)
  
  res(:) = .false.

  call assert_eq(vec_a(1), vec_b(1), "vec(1)s are not equal", res(1))
  call assert_eq(vec_a, vec_b, "vecs are not equal", res(2))
  call assert_eq(rvec_a(1), rvec_b(1), "rvec(1)s are not equal", res(3))
  call assert_eq(rvec_a, rvec_b, "rvecs are not equal", res(4))
  call assert_eq("cat", "cat", "strings are not equal", res(5))

  call assert_lt(vec_a(2), 23, "int value too large", res(6))
  call assert_lt(rvec_a(2), 23.0_ccs_real, "real value too large", res(7))
  
  call assert_gt(vec_a(3), 23, "int value too small", res(8))
  call assert_gt(rvec_a(3), 23.0_ccs_real, "real value too small", res(9))
  
  call assert_bool(1 /= 2 .and. 3 > 2, "incorrect bool eval", res(10))
  call assert_bool(vec_a == vec_b, "incorrect bool array eval", res(11))
  
  call assert_neq(vec_a(1), 12, "ints should not be equal", res(12))
  call assert_neq(rvec_a(1), 12.0_ccs_real, "reals should not be equal", res(13))
  call assert_neq("cat", "dog", "strings should not be equal", res(14))

  if (.not. all(res)) then
    outmsg = ""
    do i = 1, n_tests 
      outmsg = trim(outmsg) // str(i) // " " // str(res(i)) // " " // new_line('a')
    end do
    
    call stop_test("Some tests failed." // new_line('a') // outmsg)
  end if
  
  res(:) = .true.

  call assert_eq(vec_a(1), vec_b(2), "vec(1)s are not equal", res(1))
  call assert_eq(vec_a, vec_c, "vecs are not equal", res(2))
  call assert_eq(rvec_a(1), rvec_b(2), "rvec(1)s are not equal", res(3))
  call assert_eq(rvec_a, rvec_c, "rvecs are not equal", res(4))
  call assert_eq("cat", "dog", "strings are not equal", res(5))
  
  call assert_lt(vec_a(3), 23, "int value too large", res(6))
  call assert_lt(rvec_a(3), 23.0_ccs_real, "real value too large", res(7))
  
  call assert_gt(vec_a(2), 23, "int value too small", res(8))
  call assert_gt(rvec_a(2), 23.0_ccs_real, "real value too small", res(9))
  
  call assert_bool(1 == 2 .and. 3 < 2, "incorrect bool eval", res(10))
  call assert_bool(vec_a /= vec_b, "incorrect bool array eval", res(11))
  
  call assert_neq(vec_a(1), 11, "ints should not be equal", res(12))
  call assert_neq(rvec_a(1), 11.1_ccs_real, "reals should not be equal", res(13))
  call assert_neq("cat", "cat", "strings should not be equal", res(14))
  
  if (.not. all(.not. res)) then
    outmsg = ""
    do i = 1, n_tests 
      outmsg = trim(outmsg) // str(i) // " " // str(res(i)) // " " // new_line('a')
    end do

    call stop_test("Some tests did not fail when they should have." // new_line('a') // outmsg) 
  end if

  call fin()

end program test_asserts
