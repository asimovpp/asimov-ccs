module sorts
    !! Module including a variety of sorting implementation and tools
#include "ccs_macros.inc"

use kinds

private

public :: heapsort_int
public :: quicksort_int
public :: search_in_sorted

contains

  pure function search_in_sorted(value, a) result(found)
    !! Search for value in sorted array a (Hermann Bottenbruch binary search)
    integer(ccs_int), intent(in) :: value
    integer(ccs_int), intent(in), dimension(:) :: a
    logical :: found
    integer(ccs_int) :: n, first, last, mid

    n = size(a)
    first = 1
    last = n
    do while (first /= last)
      mid = ceiling((first + last)/2.0_ccs_real, ccs_int)
      if (a(mid) > value) then
        last = mid - 1
      else
        first = mid
      end if
    end do

    if (a(first) == value) then
      found = .true.
    else
      found = .false.
    end if
  end function

  !! Heap sort implementation from rosetta code: https://rosettacode.org/wiki/Sorting_algorithms/Heapsort#Fortran
  pure subroutine heapsort_int(a)
      !! Sort of integer array in place using heap sort
  
     integer(ccs_int), intent(inout) :: a(0:)
     integer(ccs_int) :: start, n, bottom
     integer(ccs_int) :: temp
  
     n = size(a)
     do start = (n - 2) / 2, 0, -1
       call siftdown(a, start, n);
     end do
     
     do bottom = n - 1, 1, -1
       temp = a(0)
       a(0) = a(bottom)
       a(bottom) = temp;
       call siftdown(a, 0, bottom)
     end do
  
  end subroutine heapsort_int

  pure subroutine siftdown(a, start, bottom)
    integer(ccs_int), intent(inout) :: a(0:)
    integer(ccs_int), intent(in) :: start, bottom
    integer(ccs_int) :: child, root
    integer(ccs_int) :: temp
  
    root = start
    do while(root*2 + 1 < bottom)
      child = root * 2 + 1
      
      if (child + 1 < bottom) then
        if (a(child) < a(child+1)) child = child + 1
      end if
      
      if (a(root) < a(child)) then
        temp = a(child)
        a(child) = a (root)
        a(root) = temp
        root = child
      else
        return
      end if  
    end do      
  end subroutine siftdown

  !! Quick sort implementation from rosetta code: https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
  pure recursive subroutine quicksort_int(a)
    integer, parameter  ::  changesize = 64
    integer(ccs_int), dimension(:), contiguous, intent(inout) ::  a
    integer(ccs_int) :: first
    integer(ccs_int) :: i
    integer(ccs_int) :: j
    integer(ccs_int) :: last
    logical :: stay
    integer(ccs_int) :: t
    integer(ccs_int) :: x
  
    first = 1
    last = size(a, 1)
    if ((last - first) < changesize)then
      call insertion_sort_int(a(first:last)) 
      return
    end if
  
    j = shiftr((first + last), 1) + 1
    x = a(j)
    i = first
    j = last
    stay = .true.
    do while ( stay )
      do while ( a(i)<x )
        i = i + 1
      end do
      do while ( x<a(j) )
        j = j - 1
      end do
      if ( j<i ) then
        stay = .false.
      else
        t = a(i)      ! Swap the values
        a(i) = a(j)
        a(j) = t
        i = i + 1     ! Adjust the pointers (PIVOT POINTS)
        j = j - 1
      end if
    end do
  
    if ( first<i - 1 ) then
      call quicksort_int(a(first:i - 1))   ! We still have some left to do on the lower
    end if
    if ( j + 1<last ) then 
      call quicksort_int(a(j + 1:last))     ! We still have some left to do on the upper
    end if
  
  end subroutine quicksort_int

  !! insertion sort implementation from rosetta code: https://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran
  pure subroutine insertion_sort_int(a)
      implicit none
      integer(ccs_int), dimension(:), intent(inout) :: a
      integer(ccs_int) :: x
      integer(ccs_int) :: i, j
      
      do i = 2, size(a)
          x = a(i)
          j = i - 1
          do while (j >= 1)
              if (a(j) <= x) exit
              a(j + 1) = a(j)
              j = j - 1
          end do
          a(j + 1) = x
      end do
  end subroutine

end module sorts