
! --------------------------------------------------------------------
! MODULE  Sorting:
!    This module can sort a set of numbers.  The method used is
! usually referred to as "selection" method.
! --------------------------------------------------------------------

MODULE  mod_gamma_sorting
   IMPLICIT  NONE
   PRIVATE   :: FindMinimum, Swap
    
    CONTAINS
    
    ! --------------------------------------------------------------------
    ! INTEGER FUNCTION  FindMinimum():
    !    This function returns the location of the minimum in the section
    ! between Start and End.
    ! --------------------------------------------------------------------
    INTEGER FUNCTION  FindMinimum(x, Start, Ends)
        IMPLICIT  NONE
        real, DIMENSION(1:), INTENT(IN) :: x
        INTEGER, INTENT(IN)                :: Start, Ends
        real                            :: Minimum
        INTEGER                            :: Location
        INTEGER                            :: i

        Minimum  = x(Start)          ! assume the first is the min
        Location = Start             ! record its position
        DO i = Start+1, Ends          ! start with next elements
            IF (x(i) < Minimum) THEN  !   if x(i) less than the min?
                Minimum  = x(i)        !      Yes, a new minimum found
                Location = i                !      record its position
            END IF
        END DO
        FindMinimum = Location            ! return the position
    END FUNCTION  FindMinimum

    ! --------------------------------------------------------------------
    ! SUBROUTINE  Swap():
    !    This subroutine swaps the values of its two formal arguments.
    ! --------------------------------------------------------------------
    SUBROUTINE  Swap(a, b)
        IMPLICIT  NONE
        real, INTENT(INOUT)    :: a, b
        real  :: Temp

        Temp = a
        a    = b
        b    = Temp
    END SUBROUTINE  Swap

    ! --------------------------------------------------------------------
    ! SUBROUTINE  Sort():
    !    This subroutine receives an array x() and sorts it into ascending
    ! order.
    ! --------------------------------------------------------------------
    SUBROUTINE  Sort(x)
        IMPLICIT  NONE
        REAL, DIMENSION(1:), INTENT(INOUT)    :: x
        INTEGER                               :: i
        INTEGER                               :: Location
        INTEGER                               :: n

        n=size(x)
        DO i = 1, n-1             ! except for the last
            Location = FindMinimum(x, i, n)  ! find min from this to last
            CALL  Swap(x(i), x(Location))  ! swap this and the minimum
        END DO
    END SUBROUTINE  Sort

END MODULE  mod_gamma_sorting


