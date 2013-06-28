! Author: Josh Sixsmith, joshua.sixsmith@ga.gov.au

SUBROUTINE test_bool(array, n, x)
    IMPLICIT NONE
    INTEGER*8, intent(in) :: n
    INTEGER*8, DIMENSION(n), INTENT(IN) :: array
    !f2py depend(n), array
    INTEGER*1, DIMENSION(n), INTENT(INOUT) :: x
    !f2py depend(n), x
    INTEGER*1, DIMENSION(n) :: y
    !f2py depend(n), y

    y = array .gt. 5
    x = y * y
END SUBROUTINE test_bool

SUBROUTINE test_bool2(array, n, x)
    IMPLICIT NONE
    INTEGER*8 :: n
    INTEGER*8, DIMENSION(n), INTENT(IN) :: array
    !f2py depend(n), array
    INTEGER*1, DIMENSION(n), INTENT(INOUT) :: x
    !f2py depend(n), x
    INTEGER*1, DIMENSION(n) :: y
    !f2py depend(n), y

    y = array .gt. 5
    x = y * y
END SUBROUTINE test_bool2
