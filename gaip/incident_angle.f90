SUBROUTINE incident_angle(nrow, ncol, &
    solar, sazi, theta, phit, it, azi_it)

! calculates incident angles

    integer, intent(in) :: ncol, nrow
    real, dimension(nrow, ncol), intent(in) :: solar, sazi, theta, phit
    real, dimension(nrow, ncol), intent(inout) :: it, azi_it

!   internals
    real sazi_c
    integer col, row

!f2py depend(nrow, ncol), solar, sazi, theta, phit, it, azi_it

!---------------------------------------------------------------------

    do col=1,ncol
        do row=1,nrow
            sazi_c = sazi(row, col)
            if (sazi_c .le. -180.0) sazi_c = sazi_c + 360.0
            if (sazi_c .gt.  180.0) sazi_c = sazi_c - 360.0
    !----------------------------------------------------------------
            call cal_pole(solar(row, col), sazi_c, &
                          theta(row, col), phit(row, col), &
                          it(row, col), azi_it(row, col))
            if (azi_it(row, col).le. -180.0) azi_it(row, col) = azi_it(row, col) + 360
            if (azi_it(row, col).gt.  180.0) azi_it(row, col) = azi_it(row, col) - 360
    !----------------------------------------------------------------
        enddo
    enddo
    return
END SUBROUTINE incident_angle
