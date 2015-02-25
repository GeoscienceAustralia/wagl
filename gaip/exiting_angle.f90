SUBROUTINE _exiting_angle(nrow, ncol, &
    view, azi, theta, phit, et, azi_et)

! calculates incident angles

    integer, intent(in) :: ncol, nrow
    real, dimension(nrow, ncol), intent(in) :: view, azi, theta, phit
    real, dimension(nrow, ncol), intent(inout) :: et, azi_et

!   internals
    real azi_c
    integer col, row

!f2py depend(ncol), view, azi, theta, phit, et, azi_et
!f2py depend(nrow), view, azi, theta, phit, et, azi_et

!---------------------------------------------------------------------

    do col=1,ncol
        do row=1,nrow
            azi_c = azi(row, col)
            if (azi_c .le. -180.0) azi_c = azi_c + 360.0
            if (azi_c .gt.  180.0) azi_c = azi_c - 360.0
    !----------------------------------------------------------------
            call cal_pole(view(row, col), azi_c, &
                          theta(row, col), phit(row, col), &
                          et(row, col), azi_et(row, col))
            if (azi_et(row, col).le. -180.0) azi_et(row, col) = azi_et(row, col) + 360
            if (azi_et(row, col).gt.  180.0) azi_et(row, col) = azi_et(row, col) - 360
    !----------------------------------------------------------------
        enddo
    enddo
    return
END SUBROUTINE _exiting_angle
