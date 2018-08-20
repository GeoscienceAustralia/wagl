SUBROUTINE proj_terrain(n_max, m_max, n, m, z, mask, n_off, m_off, k_max, &
     n_inc, m_inc, h_offset, zmax, htol, dem_nr, dem_nc, nlA_ori, nsA_ori)

    implicit none

!   subroutine to construct mask of shade pixels

!   z(m_max,n_max) is the main array of heights
!   A(m,n) is the (sub-)matrix of target heights in:
!   z(m_off+1,n_off+1) to z(m_off+m,n_off+n)
!   mask(m,n) is the output mask
!   on input assumed to be 1 where valid data exist 0 else
!   k_max is the number of lags in the projection
!   n_inc(k_max) real increments column number for projection
!   m_inc(k_max) real increments row number for the projection
!   h_offset(k_max) is the height of the projection hor each lag
!   zmax is the maximum altitude in the whole array
!   htol is a tolerance (m) for the test for a hit (RMS error in z)

    integer dem_nr, dem_nc, nlA_ori, nsA_ori
    integer n_max, m_max, n, m, n_off, m_off, k_max
!   NOTE: n_inc and m_inc are Floating Point arrays
    real m_inc(k_max), n_inc(k_max)
    real h_offset(k_max)
    real z(dem_nr, dem_nc)
    integer*2 mask(nlA_ori, nsA_ori)
    real zmax, t, tt, test, htol, xd, yd, zpos
    integer i, j, ii, jj, k, ipos, jpos
    integer*4 itot
    real pi,r2d,d2r,eps
    common/base/pi,r2d,d2r,eps

!   loop over object matrix A and project out into buffer


!   The main 2 loops are over the pixels in the submatrix A

!   For given A(i,j) the search for occluding terrain is done along a "line"
!   in the sun direction. The search can stop when the terrain would
!   have to be higher than the maximum value to occlude the current test
!   pixel (i,j) in the sub-matrix A

    itot=0
    do 100 i=1,m
        ii=m_off+i
        do 110 j=1,n
            jj=n_off+j
            t=z(ii,jj)
            do 120 k=1,k_max
                tt=t+h_offset(k)
                if(tt .le. zmax+htol) then
                    ipos = ifix(float(ii)+m_inc(k))
                    jpos = ifix(float(jj)+n_inc(k))
                    yd = float(ii)+m_inc(k)-float(ipos)
                    xd = float(jj)+n_inc(k)-float(jpos)
                    zpos = z(ipos,jpos)*(1.0-yd)*(1.0-xd) + &
                        z(ipos,jpos+1)*xd*(1.0-yd) + &
                        z(ipos+1,jpos+1)*xd*yd + &
                        z(ipos+1,jpos)*(1.0-xd)*yd
                    test = tt-zpos
                    if (test.le.htol) then
                        mask(i,j) = 0
                        itot = itot+1
                        goto 125
                    endif
                else
                    go to 125
                endif
120         continue
125     continue
110 continue
100 continue
    return
END SUBROUTINE proj_terrain
