SUBROUTINE cstart_cend(nrow, ncol, viewmax, view, line, ncentre, &
                       istart, iend)

       integer nrow, ncol
       real viewmax 
       real, dimension(nrow, ncol), intent(in) :: view
       integer*8, dimension(ncol), intent(in) :: line, ncentre
       integer*8, dimension(ncol), intent(inout) :: istart, iend

!f2py depend(nrow, ncol), view
!f2py depend(ncol), line, ncentre, istart, iend

!      old doco
!      program coordinator
!      program to find the line and column for 9 coordinators
!      program determines 9 coordinators from centre line
!      and view angle
!      written by Fuqin Li in 2009 and update in 2010
!      converted to a python module via f2py

!      new doco
!      find the start and end column indices for each row of data
!      takes into account whether the pixel is left or right of
!      `ncentre` (satellite track location in pixel units for the array),
!      as well as the satellite view angle

!      nrow here in fortran are ncol in python
!      ncol here in fortran are nrow in python
!      the input and input-output arrays are transposed
!      in order to align with fortran contiguous arrays


        do i=1, ncol
          kk1 = 0
          kk2 = 0
          do j=1, nrow
            if (j .lt. ncentre(i)) then
              if (view(j, i) .le. viewmax) then
                if (kk1 .eq. 0) then
                  istart(i) = j
                  kk1 = kk1 + 1
                endif
              endif
            endif
            if (j .gt. ncentre(i)) then
              if (view(j, i) .gt. viewmax) then
                if (kk2 .eq. 0) then
                  iend(i) = j - 1
                  kk2 = kk2 + 1
                endif
              endif
            endif
          enddo
          if (view(1, i) .lt. viewmax) then
            istart(i) = 1
          endif
          if (view(nrow, i) .le. viewmax) then
            iend(i) = nrow
          endif
        enddo
END SUBROUTINE cstart_cend
