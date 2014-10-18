SUBROUTINE filter_float(array, output, nrow, ncol)
    IMPLICIT NONE
!   Perform a 3x3 filtering over a two dimentional array.

    integer :: i, j, k, l
    integer :: nrow, ncol
    real*4 w(3, 3)
    real*4 array(nrow, ncol)
    real*4 output(nrow, ncol)

!f2py intent(in) array
!f2py intent(out) output
!f2py integer intent(hide),depend(array) :: nrow=shape(array,0), ncol=shape(array,1)

    data ((w(i,j),i=1,3),j=1,3)/0.009511,0.078501,0.009511,0.078501,0.647954,0.078501,0.009511,0.078501,0.009511/
!   first line and first column
    output(1,1) = &
        array(1,1)*w(1,1) + array(1,1)*w(1,2) + array(1,2)*w(1,3) + &
        array(1,1)*w(2,1) + array(1,1)*w(2,2) + array(1,2)*w(2,3) + &
        array(2,1)*w(3,1) + array(2,1)*w(3,2) + array(2,2)*w(3,3)
!   first line and last column
    output(1,ncol) = &
        array(1,ncol-1)*w(1,1) + array(1,ncol)*w(1,2) + array(1,ncol)*w(1,3) + &
        array(1,ncol-1)*w(2,1) + array(1,ncol)*w(2,2) + array(1,ncol)*w(2,3) + &
        array(2,ncol-1)*w(3,1) + array(2,ncol)*w(3,2) + array(2,ncol)*w(3,3)
!   first line for other columns
    do j=2,ncol-1
        output(1,j) = &
            array(1,j-1)*w(1,1) + array(1,j)*w(1,2) + array(1,j+1)*w(1,3) + &
            array(1,j-1)*w(2,1) + array(1,j)*w(2,2) + array(1,j+1)*w(2,3) + &
            array(2,j-1)*w(3,1) + array(2,j)*w(3,2) + array(2,j+1)*w(3,3)
    enddo
!   last line and first column
    output(nrow,1) = &
        array(nrow-1,1)*w(1,1) + array(nrow-1,1)*w(1,2) + array(nrow-1,2)*w(1,3) + &
        array(nrow  ,1)*w(2,1) + array(nrow  ,1)*w(2,2) + array(nrow  ,2)*w(2,3) + &
        array(nrow  ,1)*w(3,1) + array(nrow  ,1)*w(3,2) + array(nrow  ,2)*w(3,3)
!   last line and last column
    output(nrow,ncol) = &
        array(nrow-1,ncol-1)*w(1,1) + array(nrow-1,ncol)*w(1,2) + array(nrow-1,ncol)*w(1,3) + &
        array(nrow  ,ncol-1)*w(2,1) + array(nrow  ,ncol)*w(2,2) + array(nrow  ,ncol)*w(2,3) + &
        array(nrow  ,ncol-1)*w(3,1) + array(nrow  ,ncol)*w(3,2) + array(nrow,ncol)*w(3,3)
!   last line for other column
    do j=2,ncol-1
        output(nrow,j) = &
            array(nrow-1,j-1)*w(1,1) + array(nrow-1,j)*w(1,2) + array(nrow-1,j+1)*w(1,3) + &
            array(nrow  ,j-1)*w(2,1) + array(nrow  ,j)*w(2,2) + array(nrow  ,j+1)*w(2,3) + &
            array(nrow,  j-1)*w(3,1) + array(nrow  ,j)*w(3,2) + array(nrow  ,j+1)*w(3,3)
    enddo
    do i=2,nrow-1
!       for first column
        output(i,1) = &
            array(i-1,1)*w(1,1) + array(i-1,1)*w(1,2) + array(i-1,2)*w(1,3) + &
            array(i  ,1)*w(2,1) + array(i  ,1)*w(2,2) + array(i  ,2)*w(2,3) + &
            array(i+1,1)*w(3,1) + array(i+1,1)*w(3,2) + array(i+1,2)*w(3,3)
!       for last column
        output(i,ncol) = &
            array(i-1,ncol-1)*w(1,1) + array(i-1,ncol)*w(1,2) + array(i-1,ncol)*w(1,3) + &
            array(i  ,ncol-1)*w(2,1) + array(i  ,ncol)*w(2,2) + array(i  ,ncol)*w(2,3) + &
            array(i+1,ncol-1)*w(3,1) + array(i+1,ncol)*w(3,2) + array(i+1,ncol)*w(3,3)
!       for rest of the image
        do j=2,ncol-1
            output(i,j)=0
            do k=1,3
                do l=1,3
                    output(i,j) = output(i,j) + w(k,l)*array(i+k-2,j+l-2)
                enddo
            enddo
        enddo
    enddo
END SUBROUTINE filter_float
