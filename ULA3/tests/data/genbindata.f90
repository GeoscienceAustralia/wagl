program genbindata
    integer*1 A(10,4)
    integer*2 B(10,4)
    integer*4 C(10,4)
    real*4 D(10,4)

    do i=1,10
        do j=1,4
            A(i,j) = (i-1)*4 + j
            B(i,j) = (i-1)*4 + j
            C(i,j) = (i-1)*4 + j
            D(i,j) = (i-1)*4.0 + j
        enddo
    enddo

    open (1, file='int8.bin', form='unformatted', access='direct', recl=4)
    open (2, file='int16.bin', form='unformatted', access='direct', recl=2*4)
    open (3, file='int32.bin', form='unformatted', access='direct', recl=4*4)
    open (4, file='float32.bin', form='unformatted', access='direct', recl=4*4)

    do i=1,10
        write(1, rec=i) A(i,:)
        write(2, rec=i) B(i,:)
        write(3, rec=i) C(i,:)
        write(4, rec=i) D(i,:)
    end do

    close(1)
    close(2)
    close(3)
    close(4)
end
