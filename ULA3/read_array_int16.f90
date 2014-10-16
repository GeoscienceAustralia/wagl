SUBROUTINE read_array_int16(fname, nrow, ncol, out_data)
    character fname*150
    integer nrow, ncol
    integer*2 out_data(nrow, ncol)

!f2py intent(in) fname, nrow, ncol
!f2py intent(out) out_data

    open(2, file=fname, status='old', access='direct', recl=2*ncol)
    do i=1,nrow
        read (2,rec=i)(out_data(i,j),j=1,ncol)
    enddo
END SUBROUTINE read_array_int16
