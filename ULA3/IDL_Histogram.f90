MODULE idl_histogram
    IMPLICIT NONE
    ! Author: Josh Sixsmith, joshua.sixsmith@ga.gov.au
CONTAINS

    SUBROUTINE histogram_int(array, hist, a_sz, nbins, min_, max_, max_bin, binsz)
       IMPLICIT NONE

       INTEGER*8 :: i, a_sz, y, ind, nbins
       INTEGER*2, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array
       INTEGER*4, DIMENSION(nbins), INTENT(INOUT) :: hist
       !f2py depend(nbins), hist

       INTEGER :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ! need to check that the value of array(i) is le max
       do i = 1, a_sz
          tf = (array(i) .le. max_) .and. (array(i) .ge. min_) .and. (array(i) .lt. max_bin)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
       enddo


    END SUBROUTINE histogram_int

    SUBROUTINE histogram_long(array, hist, a_sz, nbins, min_, max_, max_bin, binsz)

       IMPLICIT NONE

       INTEGER*8 :: i, a_sz, y, ind, nbins
       INTEGER*4, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins), INTENT(INOUT) :: hist
       !f2py depend(nbins), hist

       INTEGER*4 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ! need to check that the value of array(i) is le max
       do i = 1, a_sz
          tf = (array(i) .le. max_) .and. (array(i) .ge. min_) .and. (array(i) .lt. max_bin)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
       enddo


    END SUBROUTINE histogram_long

    SUBROUTINE histogram_dlong(array, hist, a_sz, nbins, min_, max_, max_bin, binsz)

       IMPLICIT NONE

       INTEGER*8 :: i, a_sz, y, ind, nbins
       INTEGER*8, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins), INTENT(INOUT) :: hist
       !f2py depend(nbins), hist

       INTEGER*8 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ! need to check that the value of array(i) is le max
       do i = 1, a_sz
          tf = (array(i) .le. max_) .and. (array(i) .ge. min_) .and. (array(i) .lt. max_bin)
          y = tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
       enddo


    END SUBROUTINE histogram_dlong

    SUBROUTINE histogram_float(array, hist, a_sz, nbins, min_, max_, max_bin, binsz)

       IMPLICIT NONE

       INTEGER*8 :: i, a_sz, y, ind, nbins
       REAL*4, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins), INTENT(INOUT) :: hist
       !f2py depend(nbins), hist

       REAL*4 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ! need to check that the value of array(i) is le max
       do i = 1, a_sz
          tf = (array(i) .le. max_) .and. (array(i) .ge. min_) .and. (array(i) .lt. max_bin)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
       enddo


    END SUBROUTINE histogram_float

    SUBROUTINE histogram_dfloat(array, hist, a_sz, nbins, min_, max_, max_bin, binsz)

       IMPLICIT NONE

       INTEGER*8 :: i, a_sz, y, ind, nbins
       REAL*8, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins), INTENT(INOUT) :: hist
       !f2py depend(nbins), hist

       REAL*8 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ! need to check that the value of array(i) is le max
       do i = 1, a_sz
          tf = (array(i) .lt. max_bin) .and. (array(i) .ge. min_) .and. (array(i) .le. max_)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
       enddo


    END SUBROUTINE histogram_dfloat

    SUBROUTINE reverse_indices_int(array, hist, ri, nbins, a_sz, ri_sz, min_, max_, max_bin, binsz)
       IMPLICIT NONE

       INTEGER*8 :: i, n, ri_sz, a_sz, y, ind, nbins
       INTEGER*2, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins), INTENT(INOUT):: hist
       !f2py depend(nbins), hist

       INTEGER*4, DIMENSION(ri_sz), INTENT(INOUT) :: ri
       !f2py depend(ri_sz), ri

       INTEGER*4 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ri(2) = nbins
       hist(1) = 0

       !print*, 'compute ivec'

       do n = 2, nbins
          ri(n+1) = ri(n) + hist(n)
       enddo

       hist = 0

       !print*, 'compute ovec'

       do i = 1, a_sz
          tf = (array(i) .lt. max_bin) .and. (array(i) .ge. min_) .and. (array(i) .le. max_)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
          ri(ri(ind) + hist(ind) + 1) = i - 1
          ri(1) = 1
          hist(1) = 0
          ri(2) = nbins
       enddo

    END SUBROUTINE reverse_indices_int

    SUBROUTINE reverse_indices_long(array, hist, ri, nbins, a_sz, ri_sz, min_, max_, max_bin, binsz)
       IMPLICIT NONE

       INTEGER*8 :: i, n, ri_sz, a_sz, y, ind, nbins
       INTEGER*4, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins), INTENT(INOUT) :: hist
       !f2py depend(nbins), hist

       INTEGER*4, DIMENSION(ri_sz), INTENT(INOUT) :: ri
       !f2py depend(ri_sz), ri

       INTEGER*4 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       hist(1) = 0

       !print*, 'compute ivec'
       do n = 2, nbins
          ri(n+1) = ri(n) + hist(n)
       enddo

       hist = 0

       !print*, 'compute ovec'
       do i = 1, a_sz
          tf = (array(i) .lt. max_bin) .and. (array(i) .ge. min_) .and. (array(i) .le. max_)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
          ri(ri(ind) + hist(ind) + 1) = i - 1
          ri(1) = 1
          hist(1) = 0
          ri(2) = nbins
       enddo


    END SUBROUTINE reverse_indices_long

    SUBROUTINE reverse_indices_dlong(array, hist, ri, nbins, a_sz, ri_sz, min_, max_, max_bin, binsz)
       IMPLICIT NONE

       INTEGER*8 :: i, n, ri_sz, a_sz, y, ind, nbins
       INTEGER*8, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins) :: hist
       !f2py depend(nbins), hist

       INTEGER*4, DIMENSION(ri_sz), INTENT(INOUT) :: ri
       !f2py depend(ri_sz), ri

       INTEGER*8 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ri(2) = nbins
       hist(1) = 0

       !print*, 'compute ivec'
       do n = 2, nbins
          ri(n+1) = ri(n) + hist(n)
       enddo

       hist = 0

       !print*, 'compute ovec'
       do i = 1, a_sz
          tf = (array(i) .lt. max_bin) .and. (array(i) .ge. min_) .and. (array(i) .le. max_)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
          ri(ri(ind) + hist(ind) + 1) = i - 1
          ri(1) = 1
          hist(1) = 0
          ri(2) = nbins
       enddo

    END SUBROUTINE reverse_indices_dlong

    SUBROUTINE reverse_indices_float(array, hist, ri, nbins, a_sz, ri_sz, min_, max_, max_bin, binsz)
       IMPLICIT NONE

       INTEGER*8 :: i, n, ri_sz, a_sz, y, ind, nbins
       REAL*4, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins) :: hist
       !f2py depend(nbins), hist

       INTEGER*4, DIMENSION(ri_sz), INTENT(INOUT) :: ri
       !f2py depend(ri_sz), ri

       REAL*4 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ri(2) = nbins
       hist(1) = 0

       !print*, 'compute ivec'
       do n = 2, nbins
          ri(n+1) = ri(n) + hist(n)
       enddo

       hist = 0

       !print*, 'compute ovec'
       do i = 1, a_sz
          tf = (array(i) .lt. max_bin) .and. (array(i) .ge. min_) .and. (array(i) .le. max_)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
          ri(ri(ind) + hist(ind) + 1) = i - 1
          ri(1) = 1
          hist(1) = 0
          ri(2) = nbins
       enddo

    END SUBROUTINE reverse_indices_float

    SUBROUTINE reverse_indices_dfloat(array, hist, ri, nbins, a_sz, ri_sz, min_, max_, max_bin, binsz)
       IMPLICIT NONE

       INTEGER*8 :: i, n, ri_sz, a_sz, y, ind, nbins
       REAL*8, DIMENSION(a_sz), INTENT(IN) :: array
       !f2py depend(a_sz), array

       INTEGER*4, DIMENSION(nbins) :: hist
       !f2py depend(nbins), hist

       INTEGER*4, DIMENSION(ri_sz), INTENT(INOUT) :: ri
       !f2py depend(ri_sz), ri

       REAL*8 :: min_, max_
       REAL*8 :: binsz, max_bin
       INTEGER :: tf

       ri(2) = nbins
       hist(1) = 0

       !print*, 'compute ivec'
       do n = 2, nbins
          ri(n+1) = ri(n) + hist(n)
       enddo

       hist = 0

       !print*, 'compute ovec'
       do i = 1, a_sz
          tf = (array(i) .lt. max_bin) .and. (array(i) .ge. min_) .and. (array(i) .le. max_)
          y = tf * tf
          ind = 1 + ((floor((array(i) - min_) / binsz) + 1) * y)
          hist(ind) = hist(ind) + y
          ri(ri(ind) + hist(ind) + 1) = i - 1
          ri(1) = 1
          hist(1) = 0
          ri(2) = nbins
       enddo

    END SUBROUTINE reverse_indices_dfloat
END MODULE idl_histogram

