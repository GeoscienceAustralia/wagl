SUBROUTINE bilinear_interpolation(nrow, ncol, coordinator, s1, s2, s3, s4, &
                                  cstart, cend, centreline, res)

        integer, intent(in) :: nrow, ncol
        integer, dimension(9, 2), intent(in) :: coordinator
        real*4, dimension(4), intent(in) :: s1, s2, s3, s4
        integer, dimension(ncol), intent(in) :: cstart, cend, centreline
        real*4, dimension(nrow, ncol), intent(inout) :: res

        integer, parameter :: nr=4, nl=4
        real*4, dimension(4, 4) :: a1, a2, a3, a4
        real*4, dimension(4) :: xx, yy
        double precision, dimension(4) :: cof1, cof2, cof3, cof4
        logical solut
        integer i, j, jj, line, istart, iend, icenter

!f2py depend(ncol), cstart, cend, centreline
!f2py depend(nrow, ncol), res

!       the result array, res, is transposed as it is passed into
!       fortran. So nrow, and ncol labelled here in fortran are the
!       cols and rows back in python land.
!       The program was originally written by Fuqin Li,
!       and now converted to a python module via f2py

!       cof is bilinear coefficients
!       boxline (line, istart, iend); now replaced by cstart, cend.

!       program is used to interpolate modtran output from four coordinators
!       nrow image lines, ncol, image column
!       nr,nl are the four coordinators which are used to
!       generate the binear interpolation
!       bilinear equation is expressed as:
!       f(x,y)=a+b(x-x0)+c(y-y0)+d(x-x0)(y-y0)
!       where a,b,c,d are the four coefficients you need to calculate
!       from here. x,y are the pixel coordinator and x0,y0 are the
!       centre pixel coordinators of box

!       reference the 9 point coordinates
        ixul = coordinator(1, 1)
        iyul = coordinator(1, 2)
        ixum = coordinator(2, 1) 
        iyum = coordinator(2, 2)
        ixur = coordinator(3, 1)
        iyur = coordinator(3, 2)
        ixml = coordinator(4, 1)
        iyml = coordinator(4, 2)
        ixmm = coordinator(5, 1)
        iymm = coordinator(5, 2)
        ixmr = coordinator(6, 1)
        iymr = coordinator(6, 2)
        ixll = coordinator(7, 1)
        iyll = coordinator(7, 2)
        ixlm = coordinator(8, 1)
        iylm = coordinator(8, 2)
        ixlr = coordinator(9, 1)
        iylr = coordinator(9, 2)

!       for first box
        ix01=(ixul+ixml)/2-ixul+1
        iy01=(iyul+iyum)/2-iyul+1
        xx(1)=ixul-ixul+1
        yy(1)=iyul-iyul+1
        xx(2)=ixum-ixul+1
        yy(2)=iyum-iyul+1
        xx(3)=ixml-ixul+1
        yy(3)=iyml-iyml+1
        xx(4)=ixmm-ixum+1
        yy(4)=iymm-iyml+1
        do i=1,nr
          do j=1,nl
            a1(i,1)=1
            a1(i,2)=(xx(i)-ix01)
            a1(i,3)=(yy(i)-iy01)
            a1(i,4)=(xx(i)-ix01)*(yy(i)-iy01)
          enddo
        enddo

!       for the second box
        ix02=(ixum+ixmm)/2-ixum+1
        iy02=(iyur+iyum)/2-iyul+1
        xx(1)=ixum-ixum+1
        yy(1)=iyum-iyul+1
        xx(2)=ixur-ixur+1
        yy(2)=iyur-iyul+1
        xx(3)=ixmm-ixum+1
        yy(3)=iymm-iyml+1
        xx(4)=ixmr-ixur+1
        yy(4)=iymr-iyml+1
        do i=1,nr
          do j=1,nl
            a2(i,1)=1
            a2(i,2)=(xx(i)-ix02)
            a2(i,3)=(yy(i)-iy02)
            a2(i,4)=(xx(i)-ix02)*(yy(i)-iy02)
          enddo
        enddo

!       for the third box
        ix03=(ixml+ixll)/2-ixul+1
        iy03=(iyml+iymm)/2-iyml+1
        xx(1)=ixml-ixul+1
        yy(1)=iyml-iyml+1
        xx(2)=ixmm-ixum+1
        yy(2)=iymm-iyml+1
        xx(3)=ixll-ixul+1
        yy(3)=iyll-iyll+1
        xx(4)=ixlm-ixum+1
        yy(4)=iylm-iyll+1
        do i=1,nr
          do j=1,nl
            a3(i,1)=1
            a3(i,2)=(xx(i)-ix03)
            a3(i,3)=(yy(i)-iy03)
            a3(i,4)=(xx(i)-ix03)*(yy(i)-iy03)
          enddo
        enddo

!       for the fourth box
        ix04=(ixmm+ixlm)/2-ixum+1
        iy04=(iymm+iymr)/2-iyml+1
        xx(1)=ixmm-ixum+1
        yy(1)=iymm-iyml+1
        xx(2)=ixmr-ixur+1
        yy(2)=iymr-iyml+1
        xx(3)=ixlm-ixum+1
        yy(3)=iylm-iyll+1
        xx(4)=ixlr-ixur+1
        yy(4)=iylr-iyll+1
        do i=1,nr
          do j=1,nl
            a4(i,1)=1
            a4(i,2)=(xx(i)-ix04)
            a4(i,3)=(yy(i)-iy04)
            a4(i,4)=(xx(i)-ix04)*(yy(i)-iy04)
          enddo
        enddo


        call gauss(a1,nr,nl,s1,cof1,solut)

        call gauss(a2,nr,nl,s2,cof2,solut)

        call gauss(a3,nr,nl,s3,cof3,solut)

        call gauss(a4,nr,nl,s4,cof4,solut)

!       calculate bilinear coefficients for four box
        do i=1,ncol
           istart = cstart(i)
           iend = cend(i)
           icenter = centreline(i)

          do j=1,nrow
            jj=j-istart+1
            if (i .lt. ixml) then
              if (j .ge. istart .and.j .lt.icenter) then
!               box 1
                res(j,i)=cof1(1)+cof1(2)*(i-ix01)+cof1(3)*(jj-iy01)+ &
                 cof1(4)*(i-ix01)*(jj-iy01)
              elseif (j .ge. icenter .and.j .le.iend) then
!               box2
                res(j,i)=cof2(1)+cof2(2)*(i-ix02)+cof2(3)*(jj-iy02)+ &
                 cof2(4)*(i-ix02)*(jj-iy02)
              else
                res(j,i)=-999
              endif
            elseif (i .ge. ixml .and. i .le.ixll) then
              if (j .ge. istart .and.j .lt.icenter) then
!               box 3
                res(j,i)=cof3(1)+cof3(2)*(i-ix03)+cof3(3)*(jj-iy03)+ &
                  cof3(4)*(i-ix03)*(jj-iy03)
              elseif (j .ge. icenter .and.j .le.iend) then
!               box4
                res(j,i)=cof4(1)+cof4(2)*(i-ix04)+cof4(3)*(jj-iy04)+ &
                 cof4(4)*(i-ix04)*(jj-iy04)
              else
                res(j,i)=-999
              endif
            else
              res(j,i)=-999
            endif
          enddo
        enddo
        return
END SUBROUTINE bilinear_interpolation

!------------------------------------------------------

SUBROUTINE gauss(a,m,n,b,x,solut)
!     gauss elimination solve ax=b with row and column pivotion
        real a(m,n),b(m)
        double precision x(m)
        logical solut
        integer m, n, k
!       elimination
!       k----eliminate column
        do 30 k=1,n
          ik=k
!         find the best pivot element
          call spiv(a,m,n,b,ik,solut)
          if(solut) then
            do 20 i=k,m
              c=a(i,k)
              if (abs(c) .ne.0.0) then
                do 10 j=k,n
                  a(i,j)=a(i,j)/c
                  if(i.gt.k) then
                    a(i,j)=a(i,j)-a(k,j)
                  endif
10              continue
                b(i)=b(i)/c
                if (i.gt.k) then
                  b(i)=b(i)-b(k)
                endif
              endif
20          continue
          else
            return
          endif
30      continue
!       now back-substiute
        x(m)=b(m)
        do 60 i=m-1,1,-1
          sum=b(i)
          do 50 j=i+1,4
            sum=sum-a(i,j)*x(j)
50        continue
          x(i)=sum
60      continue
        return
END SUBROUTINE gauss

!------------------------------------------------------

SUBROUTINE spiv(a,m,n,b,k,solut)
        real a(m,n),b(m)
        logical solut
        integer k
        integer n, m
        solut=.false.
        do i=k,m
          if (a(i,k).ne.0.0) then
            solut=.true.
          endif
        enddo
        if (solut) then
          l=k
          do j=k,m
            if (abs(a(j,k)).gt.abs(a(l,k))) then
              l=j
            endif
          enddo
          if (l.ne.k) then
            tm=b(k)
            b(k)=b(l)
            b(l)=tm
            do j=k,m
              tm=a(k,j)
              a(k,j)=a(l,j)
              a(l,j)=tm
            enddo
          endif
        endif
        return
END SUBROUTINE spiv
