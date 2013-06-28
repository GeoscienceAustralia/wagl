      program bilinear
        parameter (nr=4,nl=4)
        character fname*150
        logical solut
        real data(20000,20000)
        real s1(nr),s2(nr),s3(nr),s4(nr)
        real a1(nr,nl),a2(nr,nl),a3(nr,nl),a4(nr,nl)
        real xx(nr),yy(nr)
        integer ix01,ix02,ix03,ix04,iy01,iy02,iy03,iy04
c       cof is bilinear coefficionts
        double precision cof1(nr),cof2(nr),cof3(nr),cof4(nr)

c       program is used to interpolate modtran output from four coordinators
c       nrow image lines, ncol, image column
c       nr,nl are the four coordinators which are used to
c       generate the binear interpolation
c       bilinear equation is expressed as:
c       f(x,y)=a+b(x-x0)+c(y-y0)+d(x-x0)(y-y0)
c       where a,b,c,d are the four coefficients you need to calculate
c       from here. x,y are the pixel coordinator and x0,y0 are the
c       centre pixel coordinators of box
        if (IARGC() .ne. 5) then
          write(*,*) 'Error: Required parameters not specified
     & properly!'
          stop 10
        endif
c
c       open coordinators file for the 9 points
        call GETARG(1, fname)
        open(1,file=fname)

c       read line and column
        read(1,*)nrow,ncol
c       read upperleft coordinator, x is line,y is column
        read(1,*)ixul,iyul
c       read uppermiddleleft coordinator
        read(1,*)ixum,iyum
c       read upperright coordinator
        read(1,*)ixur,iyur
c       read middleleft coordinator
        read(1,*)ixml,iyml
c       read middlemiddle coordinator
        read(1,*)ixmm,iymm
c       read middleright coordinator
        read(1,*)ixmr,iymr
c       read lowleft coordinator
        read(1,*)ixll,iyll
c       read lowmiddle coordinator
        read(1,*)ixlm,iylm
c       read lowright coordinator
        read(1,*)ixlr,iylr

c       for first box
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
          write(12,*)(a1(i,j),j=1,nl)
        enddo

c       for the second box
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
          write(12,*)(a2(i,j),j=1,nl)
        enddo

c       for the third box
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
          write(12,*)(a3(i,j),j=1,nl)
        enddo

c       for the fourth box
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
          write(12,*)(a4(i,j),j=1,nl)
        enddo

c       input the four coordinator data
c
c       open atmospheric parameter files, startend file
c       and centerline file
        do i=2,4
          call GETARG(i, fname)
          open(i,file=fname,status='old')
        enddo

c       open output file
        call  GETARG(5, fname)
        open(21,file=fname,form='unformatted',
     &    access ='direct',recl=4*ncol)

        read(2,*)(s1(i),i=1,nr)
        read(2,*)(s2(i),i=1,nr)
        read(2,*)(s3(i),i=1,nr)
        read(2,*)(s4(i),i=1,nr)

        call gauss(a1,nr,nl,s1,cof1,solut)
        if (solut) then
          write(11,*)(cof1(i),i=1,nr)
        endif

        call gauss(a2,nr,nl,s2,cof2,solut)

        if (solut) then
          write(11,*)(cof2(i),i=1,nr)
        endif

        call gauss(a3,nr,nl,s3,cof3,solut)

        if (solut) then
        write(11,*)(cof3(i),i=1,nr)
        endif

        call gauss(a4,nr,nl,s4,cof4,solut)
        if (solut) then
          write(11,*)(cof4(i),i=1,nr)
        endif
c       read header of center line
        read(4,*)angle
        read(4,*)mm,mn

c       calculate bilinear coefficients for four box
        do i=1,nrow
          read (3,*)line,istart,iend
          read (4,*)line,icenter,clat,clon

          do j=1,ncol
            jj=j-istart+1
            if (i .lt. ixml) then
              if (j .ge. istart .and.j .lt.icenter) then
c               box 1
                data(i,j)=cof1(1)+cof1(2)*(i-ix01)+cof1(3)*(jj-iy01)+
     &            cof1(4)*(i-ix01)*(jj-iy01)
              elseif (j .ge. icenter .and.j .le.iend) then
c               box2
                data(i,j)=cof2(1)+cof2(2)*(i-ix02)+cof2(3)*(jj-iy02)+
     &            cof2(4)*(i-ix02)*(jj-iy02)
              else
                data(i,j)=-999
              endif
            elseif (i .ge. ixml .and. i .le.ixll) then
              if (j .ge. istart .and.j .lt.icenter) then
c               box 3
                data(i,j)=cof3(1)+cof3(2)*(i-ix03)+cof3(3)*(jj-iy03)+
     &            cof3(4)*(i-ix03)*(jj-iy03)
              elseif (j .ge. icenter .and.j .le.iend) then
c               box4
                data(i,j)=cof4(1)+cof4(2)*(i-ix04)+cof4(3)*(jj-iy04)+
     &            cof4(4)*(i-ix04)*(jj-iy04)
              else
                data(i,j)=-999
              endif
            else
              data(i,j)=-999
            endif
          enddo
          write(21,rec=i)(data(i,j),j=1,ncol)
        enddo
        close(21)
        do i=1,4
          close(i)
        enddo
        stop
      end

c------------------------------------------------------

      subroutine gauss(a,m,n,b,x,solut)
c     gauss elimination solve ax=b with row and column pivotion
        real a(m,n),b(m)
        double precision x(m)
        logical solut
c       elimination
c       k----eliminate column
        do 30 k=1,n
          ik=k
c         find the best pivot element
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
c       now back-substiute
        x(m)=b(m)
        do 60 i=m-1,1,-1
          sum=b(i)
          do 50 j=i+1,4
            sum=sum-a(i,j)*x(j)
50        continue
          x(i)=sum
60      continue
        return
      end

c------------------------------------------------------

      subroutine spiv(a,m,n,b,k,solut)
        real a(m,n),b(m)
        logical solut
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
      end
