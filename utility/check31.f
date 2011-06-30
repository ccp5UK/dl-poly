c*********************************************************************
c
c     dl_poly utility to check the structure of the REVCON
c     file produced by TEST CASE 31
c
c     copyright daresbury laboratory
c     author w.smith september 2007
c
c*********************************************************************

      implicit real*8(a-h,o-z)
      dimension xxx(6),yyy(6),zzz(6),rrr(5)

c      REVCON file is standard input

      read(*,*)
      read(*,*)lev
      read(*,*)a
      read(*,*)x,b
      read(*,*)x,y,c

      do m=1,675

        do i=1,6

          read(*,*)
          read(*,*)xxx(i),yyy(i),zzz(i)
          if(lev.gt.0)read(*,*)
          if(lev.gt.1)read(*,*)

        enddo

        do i=1,5

          dx=xxx(i)-xxx(i+1)
          dy=yyy(i)-yyy(i+1)
          dz=zzz(i)-zzz(i+1)
          dx=dx-a*nint(dx/a)
          dy=dy-b*nint(dy/b)
          dz=dz-c*nint(dz/c)
          rrr(i)=sqrt(dx*dx+dy*dy+dz*dz)

        enddo

        write(*,'(i6,5f14.7)')m,(rrr(i),i=1,5)

      enddo
      end
