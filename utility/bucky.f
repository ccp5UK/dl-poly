      program bucky
c     
c*********************************************************************
c     
c     construction of buckminster fullerene molecule
c     
c     copyright daresbury laboratory  november 1994
c     author w. smith     november 1994
c     
c     itt
c     2010-10-30 17:20:49
c     1.3
c     Exp
c
c*********************************************************************
c     
      implicit real*8(a-h,o-z)
      logical same
      dimension o(9,3),c(3,61)
      pi=3.141592653589793
      dr=pi/180.d0
      write(*,*) "enter required bondlength"
      read(*,*)bnd
      c(1,1)=0.d0
      c(2,1)=0.850650808d0*bnd
      c(3,1)=2.327438436d0*bnd
      o(1,1)=cos(dr*72.d0)
      o(2,1)=sin(dr*72.d0)
      o(3,1)=0.d0
      o(4,1)=-sin(dr*72.d0)
      o(5,1)=cos(dr*72.d0)
      o(6,1)=0.d0
      o(7,1)=0.d0
      o(8,1)=0.d0
      o(9,1)=1.d0
      gam=dr*37.37736815d0
      o(1,2)=cos(dr*120.d0)
      o(2,2)=sin(dr*120.d0)*cos(gam)
      o(3,2)=sin(dr*120.d0)*sin(gam)
      o(4,2)=-sin(dr*120.d0)*cos(gam)
      o(5,2)=cos(dr*120.d0)*cos(gam)**2+sin(gam)**2
      o(6,2)=cos(gam)*sin(gam)*(cos(dr*120.d0)-1.d0)
      o(7,2)=-sin(dr*120.d0)*sin(gam)
      o(8,2)=cos(gam)*sin(gam)*(cos(dr*120.d0)-1.d0)
      o(9,2)=cos(dr*120.d0)*sin(gam)**2+cos(gam)**2
      ngp=1
      open (8,file="FULLERENE")
      write(8,'(1p,3e16.8)') c(1,1),c(2,1),c(3,1)
      do m=1,100
         if(ngp.lt.60)then
            do i=1,2
               do j=1,60
                  if(j.le.ngp)then
                     n=ngp+1
                     c(1,n)=o(1,i)*c(1,j)+o(4,i)*c(2,j)+o(7,i)*c(3,j)
                     c(2,n)=o(2,i)*c(1,j)+o(5,i)*c(2,j)+o(8,i)*c(3,j)
                     c(3,n)=o(3,i)*c(1,j)+o(6,i)*c(2,j)+o(9,i)*c(3,j)
                     same=.false.
                     do k=1,ngp
                        if(.not.same)then
                           same=.true.
                           if(abs(c(1,n)-c(1,k)).gt.1.d-4)same=.false.
                           if(abs(c(2,n)-c(2,k)).gt.1.d-4)same=.false.
                           if(abs(c(3,n)-c(3,k)).gt.1.d-4)same=.false.
                        endif
                     enddo
                     if(.not.same)then
                        ngp=ngp+1
                        write(*,'(a,i3)')"atoms created = ",ngp
                        write(8,'(1p,3e16.8)') c(1,n),c(2,n),c(3,n)
                     endif
                  endif
               enddo
            enddo
         endif
      enddo
      
      close (8)
      
      end
      
