      subroutine cgropt
     x (keyopt,natms,step,ida,hnrm,fff,grad,hhh,dxyz,xyz)

c*********************************************************************
c
c     dl_poly conjugate gradient routine
c
c     copyright - daresbury laboratory
c     author    - w.smith 2002
c
*********************************************************************

      implicit none

      integer ida(*)
      integer i,j,keyopt,natms
      real*8 ggg,stride,gam2,step
      real*8 hnrm(*),fff(*),grad(*),hhh(3,*),dxyz(3,*),xyz(3,*)

c     Magnitude of current gradient vector

      ggg=0.d0
      do i=1,natms

        ggg=ggg+(dxyz(1,i)*dxyz(1,i)+dxyz(2,i)*dxyz(2,i)+
     x    dxyz(3,i)*dxyz(3,i))

      enddo
      ggg=sqrt(ggg)

      if(keyopt .eq. 0)then

c     Set original search direction (vector hhh)

        keyopt=1
        hnrm(1)=ggg
        grad(3)=ggg
        fff(3)=fff(1)

        do i=1,natms

          j=ida(i)
          hhh(1,i)=dxyz(1,i)
          hhh(2,i)=dxyz(2,i)
          hhh(3,i)=dxyz(3,i)
          xyz(1,j)=xyz(1,j)+step*hhh(1,i)
          xyz(2,j)=xyz(2,j)+step*hhh(2,i)
          xyz(3,j)=xyz(3,j)+step*hhh(3,i)

        enddo

      else if(keyopt .eq. 1)then

c     Line search along chosen direction

        stride=step
        fff(2)=fff(3)
        fff(3)=fff(1)
        grad(2)=grad(3)

        grad(3)=0.d0
        do i=1,natms

          grad(3)=grad(3)+(hhh(1,i)*dxyz(1,i)+hhh(2,i)*dxyz(2,i)+
     x      hhh(3,i)*dxyz(3,i))

        enddo
        grad(3)=grad(3)/hnrm(1)

c     Linear extrapolation to minimum

        if(grad(3) .lt. 0)then

          stride=step*grad(3)/(grad(2)-grad(3))
          keyopt=2

        endif

        do i=1,natms

          j=ida(i)
          xyz(1,j)=xyz(1,j)+stride*hhh(1,i)
          xyz(2,j)=xyz(2,j)+stride*hhh(2,i)
          xyz(3,j)=xyz(3,j)+stride*hhh(3,i)

        enddo

      else if(keyopt .eq. 2)then

        fff(2)=fff(3)
        fff(3)=fff(1)
        grad(2)=grad(3)
        grad(1)=ggg
        grad(3)=0.d0

        do i=1,natms

          grad(3)=grad(3)+(hhh(1,i)*dxyz(1,i)+hhh(2,i)*dxyz(2,i)+
     x      hhh(3,i)*dxyz(3,i))

        enddo

        grad(3)=grad(3)/hnrm(1)
        keyopt=3
        return

      else if(keyopt .eq. 3)then

        fff(2)=fff(3)
        fff(3)=fff(1)

c     Check for global convergence

        if(abs(ggg/natms) .lt. 0.0000001d0)then

          keyopt=999
          return

        endif

c     Construct conjugate search vector

        gam2=(ggg/grad(1))**2
        hnrm(1)=0.d0
        grad(3)=0.d0

        do i=1,natms

          hhh(1,i)=dxyz(1,i)+gam2*hhh(1,i)
          hhh(2,i)=dxyz(2,i)+gam2*hhh(2,i)
          hhh(3,i)=dxyz(3,i)+gam2*hhh(3,i)
          hnrm(1)=hnrm(1)+(hhh(1,i)*hhh(1,i)+hhh(2,i)*hhh(2,i)+
     x      hhh(3,i)*hhh(3,i))
          grad(3)=grad(3)+(hhh(1,i)*dxyz(1,i)+hhh(2,i)*dxyz(2,i)+
     x      hhh(3,i)*dxyz(3,i))


        enddo
        hnrm(1)=sqrt(hnrm(1))
        grad(3)=grad(3)/hnrm(1)

        do i=1,natms

          j=ida(i)
          xyz(1,j)=xyz(1,j)+step*hhh(1,i)
          xyz(2,j)=xyz(2,j)+step*hhh(2,i)
          xyz(3,j)=xyz(3,j)+step*hhh(3,i)

        enddo

        keyopt=1
	return
      endif
      return
      end
