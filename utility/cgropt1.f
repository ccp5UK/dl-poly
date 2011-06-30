      subroutine cgropt(keyopt,natms,tol,step,hnrm,fff,grad,
     x  hhh,xyz,dxyz)

c*********************************************************************
c
c     conjugate gradient routine terminates when force per atom < tol
c
c     keyopt = 0 on first entry, keep calling until keyopt=999
c     natms = number of atoms
c     step = default step in Ansgtroms (suggest 0.05)
c     hnrm = search vector magnitude
c     fff(3) = stored values of last 3 energies
c     grad(3) = stored values of last 3 forces along search vector
c     hhh(3,natms) = search vector
c     xyz(3,natms) = atomic position vectors
c     dxyz(3,natms) = atomic force vectors
c
c     copyright - daresbury laboratory
c     author    - w.smith 2002
c
c*********************************************************************

      implicit none

      integer keyopt,i,natms
      real*8 ggg,stride,gam2,step,hnrm,tol
      real*8 dxyz(3,*),hhh(3,*),xyz(3,*),grad(3),fff(3)

c     Magnitude of current gradient vector

      ggg=0.d0
      do i=1,natms

        ggg=ggg+dxyz(1,i)**2+dxyz(2,i)**2+dxyz(3,i)**2

      enddo
      ggg=sqrt(ggg)

      if(keyopt.eq.0)then

c     Set original search direction (vector hhh)

        keyopt=1
        hnrm=ggg
        grad(3)=ggg
        fff(3)=fff(1)

        do i=1,natms

          hhh(1,i)=dxyz(1,i)
          hhh(2,i)=dxyz(2,i)
          hhh(3,i)=dxyz(3,i)
          xyz(1,i)=xyz(1,i)+step*hhh(1,i)
          xyz(2,i)=xyz(2,i)+step*hhh(2,i)
          xyz(3,i)=xyz(3,i)+step*hhh(3,i)

        enddo

      else if(keyopt.eq.1)then

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
        grad(3)=grad(3)/hnrm

c     Linear extrapolation to minimum

        if(grad(3).lt.0.d0)then

          stride=step*grad(3)/(grad(2)-grad(3))
          keyopt=2

        endif

        do i=1,natms

          xyz(1,i)=xyz(1,i)+stride*hhh(1,i)
          xyz(2,i)=xyz(2,i)+stride*hhh(2,i)
          xyz(3,i)=xyz(3,i)+stride*hhh(3,i)

        enddo

      else if(keyopt.eq.2)then

        fff(2)=fff(3)
        fff(3)=fff(1)
        grad(2)=grad(3)
        grad(1)=ggg
        grad(3)=0.d0
        do i=1,natms

          grad(3)=grad(3)+(hhh(1,i)*dxyz(1,i)+hhh(2,i)*dxyz(2,i)+
     x      hhh(3,i)*dxyz(3,i))

        enddo
        grad(3)=grad(3)/hnrm
        keyopt=3
        return

      else if(keyopt.eq.3)then

        fff(2)=fff(3)
        fff(3)=fff(1)

c     Check for global convergence

        if(abs(ggg/dble(natms)).lt.tol)then

          keyopt=999
          return

        endif

c     Construct conjugate search vector

        gam2=(ggg/grad(1))**2

        hnrm=0.d0
        grad(3)=0.d0
        do i=1,natms

          hhh(1,i)=dxyz(1,i)+gam2*hhh(1,i)
          hhh(2,i)=dxyz(2,i)+gam2*hhh(2,i)
          hhh(3,i)=dxyz(3,i)+gam2*hhh(3,i)
          hnrm=hnrm+hhh(1,i)*hhh(1,i)+hhh(2,i)*hhh(2,i)+
     x      hhh(3,i)*hhh(3,i)
          grad(3)=grad(3)+hhh(1,i)*dxyz(1,i)+hhh(2,i)*dxyz(2,i)+
     x      hhh(3,i)*dxyz(3,i)

        enddo
        hnrm=sqrt(hnrm)
        grad(3)=grad(3)/hnrm

        do i=1,natms

          xyz(1,i)=xyz(1,i)+step*hhh(1,i)
          xyz(2,i)=xyz(2,i)+step*hhh(2,i)
          xyz(3,i)=xyz(3,i)+step*hhh(3,i)

        enddo

        keyopt=1

      endif

      return
      end
