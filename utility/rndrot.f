      subroutine rndrot
     x  (lrotv,idnode,mxnode,natms,xxx,yyy,zzz,vxx,vyy,vzz)

c*********************************************************************
c
c     dl_poly routine to perform a random rotation on a 3D rigid body
c
c     author - w.smith july 2000
c     copyright daresbury laboratory 2000
c
c*********************************************************************

#include 'dl_params.inc'

      logical lrotv

      dimension rot(9),qtn(4)
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
c
c     generate random orientation

      qtn(1)=2.d0*duni()-1.d0
      qtn(2)=2.d0*duni()-1.d0
      qtn(3)=2.d0*duni()-1.d0
      qtn(4)=2.d0*duni()-1.d0
      qqq=sqrt(qtn(1)**2+qtn(2)**2+qtn(3)**2+qtn(4)**2)
      qtn(1)=qtn(1)/qqq
      qtn(2)=qtn(2)/qqq
      qtn(3)=qtn(3)/qqq
      qtn(4)=qtn(4)/qqq

      rot(1) = qtn(1)**2+qtn(2)**2-qtn(3)**2-qtn(4)**2
      rot(2) = 2.d0*(qtn(2)*qtn(3) - qtn(1)*qtn(4))
      rot(3) = 2.d0*(qtn(2)*qtn(4) + qtn(1)*qtn(3))
      rot(4) = 2.d0*(qtn(2)*qtn(3) + qtn(1)*qtn(4))
      rot(5) = qtn(1)**2-qtn(2)**2+qtn(3)**2-qtn(4)**2
      rot(6) = 2.d0*(qtn(3)*qtn(4) - qtn(1)*qtn(2))
      rot(7) = 2.d0*(qtn(2)*qtn(4) - qtn(1)*qtn(3))
      rot(8) = 2.d0*(qtn(3)*qtn(4) + qtn(1)*qtn(2))
      rot(9) = qtn(1)**2-qtn(2)**2-qtn(3)**2+qtn(4)**2
c
c     perform coordinate rotation

      do i=1,natms

        txx=xxx(i)
        tyy=yyy(i)
        tzz=zzz(i)
        xxx(i)=rot(1)*txx+rot(4)*tyy+rot(7)*tzz
        yyy(i)=rot(2)*txx+rot(5)*tyy+rot(8)*tzz
        zzz(i)=rot(3)*txx+rot(6)*tyy+rot(9)*tzz

      enddo

c
c     perform velocity rotation

      if(lrotv)then

        do i=1,natms

          txx=vxx(i)
          tyy=vyy(i)
          tzz=vzz(i)
          vxx(i)=rot(1)*txx+rot(4)*tyy+rot(7)*tzz
          vyy(i)=rot(2)*txx+rot(5)*tyy+rot(8)*tzz
          vzz(i)=rot(3)*txx+rot(6)*tyy+rot(9)*tzz

        enddo

      endif

      return
      end
