      subroutine rndrot
     x  (lrvel,natms,imcon,idnode,weight,xxx,yyy,zzz,vxx,vyy,vzz)

c*********************************************************************
c
c     dl-poly subroutine to rotate atom clusters through random angle
c     around the centre of mass
c
c     copyright daresbury laboratory 2000
c     author w.smith july 2000
c
c*********************************************************************

      include 'dl_params.inc'

      logical lrvel
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension vxx(mxatms),vyy(mxatms),vzz(mxatms)
      dimension weight(mxatms),rot(9),qtn(4)

c     check boundary condition

      if(imcon.ne.0.and.idnode.eq.0))then

        write(nprint,'(a)')'warning - system is a periodic cell'

      endif

c     construct random quaternions

      qtn(1)=(2.d0*duni()-1.d0)
      qtn(2)=(2.d0*duni()-1.d0)
      qtn(3)=(2.d0*duni()-1.d0)
      qtn(4)=(2.d0*duni()-1.d0)
      rnm=1.d0/sqrt(qtn(1)**2+qtn(2)**2+qtn(3)**2+qtn(4)**2)
      qtn(1)=rnm*qtn(1)
      qtn(2)=rnm*qtn(2)
      qtn(3)=rnm*qtn(3)
      qtn(4)=rnm*qtn(4)

c     construct rotation matrix

      rot(1) = qtn(1)**2+qtn(2)**2-qtn(3)**2-qtn(4)**2
      rot(2) = 2.d0*(qtn(2)*qtn(3) - qtn(1)*qtn(4))
      rot(3) = 2.d0*(qtn(2)*qtn(4) + qtn(1)*qtn(3))
      rot(4) = 2.d0*(qtn(2)*qtn(3) + qtn(1)*qtn(4))
      rot(5) = qtn(1)**2-qtn(2)**2+qtn(3)**2-qtn(4)**2
      rot(6) = 2.d0*(qtn(3)*qtn(4) - qtn(1)*qtn(2))
      rot(7) = 2.d0*(qtn(2)*qtn(4) - qtn(1)*qtn(3))
      rot(8) = 2.d0*(qtn(3)*qtn(4) + qtn(1)*qtn(2))
      rot(9) = qtn(1)**2-qtn(2)**2-qtn(3)**2+qtn(4)**2

c     calculate centre of mass

      cmx=0.d0
      cmy=0.d0
      cmz=0.d0
      do i=1,natms

        cmx=cmx+weight(i)*xxx(i)
        cmy=cmy+weight(i)*yyy(i)
        cmz=cmz+weight(i)*zzz(i)

      enddo
      cmx=cmx/dble(natms)
      cmy=cmy/dble(natms)
      cmz=cmz/dble(natms)

c     rotate cluster about centre of mass

      do i=1,natms

        xxs=xxx(i)-cmx
        yys=yyy(i)-cmy
        zzs=zzz(i)-cmz
        xxx(i)=rot(1)*xxs+rot(4)*yys+rot(7)*zzs+cmx
        yyy(i)=rot(2)*xxs+rot(5)*yys+rot(8)*zzs+cmy
        zzz(i)=rot(3)*xxs+rot(6)*yys+rot(9)*zzs+cmz

      enddo

c     rotate velocity vectors if requested

      if(lrvel)then
        
        do i=1,natms
          
          xxs=vxx(i)
          yys=vyy(i)
          zzs=vzz(i)
          vxx(i)=rot(1)*xxs+rot(4)*yys+rot(7)*zzs
          vyy(i)=rot(2)*xxs+rot(5)*yys+rot(8)*zzs
          vzz(i)=rot(3)*xxs+rot(6)*yys+rot(9)*zzs
          
        enddo

      endif

      return
      end
