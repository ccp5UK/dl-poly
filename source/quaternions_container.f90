!!!!!!!!!!!!!!!!!!!!! THIS IS QUATERNIONS_CONTAINER !!!!!!!!!!!!!!!!!!!!
!
! Subroutine q_setup - sets quaternions for RB dynamics
!
! Subroutine getrotmat - constructs rotation matrix from quaternions
!
! Subroutine no_squish - implements the symplectic no_squish quaternion
!                        algorithm
!
! Subroutine q_update - update the quaternions in LFV scope
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine q_setup()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for setting up RBs' quaternions
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gmax
  Use setup_module
  Use config_module,      Only : cell,xxx,yyy,zzz
  Use rigid_bodies_module

  Implicit None

  Integer           :: fail,imcon,irgd,jrgd,krgd,lrgd,rgdtyp, &
                       ill,i1,i2,i3
  Real( Kind = wp ) :: rot(1:9),aa(1:9),rsq,tol, &
                       aq,bq,cq,dq,eq,fq,gq,hq,rnorm, x,y,z

  Real( Kind = wp ), Allocatable :: gxx(:),gyy(:),gzz(:)

  fail = 0
  Allocate (gxx(1:mxlrgd*Max(mxrgd,mxtrgd)),gyy(1:mxlrgd*Max(mxrgd,mxtrgd)), &
            gzz(1:mxlrgd*Max(mxrgd,mxtrgd)), Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'q_setup allocation failure, node: ', idnode
     Call error(0)
  End If


! Recover/localise imcon

  imcon=rgdimc

! quaternions for all RB on this domain

  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then
        i1=indrgd(rgdind(1,rgdtyp),irgd)
        i2=indrgd(rgdind(2,rgdtyp),irgd)
        i3=indrgd(rgdind(3,rgdtyp),irgd)

! group basis vectors

        aa(1)=xxx(i1)-xxx(i2)
        aa(4)=yyy(i1)-yyy(i2)
        aa(7)=zzz(i1)-zzz(i2)

! minimum image convention for bond vectors

        Call images(imcon,cell,1,aa(1),aa(4),aa(7))

        ill=rgdind(0,rgdtyp)
        If (ill == 0) Then
           aa(2)=xxx(i1)-xxx(i3)
           aa(5)=yyy(i1)-yyy(i3)
           aa(8)=zzz(i1)-zzz(i3)
        Else
           rsq=Sqrt(aa(1)**2+aa(4)**2+aa(7)**2)
           If      (Abs(aa(7)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(4)**2+aa(7)**2)
              aa(2)= 0.0_wp
              aa(5)= aa(7)/rsq
              aa(8)=-aa(4)/rsq
           Else If (Abs(aa(4)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(4)**2+aa(1)**2)
              aa(2)=-aa(4)/rsq
              aa(5)= aa(1)/rsq
              aa(8)= 0.0_wp
           Else If (Abs(aa(1)/rsq) > 0.5_wp) Then
              rsq=Sqrt(aa(1)**2+aa(7)**2)
              aa(2)=-aa(7)/rsq
              aa(5)= 0.0_wp
              aa(8)= aa(1)/rsq
           End If
        End If

! minimum image convention for bond vectors

        Call images(imcon,cell,1,aa(2),aa(5),aa(8))

        aa(3)=aa(4)*aa(8)-aa(7)*aa(5)
        aa(6)=aa(7)*aa(2)-aa(1)*aa(8)
        aa(9)=aa(1)*aa(5)-aa(4)*aa(2)

! group rotational matrix

        rot(1)=rgdaxs(1,rgdtyp)*aa(1)+rgdaxs(4,rgdtyp)*aa(2)+rgdaxs(7,rgdtyp)*aa(3)
        rot(2)=rgdaxs(2,rgdtyp)*aa(1)+rgdaxs(5,rgdtyp)*aa(2)+rgdaxs(8,rgdtyp)*aa(3)
        rot(3)=rgdaxs(3,rgdtyp)*aa(1)+rgdaxs(6,rgdtyp)*aa(2)+rgdaxs(9,rgdtyp)*aa(3)
        rot(4)=rgdaxs(1,rgdtyp)*aa(4)+rgdaxs(4,rgdtyp)*aa(5)+rgdaxs(7,rgdtyp)*aa(6)
        rot(5)=rgdaxs(2,rgdtyp)*aa(4)+rgdaxs(5,rgdtyp)*aa(5)+rgdaxs(8,rgdtyp)*aa(6)
        rot(6)=rgdaxs(3,rgdtyp)*aa(4)+rgdaxs(6,rgdtyp)*aa(5)+rgdaxs(9,rgdtyp)*aa(6)
        rot(7)=rgdaxs(1,rgdtyp)*aa(7)+rgdaxs(4,rgdtyp)*aa(8)+rgdaxs(7,rgdtyp)*aa(9)
        rot(8)=rgdaxs(2,rgdtyp)*aa(7)+rgdaxs(5,rgdtyp)*aa(8)+rgdaxs(8,rgdtyp)*aa(9)
        rot(9)=rgdaxs(3,rgdtyp)*aa(7)+rgdaxs(6,rgdtyp)*aa(8)+rgdaxs(9,rgdtyp)*aa(9)

! determine quaternions from rotational matrix

        aq=rot(1)+rot(5)
        bq=rot(2)-rot(4)
        cq=rot(6)-rot(8)
        dq=rot(2)+rot(4)
        eq=rot(3)+rot(7)
        fq=rot(6)+rot(8)
        gq=rot(3)-rot(7)
        hq=rot(1)-rot(5)

        q0(irgd)=0.5_wp*Sqrt(aq+Sqrt(aq*aq+bq*bq))

        If (q0(irgd) > 1.0e-4_wp) Then
           q1(irgd)=-0.25_wp*cq/q0(irgd)
           q2(irgd)= 0.25_wp*gq/q0(irgd)
           q3(irgd)=-0.25_wp*bq/q0(irgd)
        Else
           q1(irgd)=0.5_wp*Sqrt(hq+Sqrt(hq*hq+dq*dq))

           If (q1(irgd) > 1.0e-4_wp) Then
              q2(irgd)=0.25_wp*dq/q1(irgd)
              q3(irgd)=0.25_wp*eq/q1(irgd)
           Else
              q2(irgd)=0.5_wp*Sqrt(-hq+Sqrt(hq*hq+dq*dq))

              If (q2(irgd) > 1.0e-4_wp) Then
                 q3(irgd)=0.25_wp*fq/q2(irgd)
              Else
                 q3(irgd)=1.0_wp
              End If
           End If
        End If

! normalise quaternions

        rnorm=1.0_wp/Sqrt(q0(irgd)**2+q1(irgd)**2+q2(irgd)**2+q3(irgd)**2)
        q0(irgd)=rnorm*q0(irgd)
        q1(irgd)=rnorm*q1(irgd)
        q2(irgd)=rnorm*q2(irgd)
        q3(irgd)=rnorm*q3(irgd)
     Else
        q0(irgd)=0.0_wp
        q1(irgd)=0.0_wp
        q2(irgd)=0.0_wp
        q3(irgd)=1.0_wp
     End If
  End Do

! Check of quaternion set up with atomic positions

  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then ! Not that it matters

! new rotational matrix

        Call getrotmat(q0(irgd),q1(irgd),q2(irgd),q3(irgd),rot)

        Do jrgd=1,lrgd
           krgd=krgd+1

           x=rot(1)*rgdx(jrgd,rgdtyp)+rot(2)*rgdy(jrgd,rgdtyp)+rot(3)*rgdz(jrgd,rgdtyp)+rgdxxx(irgd)
           y=rot(4)*rgdx(jrgd,rgdtyp)+rot(5)*rgdy(jrgd,rgdtyp)+rot(6)*rgdz(jrgd,rgdtyp)+rgdyyy(irgd)
           z=rot(7)*rgdx(jrgd,rgdtyp)+rot(8)*rgdy(jrgd,rgdtyp)+rot(9)*rgdz(jrgd,rgdtyp)+rgdzzz(irgd)

           gxx(krgd)=xxx(indrgd(jrgd,irgd))-x
           gyy(krgd)=yyy(indrgd(jrgd,irgd))-y
           gzz(krgd)=zzz(indrgd(jrgd,irgd))-z
        End Do

     End If
  End Do

! minimum image convention for bond vectors

  Call images(imcon,cell,krgd,gxx,gyy,gzz)

! test quaternion setup

  tol=1.0e-2_wp
  rsq=0.0_wp
  krgd=0
  Do irgd=1,ntrgd
     rgdtyp=listrgd(0,irgd)

! For all good RBs

     lrgd=listrgd(-1,irgd)
     If (rgdfrz(0,rgdtyp) < lrgd) Then ! Keep consistent with above

        Do jrgd=1,lrgd
           krgd=krgd+1

           rsq=Max(rsq,gxx(krgd)**2+gyy(krgd)**2+gzz(krgd)**2)
        End Do
     End If
  End Do
  If (mxnode > 1) Call gmax(rsq)
  If (rsq > tol) Call error(648)

  Deallocate (gxx,gyy,gzz, Stat = fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'q_setup deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine q_setup

Subroutine getrotmat(q0,q1,q2,q3,rot)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to construct rotation matrix from quaternions using
! x convention for euler angles
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: q0,q1,q2,q3
  Real( Kind = wp ), Intent(   Out ) :: rot(1:9)

  rot(1)=q0**2+q1**2-q2**2-q3**2
  rot(2)=2.0_wp*(q1*q2-q0*q3)
  rot(3)=2.0_wp*(q1*q3+q0*q2)
  rot(4)=2.0_wp*(q1*q2+q0*q3)
  rot(5)=q0**2-q1**2+q2**2-q3**2
  rot(6)=2.0_wp*(q2*q3-q0*q1)
  rot(7)=2.0_wp*(q1*q3-q0*q2)
  rot(8)=2.0_wp*(q2*q3+q0*q1)
  rot(9)=q0**2-q1**2-q2**2+q3**2

End Subroutine getrotmat

Subroutine no_squish                    &
           (tstep,rgdrix,rgdriy,rgdriz, &
           q0,q1,q2,q3,p0,p1,p2,p3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine routine to implement the symplectic no_squish
! quaternion algorithm of miller et al j.chem.phys 116 (2002) 8649
!
! copyright - daresbury laboratory
! author    - m.leslie january 2004
! adapted   - i.t.todorov october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: tstep,rgdrix,rgdriy,rgdriz
  Real( Kind = wp ), Intent( InOut ) :: q0,q1,q2,q3,p0,p1,p2,p3

  Integer, Parameter :: mrot = 10 ! Rotational timesteps
  Integer            :: m
  Real( Kind = wp )  :: rotstep,zetax,zetay,zetaz,cs,sn, &
                        qn1(0:3),pq2(0:3), &
                        qn2(0:3),pq3(0:3), &
                        qn3(0:3),pq4(0:3)

! rotation: iterate over mrot rotational time steps

  rotstep=tstep/Real(mrot,wp) ! rotational time step
  Do m=1,mrot
     zetaz=0.125_wp*rgdriz*rotstep* &
           ( -p0*q3+p1*q2-          &
              p2*q1+p3*q0 )
     cs=Cos(zetaz)
     sn=Sin(zetaz)

     qn1(0)=cs*q0-sn*q3
     qn1(1)=cs*q1+sn*q2
     qn1(2)=cs*q2-sn*q1
     qn1(3)=cs*q3+sn*q0

     pq2(0)=cs*p0-sn*p3
     pq2(1)=cs*p1+sn*p2
     pq2(2)=cs*p2-sn*p1
     pq2(3)=cs*p3+sn*p0

     zetay=0.125_wp*rgdriy*rotstep*        &
           ( -pq2(0)*qn1(2)-pq2(1)*qn1(3)+ &
              pq2(2)*qn1(0)+pq2(3)*qn1(1) )
     cs=Cos(zetay)
     sn=Sin(zetay)

     qn2(0)=cs*qn1(0)-sn*qn1(2)
     qn2(1)=cs*qn1(1)-sn*qn1(3)
     qn2(2)=cs*qn1(2)+sn*qn1(0)
     qn2(3)=cs*qn1(3)+sn*qn1(1)

     pq3(0)=cs*pq2(0)-sn*pq2(2)
     pq3(1)=cs*pq2(1)-sn*pq2(3)
     pq3(2)=cs*pq2(2)+sn*pq2(0)
     pq3(3)=cs*pq2(3)+sn*pq2(1)

     zetax=0.250_wp*rgdrix*rotstep*        &
           ( -pq3(0)*qn2(1)+pq3(1)*qn2(0)+ &
              pq3(2)*qn2(3)-pq3(3)*qn2(2) )
     cs=Cos(zetax)
     sn=Sin(zetax)

     qn3(0)=cs*qn2(0)-sn*qn2(1)
     qn3(1)=cs*qn2(1)+sn*qn2(0)
     qn3(2)=cs*qn2(2)+sn*qn2(3)
     qn3(3)=cs*qn2(3)-sn*qn2(2)

     pq4(0)=cs*pq3(0)-sn*pq3(1)
     pq4(1)=cs*pq3(1)+sn*pq3(0)
     pq4(2)=cs*pq3(2)+sn*pq3(3)
     pq4(3)=cs*pq3(3)-sn*pq3(2)

     zetay=0.125_wp*rgdriy*rotstep*        &
           ( -pq4(0)*qn3(2)-pq4(1)*qn3(3)+ &
              pq4(2)*qn3(0)+pq4(3)*qn3(1) )
     cs=Cos(zetay)
     sn=Sin(zetay)

     qn2(0)=cs*qn3(0)-sn*qn3(2)
     qn2(1)=cs*qn3(1)-sn*qn3(3)
     qn2(2)=cs*qn3(2)+sn*qn3(0)
     qn2(3)=cs*qn3(3)+sn*qn3(1)

     pq3(0)=cs*pq4(0)-sn*pq4(2)
     pq3(1)=cs*pq4(1)-sn*pq4(3)
     pq3(2)=cs*pq4(2)+sn*pq4(0)
     pq3(3)=cs*pq4(3)+sn*pq4(1)

     zetaz=0.125_wp*rgdriz*rotstep*        &
           ( -pq3(0)*qn2(3)+pq3(1)*qn2(2)- &
              pq3(2)*qn2(1)+pq3(3)*qn2(0) )
     cs=Cos(zetaz)
     sn=Sin(zetaz)

     q0=cs*qn2(0)-sn*qn2(3)
     q1=cs*qn2(1)+sn*qn2(2)
     q2=cs*qn2(2)-sn*qn2(1)
     q3=cs*qn2(3)+sn*qn2(0)

     p0=cs*pq3(0)-sn*pq3(3)
     p1=cs*pq3(1)+sn*pq3(2)
     p2=cs*pq3(2)-sn*pq3(1)
     p3=cs*pq3(3)+sn*pq3(0)
  End Do

End Subroutine no_squish

subroutine q_update                                       &
           (tstep,oxp,oyp,ozp,oxq,oyq,ozq,mxquat,quattol, &
           q0,q1,q2,q3,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to update the quaternions in LFV scope
!
! copyright - daresbury laboratory
! author    - t.forester october 1993
! adapted   - i.t.todorov october 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90

  Implicit None

  Logical,           Intent( InOut ) :: safe
  Integer,           Intent( In    ) :: mxquat
  Real( Kind = wp ), Intent( In    ) :: tstep,oxp,oyp,ozp,oxq,oyq,ozq,quattol
  Real( Kind = wp ), Intent( InOut ) :: q0,q1,q2,q3

  Integer           :: itq
  Real( Kind = wp ) :: eps,rnorm,qn0,qn1,qn2,qn3, &
                       qn0a,qn1a,qn2a,qn3a,qn0b,qn1b,qn2b,qn3b

! first estimate of new quaternions (laboratory frame)

  qn0=q0+(-q1*oxp - q2*oyp - q3*ozp)*tstep*0.5_wp
  qn1=q1+( q0*oxp - q3*oyp + q2*ozp)*tstep*0.5_wp
  qn2=q2+( q3*oxp + q0*oyp - q1*ozp)*tstep*0.5_wp
  qn3=q3+(-q2*oxp + q1*oyp + q0*ozp)*tstep*0.5_wp

! first iteration of new quaternions (lab fixed)

  qn0b=0.0_wp ; qn1b=0.0_wp ; qn2b=0.0_wp ; qn3b=0.0_wp

  itq=0 ; eps=2.0_wp*quattol
  Do While (itq < mxquat .and. eps > quattol)
     qn0a=0.5_wp*(-q1 *oxp - q2 *oyp - q3 *ozp) + &
          0.5_wp*(-qn1*oxq - qn2*oyq - qn3*ozq)
     qn1a=0.5_wp*( q0 *oxp - q3 *oyp + q2 *ozp) + &
          0.5_wp*( qn0*oxq - qn3*oyq + qn2*ozq)
     qn2a=0.5_wp*( q3 *oxp + q0 *oyp - q1 *ozp) + &
          0.5_wp*( qn3*oxq + qn0*oyq - qn1*ozq)
     qn3a=0.5_wp*(-q2 *oxp + q1 *oyp + q0 *ozp) + &
          0.5_wp*(-qn2*oxq + qn1*oyq + qn0*ozq)

     qn0=q0+0.5_wp*qn0a*tstep
     qn1=q1+0.5_wp*qn1a*tstep
     qn2=q2+0.5_wp*qn2a*tstep
     qn3=q3+0.5_wp*qn3a*tstep

     rnorm=1.0_wp/Sqrt(qn0**2+qn1**2+qn2**2+qn3**2)

     qn0=qn0*rnorm
     qn1=qn1*rnorm
     qn2=qn2*rnorm
     qn3=qn3*rnorm

! convergence test

     eps=Sqrt(((qn0a-qn0b)**2+(qn1a-qn1b)**2+(qn2a-qn2b)**2+(qn3a-qn3b)**2)*tstep**2)

     qn0b=qn0a
     qn1b=qn1a
     qn2b=qn2a
     qn3b=qn3a

     itq=itq+1
  End Do

! store new quaternions

  q0=qn0
  q1=qn1
  q2=qn2
  q3=qn3

! Check safety

  If (itq == mxquat) safe=.false.

End Subroutine q_update
