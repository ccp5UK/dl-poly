Subroutine nst_h1_scl &
           (sw,tstep,degfre,degrot,pmass,chit,volm,press,  &
           iso,ten,h_z,strext,strcon,stress,strcom,        &
           vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz,eta,strkin,strknf,strknt,engke)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to integrate and apply NsT barostat
! when singled RBs are present
!
! sw=1 coupling to NVT thermostat for nst_m ensemble and
!                                     additional scaling factor
!
! sw=0 coupling to NVT thermostat for nst_h ensemble and
!                                     no additional scaling factor
!
! iso=0 fully anisotropic barostat
! iso=1 semi-isotropic barostat to constant normal pressure & surface area
! iso=2 semi-isotropic barostat to constant normal pressure & surface tension
!                               or with orthorhombic constraints (ten=0.0_wp)
! iso=3 semi-isotropic barostat with semi-orthorhombic constraints
!
! reference: Mitsunori Ikeguchi, J. Comp. Chem. (2004), 25, p529
!
! copyright - daresbury laboratory
! author    - i.t.todorov december 2012
! contrib   - a.m.elena december 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module,       Only : nfree,lstfre
  Use rigid_bodies_module, Only : ntrgd
  Use kinetic_module,      Only : kinstresf,kinstrest

  Implicit None

  Integer,           Intent( In    ) :: sw,iso
  Integer(Kind=ip),  Intent( In    ) :: degfre,degrot

  Real( Kind = wp ), Intent( In    ) :: tstep,pmass,chit,volm,press,ten,h_z
  Real( Kind = wp ), Intent( In    ) :: strext(1:9),strcon(1:9),stress(1:9),strcom(1:9)
  Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)
  Real( Kind = wp ), Intent( InOut ) :: rgdvxx(1:mxrgd),rgdvyy(1:mxrgd),rgdvzz(1:mxrgd)
  Real( Kind = wp ), Intent( InOut ) :: eta(1:9)
  Real( Kind = wp ), Intent(   Out ) :: strkin(1:9),strknf(1:9),strknt(1:9)
  Real( Kind = wp ), Intent(   Out ) :: engke


  Logical,     Save :: newjob = .true.

  Integer           :: i,j,irgd

  Real( Kind = wp ) :: a1,a2,a3,a5,a6,a9,b1,b2,b3,b5,b6,b9, vxt,vyt,vzt

! initialise factor for Nose-Hoover ensembles

  Real( Kind = wp ), Save :: rf, factor

! uni is the diagonal unit matrix

  Real( Kind = wp ), Parameter :: &
  uni(1:9) = (/ 1.0_wp,0.0_wp,0.0_wp, 0.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp /)

  Real( Kind = wp ) :: hstep,qstep


! Initialise factor and 1/Nf for Nose-Hoover ensembles

  If (newjob) Then
     newjob = .false.

     factor = 0.0_wp
     rf = 0.0_wp
     If (sw == 1) rf=1.0_wp/Real(degfre-degrot,wp)
  End If

! timestep derivatives

  hstep=0.5_wp*tstep
  qstep=0.5_wp*hstep

! thermostat eta to 1/4*tstep

  eta=eta*Exp(-qstep*chit)

! calculate kinetic contributions to stress tensor

  Call kinstresf(vxx,vyy,vzz,strknf)
  Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

  strkin=strknf+strknt

! kinetic energy

  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! barostat eta to 1/2*tstep

  If (sw == 1) factor = 2.0_wp*engke*rf

! split anisotropic from semi-isotropic barostats (iso=0,1,2,3)

  If (iso == 0) Then
     eta=eta + hstep*(strcom+strcon+stress+strkin + factor*uni - (press*uni+strext)*volm)/pmass
  Else
     If      (iso == 2) Then
        eta(1)=eta(1) + hstep*(strcom(1)+strcon(1)+stress(1)+strkin(1) + factor - (press+strext(1)-ten/h_z)*volm)/pmass
        eta(5)=eta(5) + hstep*(strcom(5)+strcon(5)+stress(5)+strkin(5) + factor - (press+strext(5)-ten/h_z)*volm)/pmass
     Else If (iso == 3) Then
        eta(1)=0.5_wp*(eta(1)+eta(5)) + hstep*( 0.5_wp*                                            &
               (strcom(1)+strcon(1)+stress(1)+strkin(1)+strcom(5)+strcon(5)+stress(5)+strkin(5)) + &
               factor - (press+0.5_wp*(strext(1)+strext(5))-ten/h_z)*volm ) / pmass
        eta(5)=eta(1)
     End If
     eta(9)=eta(9) + hstep*(strcom(9)+strcon(9)+stress(9)+strkin(9) + factor - (press+strext(9))*volm)/pmass
  End If

! thermostat eta to 2/4*tstep

  eta=eta*Exp(-qstep*chit)

! barostat the velocities to full 1*tstep
! second order taylor expansion of Exp[-tstep*(eta+factor*I)],
! where I is the unit tensor
! factor = Tr(eta)/Nf if sw=1, where Nf is degfre-degrot,
! else if sw=0 then factor=0, by default

  If (sw == 1) factor = (eta(1)+eta(5)+eta(9))*rf

  a1 = -tstep*(eta(1)+factor)
  a2 = -tstep*eta(2)
  a3 = -tstep*eta(3)
  a5 = -tstep*(eta(5)+factor)
  a6 = -tstep*eta(6)
  a9 = -tstep*(eta(9)+factor)

  b1 = (a1*a1 + a2*a2 + a3*a3)*0.5_wp + a1 + 1.0_wp
  b2 = (a1*a2 + a2*a5 + a3*a6)*0.5_wp + a2
  b3 = (a1*a3 + a2*a6 + a3*a9)*0.5_wp + a3
  b5 = (a2*a2 + a5*a5 + a6*a6)*0.5_wp + a5 + 1.0_wp
  b6 = (a2*a3 + a5*a6 + a6*a9)*0.5_wp + a6
  b9 = (a3*a3 + a6*a6 + a9*a9)*0.5_wp + a9 + 1.0_wp

  Do j=1,nfree
     i=lstfre(j)

     vxt=vxx(i)
     vyt=vyy(i)
     vzt=vzz(i)

     vxx(i) = b1*vxt + b2*vyt + b3*vzt
     vyy(i) = b2*vxt + b5*vyt + b6*vzt
     vzz(i) = b3*vxt + b6*vyt + b9*vzt
  End Do

  Do irgd=1,ntrgd
     vxt=rgdvxx(irgd)
     vyt=rgdvyy(irgd)
     vzt=rgdvzz(irgd)

     rgdvxx(irgd) = b1*vxt + b2*vyt + b3*vzt
     rgdvyy(irgd) = b2*vxt + b5*vyt + b6*vzt
     rgdvzz(irgd) = b3*vxt + b6*vyt + b9*vzt
  End Do

! thermostat eta to 2/4*tstep

! calculate kinetic contributions to stress tensor

  Call kinstresf(vxx,vyy,vzz,strknf)
  Call kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

  strkin=strknf+strknt

! kinetic energy

  engke=0.5_wp*(strkin(1)+strkin(5)+strkin(9))

! thermostat eta to 3/4*tstep

  eta=eta*Exp(-qstep*chit)

! barostat eta to full (2/2)*tstep

  If (sw == 1) factor = 2.0_wp*engke*rf

! split anisotropic from semi-isotropic barostats (iso=0,1,2,3)

  If (iso == 0) Then
     eta=eta + hstep*(strcom+strcon+stress+strkin + factor*uni - (press*uni+strext)*volm)/pmass
  Else
     If      (iso == 2) Then
        eta(1)=eta(1) + hstep*(strcom(1)+strcon(1)+stress(1)+strkin(1) + factor - (press+strext(1)-ten/h_z)*volm)/pmass
        eta(5)=eta(5) + hstep*(strcom(5)+strcon(5)+stress(5)+strkin(5) + factor - (press+strext(5)-ten/h_z)*volm)/pmass
     Else If (iso == 3) Then
        eta(1)=0.5_wp*(eta(1)+eta(5)) + hstep*( 0.5_wp*                                            &
               (strcom(1)+strcon(1)+stress(1)+strkin(1)+strcom(5)+strcon(5)+stress(5)+strkin(5)) + &
               factor - (press+0.5_wp*(strext(1)+strext(5))-ten/h_z)*volm ) / pmass
        eta(5)=eta(1)
     End If
     eta(9)=eta(9) + hstep*(strcom(9)+strcon(9)+stress(9)+strkin(9) + factor - (press+strext(9))*volm)/pmass
  End If

! thermostat eta to full (4/4)*tstep

  eta=eta*Exp(-qstep*chit)

End Subroutine nst_h1_scl
