Subroutine coul_dddp_forces &
           (iatm,rcut,epsq,xdf,ydf,zdf,rsqdf,engcpe,vircpe,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system assuming a distance dependant dielectric
! `constant'
!
! copyright - daresbury laboratory
! author    - t.forester april 1993
! amended   - i.t.todorov may 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use config_module, Only : natms,ltg,list,chge,fxx,fyy,fzz
  Use setup_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: rcsq

  Integer           :: idi,jatm,m

  Real( Kind = wp ) :: chgea,chgprd,rsq,coul,fcoul, &
                       fix,fiy,fiz,fx,fy,fz,        &
                       strs1,strs2,strs3,strs5,strs6,strs9


! set cutoff condition for pair forces

  If (newjob) Then
     newjob = .false.

     rcsq = rcut**2
  End If

! initialise potential energy and virial

  engcpe=0.0_wp
  vircpe=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity of iatm

  idi=ltg(iatm)

! start of primary loop for forces evaluation

  chgea = chge(iatm)

! ignore interaction if the charge is zero

  If (Abs(chgea) > zero_plus) Then

     chgea = chgea*r4pie0/epsq

! load forces

     fix=fxx(iatm)
     fiy=fyy(iatm)
     fiz=fzz(iatm)

     Do m=1,list(0,iatm)

! atomic index

        jatm=list(m,iatm)
        chgprd=chge(jatm)

! ignore interaction if the charge is zero

        If (Abs(chgprd) > zero_plus) Then

! charge product

           chgprd=chgprd*chgea

! calculate interatomic distance

           rsq=rsqdf(m)

! apply truncation of potential

           If (rsq < rcsq) Then

! calculate forces

              coul = chgprd/rsq
              fcoul = 2.0_wp*coul/rsq

              fx = fcoul*xdf(m)
              fy = fcoul*ydf(m)
              fz = fcoul*zdf(m)

              fix=fix+fx
              fiy=fiy+fy
              fiz=fiz+fz

              If (jatm <= natms) Then

                 fxx(jatm)=fxx(jatm)-fx
                 fyy(jatm)=fyy(jatm)-fy
                 fzz(jatm)=fzz(jatm)-fz

              End If

              If (jatm <= natms .or. idi < ltg(jatm)) Then

! calculate potential energy

                 engcpe = engcpe + coul

! calculate stress tensor

                 strs1 = strs1 + xdf(m)*fx
                 strs2 = strs2 + xdf(m)*fy
                 strs3 = strs3 + xdf(m)*fz
                 strs5 = strs5 + ydf(m)*fy
                 strs6 = strs6 + ydf(m)*fz
                 strs9 = strs9 + zdf(m)*fz

              End If

           End If

        End If

     End Do

! load back forces

     fxx(iatm)=fix
     fyy(iatm)=fiy
     fzz(iatm)=fiz

! virial

     vircpe = -2.0_wp*engcpe

! complete stress tensor

     stress(1) = stress(1) + strs1
     stress(2) = stress(2) + strs2
     stress(3) = stress(3) + strs3
     stress(4) = stress(4) + strs2
     stress(5) = stress(5) + strs5
     stress(6) = stress(6) + strs6
     stress(7) = stress(7) + strs3
     stress(8) = stress(8) + strs6
     stress(9) = stress(9) + strs9

  End If

End Subroutine coul_dddp_forces
