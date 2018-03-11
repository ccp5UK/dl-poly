Subroutine coul_cp_forces &
           (iatm,rcut,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using 1/r potential with no truncation or damping
!
! copyright - daresbury laboratory
! author    - t.forester february 1993
! amended   - i.t.todorov november 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup_module
  Use configuration, Only : natms,ltg,list,chge,fxx,fyy,fzz

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Integer           :: idi,jatm,m

  Real( Kind = wp ) :: chgea,chgprd,rrr,coul,fcoul, &
                       fix,fiy,fiz,fx,fy,fz,        &
                       strs1,strs2,strs3,strs5,strs6,strs9

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

! ignore interaction if the charge is zero

  chgea = chge(iatm)

  If (Abs(chgea) > zero_plus) Then

     chgea = chgea*r4pie0/epsq

! load forces

     fix=fxx(iatm)
     fiy=fyy(iatm)
     fiz=fzz(iatm)

! start of primary loop for forces evaluation

     Do m=1,list(0,iatm)

! atomic index and charge

        jatm=list(m,iatm)
        chgprd=chge(jatm)

! interatomic distance

        rrr=rrt(m)

! interaction validity and truncation of potential

        If (Abs(chgprd) > zero_plus .and. rrr < rcut) Then

! charge product

           chgprd=chgprd*chgea

! calculate forces

           coul = chgprd/rrr
           fcoul = coul/rrr**2

           fx = fcoul*xxt(m)
           fy = fcoul*yyt(m)
           fz = fcoul*zzt(m)

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

              strs1 = strs1 + xxt(m)*fx
              strs2 = strs2 + xxt(m)*fy
              strs3 = strs3 + xxt(m)*fz
              strs5 = strs5 + yyt(m)*fy
              strs6 = strs6 + yyt(m)*fz
              strs9 = strs9 + zzt(m)*fz

           End If

        End If

     End Do

! load back forces

     fxx(iatm)=fix
     fyy(iatm)=fiy
     fzz(iatm)=fiz

! virial

     vircpe = -engcpe

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

End Subroutine coul_cp_forces
