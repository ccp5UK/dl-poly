Subroutine coul_chrm_forces(iatm,epsq,xxt,yyt,zzt,rrt,engcpe_ch,vircpe_ch,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! for CHARMM model intra-core-shell interactions
!
! S(r)=1-(1+u/2).exp(-u) ; u=u(r)=r_ij.(a_i+a_j)/(p_i.p_j)^(1/6)
!
! S'(r)=(1+u).u'.exp(-u)/2 ; |r|'=r/|r|
!
! Uchrm(r_ij) =  S(r_ij)*U(r_ij) ; F(r_ij) = -U'(r_ij)
! Fchrm(r_ij) = -Uchrm'(r_ij) = S(r_ij)*F(r_ij) - S'(r_ij)*U(r_ij)
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
  Use config_module, Only : natms,ltg,list,chge,fxx,fyy,fzz
  Use mpoles_module, Only : plratm,dmpatm

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: epsq
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe_ch,vircpe_ch
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Integer           :: limit,idi,jatm,m

  Real( Kind = wp ) :: chgea,chgprd,rrr,r_r,coul,fcoul,tmp, &
                       plra,plrprd,dmpa,dmpsum,u,scr,spr,   &
                       fix,fiy,fiz,fx,fy,fz,                &
                       strs1,strs2,strs3,strs5,strs6,strs9

! initialise potential energy and virial

  engcpe_ch=0.0_wp
  vircpe_ch=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity of iatm

  idi=ltg(iatm)

! charge, inverse polarisability and dumping

  chgea = chge(iatm)
  plra  = plratm(iatm)
  dmpa  = dmpatm(iatm)

! scale main charge

  chgea = chgea*r4pie0/epsq

! load forces

  fix=fxx(iatm)
  fiy=fyy(iatm)
  fiz=fzz(iatm)

! start of primary loop for forces evaluation

! Get list limit

  limit=list(-4,iatm)-list(0,iatm)

  Do m=1,limit

! interatomic distance and derivatives

     rrr=rrt(m)
     r_r=1.0_wp/rrr

! atomic index, charge & inverse polarisability products
! and total inter-atomic summed dumping

     jatm=list(list(0,iatm)+m,iatm)
     chgprd=chgea*chge(jatm)
     plrprd=plra*plratm(jatm)
     dmpsum=dmpa+dmpatm(jatm)

     u = (dmpsum*plrprd**6) * rrr         ! dimensionless

     tmp = Exp(-u)
     scr = 1.0_wp-(1.0_wp+0.5_wp*u) * tmp ! S(r)
     spr = (1.0_wp+u) * tmp * u           ! S'(r).r

! calculate forces

     coul  = scr*chgprd*r_r
     tmp   = (scr-spr)*chgprd             ! used later for the virial
     fcoul = tmp*r_r**3

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

! calculate potential energy and virial

        engcpe_ch = engcpe_ch + coul
        vircpe_ch = vircpe_ch - tmp*r_r

! calculate stress tensor

        strs1 = strs1 + xxt(m)*fx
        strs2 = strs2 + xxt(m)*fy
        strs3 = strs3 + xxt(m)*fz
        strs5 = strs5 + yyt(m)*fy
        strs6 = strs6 + yyt(m)*fz
        strs9 = strs9 + zzt(m)*fz

     End If

  End Do

! load back forces

  fxx(iatm)=fix
  fyy(iatm)=fiy
  fzz(iatm)=fiz

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

End Subroutine coul_chrm_forces
