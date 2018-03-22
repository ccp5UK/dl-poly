Module coul_spole
  Use kinds, Only : wp
  Use comms,  Only : comms_type
  Use setup_module, Only : mxlist, r4pie0, zero_plus, mxgele, nrite
  Use configuration, Only : natms,ltg,list,chge,fxx,fyy,fzz
  Implicit None

  Private

  Public :: coul_fscp_forces, coul_rfp_forces, coul_cp_forces, coul_dddp_forces

  Contains

  Subroutine coul_fscp_forces &
             (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system assuming a force shifted coulombic potential
  !
  ! U is proportional to ( 1/r + aa*r  + bb ) such that dU(rcut)/dr = 0
  ! therefore aa = 1/(rcut)**2 and U(rcut) = 0 therefore bb = -2/(rcut)
  !
  ! Note: FS potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  !
  ! copyright - daresbury laboratory
  ! author    - t.forester october 1995
  ! amended   - i.t.todorov november 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( comms_type ),                       Intent( In    ) :: comm

    Logical,           Save :: newjob = .true. , damp
    Real( Kind = wp ), Save :: drewd  = 0.0_wp , &
                               rdrewd = 0.0_wp , &
                               aa     = 0.0_wp , &
                               bb     = 0.0_wp

    Integer           :: fail,k,idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,ppp,egamma, &
                         fix,fiy,fiz,fx,fy,fz,            &
                         vk0,vk1,vk2,gk0,gk1,gk2,t1,t2,   &
                         strs1,strs2,strs3,strs5,strs6,strs9

    Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

    If (newjob) Then
       newjob = .false.

       If (alpha > zero_plus) Then
          damp = .true.
       Else
          damp = .false.
       End If

       If (damp) Then

  ! interpolation interval

          drewd = rcut/Real(mxgele-4,wp)

  ! reciprocal of interpolation interval

          rdrewd = 1.0_wp/drewd

          fail=0
          Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'coul_fscp_forces allocation failure, idnode: ', comm%idnode
             Call error(0)
          End If

  ! generate error function complement tables for ewald sum

          Call erfcgen(rcut,alpha,mxgele,erc,fer)

  ! set force and potential shifting parameters (screened terms)

          aa =   fer(mxgele-4)*rcut
          bb = -(erc(mxgele-4)+aa*rcut)

       Else

  ! set force and potential shifting parameters (screened terms)

          aa =  1.0_wp/rcut**2
          bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)

       End If
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

  ! Squared distance

             rsq=rrr**2

  ! calculate forces

             If (damp) Then
                k   = Int(rrr*rdrewd)
                ppp = rrr*rdrewd - Real(k,wp)

  ! calculate forces using 3pt interpolation

                gk0 = fer(k) ; If (k == 0) gk0 = gk0*rrr
                gk1 = fer(k+1)
                gk2 = fer(k+2)

                t1 = gk0 + (gk1 - gk0)*ppp
                t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                egamma = ((t1 + (t2-t1)*ppp*0.5_wp) - aa/rrr)*chgprd
             Else
                egamma=chgprd*(1.0_wp/rsq - aa)/rrr
             End If

             fx = egamma*xxt(m)
             fy = egamma*yyt(m)
             fz = egamma*zzt(m)

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

                If (damp) Then

  ! calculate interaction energy using 3-point interpolation

                   vk0 = erc(k)
                   vk1 = erc(k+1)
                   vk2 = erc(k+2)

                   t1 = vk0 + (vk1 - vk0)*ppp
                   t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                   engcpe = engcpe + ((t1 + (t2-t1)*ppp*0.5_wp) + aa*rrr + bb)*chgprd

                Else

                   engcpe = engcpe + chgprd*(1.0_wp/rrr + aa*rrr + bb)

                End If

                vircpe = vircpe - egamma*rsq

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

  End Subroutine coul_fscp_forces

  Subroutine coul_rfp_forces &
             (iatm,rcut,alpha,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress,comm)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system using reaction field potential (correcting for
  ! the existence of a dipole moment outside rcut)
  !
  ! Note: RF potential can be generalised (R1) by using a damping function
  ! as used for damping the real space coulombic interaction in the
  ! standard Ewald summation.  This generalisation applies when alpha > 0.
  !
  ! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
  ! R2: M. Neumann, J. Chem. Phys., 82 (12), 5663, (1985)
  !
  ! copyright - daresbury laboratory
  ! author    - t.forester february 1995
  ! amended   - i.t.todorov november 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Implicit None

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
    Type( comms_type ),                       Intent( In    ) :: comm

    Logical,           Save :: newjob = .true. , damp
    Real( Kind = wp ), Save :: b0     = 0.0_wp , &
                               rfld0  = 0.0_wp , &
                               rfld1  = 0.0_wp , &
                               rfld2  = 0.0_wp , &
                               drewd  = 0.0_wp , &
                               rdrewd = 0.0_wp , &
                               aa     = 0.0_wp , &
                               bb     = 0.0_wp , &
                               rcsq   = 0.0_wp


    Integer           :: fail,k,idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,ppp,egamma, &
                         fix,fiy,fiz,fx,fy,fz,            &
                         vk0,vk1,vk2,gk0,gk1,gk2,t1,t2,   &
                         strs1,strs2,strs3,strs5,strs6,strs9

    Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

    If (newjob) Then
       newjob = .false.

       If (alpha > zero_plus) Then
          damp = .true.
       Else
          damp = .false.
       End If

  ! reaction field terms

       b0    = 2.0_wp*(epsq - 1.0_wp)/(2.0_wp*epsq + 1.0_wp)
       rfld0 = b0/rcut**3
       rfld1 = (1.0_wp + 0.5_wp*b0)/rcut
       rfld2 = 0.5_wp*rfld0

       If (damp) Then

  ! interpolation interval

          drewd = rcut/Real(mxgele-4,wp)

  ! reciprocal of interpolation interval

          rdrewd = 1.0_wp/drewd

          fail=0
          Allocate (erc(0:mxgele),fer(0:mxgele), Stat=fail)
          If (fail > 0) Then
             Write(nrite,'(/,1x,a,i0)') 'coul_fscp_forces allocation failure, idnode: ', comm%idnode
             Call error(0)
          End If

  ! generate error function complement tables for ewald sum

          Call erfcgen(rcut,alpha,mxgele,erc,fer)

  ! set force and potential shifting parameters (screened terms)

          aa =   fer(mxgele-4)*rcut
          bb = -(erc(mxgele-4)+aa*rcut)

  ! Cutoff squared

          rcsq = rcut**2

       End If
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

  ! Squared distance

             rsq=rrr**2

  ! calculate forces

             If (damp) Then
                k   = Int(rrr*rdrewd)
                ppp = rrr*rdrewd - Real(k,wp)

  ! calculate forces using 3pt interpolation

                gk0 = fer(k) ; If (k == 0) gk0 = gk0*rrr
                gk1 = fer(k+1)
                gk2 = fer(k+2)

                t1 = gk0 + (gk1 - gk0)*ppp
                t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                egamma = ((t1 + (t2-t1)*ppp*0.5_wp) - aa/rrr - rfld0)*chgprd
             Else
                egamma=chgprd*(1.0_wp/rsq/rrr - rfld0)
             End If

             fx = egamma*xxt(m)
             fy = egamma*yyt(m)
             fz = egamma*zzt(m)

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

                If (damp) Then

  ! calculate interaction energy using 3-point interpolation

                   vk0 = erc(k)
                   vk1 = erc(k+1)
                   vk2 = erc(k+2)

                   t1 = vk0 + (vk1 - vk0)*ppp
                   t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                   engcpe = engcpe + ((t1 + (t2-t1)*ppp*0.5_wp) + aa*rrr + bb + &
                                      rfld2*(rsq-rcsq))*chgprd

                Else

                   engcpe = engcpe + chgprd*(1.0_wp/rrr + rfld2*rsq - rfld1)

                End If

                vircpe = vircpe - egamma*rsq

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

  End Subroutine coul_rfp_forces

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

  Subroutine coul_dddp_forces &
             (iatm,rcut,epsq,xxt,yyt,zzt,rrt,engcpe,vircpe,stress)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 subroutine for calculating coulombic energy and force terms
  ! in a periodic system assuming a distance dependant dielectric
  ! `constant'
  !
  ! copyright - daresbury laboratory
  ! author    - t.forester april 1993
  ! amended   - i.t.todorov november 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,                                  Intent( In    ) :: iatm
    Real( Kind = wp ),                        Intent( In    ) :: rcut,epsq
    Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xxt,yyt,zzt,rrt
    Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
    Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

    Integer           :: idi,jatm,m

    Real( Kind = wp ) :: chgea,chgprd,rrr,rsq,coul,fcoul, &
                         fix,fiy,fiz,fx,fy,fz,            &
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

  ! Squared distance

             rsq=rrr**2

  ! calculate forces

             coul = chgprd/rsq
             fcoul = 2.0_wp*coul/rsq

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
End Module coul_spole