Subroutine coul_fscp_forces &
           (iatm,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf,engcpe,vircpe,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system assuming a force shifted coulombic potential
!
! U is proportional to ( 1/r + aa*r  + bb ) such that dU(rcut)/dr = 0
! therefore aa = 1/(rcut)**2 and U(rcut) = 0 therefore bb = -2/(rcut)
!
! Note: FS potential can be generalised (R1) by using a damping function
! as used for damping the real space coulombic interation in the
! standard Ewald summation.  This generalisation applies when alpha > 0.
!
! R1: C.J. Fennell and J.D. Gezelter J. Chem. Phys. 124, 234104 (2006)
!
! copyright - daresbury laboratory
! author    - t.forester october 1995
! amended   - i.t.todorov march 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use config_module, Only : natms,ltg,list,chge,fxx,fyy,fzz
  Use setup_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical,           Save :: newjob = .true. , damp
  Real( Kind = wp ), Save :: drewd  = 0.0_wp , &
                             rdrewd = 0.0_wp , &
                             rcsq   = 0.0_wp , &
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

! set cutoff condition for pair forces

     rcsq = rcut**2

     If (damp) Then

! interpolation interval

        drewd = rcut/Real(mxgrid-4,wp)

! reciprocal of interpolation interval

        rdrewd = 1.0_wp/drewd

        fail=0
        Allocate (erc(1:mxgrid),fer(1:mxgrid), Stat=fail)
        If (fail > 0) Then
           Write(nrite,'(/,1x,a,i0)') 'coul_fscp_forces allocation failure, idnode: ', idnode
           Call error(0)
        End If

! generate error function complement tables for ewald sum

        Call erfcgen(rcut,alpha,mxgrid,erc,fer)

! set force and potential shifting parameters (screened terms)

        aa =   fer(mxgrid-4)*rcut
        bb = -(erc(mxgrid-4)+aa*rcut)

     Else

! set force and potential shifting parameters (screened terms)

        aa =  1.0_wp/rcsq
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

              rrr = Sqrt(rsq)
              If (damp) Then
                 k   = Int(rrr*rdrewd)
                 ppp = rrr*rdrewd - Real(k,wp)

! calculate forces using 3pt interpolation

                 gk0 = fer(k)
                 gk1 = fer(k+1)
                 gk2 = fer(k+2)

                 t1 = gk0 + (gk1 - gk0)*ppp
                 t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                 egamma = ((t1 + (t2-t1)*ppp*0.5_wp) - aa/rrr)*chgprd
              Else
                 egamma=chgprd*(1.0_wp/rsq - aa)/rrr
              End If

              fx = egamma*xdf(m)
              fy = egamma*ydf(m)
              fz = egamma*zdf(m)

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
