Subroutine ewald_real_forces_helper                 &
           (iatm,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf, &
           engcpe_rl,vircpe_rl,stress,erc,fer,drewd,&
           rdrewd,rcsq)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using ewald's method
!
! Note: real space terms
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov february 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Centre for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
  Use config_module, Only : natms,ltg,list,chge,fxx,fyy,fzz

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mx_two ), Intent( In    ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe_rl,vircpe_rl
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress
  Real( Kind = wp ), Dimension( 1:mxgrid),  Intent( In    ) :: erc,fer
  Real( Kind = wp ),                        Intent( In    ) :: drewd,rdrewd,rcsq
!  Logical,           Save :: newjob = .true.
!  Real( Kind = wp ), Save :: drewd,rdrewd,rcsq

  Integer           :: fail,m,idi,jatm,k, mstride
  Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,ppp,egamma, &
                       fix,fiy,fiz,fx,fy,fz,            &
                       vk0,vk1,vk2,gk0,gk1,gk2,t1,t2,   &
                       strs1,strs2,strs3,strs5,strs6,strs9


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

     fix=0.0_wp!fxx(iatm)
     fiy=0.0_wp!fyy(iatm)
     fiz=0.0_wp!fzz(iatm)

!$OMP DO
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

              rrr = Sqrt(rsq)
              k   = Int(rrr*rdrewd)
              ppp = rrr*rdrewd - Real(k,wp)

! calculate forces using 3pt interpolation

              gk0 = fer(k)
              gk1 = fer(k+1)
              gk2 = fer(k+2)

              t1 = gk0 + (gk1 - gk0)*ppp
              t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

              egamma = (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

! calculate forces

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

! calculate interaction energy using 3-point interpolation

                 vk0 = erc(k)
                 vk1 = erc(k+1)
                 vk2 = erc(k+2)

                 t1 = vk0 + (vk1 - vk0)*ppp
                 t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                 engcpe_rl = engcpe_rl + (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

! calculate virial

                 vircpe_rl = vircpe_rl - egamma*rsq

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
!        m = m + mstride
     End Do
!$OMP END DO NOWAIT

! load back forces

!$OMP ATOMIC
     fxx(iatm)=fxx(iatm)+fix
!$OMP ATOMIC
     fyy(iatm)=fyy(iatm)+fiy
!$OMP ATOMIC
     fzz(iatm)=fzz(iatm)+fiz

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
!!$OMP MASTER
!  Call stop_timing_ewald_real_forces()
!!$OMP END MASTER
End Subroutine ewald_real_forces_helper

Subroutine ewald_real_forces                        &
           (iatm,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf, &
           engcpe_rl,vircpe_rl,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating coulombic energy and force terms
! in a periodic system using ewald's method
!
! Note: real space terms
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov february 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode
  Use setup_module
  Use config_module, Only : natms,ltg,list,chge,fxx,fyy,fzz

#ifdef COMPILE_CUDA
  Use dl_poly_cuda_module
#endif

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha,epsq
  Real( Kind = wp ), Dimension( 1:mx_two ), Intent( In    ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe_rl,vircpe_rl
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: drewd,rdrewd,rcsq

  Integer           :: fail,m,idi,jatm,k
  Real( Kind = wp ) :: chgea,chgprd,rsq,rrr,ppp,egamma, &
                       fix,fiy,fiz,fx,fy,fz,            &
                       vk0,vk1,vk2,gk0,gk1,gk2,t1,t2,   &
                       strs1,strs2,strs3,strs5,strs6,strs9

  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer

#ifdef COMPILE_CUDA
  Call start_timing_ewald_real_forces()
#endif

  If (newjob) Then
     newjob = .false.

     fail=0
     Allocate (erc(1:mxgrid),fer(1:mxgrid), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_real_forces allocation failure, node: ', idnode
        Call error(0)
     End If
! set cutoff condition for pair forces

     rcsq = rcut**2

! interpolation interval

     drewd = rcut/Real(mxgrid-4,wp)

! reciprocal of interpolation interval

     rdrewd = 1.0_wp/drewd

! generate error function complement tables for ewald sum

     Call erfcgen(rcut,alpha,mxgrid,erc,fer)
  End If


! initialise potential energy and virial

  engcpe_rl=0.0_wp
  vircpe_rl=0.0_wp

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

     fix=0.0_wp!fxx(iatm)
     fiy=0.0_wp!fyy(iatm)
     fiz=0.0_wp!fzz(iatm)

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

              rrr = Sqrt(rsq)
              k   = Int(rrr*rdrewd)
              ppp = rrr*rdrewd - Real(k,wp)

! calculate forces using 3pt interpolation

              gk0 = fer(k)
              gk1 = fer(k+1)
              gk2 = fer(k+2)

              t1 = gk0 + (gk1 - gk0)*ppp
              t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

              egamma = (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

! calculate forces

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

! calculate interaction energy using 3-point interpolation

                 vk0 = erc(k)
                 vk1 = erc(k+1)
                 vk2 = erc(k+2)

                 t1 = vk0 + (vk1 - vk0)*ppp
                 t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                 engcpe_rl = engcpe_rl + (t1 + (t2-t1)*ppp*0.5_wp)*chgprd

! calculate virial

                 vircpe_rl = vircpe_rl - egamma*rsq

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

     fxx(iatm)=fxx(iatm)+fix
     fyy(iatm)=fyy(iatm)+fiy
     fzz(iatm)=fzz(iatm)+fiz

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

#ifdef COMPILE_CUDA
  Call stop_timing_ewald_real_forces()
#endif

End Subroutine ewald_real_forces
