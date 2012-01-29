! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Center for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.

Subroutine two_body_forces_cuda_helper(&
     iatmb,iters,keyfce,rho,safe,&
     imcon,rcut,rvdw,alpha,epsq,stress,&
     engmet,virmet,engvdw,virvdw,engcpe_rl,vircpe_rl)
  Use kinds_f90
  Use config_module, Only : cell,list,xxx,yyy,zzz,ltype,ltg,chge,fxx,fyy,fzz
  Use setup_module,  Only : mxatms,mxgrid,nrite
  Use vdw_module
  Use metal_module
  Use comms_module,  Only : idnode

  Implicit None
  Integer,           Intent(In   ) :: iatmb,iters,imcon,keyfce
  Real( Kind = wp ), Intent(In   ) :: rcut,rvdw,alpha,epsq
  Real( Kind = wp),  Dimension(1:mxatms), Intent(In) :: rho
  Real( Kind = wp ), Intent(InOut) :: engmet,virmet,engvdw,virvdw,engcpe_rl,vircpe_rl
  Real( Kind = wp ), Dimension( 1:9 ),&
                     Intent(InOut) :: stress
  Logical,           Intent(InOut) :: safe
  Integer :: i,j,k,limit
  Real( Kind = wp ), Dimension( : ), Allocatable :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ) :: engacc,viracc,aaa,bbb,ccc
  Real( Kind = wp ), Dimension(1:9) :: stressloc
  Real( Kind = wp) :: engmetloc,virmetloc,engvdwloc,virvdwloc,engcpe_rlloc,vircpe_rlloc
  Logical,           Save :: newjob = .true.

! vdw-specific stuff
  Real( Kind = wp ), Save :: dlrpot,rdr,vdw_rcsq
! ewlds_real-stuff
  Real( Kind = wp ), Dimension( : ), Allocatable, Save :: erc,fer
  Real( Kind = wp ), Save :: erl_rcsq,drewd,rdrewd
  Integer           :: fail,l
  Integer, Save     :: keypot

! TODO: reduction against 'safe'

! Perform the static stuff here to simplify the threaded logic (eliminates
! a number of directives)
  If (newjob) Then
     newjob = .false.
! vdw stuff
     dlrpot = rvdw/Real(mxgrid-4,wp)
     rdr    = 1.0_wp/dlrpot
     vdw_rcsq   = rvdw**2

! ewald_real stuff
     Allocate (erc(1:mxgrid),fer(1:mxgrid), Stat=fail)
     If (fail > 0) Then
        Write(nrite,'(/,1x,a,i0)') 'ewald_real_forces allocation failure, idnode', idnode
        Call error(0)
     End If
     erl_rcsq = rcut**2
     drewd = rcut/Real(mxgrid-4,wp)
     rdrewd = 1.0_wp/drewd
     Call erfcgen(rcut,alpha,mxgrid,erc,fer)

     If (ntpmet>0) Then
        keypot=0
        Do l=1,ntpmet
           keypot=ltpmet(l)
           If (l > 1) Then
              If (keypot /= ltpmet(l-1)) Call error(92)
           End If
        End Do
     End If
  End If

  Allocate (xdf(1:512),ydf(1:512),zdf(1:512),rsqdf(1:512))

!$OMP PARALLEL PRIVATE(i,limit,aaa,bbb,ccc,engmetloc,virmetloc,&
!$OMP engvdwloc,virvdwloc,engcpe_rlloc,vircpe_rlloc,stressloc,k,j)
  engmetloc=0.0_wp
  virmetloc=0.0_wp
  engvdwloc=0.0_wp
  virvdwloc=0.0_wp
  engcpe_rlloc=0.0_wp
  vircpe_rlloc=0.0_wp
  stressloc=0.0_wp

  Do i = iatmb, iatmb+iters-1
     limit=list(0,i)

!$OMP DO
     Do k=1,limit
        j=list(k,i)
        xdf(k)=xxx(i)-xxx(j)
        ydf(k)=yyy(i)-yyy(j)
        zdf(k)=zzz(i)-zzz(j)
     End Do
!$OMP END DO NOWAIT

     Call images_helper_orphan(imcon,cell,limit,xdf,ydf,zdf)

!$OMP DO
     Do k=1,limit
        rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
     End Do
!$OMP END DO NOWAIT

     If (ntpmet > 0) Then
        Call metal_forces_helper(&
             i,xdf,ydf,zdf,rsqdf,rho,keypot,&
             engmetloc,virmetloc,stressloc,safe)
     End If


     If (ntpvdw > 0) Then
        Call vdw_forces_helper(&
             i,rvdw,xdf,ydf,zdf,rsqdf,engvdwloc,&
             virvdwloc,stressloc,dlrpot,rdr,vdw_rcsq)
     End If

     If (keyfce == 2) Then
        Call ewald_real_forces_helper &
             (i,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf,engcpe_rlloc,&
             vircpe_rlloc,stressloc,erc,fer,drewd,rdrewd,erl_rcsq)
     End If
!$OMP BARRIER
  End Do

!$OMP CRITICAL
  If (ntpmet>0) Then
     engmet = engmet + engmetloc
     virmet = virmet + virmetloc
  End If

  If (ntpvdw > 0) Then
     engvdw = engvdw + engvdwloc
     virvdw = virvdw + virvdwloc
  End If

  If (keyfce>0 .and. Mod(keyfce,2)==0 .and. keyfce<=10) Then
     engcpe_rl=engcpe_rl + engcpe_rlloc
     vircpe_rl=vircpe_rl + vircpe_rlloc
  End If

  stress = stress + stressloc
!$OMP END CRITICAL

!$OMP END PARALLEL
  Deallocate (xdf,ydf,zdf,rsqdf)

End Subroutine two_body_forces_cuda_helper




Subroutine two_body_forces                        &
           (imcon,rcut,rvdw,rmet,keyens,          &
           alpha,epsq,keyfce,nstfce,lbook,megfrz, &
           lrdf,nstrdf,leql,nsteql,nstep,         &
           elrc,virlrc,elrcm,vlrcm,               &
           engcpe,vircpe,engsrp,virsrp,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating interatomic forces and rdf
! using the verlet neighbour list
!
! ntpvdw > 0 ------ switch for vdw potentials calculation
! ntpmet > 0 ------ switch for metal local density and potentials
!                   calculations
!
! ELECTROSTATICS KEYS
!
! keyfce = 0 ------ no electrostatics
! keyfce = 2 ------ Ewald sum (ewald_spme,ewald_real,ewald_excl)
! keyfce = 4 ------ distance dependent dielectric potential (coul_dddp)
! keyfce = 6 ------ coulombic 1/r potential (coul_cp)
! keyfce = 8 ------ force-shifted and damped coulombic potential (coul_fscp)
! keyfce =10 ------ reaction field and damped coulombic potential (coul_rfp)
!
! nstfce - the rate at which the k-space contributions of SPME are
!          refreshed.  Once every 1 <= nstfce <= 7 steps.
!
! copyright - daresbury laboratory
! author    - i.t.todorov june 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsum
  Use setup_module
  Use config_module,     Only : cell,volm,sumchg,natms,  &
                                lexatm,list,xxx,yyy,zzz, &
                                ltg,fxx,fyy,fzz,ltype,chge
  Use ewald_module
  Use vdw_module,        Only : ntpvdw, prmvdw
  Use metal_module,      Only : ntpmet,ltpmet,lstmet,vmet,dmet
  Use statistics_module, Only : numrdf

#ifdef COMPILE_CUDA
! Needed by the CUDA code
  Use vdw_module,        Only : ls_vdw,lstvdw,ltpvdw,vvdw,gvdw
  Use dl_poly_cuda_module
  Use iso_c_binding
#endif

  Implicit None

  Logical,                                  Intent( In    ) :: lbook,lrdf,leql
  Integer,                                  Intent( In    ) :: imcon,keyens,  &
                                                               keyfce,nstfce, &
                                                               megfrz,nstrdf, &
                                                               nsteql,nstep
  Real( Kind = wp ),                        Intent( In    ) :: rcut,rvdw,rmet,&
                                                               alpha,epsq
  Real( Kind = wp ),                        Intent( In    ) :: elrc,virlrc
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent( InOut ) :: elrcm,vlrcm
  Real( Kind = wp ),                        Intent(   Out ) :: engcpe,vircpe, &
                                                               engsrp,virsrp
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress


  Logical,           Save :: new_nz    = .true.
  Real( Kind = wp ), Save :: factor_nz = 0.0_wp

  Logical           :: safe = .true., l_do_rdf
  Integer           :: fail(1:2),i,j,k,limit
  Real( Kind = wp ) :: engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl, &
                       engcpe_ex,vircpe_ex,engcpe_fr,vircpe_fr, &
                       engcpe_nz,vircpe_nz,                     &
                       engden,virden,engmet,virmet,             &
                       engvdw,virvdw,engacc,viracc,tmp,buffer(0:14)

  Real( Kind = wp ), Dimension( : ), Allocatable :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ), Dimension( : ), Allocatable :: rho

#ifdef COMPILE_CUDA
! CUDA-specific; vdw_forces code moved here to support the initialiser:
  Logical,           Save  :: newjob = .true.
  Real( Kind = wp ), Save  :: dlrpot,rdr,rcsq

! CUDA: holds the truth of the list structure being online after its
!   construction by the accelerating function.
  Integer                 :: cuda_ilo, keypot

! Create a copy of ls_vdw for passing to C
  Logical( c_bool ) :: ls_vdw_c
  ls_vdw_c = ls_vdw
#endif

  fail=0
  Allocate (xdf(1:mxlist),ydf(1:mxlist),zdf(1:mxlist),rsqdf(1:mxlist), Stat=fail(1))
  If (ntpmet > 0) Allocate (rho(1:mxatms),                             Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'two_body_forces allocation failure, node: ', idnode
     Call error(0)
  End If

#ifdef COMPILE_CUDA
  Call start_timing_two_body_forces()
#endif

  l_do_rdf = (lrdf .and. ((.not.leql) .or. nstep >= nsteql) .and. Mod(nstep,nstrdf) == 0)

! If k-space SPME is evaluated in full infrequently
! check wheather at this timestep to evaluate or "refresh"
! with old values.  At restart allocate the "refresh"
! k-space SPME arrays and force the full force evaluation.

  If (keyfce == 2) Call ewald_check(nstep,nsteql,nstfce)

! initialise energy and virial accumulators

  engden    = 0.0_wp
  virden    = 0.0_wp

  engmet    = 0.0_wp
  virmet    = 0.0_wp

  engvdw    = 0.0_wp
  virvdw    = 0.0_wp

  engsrp    = 0.0_wp
  virsrp    = 0.0_wp


  engcpe_rc = 0.0_wp
  vircpe_rc = 0.0_wp

  engcpe_rl = 0.0_wp
  vircpe_rl = 0.0_wp

  engcpe_ex = 0.0_wp
  vircpe_ex = 0.0_wp

  engcpe_fr = 0.0_wp
  vircpe_fr = 0.0_wp

  engcpe_nz = 0.0_wp
  vircpe_nz = 0.0_wp

  engcpe    = 0.0_wp
  vircpe    = 0.0_wp

! Set up non-bonded interaction (verlet) list using link cells

  Call link_cell_pairs(imcon,rcut,lbook,megfrz)

#ifdef COMPILE_CUDA
  If (dl_poly_cuda_offload_tbforces() .and. &
      dl_poly_cuda_offload_link_cell_pairs() .and. dl_poly_cuda_is_cuda_capable() .and. dl_poly_cuda_offload_link_cell_pairs_re()) Then
     cuda_ilo = 1;
  Else
     cuda_ilo = 0;
  End If
#endif

  If (ntpmet > 0) Then

#ifdef COMPILE_CUDA
     Call metal_ld_compute_get_keypot(keypot)

! CUDA: (1) Unfortunately, due to the global stuff that metal_ld_compute
!       invokes (metal_ld_export), these functions cannot be easilly
!       integrated with the two_body_forces_cuda_invoke.
!       (2) If the list is online (cuda_ilo==1) and the following bits need
!       to be invoked, then copy the list back to the host.
     If (cuda_ilo==1 .and. keypot /= 0) Then
! CUDA: TODO: We need to signal the two body forces' accelerator that the
!             list has been pulled in.
        Call link_cell_pairs_cuda_pull_lists(natms)
     End If
#endif

! Reset metal long-range corrections (constant pressure/stress only)

     If (keyens >= 20) Call metal_lrc(imcon,rmet,elrcm,vlrcm)

! calculate local density in metals
     Call metal_ld_compute             &
          (imcon,rmet,keyfce,elrcm,vlrcm, &
          xdf,ydf,zdf,rsqdf,              &
          rho,engden,virden,stress)
  End If

! calculate coulombic forces, Ewald sum - fourier contribution

  If (keyfce == 2 .and. l_fce) Call ewald_spme_forces(alpha,epsq,engcpe_rc,vircpe_rc,stress)

#ifdef COMPILE_CUDA
  Call start_timing_two_body_forces_any()

! The CUDA port implements the features of dl_poly_3 - a subset of those found in dl_poly_4
! In particular, only tabulated calculations are available in metal_forces and vdw_forces
! The CUDA acceleration is not called if direct calculation is required
 If (dl_poly_cuda_offload_tbforces() .and. dl_poly_cuda_is_cuda_capable()) Then

     Call two_body_forces_cuda_initialise(&
          cuda_ilo, natms,&
          mxlist,mxatms,mxatdm,mxmet,ntpmet,ntpvdw,keyfce,imcon,&
          cell,xxx,yyy,zzz,list,ltype,             &
          ltpmet,lstmet,vmet,dmet,rho,             &
          mxgrid, mxvdw, ltpvdw, lstvdw, ls_vdw_c, &
          vvdw, gvdw, rvdw, ltg,                   &
          chge, rcut, alpha, epsq)
     Call start_timing_two_body_forces_any()

     Call two_body_forces_cuda_invoke(&
             fxx, fyy, fzz, stress, engmet, virmet, engvdw, virvdw, engcpe_rl, vircpe_rl, safe)

     Call two_body_forces_cuda_finalise()
  Else
#endif

  Do i=1,natms, 1

! outer loop over atoms
! Get list limit

     limit=list(0,i)

! calculate interatomic distances

     Do k=1,limit
        j=list(k,i)

        xdf(k)=xxx(i)-xxx(j)
        ydf(k)=yyy(i)-yyy(j)
        zdf(k)=zzz(i)-zzz(j)
     End Do

! periodic boundary conditions

     Call images(imcon,cell,limit,xdf,ydf,zdf)

! square of distances

     Do k=1,limit
        rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
     End Do

! calculate metal forces and potential

     If (ntpmet > 0) Then
        Call metal_forces &
       (i,rmet,xdf,ydf,zdf,rsqdf,rho,engacc,viracc,stress,safe)

        engmet=engmet+engacc
        virmet=virmet+viracc
     End If

! calculate short-range force and potential terms

     If (ntpvdw > 0) Then
        Call vdw_forces &
       (i,rvdw,xdf,ydf,zdf,rsqdf,engacc,viracc,stress)

        engvdw=engvdw+engacc
        virvdw=virvdw+viracc
     End If

!!!!!!!!!!!!!!!!!!!!!!!!!
! COULOMBIC CONTRIBUTIONS
!!!!!!!!!!!!!!!!!!!!!!!!1

     If (keyfce == 2) Then

! calculate coulombic forces, Ewald sum - real space contribution

        Call ewald_real_forces &
       (i,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf,engacc,viracc,stress)

        engcpe_rl=engcpe_rl+engacc
        vircpe_rl=vircpe_rl+viracc

     Else If (keyfce == 4) Then

! distance dependant dielectric potential

        Call coul_dddp_forces &
       (i,rcut,epsq,xdf,ydf,zdf,rsqdf,engacc,viracc,stress)

        engcpe_rl=engcpe_rl+engacc
        vircpe_rl=vircpe_rl+viracc

     Else If (keyfce == 6) Then

! coulombic 1/r potential with no truncation or damping

        Call coul_cp_forces &
       (i,rcut,epsq,xdf,ydf,zdf,rsqdf,engacc,viracc,stress)

        engcpe_rl=engcpe_rl+engacc
        vircpe_rl=vircpe_rl+viracc

     Else If (keyfce == 8) Then

! force-shifted coulomb potentials

        Call coul_fscp_forces &
       (i,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf,engacc,viracc,stress)

        engcpe_rl=engcpe_rl+engacc
        vircpe_rl=vircpe_rl+viracc

     Else If (keyfce == 10) Then

! reaction field potential

        Call coul_rfp_forces &
       (i,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf,engacc,viracc,stress)

        engcpe_rl=engcpe_rl+engacc
        vircpe_rl=vircpe_rl+viracc

     End If
! accumulate radial distribution functions


     If (l_do_rdf) Call rdf_collect(i,rcut,rsqdf)  


  End Do

  If (safe) Then
     tmp=0.0_wp
  Else
     tmp=1.0_wp
  End If

! In the case of excluded interactions
! accumulate furhter radial distribution functions and/or
! calculate Ewald corrections due to long-range exclusions
! if (keyfce == 2 .and. l_fce)  

#ifdef COMPILE_CUDA
  End If
  Call stop_timing_two_body_forces_any()

! 20100129/ck:
!  see the note in link_cell_pairs part where the initialiser is
!  invoked -- we keep the list structure online so that we don't
!  have to copy it over again and again
  If (cuda_ilo==1) Then
     Call link_cell_pairs_cuda_finalise()
  End If
#endif


! calculate coulombic forces, Ewald sum - exclusion corrections
! from intra-like molecular interactions

      If ( lbook .and. (l_do_rdf .or. (keyfce == 2 .and. l_fce)) ) Then

! outer loop over atoms

         Do i=1,natms

!Get list limit

            limit=list(-1,i)-list(0,i)
            If (limit > 0) Then

! calculate interatomic distances

               Do k=1,limit
                  j=list(list(0,i)+k,i)
                  
                  xdf(k)=xxx(i)-xxx(j)
                  ydf(k)=yyy(i)-yyy(j)
                  zdf(k)=zzz(i)-zzz(j) 
               End Do

! periodic boundary condition

               Call images(imcon,cell,limit,xdf,ydf,zdf)

! square of distance

               Do k=1,limit
                    rsqdf(k)=xdf(k)**2+ydf(k)**2+zdf(k)**2
               End Do

!accumulate radial distribution functions

               If (l_do_rdf) Call rdf_excl_collect(i,rcut,rsqdf)

! calculate correction terms

               If (keyfce == 2 .and. l_fce) Then
                  Call ewald_excl_forces &
                (i,rcut,alpha,epsq,xdf,ydf,zdf,rsqdf,engacc,viracc,stress)
             

                  engcpe_ex=engcpe_ex+engacc
                  vircpe_ex=vircpe_ex+viracc
               End If

            End If

         End Do

      End If

!counter for rdf statistics outside loop structures

      If (l_do_rdf) numrdf = numrdf+1
      
      Deallocate (xdf,ydf,zdf,rsqdf,   Stat=fail(1))
      If (ntpmet > 0) Deallocate (rho, Stat=fail(2))
      If (Any(fail > 0)) Then
         Write(nrite,'(/,1x,a,i0)') 'two_body forces deallocation failure, node: ', idnode
         Call error(0)
      End If

! Furhter Ewald corrections or an infrequent path

      If (keyfce == 2) Then
         If (l_fce) Then

! frozen pairs corrections to coulombic forces

             If (megfrz /= 0) Call ewald_frozen_forces &
                (imcon,rcut,alpha,epsq,keyens,engcpe_fr,vircpe_fr,stress)

         Else 

 ! Refresh all Ewald k-space contributions but the non-zero system charge

             Call ewald_refresh(engcpe_rc,vircpe_rc,engcpe_ex,vircpe_ex,engcpe_fr,vircpe_fr,stress)

         End If


! non-zero total system charge correction (for the whole system)
! ( Fuchs, Proc. R. Soc., A, 151, (585),1935 )

      If (Abs(sumchg) > 1.0e-6_wp) Then
         If (new_nz) Then
            new_nz = .false.
            factor_nz = -0.5_wp * (pi*r4pie0/epsq) * (sumchg/alpha)**2
         End If

         engcpe_nz=factor_nz/volm
         vircpe_nz=-3.0_wp*engcpe_nz
      End If
  End If

! sum up contributions to potentials

  If (mxnode > 1) Then
     buffer( 0) = tmp
     buffer( 1) = engden
     buffer( 2) = virden
     buffer( 3) = engmet
     buffer( 4) = virmet
     buffer( 5) = engvdw
     buffer( 6) = virvdw
     buffer( 7) = engcpe_rc
     buffer( 8) = vircpe_rc
     buffer( 9) = engcpe_rl
     buffer(10) = vircpe_rl
     buffer(11) = engcpe_ex
     buffer(12) = vircpe_ex
     buffer(13) = engcpe_fr
     buffer(14) = vircpe_fr

     Call gsum(buffer(0:14))

     tmp       = buffer( 0)
     engden    = buffer( 1)
     virden    = buffer( 2)
     engmet    = buffer( 3)
     virmet    = buffer( 4)
     engvdw    = buffer( 5)
     virvdw    = buffer( 6)
     engcpe_rc = buffer( 7)
     vircpe_rc = buffer( 8)
     engcpe_rl = buffer( 9)
     vircpe_rl = buffer(10)
     engcpe_ex = buffer(11)
     vircpe_ex = buffer(12)
     engcpe_fr = buffer(13)
     vircpe_fr = buffer(14)
  End If

  safe=(tmp < 0.5_wp)
  If (.not.safe) Call error(505)

! Globalise coulombic contributions: cpe

  engcpe = engcpe_rc + engcpe_rl + engcpe_ex + engcpe_fr + engcpe_nz
  vircpe = vircpe_rc + vircpe_rl + vircpe_ex + vircpe_fr + vircpe_nz

! Add non-zero total system charge correction to
! diagonal terms of stress tensor (per node)

  tmp = - vircpe_nz/(3.0_wp*Real(mxnode,wp))
  stress(1) = stress(1) + tmp
  stress(5) = stress(5) + tmp
  stress(9) = stress(9) + tmp

! Globalise short-range & metal interactions with
! their long-range corrections contributions: srp

  engsrp = (engden + engmet + elrcm(0)) + (engvdw + elrc)
  virsrp = (virden + virmet + vlrcm(0)) + (virvdw + virlrc)

! Add long-range corrections to diagonal terms of stress tensor (per node)

  tmp = - (virlrc+vlrcm(0))/(3.0_wp*Real(mxnode,wp))
  stress(1) = stress(1) + tmp
  stress(5) = stress(5) + tmp
  stress(9) = stress(9) + tmp

#if COMPILE_CUDA
  Call stop_timing_two_body_forces()
#endif

End Subroutine two_body_forces
