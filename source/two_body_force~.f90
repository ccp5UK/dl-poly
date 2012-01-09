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
! author    - i.t.todorov december 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,      Only : idnode,mxnode,gsum
  Use setup_module
  Use config_module,     Only : cell,volm,sumchg,natms,nlast,lsi,lsa, &
                                lexatm,list,xxx,yyy,zzz
  Use ewald_module
  Use vdw_module,        Only : ntpvdw
  Use metal_module,      Only : ntpmet
  Use statistics_module, Only : numrdf

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

  Logical           :: safe = .true.
  Integer           :: fail(1:2),i,j,k,limit,local_index
  Real( Kind = wp ) :: engcpe_rc,vircpe_rc,engcpe_rl,vircpe_rl, &
                       engcpe_ex,vircpe_ex,engcpe_fr,vircpe_fr, &
                       engcpe_nz,vircpe_nz,                     &
                       engden,virden,engmet,virmet,             &
                       engvdw,virvdw,engacc,viracc,tmp,buffer(0:14)

  Real( Kind = wp ), Dimension( : ), Allocatable :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ), Dimension( : ), Allocatable :: rho

  fail=0
  Allocate (xdf(1:mxlist),ydf(1:mxlist),zdf(1:mxlist),rsqdf(1:mxlist), Stat=fail(1))
  If (ntpmet > 0) Allocate (rho(1:mxatms),                             Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'two_body_forces allocation failure, node: ', idnode
     Call error(0)
  End If

! If k-space SPME is evaluated in full infrequently
! check whether at this timestep to evaluate or "refresh"
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

  If (ntpmet > 0) Then

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

! outer loop over atoms

  Do i=1,natms

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

     If ( lrdf .and. ((.not.leql) .or. nstep >= nsteql) .and. &
          Mod(nstep,nstrdf) == 0 ) Call rdf_collect(i,rcut,rsqdf)

  End Do

! counter for rdf statistics outside loop structure

  If ( lrdf .and. ((.not.leql) .or. nstep >= nsteql) .and. &
       Mod(nstep,nstrdf) == 0 ) numrdf = numrdf + 1

! Ewald corrections

  If (keyfce == 2) Then

! calculate coulombic forces, Ewald sum - exclusion corrections
! from intra-like molecular interactions

     If (lbook .and. l_fce) Then

! outer loop over atoms

        Do i=1,natms

! calculate interatomic distances

           Do k=1,lexatm(0,i)

              j=local_index(lexatm(k,i),nlast,lsi,lsa)

              If (j > 0) Then
                 xdf(k)=xxx(i)-xxx(j)
                 ydf(k)=yyy(i)-yyy(j)
                 zdf(k)=zzz(i)-zzz(j)
              Else
                 xdf(k)=0.0_wp
                 ydf(k)=0.0_wp
                 zdf(k)=0.0_wp
              End If
           End Do

! periodic boundary condition

           Call images(imcon,cell,lexatm(0,i),xdf,ydf,zdf)

! calculate correction terms

           Call ewald_excl_forces &
                (i,alpha,epsq,xdf,ydf,zdf,engacc,viracc,stress)

           engcpe_ex=engcpe_ex+engacc
           vircpe_ex=vircpe_ex+viracc

        End Do

     End If

! frozen pairs corrections to coulombic forces

     If (megfrz /= 0 .and. l_fce) Call ewald_frozen_forces &
        (imcon,rcut,alpha,epsq,megfrz,engcpe_fr,vircpe_fr,stress)

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

! Refresh ewald k-space contributions

  If (keyfce == 2 .and. (.not.l_fce)) &
     Call ewald_refresh(engcpe_rc,vircpe_rc,engcpe_ex,vircpe_ex,engcpe_fr,vircpe_fr,stress)

  If (safe) Then
     tmp=0.0_wp
  Else
     tmp=1.0_wp
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

  Deallocate (xdf,ydf,zdf,rsqdf,   Stat=fail(1))
  If (ntpmet > 0) Deallocate (rho, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'two_body_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine two_body_forces
