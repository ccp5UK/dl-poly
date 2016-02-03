Subroutine relocate_particles        &
           (imcon,rlnk,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange relocation of data between neighbouring
! domains/nodes after positions updates
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov december 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : mxnode,gsum,gmax,gcheck
  Use setup_module
  Use domains_module

  Use site_module
  Use config_module

  Use core_shell_module,   Only : ntshl, listshl,legshl,lshmv_shl,lishp_shl,lashp_shl

  Use constraints_module,  Only : ntcons,listcon,legcon,lshmv_con,lishp_con,lashp_con
  Use pmf_module,          Only : ntpmf

  Use rigid_bodies_module, Only : ntrgd, listrgd,legrgd,lshmv_rgd,lishp_rgd,lashp_rgd

  Use tethers_module,      Only : ntteth,listtet,legtet

  Use bonds_module,        Only : ntbond,listbnd,legbnd
  Use angles_module,       Only : ntangl,listang,legang
  Use dihedrals_module,    Only : ntdihd,listdih,legdih
  Use inversions_module,   Only : ntinv, listinv,leginv

  Implicit None

  Logical,           Intent( In    ) :: lbook
  Integer,           Intent( In    ) :: imcon,megatm,        &
                                        megshl,m_con,megpmf, &
                                        m_rgd,megtet,        &
                                        megbnd,megang,megdih,meginv
  Real( Kind = wp ), Intent( In    ) :: rlnk

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: cut

  Logical           :: safe(1:9)
  Integer           :: i,nlimit,ipx,ipy,ipz
  Real( Kind = wp ) :: big(1:3),det,celprp(1:10),rcell(1:9),x,y,z

  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Define cut

     cut=rlnk+1.0e-6_wp
  End If

! rescale mock cell vectors for non-periodic system

  If (imcon == 0 .or. imcon == 6) Then

! find maximum x,y,z positions

     big=0.0_wp

     Do i =1,natms
        big(1)=Max(big(1),Abs(xxx(i)))
        big(2)=Max(big(2),Abs(yyy(i)))
        big(3)=Max(big(3),Abs(zzz(i)))
     End Do

     If (mxnode > 1) Call gmax(big)

     If (imcon == 0) Then

        cell(1)=Max(2.0_wp*big(1)+cut,3.0_wp*cut,cell(1))
        cell(5)=Max(2.0_wp*big(2)+cut,3.0_wp*cut,cell(5))
        cell(9)=Max(2.0_wp*big(3)+cut,3.0_wp*cut,cell(9))

        cell(2)=0.0_wp
        cell(3)=0.0_wp
        cell(4)=0.0_wp
        cell(6)=0.0_wp
        cell(7)=0.0_wp
        cell(8)=0.0_wp

     Else If (imcon == 6) Then

        cell(9)=Max(2.0_wp*big(3)+cut,3.0_wp*cut,cell(9))

     End If

  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Only when we have more than one domain we need real relocation

  If (mxnode > 1) Then

     Call invert(cell,rcell,det)

! Convert atomic positions from MD centred
! Cartesian coordinates to reduced space ones.
! Populate the move (former halo) indicator array.
! Here we assume no particle has moved more than
! a link-cell width (> rlnk) in any direction
! since last call to relocate as we don't test this!!!

     ixyz(1:natms)=0 ! Initialise move (former halo) indicator
     Do i=1,natms
        x=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        y=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        z=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

! assign domain coordinates (call for errors)

        ipx=Int((x+0.5_wp)*nprx_r)
        ipy=Int((y+0.5_wp)*npry_r)
        ipz=Int((z+0.5_wp)*nprz_r)

        If (idx == 0) Then
           If (x < -half_plus) ixyz(i)=ixyz(i)+1
        Else
           If (ipx < idx) ixyz(i)=ixyz(i)+1
        End If
        If (idx == nprx-1) Then
           If (x >= half_minus) ixyz(i)=ixyz(i)+2
        Else
           If (ipx > idx) ixyz(i)=ixyz(i)+2
        End If

        If (idy == 0) Then
           If (y < -half_plus) ixyz(i)=ixyz(i)+10
        Else
           If (ipy < idy) ixyz(i)=ixyz(i)+10
        End If
        If (idy == npry-1) Then
           If (y >= half_minus) ixyz(i)=ixyz(i)+20
        Else
           If (ipy > idy) ixyz(i)=ixyz(i)+20
        End If

        If (idz == 0) Then
           If (z < -half_plus) ixyz(i)=ixyz(i)+100
        Else
           If (ipz < idz) ixyz(i)=ixyz(i)+100
        End If
        If (idz == nprz-1) Then
           If (z >= half_minus) ixyz(i)=ixyz(i)+200
        Else
           If (ipz > idz) ixyz(i)=ixyz(i)+200
        End If
     End Do

! exchange atom data in -/+ x directions

     Call deport_atomic_data(-1,lbook)
     Call deport_atomic_data( 1,lbook)

! exchange atom data in -/+ y directions

     Call deport_atomic_data(-2,lbook)
     Call deport_atomic_data( 2,lbook)

! exchange atom data in -/+ z directions

     Call deport_atomic_data(-3,lbook)
     Call deport_atomic_data( 3,lbook)

! check system for loss of atoms

     safe(1)=(All(ixyz(1:natms) == 0)) ; Call gcheck(safe(1),"enforce")
     nlimit=natms ; Call gsum(nlimit)
     If ((.not.safe(1)) .or. nlimit /= megatm) Call error(58)

! reassign atom properties

     Do i=1,natms
        atmnam(i)=sitnam(lsite(i))
        ltype(i)=typsit(lsite(i))
        chge(i)=chgsit(lsite(i))
        weight(i)=wgtsit(lsite(i))
        lfrzn(i)=frzsit(lsite(i))
        lfree(i)=fresit(lsite(i))
     End Do

     If (lbook) Then
        safe=.true. ! Initialise safety flag

! Change nlast and refresh record global atom indices for local sorting
! since particles may have been relocated across domains as the
! the old halo is invalid and a new one is not set yet.  This allows for
! local_index search over natms in pmf_units_set and compress_book_intra.
! Otherwise, everywhere else in the code, the search is over nlast as
! domain only indices are caught by the condition (1 >= index <= natms)!!!

        nlast=natms
        Do i=1,nlast
           lsi(i)=i
           lsa(i)=ltg(i)
        End Do
        Call shellsort2(nlast,lsi,lsa)

! Check safety of working arrays for all active bookkeeping arrays

        If (megshl > 0) safe(1)=(ntshl  <= mxshl )
        If (m_con  > 0) safe(2)=(ntcons <= mxcons)
        If (megpmf > 0) safe(3)=(ntpmf  <= mxpmf )
        If (m_rgd  > 0) safe(4)=(ntrgd  <= mxrgd )
        If (megtet > 0) safe(5)=(ntteth <= mxteth)
        If (megbnd > 0) safe(6)=(ntbond <= mxbond)
        If (megang > 0) safe(7)=(ntangl <= mxangl)
        If (megdih > 0) safe(8)=(ntdihd <= mxdihd)
        If (meginv > 0) safe(9)=(ntinv  <= mxinv )

        If (mxnode > 1) Call gcheck(safe)
        If (.not.safe(1)) Call error( 59)
        If (.not.safe(2)) Call error( 41)
        If (.not.safe(3)) Call error(488)
        If (.not.safe(4)) Call error(640)
        If (.not.safe(5)) Call error( 63)
        If (.not.safe(6)) Call error( 31)
        If (.not.safe(7)) Call error( 51)
        If (.not.safe(8)) Call error( 61)
        If (.not.safe(9)) Call error( 77)

! Update shared core-shell, constraint, PMF and RB units

        If (megshl > 0) Call pass_shared_units &
     (mxshl, Lbound(listshl,Dim=1),Ubound(listshl,Dim=1),ntshl, listshl,mxfshl,legshl,lshmv_shl,lishp_shl,lashp_shl)

        If (m_con  > 0) Call pass_shared_units &
     (mxcons,Lbound(listcon,Dim=1),Ubound(listcon,Dim=1),ntcons,listcon,mxfcon,legcon,lshmv_con,lishp_con,lashp_con)

        If (megpmf > 0) Call pmf_units_set()

        If (m_rgd  > 0) Call pass_shared_units &
     (mxrgd, Lbound(listrgd,Dim=1),Ubound(listrgd,Dim=1),ntrgd, listrgd,mxfrgd,legrgd,lshmv_rgd,lishp_rgd,lashp_rgd)

! Compress the rest of the bookkeeping arrays if needed

        If (megtet > 0) Call compress_book_intra &
           (mxteth,ntteth,Ubound(listtet,Dim=1),listtet,mxftet,legtet)

        If (megbnd > 0) Call compress_book_intra &
           (mxbond,ntbond,Ubound(listbnd,Dim=1),listbnd,mxfbnd,legbnd)
        If (megang > 0) Call compress_book_intra &
           (mxangl,ntangl,Ubound(listang,Dim=1),listang,mxfang,legang)
        If (megdih > 0) Call compress_book_intra &
           (mxdihd,ntdihd,Ubound(listdih,Dim=1),listdih,mxfdih,legdih)
        If (meginv > 0) Call compress_book_intra &
           (mxinv,ntinv,  Ubound(listinv,Dim=1),listinv,mxfinv,leginv)

     End If

  Else

! Restore periodic boundaries

     Call pbcshift(imcon,cell,natms,xxx,yyy,zzz)

  End If

! Halt program if potential cutoff exceeds cell width

  If (rlnk > Min(celprp(7),celprp(8),celprp(9))/2.0_wp) Call error(95)

End Subroutine relocate_particles
