Subroutine relocate_particles        &
           (imcon,rcut,lbook,megatm, &
           megshl,m_con,megpmf,      &
           m_rgd,megtet,             &
           megbnd,megang,megdih,meginv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange relocation of data between neighbouring
! domains/nodes after position update
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov march 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,       Only : idnode,mxnode,gsum,gmax,gcheck
  Use setup_module
  Use domains_module,     Only : nprx,npry,nprz

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
  Real( Kind = wp ), Intent( In    ) :: rcut

  Logical,           Save :: newjob = .true.
  Integer,           Save :: idx,idy,idz
  Real( Kind = wp ), Save :: sidex,sidey,sidez,cut

  Logical           :: safe(1:9)
  Integer           :: nlx,nly,nlz,i,nlimit
  Real( Kind = wp ) :: big(1:3),celprp(1:10),rcell(1:9),det, &
                       cwx,cwy,cwz, dispx,dispy,dispz, uuu,vvv,www


  If (newjob) Then
     newjob = .false.

! image conditions not compliant with DD and link-cell

     If (imcon == 4 .or. imcon == 5 .or. imcon == 7) Call error(300)

! Get this node's (domain's) coordinates

     idz=idnode/(nprx*npry)
     idy=idnode/nprx-idz*npry
     idx=Mod(idnode,nprx)

! Get the domains' dimensions in reduced space
! (domains are geometrically equivalent)

     sidex=1.0_wp/Real(nprx,wp)
     sidey=1.0_wp/Real(npry,wp)
     sidez=1.0_wp/Real(nprz,wp)

! Define cut

     cut=rcut+1.0e-6_wp
  End If

! rescale mock cell vectors for non-periodic system

  If (imcon == 0 .or. imcon == 6) Then

! find maximum x,y,z positions

     big=0.0_wp

     Do i =1,natms
        big(1)=Max(big(1),abs(xxx(i)))
        big(2)=Max(big(2),abs(yyy(i)))
        big(3)=Max(big(3),abs(zzz(i)))
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

! Only when we use more than one CPU we need real relocation

  If (mxnode > 1) Then

! Calculate the displacements from the origin of the MD cell
! to the origin of this domain in reduced space

! First term (0.5_wp) = move to the bottom left corner of MD cell
! Second term, first term (side) = scale by the number of domains
! in the given direction
! Second term, second term, first term (id) = move to the bottom
! left corner of this domain in the given direction
! Second term, second term, second term (0.5_wp) = move to the
! middle of this domain

     dispx=0.5_wp-sidex*(Real(idx,wp)+0.5_wp)
     dispy=0.5_wp-sidey*(Real(idy,wp)+0.5_wp)
     dispz=0.5_wp-sidez*(Real(idz,wp)+0.5_wp)

! calculate link cell dimensions per node

     nlx=Int(sidex*celprp(7)/cut)
     nly=Int(sidey*celprp(8)/cut)
     nlz=Int(sidez*celprp(9)/cut)

! Calculate a link-cell width in every direction in the
! reduced space of the domain
! First term = the width of the domain in reduced space
! Second term = number of link-cells per domain per direction

     cwx=sidex/Real(nlx,wp)
     cwy=sidey/Real(nly,wp)
     cwz=sidez/Real(nlz,wp)

! Get the inverse cell matrix

     Call invert(cell,rcell,det)

! Convert atomic positions from centred Cartesian coordinates
! to reduced space coordinates of this node

     Do i=1,natms
        uuu=xxx(i)
        vvv=yyy(i)
        www=zzz(i)

        xxx(i)=rcell(1)*uuu+rcell(4)*vvv+rcell(7)*www+dispx
        yyy(i)=rcell(2)*uuu+rcell(5)*vvv+rcell(8)*www+dispy
        zzz(i)=rcell(3)*uuu+rcell(6)*vvv+rcell(9)*www+dispz
     End Do

! exchange atom data in -/+ x directions

     Call deport_atomic_data(-1,lbook,sidex,sidey,sidez,cwx,cwy,cwz)
     Call deport_atomic_data( 1,lbook,sidex,sidey,sidez,cwx,cwy,cwz)

! exchange atom data in -/+ y directions

     Call deport_atomic_data(-2,lbook,sidex,sidey,sidez,cwx,cwy,cwz)
     Call deport_atomic_data( 2,lbook,sidex,sidey,sidez,cwx,cwy,cwz)

! exchange atom data in -/+ z directions

     Call deport_atomic_data(-3,lbook,sidex,sidey,sidez,cwx,cwy,cwz)
     Call deport_atomic_data( 3,lbook,sidex,sidey,sidez,cwx,cwy,cwz)

! check system for loss of atoms

     nlimit=natms
     Call gsum(nlimit)
     If (nlimit /= megatm) Call error(58)

! restore atomic coordinates to real coordinates
! and reassign atom properties

     Do i=1,natms
        uuu=xxx(i)-dispx
        vvv=yyy(i)-dispy
        www=zzz(i)-dispz

        xxx(i)=cell(1)*uuu+cell(4)*vvv+cell(7)*www
        yyy(i)=cell(2)*uuu+cell(5)*vvv+cell(8)*www
        zzz(i)=cell(3)*uuu+cell(6)*vvv+cell(9)*www

        atmnam(i)=sitnam(lsite(i))
        ltype(i)=typsit(lsite(i))
        chge(i)=chgsit(lsite(i))
        weight(i)=wgtsit(lsite(i))
        lfrzn(i)=frzsit(lsite(i))
        lfree(i)=fresit(lsite(i))
     End Do

     If (lbook) Then
        safe=.true. ! Initialise safety flag

! Change nlast and refresh gtl or record global atom indices for local sorting
! since particles may have been relocated across domains as the
! the old halo is invalid and a new one is not set yet.  This allows for
! local_index search over natms in pmf_units_set and compress_book_intra.
! Otherwise, everywhere else in the code, the search is over nlast as
! domain only indices are caught by the condition (1 >= index <= natms)!!!

        nlast=natms
        If (gtl_b > 0) Then
           Call get_gtl(lbook)
        Else
           Do i=1,nlast
              lsi(i)=i
              lsa(i)=ltg(i)
           End Do
           Call shellsort2(nlast,lsi,lsa)
        End If

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

  End If

! Restore periodic boundaries

  Call pbcshift(imcon,cell,natms,xxx,yyy,zzz)

! Halt program if potential cutoff exceeds cell width

  If (rcut > Min(celprp(7),celprp(8),celprp(9))/2.0_wp) Call error(95)

End Subroutine relocate_particles
