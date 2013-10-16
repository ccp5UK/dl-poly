Subroutine build_book_intra              &
           (megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,         &
           megrgd,degrot,degtra,         &
           megtet,megbnd,megang,megdih,meginv)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for initialising global bookkeeping of intra-like
! interactions: core-shell, bond constraints, PMF constraints, RBs,
! tethered atoms, chemical bonds, valence angles, torsion and improper
! torsion angles, and inversion angles
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SETUP MODULES

  Use comms_module,  Only : idnode,mxnode,gcheck
  Use setup_module

! SITE MODULE

  Use site_module

! CONFIG MODULE

  Use config_module, Only : natms,nlast,lsi,lsa

! INTERACTION MODULES

  Use core_shell_module

  Use constraints_module
  Use pmf_module

  Use rigid_bodies_module

  Use tethers_module

  Use bonds_module
  Use angles_module
  Use dihedrals_module
  Use inversions_module

  Implicit None

  Integer,           Intent( In    ) :: megatm,atmfre,atmfrz, &
                                        megshl,megpmf,        &
                                        megtet,megbnd,megang, &
                                        megdih,meginv
  Integer,           Intent( InOut ) :: megfrz,megcon,megrgd
  Integer(Kind=ip),  Intent( InOut ) :: degrot,degtra

  Logical :: safe(1:11),go
  Integer :: fail(1:2),i,j,isite,itmols,imols,   &
             nsatm,neatm,nlapm,local_index,      &
             iat0,jat0,kat0,lat0,                &
             iatm,jatm,katm,latm,matm,natm,      &
             jshels,kshels,lshels,mshels,        &
             jconst,kconst,lconst,lpmf,          &
             jrigid,krigid,lrigid,mrigid,irigid, &
             jteths,kteths,lteths,               &
             jbonds,kbonds,lbonds,               &
             jangle,kangle,langle,               &
             jdihed,kdihed,ldihed,               &
             jinver,kinver,linver

  Integer, Dimension( : ), Allocatable :: iwrk,irgd,irgd0, &
                                          i1pmf,i1pmf0,i2pmf,i2pmf0

  fail=0
  Allocate (iwrk(1:mxatms),                                Stat=fail(1))
  If (m_rgd > 0) Allocate (irgd(1:mxlrgd),irgd0(1:mxlrgd), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'build_book_intra allocation failure, node: ', idnode
     Call error(0)
  End If

! Initialise safety flags

  safe=.true.

! Initialise bookkeeping indices

  isite=0

! "j" stands for running index of intra-like unit
!
! "k" stands for last index of intra-like unit
!     on last molecule of previous molecule type

  jshels=0
  kshels=0

  jconst=0
  kconst=0

! No 'jpmf' and 'kpmf' needed since PMF is defined on one and only one
! molecular type and therefore 'ntpmf' is enough and is used locally as 'jpmf'

  jrigid=0
  krigid=0

  jteths=0
  kteths=0

  jbonds=0
  kbonds=0

  jangle=0
  kangle=0

  jdihed=0
  kdihed=0

  jinver=0
  kinver=0

! global atom counter

  nsatm=0

! loop over molecule types in the system

  Do itmols=1,ntpmls

! loop over molecules of this type

     Do imols=1,nummols(itmols)

! last atom in the molecule

        neatm=nsatm+numsit(itmols)

! number of local atoms on this molecule

        nlapm=0

! From the first till the last atom of this molecule, get the number
! of localised atoms on this node and save their local_index in iwrk

        Do iatm=nsatm+1,neatm
           iat0=local_index(iatm,nlast,lsi,lsa)
           If (iat0 > natms) iat0=0

           If (iat0 > 0) Then
              nlapm=nlapm+1
              iwrk(nlapm)=iat0
           End If
        End Do

! If there are atoms of this molecule on this node, get global indices and
! corresponding local_indices of the specific intra-like unit.  If any of
! the local indices exists increase the local number of these intra-like
! units (local = belongs to this node).  If it's safe to proceed (array
! bounds check) record the global number of the intra-like unit and the
! global atomic indices of its atoms.  Tag the local numbers of the
! intra-like units and count them for each local atom (on domain) that is
! a member of them in the intra-like unit's legend array.

        If (nlapm > 0) Then

! Construct core-shell list

           Do lshels=1,numshl(itmols)
              iatm=lstshl(1,lshels+kshels)+isite
              jatm=lstshl(2,lshels+kshels)+isite

              iat0=local_index(iatm,nlast,lsi,lsa)
              jat0=local_index(jatm,nlast,lsi,lsa)

              If (iat0 > natms) iat0=0
              If (jat0 > natms) jat0=0

              If (iat0 > 0 .or. jat0 > 0) Then
                 jshels=jshels+1

                 If (jshels <= mxshl) Then
                    listshl(0,jshels)=lshels+kshels
                    listshl(1,jshels)=iatm
                    listshl(2,jshels)=jatm

                    If (iat0 > 0) Then
                       Call tag_legend(safe(1),iat0,jshels,legshl,mxfshl)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many core-shell type neigbours !!! ***",                 &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legshl(0,iat0),               &
  "***           but maximum length allowed: ", mxfshl,                       &
  "***           for particle (global ID #): ", iatm,                         &
  "***           on mol. site (local  ID #): ", lstshl(1,lshels+kshels),      &
  "***           of unit      (local  ID #): ", lshels,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (jat0 > 0) Then
                       Call tag_legend(safe(1),jat0,jshels,legshl,mxfshl)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many core-shell type neigbours !!! ***",                 &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legshl(0,jat0),               &
  "***           but maximum length allowed: ", mxfshl,                       &
  "***           for particle (global ID #): ", jatm,                         &
  "***           on mol. site (local  ID #): ", lstshl(2,lshels+kshels),      &
  "***           of unit      (local  ID #): ", lshels,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If
                 Else
                    safe(2)=.false.
                 End If
              End If
           End Do

! Construct constraint bond list

           Do lconst=1,numcon(itmols)
              iatm=lstcon(1,lconst+kconst)+isite
              jatm=lstcon(2,lconst+kconst)+isite

              iat0=local_index(iatm,nlast,lsi,lsa)
              jat0=local_index(jatm,nlast,lsi,lsa)

              If (iat0 > natms) iat0=0
              If (jat0 > natms) jat0=0

              If (iat0 > 0 .or. jat0 > 0) Then
                 jconst=jconst+1
                 If (jconst <= mxcons) Then
                    listcon(0,jconst)=lconst+kconst
                    listcon(1,jconst)=iatm
                    listcon(2,jconst)=jatm

                    If (iat0 > 0) Then
                       Call tag_legend(safe(1),iat0,jconst,legcon,mxfcon)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many constraint type neigbours !!! ***",                 &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legcon(0,iat0),               &
  "***           but maximum length allowed: ", mxfcon,                       &
  "***           for particle (global ID #): ", iatm,                         &
  "***           on mol. site (local  ID #): ", lstcon(1,lconst+kconst),      &
  "***           of unit      (local  ID #): ", lconst,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (jat0 > 0) Then
                       Call tag_legend(safe(1),jat0,jconst,legcon,mxfcon)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many constraint type neigbours !!! ***",                 &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legcon(0,jat0),               &
  "***           but maximum length allowed: ", mxfcon,                       &
  "***           for particle (global ID #): ", jatm,                         &
  "***           on mol. site (local  ID #): ", lstcon(2,lconst+kconst),      &
  "***           of unit      (local  ID #): ", lconst,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If
                 Else
                    safe(3)=.false.
                 End If
              End If
           End Do

! Construct PMF bond list
! Note: This is executed for only one given molecular
!       type as only one PMF type per MD system is allowed.

           Do lpmf=1,numpmf(itmols) ! numpmf can only be 1 or 0, so the 'Do' loop is used as an 'If' condition
              Allocate (i1pmf(1:mxtpmf(1)),i1pmf0(1:mxtpmf(1)),i2pmf(1:mxtpmf(2)),i2pmf0(1:mxtpmf(2)), Stat=fail(1))
              If (fail(1) > 0) Then
                 Write(nrite,'(/,1x,a,i0)') 'build_book_intra PMF allocation failure, node: ', idnode
                 Call error(0)
              End If

              i1pmf=0 ; i1pmf0=0
              Do i=1,mxtpmf(1)
                 i1pmf(i) =lstpmf(i,1)+isite
                 i1pmf0(i)=local_index(i1pmf(i),nlast,lsi,lsa)
                 If (i1pmf0(i) > natms) i1pmf0(i)=0
              End Do

              i2pmf=0 ; i2pmf0=0
              Do i=1,mxtpmf(2)
                 i2pmf(i) =lstpmf(i,2)+isite
                 i2pmf0(i)=local_index(i2pmf(i),nlast,lsi,lsa)
                 If (i2pmf0(i) > natms) i2pmf0(i)=0
              End Do

              If (Any(i1pmf0 > 0) .or. Any(i2pmf0 > 0)) Then
                 ntpmf=ntpmf+1
                 If (ntpmf <= mxpmf) Then

! This holds the global PMF index

                    listpmf(0,1,ntpmf)=imols

! For presence of : PMF unit 1 only - this holds 1
!                   PMF unit 2 only - this holds 2
!                   both units 1&2  - this holds 3
! It CANNOT and MUST NOT hold ZERO

                    listpmf(0,2,ntpmf)=0

                    Do i=1,mxtpmf(1)
                       listpmf(i,1,ntpmf)=i1pmf(i)
                       If (i1pmf0(i) > 0) Then
                          Call tag_legend(safe(1),i1pmf0(i),ntpmf,legpmf,mxfpmf)
                          If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many PMF type neigbours !!! ***",                           &
  "***           on node      (MPI  rank #): ", idnode,                          &
  "***           requiring a list length of: ", legpmf(0,i1pmf0(i)),             &
  "***           but maximum length allowed: ", mxfpmf,                          &
  "***           on mol. site (local  ID #): ", lstpmf(i,1),                     &
  "***           for particle (global ID #): ", i1pmf(i),                        &
  "***           of PMF unit  (1 or 2 only): ", 1,                               &
  "***           in molecule  (local  ID #): ", imols,                           &
  "***           of type      (       ID #): ", itmols
                       End If
                    End Do
                    If (Any(i1pmf0 > 0)) listpmf(0,2,ntpmf)=listpmf(0,2,ntpmf)+1

                    Do i=1,mxtpmf(2)
                       listpmf(i,2,ntpmf)=i2pmf(i)
                       If (i2pmf0(i) > 0) Then
                          Call tag_legend(safe(1),i2pmf0(i),ntpmf,legpmf,mxfpmf)
                          If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many PMF type neigbours !!! ***",                           &
  "***           on node      (MPI  rank #): ", idnode,                          &
  "***           requiring a list length of: ", legpmf(0,i2pmf0(i)),             &
  "***           but maximum length allowed: ", mxfpmf,                          &
  "***           for particle (global ID #): ", i2pmf(i),                        &
  "***           on mol. site (local  ID #): ", lstpmf(i,2),                     &
  "***           of PMF unit  (1 or 2 only): ", 2,                               &
  "***           in molecule  (local  ID #): ", imols,                           &
  "***           of type      (       ID #): ", itmols
                       End If
                    End Do
                    If (Any(i2pmf0 > 0)) listpmf(0,2,ntpmf)=listpmf(0,2,ntpmf)+2

                 Else
                    safe(4)=.false.
                 End If
              End If

              Deallocate (i1pmf,i1pmf0,i2pmf,i2pmf0, Stat=fail(1))
              If (fail(1) > 0) Then
                 Write(nrite,'(/,1x,a,i0)') 'build_book_intra PMF deallocation failure, node: ', idnode
                 Call error(0)
              End If
           End Do

! Construct RBs list

           Do lrigid=1,numrgd(itmols)
              mrigid=lstrgd(0,lrigid+krigid)

              irgd=0 ; irgd0=0
              Do irigid=1,mrigid
                 irgd(irigid)=lstrgd(irigid,lrigid+krigid)+isite
                 irgd0(irigid)=local_index(irgd(irigid),nlast,lsi,lsa)
                 If (irgd0(irigid) > natms) irgd0(irigid)=0
              End Do

              If (Any(irgd0 > 0)) Then
                 jrigid=jrigid+1
                 If (jrigid <= mxrgd) Then
                    listrgd(-1,jrigid)=mrigid
                    listrgd( 0,jrigid)=lrigid+krigid
                    Do irigid=1,mrigid
                       listrgd(irigid,jrigid)=irgd(irigid)
                       If (irgd0(irigid) > 0) Then
                          Call tag_legend(safe(1),irgd0(irigid),jrigid,legrgd,mxfrgd)
                          If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many PMF type neigbours !!! ***",                           &
  "***           on node      (MPI  rank #): ", idnode,                          &
  "***           requiring a list length of: ", legrgd(0,irgd0(irigid)),         &
  "***           but maximum length allowed: ", mxfrgd,                          &
  "***           for particle (global ID #): ", irgd(irigid),                    &
  "***           on mol. site (local  ID #): ", lstrgd(irigid,lrigid+krigid),    &
  "***           of unit      (local  ID #): ", lrigid,                          &
  "***           in molecule  (local  ID #): ", imols,                           &
  "***           of type      (       ID #): ", itmols
                       End If
                    End Do
                 Else
                    safe(5)=.false.
                 End If
              End If
           End Do

! Construct tethered atoms interaction list

           Do lteths=1,numteth(itmols)
              iatm=lsttet(lteths+kteths)+isite
              iat0=local_index(iatm,nlast,lsi,lsa)
              If (iat0 > natms) iat0=0

              If (iat0 > 0) Then
                 jteths=jteths+1
                 If (jteths <= mxteth) Then
                    listtet(0,jteths)=lteths+kteths
                    listtet(1,jteths)=iatm

                    Call tag_legend(safe(1),iat0,jteths,legtet,mxftet)
                    If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many tether type neigbours !!! ***",                  &
  "***           on node      (MPI  rank #): ", idnode,                    &
  "***           requiring a list length of: ", legtet(0,iat0),            &
  "***           but maximum length allowed: ", mxftet,                    &
  "***           for particle (global ID #): ", iatm,                      &
  "***           on mol. site (local  ID #): ", lsttet(lteths+kteths),     &
  "***           of unit      (local  ID #): ", lteths,                    &
  "***           in molecule  (local  ID #): ", imols,                     &
  "***           of type      (       ID #): ", itmols
                 Else
                    safe(6)=.false.
                 End If
              End If
           End Do

! Construct chemical bond interaction list

           Do lbonds=1,numbonds(itmols)
              iatm=lstbnd(1,lbonds+kbonds)+isite
              jatm=lstbnd(2,lbonds+kbonds)+isite

              iat0=local_index(iatm,nlast,lsi,lsa)
              jat0=local_index(jatm,nlast,lsi,lsa)

              If (iat0 > natms) iat0=0
              If (jat0 > natms) jat0=0

              If (iat0 > 0 .or. jat0 > 0) Then
                 jbonds=jbonds+1
                 If (jbonds <= mxbond) Then
                    listbnd(0,jbonds)=lbonds+kbonds
                    listbnd(1,jbonds)=iatm
                    listbnd(2,jbonds)=jatm

                    If (iat0 > 0) Then
                       Call tag_legend(safe(1),iat0,jbonds,legbnd,mxfbnd)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many bond type neigbours !!! ***",                       &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legbnd(0,iat0),               &
  "***           but maximum length allowed: ", mxfbnd,                       &
  "***           for particle (global ID #): ", iatm,                         &
  "***           on mol. site (local  ID #): ", lstbnd(1,lbonds+kbonds),      &
  "***           of unit      (local  ID #): ", lbonds,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (jat0 > 0) Then
                       Call tag_legend(safe(1),jat0,jbonds,legbnd,mxfbnd)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many bond type neigbours !!! ***",                       &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legbnd(0,jat0),               &
  "***           but maximum length allowed: ", mxfbnd,                       &
  "***           for particle (global ID #): ", jatm,                         &
  "***           on mol. site (local  ID #): ", lstbnd(2,lbonds+kbonds),      &
  "***           of unit      (local  ID #): ", lbonds,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If
                 Else
                    safe(7)=.false.
                 End If
              End If

           End Do

! Construct valence angle interaction list

           Do langle=1,numang(itmols)
              iatm=lstang(1,langle+kangle)+isite
              jatm=lstang(2,langle+kangle)+isite
              katm=lstang(3,langle+kangle)+isite

              iat0=local_index(iatm,nlast,lsi,lsa)
              jat0=local_index(jatm,nlast,lsi,lsa)
              kat0=local_index(katm,nlast,lsi,lsa)

              If (iat0 > natms) iat0=0
              If (jat0 > natms) jat0=0
              If (kat0 > natms) kat0=0

              If (iat0 > 0 .or. jat0 > 0 .or. kat0 > 0) Then
                 jangle=jangle+1
                 If (jangle <= mxangl) Then
                    listang(0,jangle)=langle+kangle
                    listang(1,jangle)=iatm
                    listang(2,jangle)=jatm
                    listang(3,jangle)=katm

                    If (iat0 > 0) Then
                       Call tag_legend(safe(1),iat0,jangle,legang,mxfang)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many angle type neigbours !!! ***",                      &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legang(0,iat0),               &
  "***           but maximum length allowed: ", mxfang,                       &
  "***           for particle (global ID #): ", iatm,                         &
  "***           on mol. site (local  ID #): ", lstang(1,langle+kangle),      &
  "***           of unit      (local  ID #): ", langle,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (jat0 > 0) Then
                       Call tag_legend(safe(1),jat0,jangle,legang,mxfang)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many angle type neigbours !!! ***",                      &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legang(0,jat0),               &
  "***           but maximum length allowed: ", mxfang,                       &
  "***           for particle (global ID #): ", jatm,                         &
  "***           on mol. site (local  ID #): ", lstang(2,langle+kangle),      &
  "***           of unit      (local  ID #): ", langle,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (kat0 > 0) Then
                       Call tag_legend(safe(1),kat0,jangle,legang,mxfang)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many angle type neigbours !!! ***",                      &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legang(0,kat0),               &
  "***           but maximum length allowed: ", mxfang,                       &
  "***           for particle (global ID #): ", katm,                         &
  "***           on mol. site (local  ID #): ", lstang(3,langle+kangle),      &
  "***           of unit      (local  ID #): ", langle,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If
                 Else
                    safe(8)=.false.
                 End If
              End If
           End Do

! Construct dihedral angle interaction list

           Do ldihed=1,numdih(itmols)
              iatm=lstdih(1,ldihed+kdihed)+isite
              jatm=lstdih(2,ldihed+kdihed)+isite
              katm=lstdih(3,ldihed+kdihed)+isite
              latm=lstdih(4,ldihed+kdihed)+isite
              If (lx_dih) Then
                 matm=lstdih(5,ldihed+kdihed)+isite
                 natm=lstdih(6,ldihed+kdihed)+isite
              End If

              iat0=local_index(iatm,nlast,lsi,lsa)
              jat0=local_index(jatm,nlast,lsi,lsa)
              kat0=local_index(katm,nlast,lsi,lsa)
              lat0=local_index(latm,nlast,lsi,lsa)

              If (iat0 > natms) iat0=0
              If (jat0 > natms) jat0=0
              If (kat0 > natms) kat0=0
              If (lat0 > natms) lat0=0

              If (iat0 > 0 .or. jat0 > 0 .or. kat0 > 0 .or. lat0 > 0) Then
                 jdihed=jdihed+1
                 If (jdihed <= mxdihd) Then
                    listdih(0,jdihed)=ldihed+kdihed
                    listdih(1,jdihed)=iatm
                    listdih(2,jdihed)=jatm
                    listdih(3,jdihed)=katm
                    listdih(4,jdihed)=latm
                    If (lx_dih) Then
                       listdih(5,jdihed)=matm
                       listdih(6,jdihed)=natm
                    End If

                    If (iat0 > 0) Then
                       Call tag_legend(safe(1),iat0,jdihed,legdih,mxfdih)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many dihedral type neigbours !!! ***",                   &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legdih(0,iat0),               &
  "***           but maximum length allowed: ", mxfdih,                       &
  "***           for particle (global ID #): ", iatm,                         &
  "***           on mol. site (local  ID #): ", lstdih(1,ldihed+kdihed),      &
  "***           of unit      (local  ID #): ", ldihed,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (jat0 > 0) Then
                       Call tag_legend(safe(1),jat0,jdihed,legdih,mxfdih)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many dihedral type neigbours !!! ***",                   &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legdih(0,jat0),               &
  "***           but maximum length allowed: ", mxfdih,                       &
  "***           for particle (global ID #): ", jatm,                         &
  "***           on mol. site (local  ID #): ", lstdih(2,ldihed+kdihed),      &
  "***           of unit      (local  ID #): ", ldihed,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (kat0 > 0) Then
                       Call tag_legend(safe(1),kat0,jdihed,legdih,mxfdih)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many dihedral type neigbours !!! ***",                   &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legdih(0,kat0),               &
  "***           but maximum length allowed: ", mxfdih,                       &
  "***           for particle (global ID #): ", katm,                         &
  "***           on mol. site (local  ID #): ", lstdih(3,ldihed+kdihed),      &
  "***           of unit      (local  ID #): ", ldihed,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (lat0 > 0) Then
                       Call tag_legend(safe(1),lat0,jdihed,legdih,mxfdih)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many dihedral type neigbours !!! ***",                   &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", legdih(0,lat0),               &
  "***           but maximum length allowed: ", mxfdih,                       &
  "***           for particle (global ID #): ", latm,                         &
  "***           on mol. site (local  ID #): ", lstdih(4,ldihed+kdihed),      &
  "***           of unit      (local  ID #): ", ldihed,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If
                 Else
                    safe(9)=.false.
                 End If
              End If
           End Do

! Construct inversion potential interaction list

           Do linver=1,numinv(itmols)
              iatm=lstinv(1,linver+kinver)+isite
              jatm=lstinv(2,linver+kinver)+isite
              katm=lstinv(3,linver+kinver)+isite
              latm=lstinv(4,linver+kinver)+isite

              iat0=local_index(iatm,nlast,lsi,lsa)
              jat0=local_index(jatm,nlast,lsi,lsa)
              kat0=local_index(katm,nlast,lsi,lsa)
              lat0=local_index(latm,nlast,lsi,lsa)

              If (iat0 > natms) iat0=0
              If (jat0 > natms) jat0=0
              If (kat0 > natms) kat0=0
              If (lat0 > natms) lat0=0

              If (iat0 > 0 .or. jat0 > 0 .or. kat0 > 0 .or. lat0 > 0) Then
                 jinver=jinver+1
                 If (jinver <= mxinv) Then
                    listinv(0,jinver)=linver+kinver
                    listinv(1,jinver)=iatm
                    listinv(2,jinver)=jatm
                    listinv(3,jinver)=katm
                    listinv(4,jinver)=latm

                    If (iat0 > 0) Then
                       Call tag_legend(safe(1),iat0,jinver,leginv,mxfinv)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many inversion type neigbours !!! ***",                  &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", leginv(0,iat0),               &
  "***           but maximum length allowed: ", mxfinv,                       &
  "***           for particle (global ID #): ", iatm,                         &
  "***           on mol. site (local  ID #): ", lstinv(1,linver+kinver),      &
  "***           of unit      (local  ID #): ", linver,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (jat0 > 0) Then
                       Call tag_legend(safe(1),jat0,jinver,leginv,mxfinv)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many inversion type neigbours !!! ***",                  &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", leginv(0,jat0),               &
  "***           but maximum length allowed: ", mxfinv,                       &
  "***           for particle (global ID #): ", jatm,                         &
  "***           on mol. site (local  ID #): ", lstinv(2,linver+kinver),      &
  "***           of unit      (local  ID #): ", linver,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (kat0 > 0) Then
                       Call tag_legend(safe(1),kat0,jinver,leginv,mxfinv)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many inversion type neigbours !!! ***",                  &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", leginv(0,kat0),               &
  "***           but maximum length allowed: ", mxfinv,                       &
  "***           for particle (global ID #): ", katm,                         &
  "***           on mol. site (local  ID #): ", lstinv(3,linver+kinver),      &
  "***           of unit      (local  ID #): ", linver,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If

                    If (lat0 > 0) Then
                       Call tag_legend(safe(1),lat0,jinver,leginv,mxfinv)
                       If (.not.safe(1)) Write(nrite,'(/,1x,a,8(/,1x,a,i0))') &
  "*** warning - too many inversion type neigbours !!! ***",                  &
  "***           on node      (MPI  rank #): ", idnode,                       &
  "***           requiring a list length of: ", leginv(0,lat0),               &
  "***           but maximum length allowed: ", mxfinv,                       &
  "***           for particle (global ID #): ", latm,                         &
  "***           on mol. site (local  ID #): ", lstinv(4,linver+kinver),      &
  "***           of unit      (local  ID #): ", linver,                       &
  "***           in molecule  (local  ID #): ", imols,                        &
  "***           of type      (       ID #): ", itmols
                    End If
                 Else
                    safe(10)=.false.
                 End If
              End If
           End Do

        End If

        isite=isite+numsit(itmols)
        nsatm=neatm

     End Do

! Update global unit numbers for all passed molecules so far

     kshels=kshels+numshl(itmols)

     kconst=kconst+numcon(itmols)
! No 'kpmf' needed since PMF is defined on one and only one molecular type

     krigid=krigid+numrgd(itmols)

     kteths=kteths+numteth(itmols)

     kbonds=kbonds+numbonds(itmols)
     kangle=kangle+numang(itmols)
     kdihed=kdihed+numdih(itmols)
     kinver=kinver+numinv(itmols)

  End Do

! Store array counters for bookkeeping

  ntshl =jshels

  ntcons=jconst
! 'ntpmf' is updated locally as PMFs are global and one type only

  ntrgd =jrigid

  ntteth=jteths

  ntbond=jbonds
  ntangl=jangle
  ntdihd=jdihed
  ntinv =jinver

  If (megshl == 0) Then
     ntshl1 =ntshl

     ntcons1=ntcons

     ntrgd1 =ntrgd

     ntbond1=ntbond
     ntangl1=ntangl
     ntdihd1=ntdihd
     ntinv1 =ntinv

     ntshl2 =ntshl1

     Go To 400
  End If

! Cycle through all constraint, RB, bond, angle, dihedral and inversion
! units on this node and record the non-local index particles

  iwrk=0
  mshels=0
  Do i=1,ntcons
     iatm=listcon(1,i)
     jatm=listcon(2,i)

     iat0=local_index(iatm,nlast,lsi,lsa)
     jat0=local_index(jatm,nlast,lsi,lsa)

     If (iat0 > natms .and. (.not.Any(iwrk(1:mshels) == iatm))) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (jat0 > natms .and. (.not.Any(iwrk(1:mshels) == jatm))) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If
  End Do
  Do i=1,ntrgd
     mrigid=listrgd(-1,i)

     Do lrigid=1,mrigid
        iatm=listrgd(lrigid,i)
        iat0=local_index(iatm,nlast,lsi,lsa)

        If (iat0 > natms .and. (.not.Any(iwrk(1:mshels) == iatm))) Then
           mshels=mshels+1
           iwrk(mshels)=iatm
        End If
     End Do
  End Do
  Do i=1,ntbond
     iatm=listbnd(1,i)
     jatm=listbnd(2,i)

     iat0=local_index(iatm,nlast,lsi,lsa)
     jat0=local_index(jatm,nlast,lsi,lsa)

     If (iat0 > natms .and. (.not.Any(iwrk(1:mshels) == iatm))) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (jat0 > natms .and. (.not.Any(iwrk(1:mshels) == jatm))) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If
  End Do
  Do i=1,ntangl
     iatm=listang(1,i)
     jatm=listang(2,i)
     katm=listang(3,i)

     iat0=local_index(iatm,nlast,lsi,lsa)
     jat0=local_index(jatm,nlast,lsi,lsa)
     kat0=local_index(katm,nlast,lsi,lsa)

     If (iat0 > natms .and. (.not.Any(iwrk(1:mshels) == iatm))) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (jat0 > natms .and. (.not.Any(iwrk(1:mshels) == jatm))) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If

     If (kat0 > natms .and. (.not.Any(iwrk(1:mshels) == katm))) Then
        mshels=mshels+1
        iwrk(mshels)=katm
     End If
  End Do
  Do i=1,ntdihd
     iatm=listdih(1,i)
     jatm=listdih(2,i)
     katm=listdih(3,i)
     latm=listdih(4,i)

     iat0=local_index(iatm,nlast,lsi,lsa)
     jat0=local_index(jatm,nlast,lsi,lsa)
     kat0=local_index(katm,nlast,lsi,lsa)
     lat0=local_index(latm,nlast,lsi,lsa)

     If (iat0 > natms .and. (.not.Any(iwrk(1:mshels) == iatm))) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (jat0 > natms .and. (.not.Any(iwrk(1:mshels) == jatm))) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If

     If (kat0 > natms .and. (.not.Any(iwrk(1:mshels) == katm))) Then
        mshels=mshels+1
        iwrk(mshels)=katm
     End If
     If (lat0 > natms .and. (.not.Any(iwrk(1:mshels) == latm))) Then
        mshels=mshels+1
        iwrk(mshels)=latm
     End If
  End Do
  Do i=1,ntinv
     iatm=listinv(1,i)
     jatm=listinv(2,i)
     katm=listinv(3,i)
     latm=listinv(4,i)

     iat0=local_index(iatm,nlast,lsi,lsa)
     jat0=local_index(jatm,nlast,lsi,lsa)
     kat0=local_index(katm,nlast,lsi,lsa)
     lat0=local_index(latm,nlast,lsi,lsa)

     If (iat0 > natms .and. (.not.Any(iwrk(1:mshels) == iatm))) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (jat0 > natms .and. (.not.Any(iwrk(1:mshels) == jatm))) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If

     If (kat0 > natms .and. (.not.Any(iwrk(1:mshels) == katm))) Then
        mshels=mshels+1
        iwrk(mshels)=katm
     End If
     If (lat0 > natms .and. (.not.Any(iwrk(1:mshels) == latm))) Then
        mshels=mshels+1
        iwrk(mshels)=latm
     End If
  End Do

  If (mshels == 0) Go To 100

! Include in local core-shell units' description
! non-local units that are connected to partly shared
! constraints, RBs, bonds, angles, dihedrals and inversions

  isite=0                            ! initialise bookkeeping indices
  kshels=0                           ! last index of core-shell unit
  nsatm =0                           ! global atom counter
  Do itmols=1,ntpmls                 ! loop over molecule types in the system
     Do imols=1,nummols(itmols)      ! loop over molecules of this type
        neatm=nsatm+numsit(itmols)   ! last atom in the molecule

! From the first till the last atom of this molecule, is there a
! non-local atom iwrk(1:mshels)

        i=0                          ! There is not
        Do iatm=nsatm+1,neatm
           If (Any(iwrk(1:mshels) == iatm)) i=i+1
        End Do

! If there is a non-local atom iwrk(1:mshels) on this
! molecule on this node, extend listshl

        If (i > 0) Then

! Extend core-shell units interaction list

           Do lshels=1,numshl(itmols)
              iatm=lstshl(1,lshels+kshels)+isite
              jatm=lstshl(2,lshels+kshels)+isite

              If ( Any(iwrk(1:mshels) == iatm) .or. &
                   Any(iwrk(1:mshels) == jatm) ) Then
                 jshels=jshels+1
                 If (jshels <= mxshl) Then
                    listshl(0,jshels)=lshels+kshels
                    listshl(1,jshels)=iatm
                    listshl(2,jshels)=jatm
                 Else
                    safe(2)=.false.
                 End If
              End If
           End Do

        End If

        isite=isite+numsit(itmols)
        nsatm=neatm
     End Do

! Update core-shell units number for all passed molecules so far

     kshels=kshels+numshl(itmols)
  End Do

100 Continue

! Store first (local+non-local) extended array counter
! for bookkeeping and exclusion of core-shell units

  ntshl1 =jshels

! Cycle through all local and partly shared core-shell units on
! this node and record the non-local indices of cross-domained
! core-shell unit particles

  iwrk=0
  mshels=0
  Do i=1,ntshl
     iatm=listshl(1,i)
     jatm=listshl(2,i)

     iat0=local_index(iatm,nlast,lsi,lsa) ! This is a core
     jat0=local_index(jatm,nlast,lsi,lsa) ! This is a shell

     If (iat0 == 0 .and. jat0 == 0) safe(11)=.false.

! Core is out

     If (iat0 > natms) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

! Shell is out

     If (jat0 > natms) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If
  End Do

  If (mshels == 0) Go To 200

! Include in (local constraint, RB, bond, angle, dihedral
! and inversion) units' descriptions non-local units that
! are connected to partly shared core-shell units

  isite=0                            ! initialise bookkeeping indices
  kconst=0                           ! last index of constraint unit
  krigid=0                           ! last index of RB unit
  kbonds=0                           ! last index of bond unit
  kangle=0                           ! last index of angle unit
  kdihed=0                           ! last index of dihedral unit
  kinver=0                           ! last index of inversion unit
  nsatm =0                           ! global atom counter
  Do itmols=1,ntpmls                 ! loop over molecule types in the system
     Do imols=1,nummols(itmols)      ! loop over molecules of this type
        neatm=nsatm+numsit(itmols)   ! last atom in the molecule

! From the first till the last atom of this molecule, is there a
! non-local, cross-domained core-shell unit atom

        i=0                          ! There is not
        Do iatm=nsatm+1,neatm
           If (Any(iwrk(1:mshels) == iatm)) i=i+1
        End Do

! If there is a non-local, cross-domained core-shell unit atom on this
! molecule on this node, extend listcon, listrgd, listbnd and listang

        If (i > 0) Then

! Extend constraint bond list

           Do lconst=1,numcon(itmols)
              iatm=lstcon(1,lconst+kconst)+isite
              jatm=lstcon(2,lconst+kconst)+isite

              If ( Any(iwrk(1:mshels) == iatm) .or. &
                   Any(iwrk(1:mshels) == jatm) ) Then
                 jconst=jconst+1
                 If (jconst <= mxcons) Then
                    listcon(0,jconst)=lconst+kconst
                    listcon(1,jconst)=iatm
                    listcon(2,jconst)=jatm
                 Else
                    safe(3)=.false.
                 End If
              End If
           End Do

! Extend RB list

           Do lrigid=1,numrgd(itmols)
              mrigid=lstrgd(0,lrigid+krigid)

              irgd=0 ; irgd0=0 ; go=.false.
              Do irigid=1,mrigid
                 irgd(irigid)=lstrgd(irigid,lrigid+krigid)+isite
                 go=(go .or. Any(iwrk(1:mshels) == irgd(irigid)))
              End Do

              If (go) Then
                 jrigid=jrigid+1
                 If (jrigid <= mxrgd) Then
                    listrgd(-1,jrigid)=lstrgd(0,lrigid+krigid)
                    listrgd( 0,jrigid)=lrigid+krigid
                    Do irigid=1,mrigid
                       listrgd(irigid,jrigid)=irgd(irigid)
                    End Do
                 Else
                    safe(5)=.false.
                 End If
              End If
           End Do

! Extend chemical bond interaction list

           Do lbonds=1,numbonds(itmols)
              iatm=lstbnd(1,lbonds+kbonds)+isite
              jatm=lstbnd(2,lbonds+kbonds)+isite

              If ( Any(iwrk(1:mshels) == iatm) .or. &
                   Any(iwrk(1:mshels) == jatm) ) Then
                 jbonds=jbonds+1
                 If (jbonds <= mxbond) Then
                    listbnd(0,jbonds)=lbonds+kbonds
                    listbnd(1,jbonds)=iatm
                    listbnd(2,jbonds)=jatm
                 Else
                    safe(7)=.false.
                 End If
              End If
           End Do

! Extend valence angle interaction list

           Do langle=1,numang(itmols)
              iatm=lstang(1,langle+kangle)+isite
              jatm=lstang(2,langle+kangle)+isite
              katm=lstang(3,langle+kangle)+isite

              If ( Any(iwrk(1:mshels) == iatm) .or. &
                   Any(iwrk(1:mshels) == jatm) .or. &
                   Any(iwrk(1:mshels) == katm) ) Then
                 jangle=jangle+1
                 If (jangle <= mxangl) Then
                    listang(0,jangle)=langle+kangle
                    listang(1,jangle)=iatm
                    listang(2,jangle)=jatm
                    listang(3,jangle)=katm
                 Else
                    safe(8)=.false.
                 End If
              End If
           End Do

! Extend dihedral angle interaction list

           Do ldihed=1,numdih(itmols)
              iatm=lstdih(1,ldihed+kdihed)+isite
              jatm=lstdih(2,ldihed+kdihed)+isite
              katm=lstdih(3,ldihed+kdihed)+isite
              latm=lstdih(4,ldihed+kdihed)+isite
              If (lx_dih) Then
                 matm=lstdih(5,ldihed+kdihed)+isite
                 natm=lstdih(6,ldihed+kdihed)+isite
              End If

              If ( Any(iwrk(1:mshels) == iatm) .or. &
                   Any(iwrk(1:mshels) == jatm) .or. &
                   Any(iwrk(1:mshels) == katm) .or. &
                   Any(iwrk(1:mshels) == latm) ) Then
                 jdihed=jdihed+1
                 If (jdihed <= mxdihd) Then
                    listdih(0,jdihed)=ldihed+kdihed
                    listdih(1,jdihed)=iatm
                    listdih(2,jdihed)=jatm
                    listdih(3,jdihed)=katm
                    listdih(4,jdihed)=latm
                    If (lx_dih) Then
                       listdih(5,jdihed)=matm
                       listdih(6,jdihed)=natm
                    End If
                 Else
                    safe(9)=.false.
                 End If
              End If
           End Do

! Extend inversion potential interaction list

           Do linver=1,numinv(itmols)
              iatm=lstinv(1,linver+kinver)+isite
              jatm=lstinv(2,linver+kinver)+isite
              katm=lstinv(3,linver+kinver)+isite
              latm=lstinv(4,linver+kinver)+isite

              If ( Any(iwrk(1:mshels) == iatm) .or. &
                   Any(iwrk(1:mshels) == jatm) .or. &
                   Any(iwrk(1:mshels) == katm) .or. &
                   Any(iwrk(1:mshels) == latm) ) Then
                 jinver=jinver+1
                 If (jinver <= mxinv) Then
                    listinv(0,jinver)=linver+kinver
                    listinv(1,jinver)=iatm
                    listinv(2,jinver)=jatm
                    listinv(3,jinver)=katm
                    listinv(4,jinver)=latm
                 Else
                    safe(10)=.false.
                 End If
              End If
           End Do

        End If

        isite=isite+numsit(itmols)
        nsatm=neatm
     End Do

! Update constraint, RB, bond, angle, dihedral and inversion
! units numbers for all passed molecules so far

     kconst=kconst+numcon(itmols)

     krigid=krigid+numrgd(itmols)

     kbonds=kbonds+numbonds(itmols)
     kangle=kangle+numang(itmols)
     kdihed=kdihed+numdih(itmols)
     kinver=kinver+numinv(itmols)
  End Do

200 Continue

! Store first extended array counters for bookkeeping and exclusion
! of constraint, RB, bond, angle, dihedral and inversion units

  ntcons1=jconst

  ntrgd1 =jrigid

  ntbond1=jbonds
  ntangl1=jangle
  ntdihd1=jdihed
  ntinv1 =jinver

! Cycle through the extended -
! constraint, RB, bond, angle, dihedral and inversion units
! on this node and record the non-local index particles

  iwrk=0
  mshels=0
  Do i=ntcons+1,ntcons1
     iatm=listcon(1,i)
     jatm=listcon(2,i)

     If (.not.Any(iwrk(1:mshels) == iatm)) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (.not.Any(iwrk(1:mshels) == jatm)) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If
  End Do
  Do i=ntrgd+1,ntrgd1
     mrigid=listrgd(-1,i)

     Do j=1,mrigid
        iatm=listrgd(j,i)
        If (.not.Any(iwrk(1:mshels) == iatm)) Then
           mshels=mshels+1
           iwrk(mshels)=iatm
        End If
     End Do
  End Do
  Do i=ntbond+1,ntbond1
     iatm=listbnd(1,i)
     jatm=listbnd(2,i)

     If (.not.Any(iwrk(1:mshels) == iatm)) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (.not.Any(iwrk(1:mshels) == jatm)) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If
  End Do
  Do i=ntangl+1,ntangl1
     iatm=listang(1,i)
     jatm=listang(2,i)
     katm=listang(3,i)

     If (.not.Any(iwrk(1:mshels) == iatm)) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (.not.Any(iwrk(1:mshels) == jatm)) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If

     If (.not.Any(iwrk(1:mshels) == katm)) Then
        mshels=mshels+1
        iwrk(mshels)=katm
     End If
  End Do
  Do i=ntdihd+1,ntdihd1
     iatm=listdih(1,i)
     jatm=listdih(2,i)
     katm=listdih(3,i)
     latm=listdih(4,i)

     If (.not.Any(iwrk(1:mshels) == iatm)) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (.not.Any(iwrk(1:mshels) == jatm)) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If

     If (.not.Any(iwrk(1:mshels) == katm)) Then
        mshels=mshels+1
        iwrk(mshels)=katm
     End If

     If (.not.Any(iwrk(1:mshels) == latm)) Then
        mshels=mshels+1
        iwrk(mshels)=latm
     End If
  End Do
  Do i=ntinv+1,ntinv1
     iatm=listinv(1,i)
     jatm=listinv(2,i)
     katm=listinv(3,i)
     latm=listinv(4,i)

     If (.not.Any(iwrk(1:mshels) == iatm)) Then
        mshels=mshels+1
        iwrk(mshels)=iatm
     End If

     If (.not.Any(iwrk(1:mshels) == jatm)) Then
        mshels=mshels+1
        iwrk(mshels)=jatm
     End If

     If (.not.Any(iwrk(1:mshels) == katm)) Then
        mshels=mshels+1
        iwrk(mshels)=katm
     End If

     If (.not.Any(iwrk(1:mshels) == latm)) Then
        mshels=mshels+1
        iwrk(mshels)=latm
     End If
  End Do

  If (mshels == 0) Go To 300

! Include in local and non-local - connected to partly shared
! constraints, RBs, bonds, angles, dihedrals or inversions -
! core-shell unit description of foreign units that are connected to
! non-local constraints, RBs, bonds, angles, dihedrals or inversions
! by partly shared core-shell units

  isite=0                            ! initialise bookkeeping indices
  kshels=0                           ! last index of core-shell unit
  nsatm =0                           ! global atom counter
  Do itmols=1,ntpmls                 ! loop over molecule types in the system
     Do imols=1,nummols(itmols)      ! loop over molecules of this type
        neatm=nsatm+numsit(itmols)   ! last atom in the molecule

! From the first till the last atom of this molecule, is there a
! non-local atom iwrk(1:mshels)

        i=0                          ! There is not
        Do iatm=nsatm+1,neatm
           If (Any(iwrk(1:mshels) == iatm)) i=i+1
        End Do

! If there is a non-local atom iwrk(1:mshels) on this
! molecule on this node, extend listshl

        If (i > 0) Then

! Extend core-shell units interaction list

           Do lshels=1,numshl(itmols)
              iatm=lstshl(1,lshels+kshels)+isite
              jatm=lstshl(2,lshels+kshels)+isite

              If ( Any(iwrk(1:mshels) == iatm) .or. &
                   Any(iwrk(1:mshels) == jatm) ) Then
                 jshels=jshels+1
                 If (jshels <= mxshl) Then
                    listshl(0,jshels)=lshels+kshels
                    listshl(1,jshels)=iatm
                    listshl(2,jshels)=jatm
                 Else
                    safe(2)=.false.
                 End If
              End If
           End Do

        End If

        isite=isite+numsit(itmols)
        nsatm=neatm
     End Do

! Update core-shell units number for all passed molecules so far

     kshels=kshels+numshl(itmols)
  End Do

300 Continue

! Store second (local+non-local+foreign) extended array counter
! for bookkeeping and exclusion of core-shell units

  ntshl2 =jshels

  If (lx_dih) ntdihd=ntdihd1 ! extend the dihedrals' set

400 Continue

! error exit for all error conditions (size of work arrays)

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe( 1)) Call error( 88)
  If (.not.safe( 2)) Call error( 59)
  If (.not.safe( 3)) Call error( 41)
  If (.not.safe( 4)) Call error(488)
  If (.not.safe( 5)) Call error(640)
  If (.not.safe( 6)) Call error( 63)
  If (.not.safe( 7)) Call error( 31)
  If (.not.safe( 8)) Call error( 51)
  If (.not.safe( 9)) Call error( 61)
  If (.not.safe(10)) Call error( 77)
  If (.not.safe(11)) Call error( 64)

  Deallocate (iwrk,                      Stat=fail(1))
  If (m_rgd > 0) Deallocate (irgd,irgd0, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'build_book_intra deallocation failure, node: ', idnode
     Call error(0)
  End If

! Set RB particulars and quaternions

  If (m_rgd > 0) Call rigid_bodies_setup(megatm,megfrz,megrgd,degtra,degrot)

  Call report_topology                   &
           (megatm,megfrz,atmfre,atmfrz, &
           megshl,megcon,megpmf,megrgd,  &
           megtet,megbnd,megang,megdih,meginv)

! DEALLOCATE INTER-LIKE SITE INTERACTION ARRAYS since no longer needed

  Call deallocate_core_shell_arrays()

  Call deallocate_constraints_arrays()
  Call deallocate_pmf_arrays()

  Call deallocate_rigid_bodies_arrays()

  Call deallocate_tethers_arrays()

  Call deallocate_bonds_arrays()
  Call deallocate_angles_arrays()
  Call deallocate_dihedrals_arrays()
  Call deallocate_inversions_arrays()

! Update shared core-shell, constraint and RB units
! (pmf data updated by construction)

  If (megshl > 0 .and. mxnode > 1) Call pass_shared_units &
     (mxshl, Lbound(listshl,Dim=1),Ubound(listshl,Dim=1),ntshl, listshl,mxfshl,legshl,lshmv_shl,lishp_shl,lashp_shl)

  If (m_con > 0 .and. mxnode > 1) Call pass_shared_units &
     (mxcons,Lbound(listcon,Dim=1),Ubound(listcon,Dim=1),ntcons,listcon,mxfcon,legcon,lshmv_con,lishp_con,lashp_con)

  If (m_rgd > 0 .and. mxnode > 1) Call pass_shared_units &
     (mxrgd, Lbound(listrgd,Dim=1),Ubound(listrgd,Dim=1),ntrgd, listrgd,mxfrgd,legrgd,lshmv_rgd,lishp_rgd,lashp_rgd)

End Subroutine build_book_intra
