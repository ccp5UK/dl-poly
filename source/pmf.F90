Module pmf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for defining global PMF constraints' variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use comms,           Only : comms_type,gcheck,gsum,gsync
  Use setup
  Use configuration,   Only : imcon,cell,natms,xxx,yyy,zzz,lfrzn, &
                              vxx,vyy,vzz,nlast,lsi,lsa
  Use errors_warnings, Only : error,warning,info
  Use numerics,        Only : images,local_index,dcell
  Implicit None

  Integer,                        Save :: ntpmf  = 0

  Real( Kind = wp ),              Save :: prmpmf = 0.0_wp
  Real( Kind = wp ),              Save :: passpmq(1:5) = (/ & ! QUENCHING per call
                                          0.0_wp         ,  & ! cycles counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average cycles
                                          999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
                                          0.0_wp /)           ! maximum cycles
  Real( Kind = wp ),              Save :: passpmf(1:5,1:2,1:2) = Reshape( (/ & ! dim::1-shake, dim:1:-per-call
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::1-shake, dim:2:-per-tst
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:1:-per-call
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:2:-per-tst
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp /) , (/5,2,2/) )


  Integer,           Allocatable, Save :: numpmf(:),pmffrz(:)
  Integer,           Allocatable, Save :: lstpmf(:,:),listpmf(:,:,:),legpmf(:,:)

  Real( Kind = wp ), Allocatable, Save :: pmfwgt(:,:),pmfwg1(:,:)

  Public :: allocate_pmf_arrays , deallocate_pmf_arrays

Contains

  Subroutine allocate_pmf_arrays()

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (numpmf(1:mxtmls),pmffrz(1:2),                    Stat = fail(1))
    Allocate (lstpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2),          Stat = fail(2))
    Allocate (listpmf(0:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf), Stat = fail(3))
    Allocate (legpmf(0:mxfpmf,1:mxatdm),                       Stat = fail(4))
    Allocate (pmfwgt(0:Max(mxtpmf(1),mxtpmf(2)),1:2),          Stat = fail(5))
    Allocate (pmfwg1(0:Max(mxtpmf(1),mxtpmf(2)),1:2),          Stat = fail(6))

    If (Any(fail > 0)) Call error(1036)

    numpmf  = 0
    pmffrz  = 0
    lstpmf  = 0
    listpmf = 0
    legpmf  = 0

    pmfwgt = 0.0_wp
    pmfwg1 = 0.0_wp

  End Subroutine allocate_pmf_arrays

  Subroutine deallocate_pmf_arrays()

    Integer :: fail

    fail = 0

    Deallocate (numpmf,lstpmf, Stat = fail)

    If (fail > 0) Call error(1037)

  End Subroutine deallocate_pmf_arrays
  
  Subroutine pmf_coms(indpmf,pxx,pyy,pzz,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing PMF units' COM vectors for
! iterative (pmf_quench,pmf_shake,pmf_rattle) constraint algorithms
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)

  Real( Kind = wp ), Intent(   Out ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical                 :: safe(1:2)

  Real( Kind = wp )       :: celprp(1:10),width,xmin,xmax,ymin,ymax,zmin,zmax

  Integer                 :: fail(1:3),gpmf,gpmf1,gpmf2,ipmf,jpmf,iadd, &
                             j,k,l,m

  Real( Kind = wp ), Dimension( : ),    Allocatable :: xxt,yyt,zzt
  Real( Kind = wp ), Dimension( :, : ), Allocatable :: xpmf,ypmf,zpmf
  Real( Kind = wp ), Dimension( : ),    Allocatable :: buffer

  Character( Len = 256 ) :: message
  fail=0
  Allocate (xxt(1:Max(mxtpmf(1),mxtpmf(2))),yyt(1:Max(mxtpmf(1),mxtpmf(2))),zzt(1:Max(mxtpmf(1),mxtpmf(2))), Stat=fail(1))
  Allocate (xpmf(1:2,1:mxpmf),ypmf(1:2,1:mxpmf),zpmf(1:2,1:mxpmf),                                           Stat=fail(2))
  Allocate (buffer(1:(mxtpmf(1)+mxtpmf(2))*(mxpmf+2)),                                                       Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_coms allocation failure'
     Call error(0,message)
  End If


! Initialise safety flags

  safe=.true.

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)
  width=Min(celprp(7),celprp(8),celprp(9))

! Initialise PMF COMs' and inter-COMs' vector arrays

  xpmf = 0.0_wp ; ypmf = 0.0_wp ; zpmf = 0.0_wp
  pxx  = 0.0_wp ; pyy  = 0.0_wp ; pzz  = 0.0_wp

! Loop over all global PMF constraints

  gpmf1=1          ! number of passed global PMFs
  buffer=0.0_wp    ! coordinates buffer
  iadd=3           ! adding three coordinates only per particle
  Do gpmf=1,mxpmf

! Loop over all local to this node PMF constraints matching the global one

     Do ipmf=1,ntpmf
        If (listpmf(0,1,ipmf) == gpmf) Then

! Loop over all PMF units present on my domain (no halo)

           Do jpmf=1,2
              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then

! Loop over their members present on my domain (no halo)

                 Do k=1,mxtpmf(jpmf)
                    j=indpmf(k,jpmf,ipmf)

! Copy particles' coordinates to buffer in an orderly manner

                    If (j > 0 .and. j <= natms) Then ! j is a domain particle
                       l=((gpmf-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1)+(k-1))*iadd
                       buffer(l+1)=xxx(j)
                       buffer(l+2)=yyy(j)
                       buffer(l+3)=zzz(j)
                    End If
                 End Do
              End If
           End Do
        End If
     End Do

! Check if it safe to fill up the buffer

     safe(1)=(gpmf-gpmf1+2 < (mxpmf+2)/3)

! If not safe or we've finished looping over all global PMFs

     If ((.not.safe(1)) .or. gpmf == mxpmf) Then
        Call gsum(comm,buffer)

        Do gpmf2=gpmf1,gpmf
           Do ipmf=1,ntpmf
              If (listpmf(0,1,ipmf) == gpmf2) Then
                 Do jpmf=1,2

! Get in the local scope of the unit

                    m=((gpmf2-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1))*iadd
                    Do k=1,mxtpmf(jpmf)
                       l=m+(k-1)*iadd

                       xxt(k)=buffer(l+1)-buffer(m+1)
                       yyt(k)=buffer(l+2)-buffer(m+2)
                       zzt(k)=buffer(l+3)-buffer(m+3)
                    End Do

                    Call images(imcon,cell,mxtpmf(jpmf),xxt,yyt,zzt)

                    xmin=0.0_wp ; xmax = 0.0_wp
                    ymin=0.0_wp ; ymax = 0.0_wp
                    zmin=0.0_wp ; zmax = 0.0_wp
                    Do k=1,mxtpmf(jpmf)
                       xmin=Min(xmin,xxt(k)) ; xmax=Max(xmax,xxt(k))
                       ymin=Min(ymin,yyt(k)) ; ymax=Max(ymax,yyt(k))
                       zmin=Min(zmin,zzt(k)) ; zmax=Max(zmax,zzt(k))
                       safe(2)=safe(2) .and. (Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2) < width/2.0_wp)
                    End Do
                    safe(2)=safe(2) .and. (xmax-xmin < width/2.0_wp) &
                                    .and. (ymax-ymin < width/2.0_wp) &
                                    .and. (zmax-zmin < width/2.0_wp)

! Get the COM of this unit

                    Do k=1,mxtpmf(jpmf)
                       xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf) + pmfwgt(k,jpmf)*xxt(k)
                       ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf) + pmfwgt(k,jpmf)*yyt(k)
                       zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf) + pmfwgt(k,jpmf)*zzt(k)
                    End Do

! Get out of local frame

                    xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf)*pmfwgt(0,jpmf) + buffer(m+1)
                    ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf)*pmfwgt(0,jpmf) + buffer(m+2)
                    zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf)*pmfwgt(0,jpmf) + buffer(m+3)

                 End Do
              End If
           End Do
        End Do

        gpmf1=gpmf2
        buffer=0.0_wp
     End If

! Check if a PMF unit has a diameter > of the minimum of all half-cell width

     Call gcheck(comm,safe(2))
     If (.not.safe(2)) Call error(492)

   End Do

! Loop over all PMF constraints on this node and calculate PMF constraint vectors

  Do ipmf=1,ntpmf
     pxx(ipmf) = xpmf(2,ipmf) - xpmf(1,ipmf)
     pyy(ipmf) = ypmf(2,ipmf) - ypmf(1,ipmf)
     pzz(ipmf) = zpmf(2,ipmf) - zpmf(1,ipmf)
  End Do

! Minimum image convention for bond vectors

  Call images(imcon,cell,ntpmf,pxx,pyy,pzz)

  Deallocate (xxt,yyt,zzt,    Stat=fail(1))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(2))
  Deallocate (buffer,         Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_coms deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_coms

Subroutine pmf_pseudo_bonds(indpmf,pxx,pyy,pzz,gxx,gyy,gzz,engpmf,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for treating PMF constraints as stiff harmonic
! springs for use with the conjugate gradient method (minimise_relax.f90)
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov august 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent( In    ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: engpmf
  Type( comms_type), Intent( InOut ) :: comm

  Real( Kind = wp ), Parameter :: rigid=1.0e6_wp

  Integer           :: ipmf,jpmf,k,l
  Real( Kind = wp ) :: r,r0,ebond,gamma,tmp

  r0=prmpmf

  engpmf=0.0_wp
  Do ipmf=1,ntpmf

     r=Sqrt(pxx(ipmf)**2+pyy(ipmf)**2+pzz(ipmf)**2)

     gamma=rigid*(r-r0)/Real(mxtpmf(1)+mxtpmf(2),wp)
     ebond=gamma*0.5_wp*(r-r0)
     gamma=gamma/r

     Do jpmf=1,2
        tmp=Real(1-2*Mod(jpmf,2),wp)*gamma

! If this unit is present on my domain

        If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then

           Do k=1,mxtpmf(jpmf)
              l=indpmf(k,jpmf,ipmf)

! For domain particles

              If (l > 0 .and. l <= natms) Then

! Accumulate energy

                 engpmf=engpmf+ebond

! Add forces

                 If (lfrzn(l) == 0) Then
                    gxx(l)=gxx(l)+pxx(ipmf)*tmp
                    gyy(l)=gyy(l)+pyy(ipmf)*tmp
                    gzz(l)=gzz(l)+pzz(ipmf)*tmp
                 End If

              End If

           End Do

        End If

     End Do

  End Do

! Global sum of energy

  Call gsum(comm,engpmf)

End Subroutine pmf_pseudo_bonds

Subroutine pmf_quench(mxshak,tolnce,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for applying PMF constraint quench
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce
  Type( comms_type), intent( InOut ) :: comm

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: amt(1:2)

  Logical                 :: safe
  Integer                 :: fail(1:5),ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: dis,esig,gamma,gamm(1:2)

  Logical,           Allocatable :: lstitr(:)
  Integer,           Allocatable :: indpmf(:,:,:)
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: xpmf(:,:),ypmf(:,:),zpmf(:,:)
  Character( Len = 256 ) :: message

  fail=0
  Allocate (lstitr(1:mxatms),                                      Stat=fail(1))
  Allocate (indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf),        Stat=fail(2))
  Allocate (pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf),                Stat=fail(3))
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),             Stat=fail(4))
  Allocate (xpmf(1:2,1:mxpmf),ypmf(1:2,1:mxpmf),zpmf(1:2,1:mxpmf), Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_quench allocation failure'
     Call error(0,message)
  End If

! Get PMF units' reciprocal masses

  If (newjob) Then
     newjob = .false.

     Do jpmf=1,2
        If (pmffrz(jpmf) == mxtpmf(jpmf)) Then
           amt(jpmf)=0.0_wp
        Else
           amt(jpmf)=pmfwg1(0,jpmf)
        End If
     End Do
  End If

  lstitr(1:natms)=.false. ! initialise lstitr
  Call pmf_tags(lstitr,indpmf,pxx,pyy,pzz,comm)

! normalise PMF constraint vectors

  Do ipmf=1,ntpmf
     dis=1.0_wp/Sqrt(pxx(ipmf)**2+pyy(ipmf)**2+pzz(ipmf)**2)
     pxx(ipmf)=pxx(ipmf)*dis
     pyy(ipmf)=pyy(ipmf)*dis
     pzz(ipmf)=pzz(ipmf)*dis
  End Do

! application of PMF constraint (rattle) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  icyc=0
  safe=.false.

  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! initialise velocity correction arrays

     Do l=1,natms
        vxt(l)=0.0_wp
        vyt(l)=0.0_wp
        vzt(l)=0.0_wp
     End Do

! calculate temporary COM velocity of each unit

     Call pmf_vcoms(indpmf,xpmf,ypmf,zpmf,comm)

! calculate PMF velocity corrections

     esig=0.0_wp
     Do ipmf=1,ntpmf

! calculate constraint force parameter - gamma

        gamma = pxx(ipmf)*(xpmf(1,ipmf)-xpmf(2,ipmf)) + &
                pyy(ipmf)*(ypmf(1,ipmf)-ypmf(2,ipmf)) + &
                pzz(ipmf)*(zpmf(1,ipmf)-zpmf(2,ipmf))

        esig=Max(esig,0.5_wp*Abs(gamma))

        gamma = gamma / (amt(1)+amt(2))

        Do jpmf=1,2

! If this unit is present on my domain

           If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
              gamm(jpmf) = Real(1-2*Mod(jpmf,2),wp)*gamma*amt(jpmf)
              Do k=1,mxtpmf(jpmf)
                 l=indpmf(k,jpmf,ipmf)

! improve approximate PMF particles velocity and force (if non-frozen)

                 If (l > 0 .and. l <= natms) Then ! l is a domain particle
                    If (lfrzn(l) == 0) Then
                       vxt(l)=vxt(l)+pxx(ipmf)*gamm(jpmf)
                       vyt(l)=vyt(l)+pyy(ipmf)*gamm(jpmf)
                       vzt(l)=vzt(l)+pzz(ipmf)*gamm(jpmf)
                    End If
                 End If
              End Do
           End If

        End Do

     End Do

! global verification of convergence

     safe=(esig < tolnce)
     Call gcheck(comm,safe,"enforce")

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then
        Do ipmf=1,ntpmf
           Do jpmf=1,2
              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
                 Do k=1,mxtpmf(jpmf)
                    l=indpmf(k,jpmf,ipmf)
                    If (l > 0 .and. l <= natms) Then ! l is a domain particle
                       If (lfrzn(l) == 0) Then
                          vxx(l)=vxx(l)+vxt(l)
                          vyy(l)=vyy(l)+vyt(l)
                          vzz(l)=vzz(l)+vzt(l)
                       End If
                    End If
                 End Do
              End If
           End Do
        End Do
     End If
  End Do

  If (.not.safe) Then ! error exit if quenching fails
     Call error(497)
  Else ! Collect per call passage statistics
     passpmq(1)=Real(icyc-1,wp)
     passpmq(3)=passpmq(2)*passpmq(3)
     passpmq(2)=passpmq(2)+1.0_wp
     passpmq(3)=passpmq(3)/passpmq(2)+passpmq(1)/passpmq(2)
     passpmq(4)=Min(passpmq(1),passpmq(4))
     passpmq(5)=Max(passpmq(1),passpmq(5))
     passpmq(1)=0.0_wp ! Reset
  End If

  Deallocate (lstitr,         Stat=fail(1))
  Deallocate (indpmf,         Stat=fail(2))
  Deallocate (pxx,pyy,pzz,    Stat=fail(3))
  Deallocate (vxt,vyt,vzt,    Stat=fail(4))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_quench deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_quench

Subroutine pmf_tags(lstitr,indpmf,pxx,pyy,pzz,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for identifying, constructing and indexing PMF
! constraints' vectors for iterative (constraints) algorithms
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( InOut ) :: lstitr(1:mxatms)
  Integer,           Intent(   Out ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent(   Out ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Type( comms_type ), Intent( InOut ) :: comm

  Integer :: ipmf,jpmf,j,k

! Loop over all local to this node PMF constraints, their units and
! save the indices of the members that are present on my domain (no halo)
! update lstitr

  Do ipmf=1,ntpmf

! Initialise indices

     indpmf(:,:,ipmf) = 0

     Do jpmf=1,2
        If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
           Do k=1,mxtpmf(jpmf)
              j=local_index(listpmf(k,jpmf,ipmf),nlast,lsi,lsa)
              indpmf(k,jpmf,ipmf)=j
              If (j > 0 .and. j <= natms) lstitr(j)=(lfrzn(j) == 0)
           End Do
        End If
     End Do

  End Do

! Get PMF units' COM vectors

  Call pmf_coms(indpmf,pxx,pyy,pzz,comm)

End Subroutine pmf_tags

Subroutine pmf_units_set(comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for setting the existence of the PMF units on the
! node where supposedly PMF constraints exist
!
! Note: (1) Deals with listpmf, legpmf and ntpmf
!       (2) Applies only at the end of relocate_particles
!           if megpmf>0 and mxnode>1
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type( comms_type ), Intent( InOut) :: comm
  Logical :: safe,ok
  Integer :: fail,ipmf,jpmf,gpmf

  Integer, Dimension( : ), Allocatable :: i1pmf0,i2pmf0
  Character( Len = 256 ) :: message

  fail=0
  Allocate (i1pmf0(1:mxtpmf(1)),i2pmf0(1:mxtpmf(2)), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_units_set allocation failure'
     Call error(0,message)
  End If

! Initialise safety flag

  safe=.true.

! is it ok not to compress the bookkeeping list arrays
! since it's safe - there's enough buffering space

  ok=.true.
  If (mxpmf > 0) ok=.not.(Real(ntpmf,wp)/Real(mxpmf,wp) > 0.85_wp)

! Sort out local PMF units presence and legend array

  ipmf=0
  Do While (ipmf < ntpmf)
     ipmf=ipmf+1
10   Continue

! This holds the global PMF index

     gpmf=listpmf(0,1,ipmf)
     If (gpmf > 0) Then

! For presence of : PMF unit 1 only - this holds 1
!                   PMF unit 2 only - this holds 2
!                   both units 1&2  - this holds 3
! It CANNOT and MUST NOT hold ZERO

        listpmf(0,2,ipmf)=0

        Do jpmf=1,mxtpmf(1)
           i1pmf0(jpmf)=local_index(listpmf(jpmf,1,ipmf),natms,lsi,lsa)

           If (i1pmf0(jpmf) > 0) Then

! This identifies which local PMF constraint the particle belongs to

              legpmf(0,i1pmf0(jpmf))=1
              legpmf(1,i1pmf0(jpmf))=ipmf

           End If
        End Do
        If (Any(i1pmf0 > 0)) listpmf(0,2,ipmf)=listpmf(0,2,ipmf)+1

        Do jpmf=1,mxtpmf(2)
           i2pmf0(jpmf)=local_index(listpmf(jpmf,2,ipmf),natms,lsi,lsa)

           If (i2pmf0(jpmf) > 0) Then

! This identifies which local PMF constraint the particle belongs to

              legpmf(0,i2pmf0(jpmf))=1
              legpmf(1,i2pmf0(jpmf))=ipmf

           End If
        End Do
        If (Any(i2pmf0 > 0)) listpmf(0,2,ipmf)=listpmf(0,2,ipmf)+2

! If the PMF has moved to another node, listpmf(0,2,ipmf)=0,
! compress listpmf

        If (listpmf(0,2,ipmf) == 0 .and. (.not.ok)) Then
           If      (ipmf  < ntpmf) Then
              listpmf(:,:,ipmf)=listpmf(:,:,ntpmf) ! Copy list content from 'ipmf' to 'ntpmf'
              listpmf(:,:,ntpmf)=0                 ! Remove list content in 'ntpmf'
              ntpmf=ntpmf-1                        ! Reduce 'ntpmf' pointer

              Go To 10 ! Go back and check it all again for the new list content in 'ipmf'
           Else If (ipmf == ntpmf) Then
              listpmf(:,:,ntpmf)=0                 ! Remove list content in 'ntpmf=ipmf'
              ntpmf=ntpmf-1                        ! Reduce 'ntpmf' pointer
           End If
        End If

     Else

        safe=.false.

     End If
  End Do

! check for inconsistently built local list

  Call gcheck(comm,safe)
  If (.not.safe) Call error(490)

  Deallocate (i1pmf0,i2pmf0, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_units_set deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_units_set
Subroutine pmf_vcoms(indpmf,xpmf,ypmf,zpmf,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for constructing PMF units' c.o.m. velocity
! vectors for iterative constraint (pmf_quench,pmf_rattle) algorithms
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent(   Out ) :: xpmf(1:2,1:mxpmf),ypmf(1:2,1:mxpmf),zpmf(1:2,1:mxpmf)
  Type( comms_type), Intent( InOut ) :: comm

  Logical                 :: safe

  Integer                 :: fail,gpmf,gpmf1,gpmf2,ipmf,jpmf,iadd, &
                             j,k,l

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message

  fail=0
  Allocate (buffer(1:(mxtpmf(1)+mxtpmf(2))*(mxpmf+2)), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_vcoms allocation failure'
     Call error(0,message)
  End If


! Initialise safety flag

  safe=.true.

! Initialise PMF c.o.m. momentum arrays

  xpmf = 0.0_wp ; ypmf = 0.0_wp ; zpmf = 0.0_wp

! Loop over all global PMF constraints

  gpmf1=1          ! number of passed global PMFs
  buffer=0.0_wp    ! velocities buffer
  iadd=3           ! adding three velocities only per particle
  Do gpmf=1,mxpmf

! Loop over all local to this node PMF constraints matching the global one

     Do ipmf=1,ntpmf
        If (listpmf(0,1,ipmf) == gpmf) Then

! Loop over all PMF units present on my domain (no halo)

           Do jpmf=1,2
              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then

! Loop over their members present on my domain (no halo)

                 Do k=1,mxtpmf(jpmf)
                    j=indpmf(k,jpmf,ipmf)

! Copy particles' velocities to buffer in an orderly manner

                    If (j > 0 .and. j <= natms) Then ! j is a domain particle
                       l=((gpmf-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1)+(k-1))*iadd
                       buffer(l+1)=vxx(j)
                       buffer(l+2)=vyy(j)
                       buffer(l+3)=vzz(j)
                    End If
                 End Do
              End If
           End Do
        End If
     End Do

! Check if it safe to fill up the buffer

     safe=(gpmf-gpmf1+2 < (mxpmf+2)/3)

! If not safe or we've finished looping over all global PMFs

     If ((.not.safe) .or. gpmf == mxpmf) Then
        Call gsum(comm,buffer)

        Do gpmf2=gpmf1,gpmf
           Do ipmf=1,ntpmf
              If (listpmf(0,1,ipmf) == gpmf2) Then
                 Do jpmf=1,2

! Get the COM momentum of this unit

                    Do k=1,mxtpmf(jpmf)
                       l=((gpmf2-gpmf1)*(mxtpmf(1)+mxtpmf(2))+(jpmf-1)*mxtpmf(1)+(k-1))*iadd

                       xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf) + pmfwg1(k,jpmf)*buffer(l+1)
                       ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf) + pmfwg1(k,jpmf)*buffer(l+2)
                       zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf) + pmfwg1(k,jpmf)*buffer(l+3)
                    End Do

                    xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf)*pmfwg1(0,jpmf)
                    ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf)*pmfwg1(0,jpmf)
                    zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf)*pmfwg1(0,jpmf)

                 End Do
              End If
           End Do
        End Do

        gpmf1=gpmf2
        buffer=0.0_wp
     End If

  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_vcoms deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_vcoms

Subroutine pmf_shake_vv          &
           (mxshak,tolnce,tstep, &
           indpmf,pxx,pyy,pzz,   &
           xxx,yyy,zzz,strpmf,   &
           virpmf,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for applying PMF constraint corrections after
! possible constrained motion
!
! Note: must be used in conjunction with integration algorithms
!       VV compliant
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce,tstep
  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent( In    ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: xxx(1:mxatms),yyy(1:mxatms),zzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: strpmf(1:9),virpmf
  Type( comms_type ), Intent( InOut ) :: comm

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: rmass_pmf_unit(1:2),dis,dis2

  Logical                 :: safe
  Integer                 :: fail,ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: amt(1:2),gamma,gamm(1:2),tstep2,tmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: pxt,pyt,pzt,pt2,esig
  Character( Len = 256 ) :: message

  fail=0
  Allocate (pxt(1:mxpmf),pyt(1:mxpmf),pzt(1:mxpmf),pt2(1:mxpmf),esig(1:mxpmf), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_shake allocation failure'
     Call error(0,message)
  End If


  If (newjob) Then
     newjob = .false.

! Get reciprocal PMF units' masses

     Do jpmf=1,2
        If (pmffrz(jpmf) == mxtpmf(jpmf)) Then
           rmass_pmf_unit(jpmf)=0.0_wp
        Else
           rmass_pmf_unit(jpmf)=pmfwg1(0,jpmf)
        End If
     End Do

! set PMF constraint parameters

     dis=prmpmf
     dis2=dis**2
  End If

! squared timestep and reciprocal masses

  tstep2 = tstep*tstep
  amt = tstep2*rmass_pmf_unit

! Initialise constraint virial and stress

  virpmf=0.0_wp
  strpmf=0.0_wp

! application of PMF constraint (shake) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  safe=.false.
  icyc=0
  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! calculate temporary PMF units' COM vectors

     Call pmf_coms(indpmf,pxt,pyt,pzt,comm)

! calculate maximum error in bondlength

     Do ipmf=1,ntpmf
        pt2(ipmf) =pxt(ipmf)**2+pyt(ipmf)**2+pzt(ipmf)**2 - dis2
        esig(ipmf)=0.5_wp*Abs(pt2(ipmf))/dis
     End Do

! global verification of convergence

     safe=Merge(Maxval(esig(1:ntpmf)) < tolnce,.true.,ntpmf > 0)
     Call gcheck(comm,safe,"enforce")

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then

! calculate PMF constraint forces

        Do ipmf=1,ntpmf

! calculate PMF constraint force parameter

           gamma = -pt2(ipmf) / &
                   ((amt(1)+amt(2))*(pxx(ipmf)*pxt(ipmf)+pyy(ipmf)*pyt(ipmf)+pzz(ipmf)*pzt(ipmf)))
           tmp   = gamma / Real(mxtpmf(1)+mxtpmf(2),wp)

           Do jpmf=1,2

! If this unit is present on my domain

              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
                 gamm(jpmf) = 0.5_wp*Real(1-2*Mod(jpmf,2),wp)*gamma*amt(jpmf)

                 Do k=1,mxtpmf(jpmf)
                    l=indpmf(k,jpmf,ipmf)

! for non-frozen domain particles accumulate PMF constraint stress
! and atomic position corrections

                    If (l > 0 .and. l <= natms) Then ! l is a domain particle
                       If (lfrzn(l) == 0) Then
                          strpmf(1) = strpmf(1) - tmp*pxx(ipmf)*pxx(ipmf)
                          strpmf(2) = strpmf(2) - tmp*pxx(ipmf)*pyy(ipmf)
                          strpmf(3) = strpmf(3) - tmp*pxx(ipmf)*pzz(ipmf)
                          strpmf(5) = strpmf(5) - tmp*pyy(ipmf)*pyy(ipmf)
                          strpmf(6) = strpmf(6) - tmp*pyy(ipmf)*pzz(ipmf)
                          strpmf(9) = strpmf(9) - tmp*pzz(ipmf)*pzz(ipmf)

                          xxx(l)=xxx(l)+pxx(ipmf)*gamm(jpmf)
                          yyy(l)=yyy(l)+pyy(ipmf)*gamm(jpmf)
                          zzz(l)=zzz(l)+pzz(ipmf)*gamm(jpmf)
                       End If
                    End If
                 End Do
              End If

           End Do

        End Do

     End If
  End Do

  If (.not.safe) Then ! error exit for non-convergence
     Do k=0,comm%mxnode-1
        If (comm%idnode == k) Then
           Do ipmf=1,ntpmf
             If (esig(ipmf) >= tolnce) Then
                 Write(message,'(3(a,i10),a)') &
                 'global PMF constraint number', ipmf, &
                 ' , with head particle numbers, U1:', listpmf(1,1,ipmf), &
                 ' & U2:', listpmf(1,2,ipmf)
               Call info(message)
               Write(message,'(a,f8.2,a,1p,e12.4)') &
                 'converges to a length of', Sqrt(pt2(ipmf)+dis2), &
                 ' Angstroms with factor', esig(ipmf)
               Call info(message)
               Call warning('Contributes towards next error',.true.)
             End If
           End Do
        End If
        Call gsync(comm)
     End Do
     Call error(498)
  Else ! Collect per call and per step passage statistics
     passpmf(1,1,1)=Real(icyc-1,wp)
     passpmf(3,1,1)=passpmf(2,1,1)*passpmf(3,1,1)
     passpmf(2,1,1)=passpmf(2,1,1)+1.0_wp
     passpmf(3,1,1)=passpmf(3,1,1)/passpmf(2,1,1)+passpmf(1,1,1)/passpmf(2,1,1)
     passpmf(4,1,1)=Min(passpmf(1,1,1),passpmf(4,1,1))
     passpmf(5,1,1)=Max(passpmf(1,1,1),passpmf(5,1,1))

     passpmf(1,2,1)=passpmf(1,2,1)+passpmf(1,1,1)
     passpmf(1,1,1)=0.0_wp ! Reset
  End If

! global sum of stress tensor

  Call gsum(comm,strpmf)

! complete stress tensor (symmetrise)

  strpmf(4) = strpmf(2)
  strpmf(7) = strpmf(3)
  strpmf(8) = strpmf(6)

! total PMF constraint virial

  virpmf=-(strpmf(1)+strpmf(5)+strpmf(9))

  Deallocate (pxt,pyt,pzt,pt2,esig, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_shake deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_shake_vv

Subroutine pmf_rattle                      &
           (mxshak,tolnce,tstep,lfst,lcol, &
           indpmf,pxx,pyy,pzz,             &
           vxx,vyy,vzz,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for applying PMF constraint corrections after
! possible constrained motion
!
! Note: must be used in conjunction with integration algorithms
!       VV compliant
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: mxshak
  Real( Kind = wp ), Intent( In    ) :: tolnce,tstep
  Logical,           Intent( In    ) :: lfst,lcol
  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: amt(1:2)

  Logical                 :: safe
  Integer                 :: fail(1:2),ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: dis,esig,gamma,gamm(1:2)

  Real( Kind = wp ), Dimension( : )   , Allocatable :: vxt,vyt,vzt
  Real( Kind = wp ), Dimension( :, : ), Allocatable :: xpmf,ypmf,zpmf
  Character( Len = 256 ) :: message

  fail=0
  Allocate (vxt(1:mxatms),vyt(1:mxatms),vzt(1:mxatms),             Stat=fail(1))
  Allocate (xpmf(1:2,1:mxpmf),ypmf(1:2,1:mxpmf),zpmf(1:2,1:mxpmf), Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_rattle allocation failure'
     Call error(0,message)
  End If

! Get PMF units' reciprocal masses

  If (newjob) Then
     newjob = .false.

     Do jpmf=1,2
        If (pmffrz(jpmf) == mxtpmf(jpmf)) Then
           amt(jpmf)=0.0_wp
        Else
           amt(jpmf)=pmfwg1(0,jpmf)
        End If
     End Do
  End If

! normalise PMF constraint vectors on first pass outside

  If (lfst) Then
     Do ipmf=1,ntpmf
        dis=1.0_wp/Sqrt(pxx(ipmf)**2+pyy(ipmf)**2+pzz(ipmf)**2)
        pxx(ipmf)=pxx(ipmf)*dis
        pyy(ipmf)=pyy(ipmf)*dis
        pzz(ipmf)=pzz(ipmf)*dis
     End Do
  End If

! application of PMF constraint (rattle) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  safe=.false.
  icyc=0
  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! initialise velocity correction arrays

     Do l=1,natms
        vxt(l)=0.0_wp
        vyt(l)=0.0_wp
        vzt(l)=0.0_wp
     End Do

! calculate temporary COM velocity of each unit

     Call pmf_vcoms(indpmf,xpmf,ypmf,zpmf,comm)

! calculate PMF velocity corrections

     esig=0.0_wp
     Do ipmf=1,ntpmf

! calculate constraint force parameter - gamma

        gamma = pxx(ipmf)*(xpmf(1,ipmf)-xpmf(2,ipmf)) + &
                pyy(ipmf)*(ypmf(1,ipmf)-ypmf(2,ipmf)) + &
                pzz(ipmf)*(zpmf(1,ipmf)-zpmf(2,ipmf))

        esig=Max(esig,0.5_wp*tstep*Abs(gamma))

        gamma = gamma / (amt(1)+amt(2))

        Do jpmf=1,2

! If this unit is present on my domain

           If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
              gamm(jpmf) = Real(1-2*Mod(jpmf,2),wp)*gamma*amt(jpmf)
              Do k=1,mxtpmf(jpmf)
                 l=indpmf(k,jpmf,ipmf)

! improve approximate PMF particles velocity and force (if non-frozen)

                 If (l > 0 .and. l <= natms) Then ! l is a domain particle
                    If (lfrzn(l) == 0) Then
                       vxt(l)=vxt(l)+pxx(ipmf)*gamm(jpmf)
                       vyt(l)=vyt(l)+pyy(ipmf)*gamm(jpmf)
                       vzt(l)=vzt(l)+pzz(ipmf)*gamm(jpmf)
                    End If
                 End If
              End Do
           End If

        End Do

     End Do

! global verification of convergence

     safe=(esig < tolnce)
     Call gcheck(comm,safe,"enforce")

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then
        Do ipmf=1,ntpmf
           Do jpmf=1,2
              If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
                 Do k=1,mxtpmf(jpmf)
                    l=indpmf(k,jpmf,ipmf)
                    If (l > 0 .and. l <= natms) Then ! l is a domain particle
                       If (lfrzn(l) == 0) Then
                          vxx(l)=vxx(l)+vxt(l)
                          vyy(l)=vyy(l)+vyt(l)
                          vzz(l)=vzz(l)+vzt(l)
                       End If
                    End If
                 End Do
              End If
           End Do
        End Do
     End If
  End Do

  If (.not.safe) Then ! error exit for non-convergence
     Call error(499)
  Else ! Collect per call and per step passage statistics
     passpmf(1,1,2)=Real(icyc-1,wp)
     passpmf(3,1,2)=passpmf(2,1,2)*passpmf(3,1,2)
     passpmf(2,1,2)=passpmf(2,1,2)+1.0_wp
     passpmf(3,1,2)=passpmf(3,1,2)/passpmf(2,1,2)+passpmf(1,1,2)/passpmf(2,1,2)
     passpmf(4,1,2)=Min(passpmf(1,1,2),passpmf(4,1,2))
     passpmf(5,1,2)=Max(passpmf(1,1,2),passpmf(5,1,2))

     passpmf(1,2,2)=passpmf(1,2,2)+passpmf(1,1,2)
     If (lcol) Then ! Collect
        passpmf(3,2,2)=passpmf(2,2,2)*passpmf(3,2,2)
        passpmf(2,2,2)=passpmf(2,2,2)+1.0_wp
        passpmf(3,2,2)=passpmf(3,2,2)/passpmf(2,2,2)+passpmf(1,2,2)/passpmf(2,2,2)
        passpmf(4,2,2)=Min(passpmf(1,2,2),passpmf(4,2,2))
        passpmf(5,2,2)=Max(passpmf(1,2,2),passpmf(5,2,2))
        passpmf(1,2,2)=0.0_wp ! Reset
     End If
     passpmf(1,1,2)=0.0_wp ! Reset
  End If

  Deallocate (vxt,vyt,vzt,    Stat=fail(1))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_rattle deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_rattle


End module pmf
