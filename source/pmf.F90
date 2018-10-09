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
  Use configuration,   Only : configuration_type
  Use errors_warnings, Only : error,warning,info
  Use numerics,        Only : images,local_index,dcell
  Use statistics, Only : stats_type
  Use particle,        Only : corePart
  Implicit None
  Private
  
  Type, Public ::  pmf_type
  Private

  Integer,                        Public :: ntpmf  = 0
  Integer, Public :: mxtpmf(1:2),mxpmf,mxfpmf,megpmf
  Logical            :: newjob_quench = .true.
  Real( Kind = wp )  :: amt(1:2)
  Logical           :: newjob_shake = .true.
  Real( Kind = wp ) :: rmass_pmf_unit(1:2),dis,dis2
  Logical            :: newjob_rattle = .true.
  Real( Kind = wp )  :: amtr(1:2)

  Real( Kind = wp ),              Public :: prmpmf = 0.0_wp
  

  Integer,           Allocatable, Public :: numpmf(:),pmffrz(:)
  Integer,           Allocatable, Public :: lstpmf(:,:),listpmf(:,:,:),legpmf(:,:)
  Integer,           Allocatable :: indpmf(:,:,:)
  
  Real( Kind = wp ), Allocatable :: pxx(:),pyy(:),pzz(:)
  Real( Kind = wp ), Allocatable, Public :: pmfwgt(:,:),pmfwg1(:,:)

  Contains
  Private
  Procedure, Public :: init => allocate_pmf_arrays
  Procedure, Public :: deallocate_pmf_tmp_arrays
  Procedure, Public :: allocate_work
  Procedure, Public :: deallocate_work
  Procedure, Public :: deallocate_pmf_arrays
  End Type pmf_type
  Public :: pmf_shake_vv
  Public :: pmf_rattle
  Public :: pmf_tags
  Public :: pmf_pseudo_bonds
  Public :: pmf_units_set
  Public :: pmf_quench

Contains
Subroutine allocate_work(T)
  Class(pmf_type) :: T
  Integer :: fail(2)
  Character(Len=100) :: message
  
  fail=0
        If (T%megpmf > 0) Then
        Allocate (T%indpmf(1:Max(T%mxtpmf(1),T%mxtpmf(2)),1:2,1:T%mxpmf), Stat=fail( 1))
        Allocate (T%pxx(1:T%mxpmf),T%pyy(1:T%mxpmf),T%pzz(1:T%mxpmf),         Stat=fail( 2))
      End If
         If (Any(fail > 0)) Then
      Write(message,'(a)') 'failed to allocate work arrays for pmf'
      Call error(0,message)
    End If
  End Subroutine allocate_work
  
  Subroutine deallocate_work(T)
  Class(pmf_type) :: T
  Integer :: fail(4)
  Character(Len=100) :: message
  
  fail=0
  If (T%megpmf >0) Then
  
  If (Allocated(T%indpmf)) Deallocate(T%indpmf, stat=fail(1))
  If (Allocated(T%pxx)) Deallocate(T%pxx, stat=fail(2))
  If (Allocated(T%pyy)) Deallocate(T%pyy, stat=fail(3))
  If (Allocated(T%pzz)) Deallocate(T%pzz, stat=fail(4))
  
  If (Any(fail > 0)) Then
      Write(message,'(a)') 'failed to allocate work arrays for pmf'
      Call error(0,message)
  End If 
  End If
  
  End Subroutine deallocate_work
  
  Subroutine allocate_pmf_arrays(T,mxtmls,mxatdm)
  Class(pmf_type) :: T
  Integer, Intent ( In ) :: mxtmls, mxatdm
    Integer, Dimension( 1:6 ) :: fail
    

    fail = 0

    Allocate (T%numpmf(1:mxtmls),T%pmffrz(1:2),                    Stat = fail(1))
    Allocate (T%lstpmf(1:Max(T%mxtpmf(1),T%mxtpmf(2)),1:2),          Stat = fail(2))
    Allocate (T%listpmf(0:Max(T%mxtpmf(1),T%mxtpmf(2)),1:2,1:T%mxpmf), Stat = fail(3))
    Allocate (T%legpmf(0:T%mxfpmf,1:mxatdm),                       Stat = fail(4))
    Allocate (T%pmfwgt(0:Max(T%mxtpmf(1),T%mxtpmf(2)),1:2),          Stat = fail(5))
    Allocate (T%pmfwg1(0:Max(T%mxtpmf(1),T%mxtpmf(2)),1:2),          Stat = fail(6))

    If (Any(fail > 0)) Call error(1036)

    T%numpmf  = 0
    T%pmffrz  = 0
    T%lstpmf  = 0
    T%listpmf = 0
    T%legpmf  = 0

    T%pmfwgt = 0.0_wp
    T%pmfwg1 = 0.0_wp

  End Subroutine allocate_pmf_arrays

  Subroutine deallocate_pmf_tmp_arrays(T)
    Class(pmf_type) :: T
    Integer :: fail(2)

    fail = 0

    If (Allocated(T%numpmf)) Deallocate (T%numpmf, Stat = fail(1))
    If (Allocated(T%lstpmf)) Deallocate (T%lstpmf, Stat = fail(2))

    If (Any(fail > 0)) Call error(1037)

  End Subroutine deallocate_pmf_tmp_arrays
  
  Subroutine deallocate_pmf_arrays(T)
    Class(pmf_type) :: T

  Integer :: fail(7)

  fail = 0

    If (Allocated(T%numpmf)) Deallocate (T%numpmf, Stat = fail(1))
    If (Allocated(T%lstpmf)) Deallocate (T%lstpmf, Stat = fail(2))
    If (Allocated(T%pmffrz)) Deallocate (T%pmffrz, Stat = fail(3))
    If (Allocated(T%listpmf)) Deallocate (T%listpmf, Stat = fail(4))
    If (Allocated(T%legpmf)) Deallocate (T%legpmf, Stat = fail(5))
    If (Allocated(T%pmfwgt)) Deallocate (T%pmfwgt, Stat = fail(6))
    If (Allocated(T%pmfwg1)) Deallocate (T%pmfwg1, Stat = fail(7))

  If (Any(fail > 0)) Call error(1037)

End Subroutine deallocate_pmf_arrays
  
  
Subroutine pmf_coms(pmf,pxx,pyy,pzz,config,comm)

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


  Type( pmf_type), Intent( InOut ) :: pmf
  Real( Kind = wp ), Intent(   Out ) :: pxx(1:),pyy(1:),pzz(1:)
  Type( configuration_type ),   Intent( InOut ) :: config
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
  Allocate (xxt(1:Max(pmf%mxtpmf(1),pmf%mxtpmf(2))),yyt(1:Max(pmf%mxtpmf(1),pmf%mxtpmf(2))),&
  zzt(1:Max(pmf%mxtpmf(1),pmf%mxtpmf(2))), Stat=fail(1))
  Allocate (xpmf(1:2,1:pmf%mxpmf),ypmf(1:2,1:pmf%mxpmf),zpmf(1:2,1:pmf%mxpmf),  Stat=fail(2))
  Allocate (buffer(1:(pmf%mxtpmf(1)+pmf%mxtpmf(2))*(pmf%mxpmf+2)),              Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_coms allocation failure'
     Call error(0,message)
  End If


! Initialise safety flags

  safe=.true.

! Get the dimensional properties of the MD cell

  Call dcell(config%cell,celprp)
  width=Min(celprp(7),celprp(8),celprp(9))

! Initialise PMF COMs' and inter-COMs' vector arrays

  xpmf = 0.0_wp ; ypmf = 0.0_wp ; zpmf = 0.0_wp
  pxx  = 0.0_wp ; pyy  = 0.0_wp ; pzz  = 0.0_wp

! Loop over all global PMF constraints

  gpmf1=1          ! number of passed global PMFs
  buffer=0.0_wp    ! coordinates buffer
  iadd=3           ! adding three coordinates only per particle
  Do gpmf=1,pmf%mxpmf

! Loop over all local to this node PMF constraints matching the global one

     Do ipmf=1,pmf%ntpmf
        If (pmf%listpmf(0,1,ipmf) == gpmf) Then

! Loop over all PMF units present on my domain (no halo)

           Do jpmf=1,2
              If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then

! Loop over their members present on my domain (no halo)

                 Do k=1,pmf%mxtpmf(jpmf)
                    j=pmf%indpmf(k,jpmf,ipmf)

! Copy particles' coordinates to buffer in an orderly manner

                    If (j > 0 .and. j <= config%natms) Then ! j is a domain particle
                       l=((gpmf-gpmf1)*(pmf%mxtpmf(1)+pmf%mxtpmf(2))+(jpmf-1)*pmf%mxtpmf(1)+(k-1))*iadd
                       buffer(l+1)=config%parts(j)%xxx
                       buffer(l+2)=config%parts(j)%yyy
                       buffer(l+3)=config%parts(j)%zzz
                    End If
                 End Do
              End If
           End Do
        End If
     End Do

! Check if it safe to fill up the buffer

     safe(1)=(gpmf-gpmf1+2 < (pmf%mxpmf+2)/3)

! If not safe or we've finished looping over all global PMFs

     If ((.not.safe(1)) .or. gpmf == pmf%mxpmf) Then
        Call gsum(comm,buffer)

        Do gpmf2=gpmf1,gpmf
           Do ipmf=1,pmf%ntpmf
              If (pmf%listpmf(0,1,ipmf) == gpmf2) Then
                 Do jpmf=1,2

! Get in the local scope of the unit

                    m=((gpmf2-gpmf1)*(pmf%mxtpmf(1)+pmf%mxtpmf(2))+(jpmf-1)*pmf%mxtpmf(1))*iadd
                    Do k=1,pmf%mxtpmf(jpmf)
                       l=m+(k-1)*iadd

                       xxt(k)=buffer(l+1)-buffer(m+1)
                       yyt(k)=buffer(l+2)-buffer(m+2)
                       zzt(k)=buffer(l+3)-buffer(m+3)
                    End Do

                    Call images(config%imcon,config%cell,pmf%mxtpmf(jpmf),xxt,yyt,zzt)

                    xmin=0.0_wp ; xmax = 0.0_wp
                    ymin=0.0_wp ; ymax = 0.0_wp
                    zmin=0.0_wp ; zmax = 0.0_wp
                    Do k=1,pmf%mxtpmf(jpmf)
                       xmin=Min(xmin,xxt(k)) ; xmax=Max(xmax,xxt(k))
                       ymin=Min(ymin,yyt(k)) ; ymax=Max(ymax,yyt(k))
                       zmin=Min(zmin,zzt(k)) ; zmax=Max(zmax,zzt(k))
                       safe(2)=safe(2) .and. (Sqrt(xxt(k)**2+yyt(k)**2+zzt(k)**2) < width/2.0_wp)
                    End Do
                    safe(2)=safe(2) .and. (xmax-xmin < width/2.0_wp) &
                                    .and. (ymax-ymin < width/2.0_wp) &
                                    .and. (zmax-zmin < width/2.0_wp)

! Get the COM of this unit

                    Do k=1,pmf%mxtpmf(jpmf)
                       xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf) + pmf%pmfwgt(k,jpmf)*xxt(k)
                       ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf) + pmf%pmfwgt(k,jpmf)*yyt(k)
                       zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf) + pmf%pmfwgt(k,jpmf)*zzt(k)
                    End Do

! Get out of local frame

                    xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf)*pmf%pmfwgt(0,jpmf) + buffer(m+1)
                    ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf)*pmf%pmfwgt(0,jpmf) + buffer(m+2)
                    zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf)*pmf%pmfwgt(0,jpmf) + buffer(m+3)

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

  Do ipmf=1,pmf%ntpmf
     pxx(ipmf) = xpmf(2,ipmf) - xpmf(1,ipmf)
     pyy(ipmf) = ypmf(2,ipmf) - ypmf(1,ipmf)
     pzz(ipmf) = zpmf(2,ipmf) - zpmf(1,ipmf)
  End Do

! Minimum image convention for bond vectors

  Call images(config%imcon,config%cell,pmf%ntpmf,pxx,pyy,pzz)

  Deallocate (xxt,yyt,zzt,    Stat=fail(1))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(2))
  Deallocate (buffer,         Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_coms deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_coms

Subroutine pmf_pseudo_bonds(gxx,gyy,gzz,stat,pmf,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for treating PMF constraints as stiff harmonic
! springs for use with the conjugate gradient method (minimise_relax.f90)
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov august 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Real( Kind = wp ), Intent( InOut ) :: gxx(1:),gyy(1:),gzz(1:)
  Type( stats_type ), Intent( InOut ) :: stat
  Type( pmf_type), Intent( InOut ) :: pmf
  Type( comms_type), Intent( InOut ) :: comm
  Type( configuration_type ), Intent( InOut ) :: config

  Real( Kind = wp ), Parameter :: rigid=1.0e6_wp

  Integer           :: ipmf,jpmf,k,l
  Real( Kind = wp ) :: r,r0,ebond,gamma,tmp

  r0=pmf%prmpmf

  stat%engpmf=0.0_wp
  Do ipmf=1,pmf%ntpmf

     r=Sqrt(pmf%pxx(ipmf)**2+pmf%pyy(ipmf)**2+pmf%pzz(ipmf)**2)

     gamma=rigid*(r-r0)/Real(pmf%mxtpmf(1)+pmf%mxtpmf(2),wp)
     ebond=gamma*0.5_wp*(r-r0)
     gamma=gamma/r

     Do jpmf=1,2
        tmp=Real(1-2*Mod(jpmf,2),wp)*gamma

! If this unit is present on my domain

        If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then

           Do k=1,pmf%mxtpmf(jpmf)
              l=pmf%indpmf(k,jpmf,ipmf)

! For domain particles

              If (l > 0 .and. l <= config%natms) Then

! Accumulate energy

                 stat%engpmf=stat%engpmf+ebond

! Add forces

                 If (config%lfrzn(l) == 0) Then
                    gxx(l)=gxx(l)+pmf%pxx(ipmf)*tmp
                    gyy(l)=gyy(l)+pmf%pyy(ipmf)*tmp
                    gzz(l)=gzz(l)+pmf%pzz(ipmf)*tmp
                 End If

              End If

           End Do

        End If

     End Do

  End Do

! Global sum of energy

  Call gsum(comm,stat%engpmf)

End Subroutine pmf_pseudo_bonds

Subroutine pmf_quench(mxshak,tolnce,stat,pmf,config,comm)

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
  Type( stats_type ), Intent( InOut ) :: stat
  Type( pmf_type), Intent( InOut ) :: pmf
  Type( configuration_type ),  Intent( InOut ) :: config
  Type( comms_type), intent( InOut ) :: comm


  Logical                 :: safe
  Integer                 :: fail(1:5),ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: dis,esig,gamma,gamm(1:2)

  Logical,           Allocatable :: lstitr(:)
  Real( Kind = wp ), Allocatable :: vxt(:),vyt(:),vzt(:)
  Real( Kind = wp ), Allocatable :: xpmf(:,:),ypmf(:,:),zpmf(:,:)
  Character( Len = 256 ) :: message

  fail=0
  Allocate (lstitr(1:config%mxatms),                                      Stat=fail(1))
  Call pmf%allocate_work()
  Allocate (vxt(1:config%mxatms),vyt(1:config%mxatms),vzt(1:config%mxatms),             Stat=fail(4))
  Allocate (xpmf(1:2,1:pmf%mxpmf),ypmf(1:2,1:pmf%mxpmf),zpmf(1:2,1:pmf%mxpmf), Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_quench allocation failure'
     Call error(0,message)
  End If

! Get PMF units' reciprocal masses

  If (pmf%newjob_quench) Then
    pmf%newjob_quench = .false.

     Do jpmf=1,2
        If (pmf%pmffrz(jpmf) == pmf%mxtpmf(jpmf)) Then
          pmf%amt(jpmf)=0.0_wp
        Else
          pmf%amt(jpmf)=pmf%pmfwg1(0,jpmf)
        End If
     End Do
  End If

  lstitr(1:config%natms)=.false. ! initialise lstitr
  Call pmf_tags(lstitr,pmf,config,comm)

! normalise PMF constraint vectors

  Do ipmf=1,pmf%ntpmf
     dis=1.0_wp/Sqrt(pmf%pxx(ipmf)**2+pmf%pyy(ipmf)**2+pmf%pzz(ipmf)**2)
     pmf%pxx(ipmf)=pmf%pxx(ipmf)*dis
     pmf%pyy(ipmf)=pmf%pyy(ipmf)*dis
     pmf%pzz(ipmf)=pmf%pzz(ipmf)*dis
  End Do

! application of PMF constraint (rattle) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

  icyc=0
  safe=.false.

  Do While ((.not.safe) .and. icyc < mxshak)
     icyc=icyc+1

! initialise velocity correction arrays

     Do l=1,config%natms
        vxt(l)=0.0_wp
        vyt(l)=0.0_wp
        vzt(l)=0.0_wp
     End Do

! calculate temporary COM velocity of each unit

     Call pmf_vcoms(xpmf,ypmf,zpmf,pmf,config,comm)

! calculate PMF velocity corrections

     esig=0.0_wp
     Do ipmf=1,pmf%ntpmf

! calculate constraint force parameter - gamma

        gamma = pmf%pxx(ipmf)*(xpmf(1,ipmf)-xpmf(2,ipmf)) + &
                pmf%pyy(ipmf)*(ypmf(1,ipmf)-ypmf(2,ipmf)) + &
                pmf%pzz(ipmf)*(zpmf(1,ipmf)-zpmf(2,ipmf))

        esig=Max(esig,0.5_wp*Abs(gamma))

        gamma = gamma / (pmf%amt(1)+pmf%amt(2))

        Do jpmf=1,2

! If this unit is present on my domain

           If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then
             gamm(jpmf) = Real(1-2*Mod(jpmf,2),wp)*gamma*pmf%amt(jpmf)
              Do k=1,pmf%mxtpmf(jpmf)
                 l=pmf%indpmf(k,jpmf,ipmf)

! improve approximate PMF particles velocity and force (if non-frozen)

                 If (l > 0 .and. l <= config%natms) Then ! l is a domain particle
                    If (config%lfrzn(l) == 0) Then
                       vxt(l)=vxt(l)+pmf%pxx(ipmf)*gamm(jpmf)
                       vyt(l)=vyt(l)+pmf%pyy(ipmf)*gamm(jpmf)
                       vzt(l)=vzt(l)+pmf%pzz(ipmf)*gamm(jpmf)
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
        Do ipmf=1,pmf%ntpmf
           Do jpmf=1,2
              If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then
                 Do k=1,pmf%mxtpmf(jpmf)
                    l=pmf%indpmf(k,jpmf,ipmf)
                    If (l > 0 .and. l <= config%natms) Then ! l is a domain particle
                       If (config%lfrzn(l) == 0) Then
                          config%vxx(l)=config%vxx(l)+vxt(l)
                          config%vyy(l)=config%vyy(l)+vyt(l)
                          config%vzz(l)=config%vzz(l)+vzt(l)
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
     stat%passpmq(1)=Real(icyc-1,wp)
     stat%passpmq(3)=stat%passpmq(2)*stat%passpmq(3)
     stat%passpmq(2)=stat%passpmq(2)+1.0_wp
     stat%passpmq(3)=stat%passpmq(3)/stat%passpmq(2)+stat%passpmq(1)/stat%passpmq(2)
     stat%passpmq(4)=Min(stat%passpmq(1),stat%passpmq(4))
     stat%passpmq(5)=Max(stat%passpmq(1),stat%passpmq(5))
     stat%passpmq(1)=0.0_wp ! Reset
  End If

  Deallocate (lstitr,         Stat=fail(1))
  Call pmf%deallocate_work()
  Deallocate (vxt,vyt,vzt,    Stat=fail(4))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_quench deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_quench

Subroutine pmf_tags(lstitr,pmf,config,comm)

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

  Logical,           Intent( InOut ) :: lstitr(1:)
  Type(pmf_type), Intent( Inout ) :: pmf
  Type( configuration_type ),   Intent( InOut ) :: config
  Type( comms_type ), Intent( InOut ) :: comm

  Integer :: ipmf,jpmf,j,k

! Loop over all local to this node PMF constraints, their units and
! save the indices of the members that are present on my domain (no halo)
! update lstitr

  Do ipmf=1,pmf%ntpmf

! Initialise indices

     pmf%indpmf(:,:,ipmf) = 0

     Do jpmf=1,2
        If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then
           Do k=1,pmf%mxtpmf(jpmf)
              j=local_index(pmf%listpmf(k,jpmf,ipmf),config%nlast,config%lsi,config%lsa)
              pmf%indpmf(k,jpmf,ipmf)=j
              If (j > 0 .and. j <= config%natms) lstitr(j)=(config%lfrzn(j) == 0)
           End Do
        End If
     End Do

  End Do

! Get PMF units' COM vectors

  Call pmf_coms(pmf,pmf%pxx,pmf%pyy,pmf%pzz,config,comm)

End Subroutine pmf_tags

Subroutine pmf_units_set(pmf,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for setting the existence of the PMF units on the
! node where supposedly PMF constraints exist
!
! Note: (1) Deals with pmf%listpmf, pmf%legpmf and pmf%ntpmf
!       (2) Applies only at the end of relocate_particles
!           if megpmf>0 and mxnode>1
!
! copyright - daresbury laboratory
! author    - i.t.todorov january 2009
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Type (pmf_type ), Intent( InOut ) :: pmf
  Type( comms_type ), Intent( InOut) :: comm
  Type( configuration_type ), Intent( InOut ) :: config
  Logical :: safe,ok
  Integer :: fail,ipmf,jpmf,gpmf

  Integer, Dimension( : ), Allocatable :: i1pmf0,i2pmf0
  Character( Len = 256 ) :: message

  fail=0
  Allocate (i1pmf0(1:pmf%mxtpmf(1)),i2pmf0(1:pmf%mxtpmf(2)), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_units_set allocation failure'
     Call error(0,message)
  End If

! Initialise safety flag

  safe=.true.

! is it ok not to compress the bookkeeping neigh%list arrays
! since it's safe - there's enough buffering space

  ok=.true.
  If (pmf%mxpmf > 0) ok=.not.(Real(pmf%ntpmf,wp)/Real(pmf%mxpmf,wp) > 0.85_wp)

! Sort out local PMF units presence and legend array

  ipmf=0
  Do While (ipmf < pmf%ntpmf)
     ipmf=ipmf+1
10   Continue

! This holds the global PMF index

     gpmf=pmf%listpmf(0,1,ipmf)
     If (gpmf > 0) Then

! For presence of : PMF unit 1 only - this holds 1
!                   PMF unit 2 only - this holds 2
!                   both units 1&2  - this holds 3
! It CANNOT and MUST NOT hold ZERO

        pmf%listpmf(0,2,ipmf)=0

        Do jpmf=1,pmf%mxtpmf(1)
           i1pmf0(jpmf)=local_index(pmf%listpmf(jpmf,1,ipmf),config%natms,config%lsi,config%lsa)

           If (i1pmf0(jpmf) > 0) Then

! This identifies which local PMF constraint the particle belongs to

              pmf%legpmf(0,i1pmf0(jpmf))=1
              pmf%legpmf(1,i1pmf0(jpmf))=ipmf

           End If
        End Do
        If (Any(i1pmf0 > 0)) pmf%listpmf(0,2,ipmf)=pmf%listpmf(0,2,ipmf)+1

        Do jpmf=1,pmf%mxtpmf(2)
           i2pmf0(jpmf)=local_index(pmf%listpmf(jpmf,2,ipmf),config%natms,config%lsi,config%lsa)

           If (i2pmf0(jpmf) > 0) Then

! This identifies which local PMF constraint the particle belongs to

              pmf%legpmf(0,i2pmf0(jpmf))=1
              pmf%legpmf(1,i2pmf0(jpmf))=ipmf

           End If
        End Do
        If (Any(i2pmf0 > 0)) pmf%listpmf(0,2,ipmf)=pmf%listpmf(0,2,ipmf)+2

! If the PMF has moved to another node, pmf%listpmf(0,2,ipmf)=0,
! compress pmf%listpmf

        If (pmf%listpmf(0,2,ipmf) == 0 .and. (.not.ok)) Then
           If      (ipmf  < pmf%ntpmf) Then
              pmf%listpmf(:,:,ipmf)=pmf%listpmf(:,:,pmf%ntpmf) ! Copy neigh%list content from 'ipmf' to 'pmf%ntpmf'
              pmf%listpmf(:,:,pmf%ntpmf)=0                 ! Remove neigh%list content in 'pmf%ntpmf'
              pmf%ntpmf=pmf%ntpmf-1                        ! Reduce 'pmf%ntpmf' pointer

              Go To 10 ! Go back and check it all again for the new neigh%list content in 'ipmf'
           Else If (ipmf == pmf%ntpmf) Then
              pmf%listpmf(:,:,pmf%ntpmf)=0                 ! Remove neigh%list content in 'pmf%ntpmf=ipmf'
              pmf%ntpmf=pmf%ntpmf-1                        ! Reduce 'pmf%ntpmf' pointer
           End If
        End If

     Else

        safe=.false.

     End If
  End Do

! check for inconsistently built local neigh%list

  Call gcheck(comm,safe)
  If (.not.safe) Call error(490)

  Deallocate (i1pmf0,i2pmf0, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_units_set deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_units_set
Subroutine pmf_vcoms(xpmf,ypmf,zpmf,pmf,config,comm)

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


  Real( Kind = wp ), Intent(   Out ) :: xpmf(1:,1:),ypmf(1:,1:),zpmf(1:,1:)
  Type(pmf_type), Intent( Inout ) :: pmf
  Type( configuration_type ), Intent( InOut ) :: config
  Type( comms_type), Intent( InOut ) :: comm

  Logical                 :: safe

  Integer                 :: fail,gpmf,gpmf1,gpmf2,ipmf,jpmf,iadd, &
                             j,k,l

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message

  fail=0
  Allocate (buffer(1:(pmf%mxtpmf(1)+pmf%mxtpmf(2))*(pmf%mxpmf+2)), Stat=fail)
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
  Do gpmf=1,pmf%mxpmf

! Loop over all local to this node PMF constraints matching the global one

     Do ipmf=1,pmf%ntpmf
        If (pmf%listpmf(0,1,ipmf) == gpmf) Then

! Loop over all PMF units present on my domain (no halo)

           Do jpmf=1,2
              If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then

! Loop over their members present on my domain (no halo)

                 Do k=1,pmf%mxtpmf(jpmf)
                    j=pmf%indpmf(k,jpmf,ipmf)

! Copy particles' velocities to buffer in an orderly manner

                    If (j > 0 .and. j <= config%natms) Then ! j is a domain particle
                       l=((gpmf-gpmf1)*(pmf%mxtpmf(1)+pmf%mxtpmf(2))+(jpmf-1)*pmf%mxtpmf(1)+(k-1))*iadd
                       buffer(l+1)=config%vxx(j)
                       buffer(l+2)=config%vyy(j)
                       buffer(l+3)=config%vzz(j)
                    End If
                 End Do
              End If
           End Do
        End If
     End Do

! Check if it safe to fill up the buffer

     safe=(gpmf-gpmf1+2 < (pmf%mxpmf+2)/3)

! If not safe or we've finished looping over all global PMFs

     If ((.not.safe) .or. gpmf == pmf%mxpmf) Then
        Call gsum(comm,buffer)

        Do gpmf2=gpmf1,gpmf
           Do ipmf=1,pmf%ntpmf
              If (pmf%listpmf(0,1,ipmf) == gpmf2) Then
                 Do jpmf=1,2

! Get the COM momentum of this unit

                    Do k=1,pmf%mxtpmf(jpmf)
                       l=((gpmf2-gpmf1)*(pmf%mxtpmf(1)+pmf%mxtpmf(2))+(jpmf-1)*pmf%mxtpmf(1)+(k-1))*iadd

                       xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf) + pmf%pmfwg1(k,jpmf)*buffer(l+1)
                       ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf) + pmf%pmfwg1(k,jpmf)*buffer(l+2)
                       zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf) + pmf%pmfwg1(k,jpmf)*buffer(l+3)
                    End Do

                    xpmf(jpmf,ipmf) = xpmf(jpmf,ipmf)*pmf%pmfwg1(0,jpmf)
                    ypmf(jpmf,ipmf) = ypmf(jpmf,ipmf)*pmf%pmfwg1(0,jpmf)
                    zpmf(jpmf,ipmf) = zpmf(jpmf,ipmf)*pmf%pmfwg1(0,jpmf)

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
           config,strpmf,   &
           virpmf,stat,pmf,comm)

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
  Type( configuration_type  ), Intent( InOut ) :: config
  Real( Kind = wp ), Intent(   Out ) :: strpmf(1:9),virpmf
  Type( stats_type ), Intent( InOut ) :: stat
  Type(pmf_type), Intent( Inout ) :: pmf
  Type( comms_type ), Intent( InOut ) :: comm


  Logical                 :: safe
  Integer                 :: fail,ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: amt(1:2),gamma,gamm(1:2),tstep2,tmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: pxt,pyt,pzt,pt2,esig
  Character( Len = 256 ) :: message

  fail=0
  Allocate (pxt(1:pmf%mxpmf),pyt(1:pmf%mxpmf),pzt(1:pmf%mxpmf),pt2(1:pmf%mxpmf),esig(1:pmf%mxpmf), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'pmf_shake allocation failure'
     Call error(0,message)
  End If


  If (pmf%newjob_shake) Then
    pmf%newjob_shake = .false.

! Get reciprocal PMF units' masses

     Do jpmf=1,2
       If (pmf%pmffrz(jpmf) == pmf%mxtpmf(jpmf)) Then
         pmf%rmass_pmf_unit(jpmf)=0.0_wp
         Else
           pmf%rmass_pmf_unit(jpmf)=pmf%pmfwg1(0,jpmf)
         End If
      End Do

! set PMF constraint parameters

     pmf%dis=pmf%prmpmf
     pmf%dis2=pmf%dis**2
  End If

! squared timestep and reciprocal masses

  tstep2 = tstep*tstep
  amt = tstep2*pmf%rmass_pmf_unit

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

     Call pmf_coms(pmf,pxt,pyt,pzt,config,comm)

! calculate maximum error in bondlength

     Do ipmf=1,pmf%ntpmf
       pt2(ipmf) =pxt(ipmf)**2+pyt(ipmf)**2+pzt(ipmf)**2 - pmf%dis2
        esig(ipmf)=0.5_wp*Abs(pt2(ipmf))/pmf%dis
     End Do

! global verification of convergence

     safe=Merge(Maxval(esig(1:pmf%ntpmf)) < tolnce,.true.,pmf%ntpmf > 0)
     Call gcheck(comm,safe,"enforce")

! bypass next section and terminate iteration if all tolerances ok

     If (.not.safe) Then

! calculate PMF constraint forces

        Do ipmf=1,pmf%ntpmf

! calculate PMF constraint force parameter

           gamma = -pt2(ipmf) / &
                   ((amt(1)+amt(2))*(pmf%pxx(ipmf)*pxt(ipmf)+pmf%pyy(ipmf)*pyt(ipmf)+pmf%pzz(ipmf)*pzt(ipmf)))
           tmp   = gamma / Real(pmf%mxtpmf(1)+pmf%mxtpmf(2),wp)

           Do jpmf=1,2

! If this unit is present on my domain

              If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then
                 gamm(jpmf) = 0.5_wp*Real(1-2*Mod(jpmf,2),wp)*gamma*amt(jpmf)

                 Do k=1,pmf%mxtpmf(jpmf)
                    l=pmf%indpmf(k,jpmf,ipmf)

! for non-frozen domain particles accumulate PMF constraint stress
! and atomic position corrections

                    If (l > 0 .and. l <= config%natms) Then ! l is a domain particle
                       If (config%lfrzn(l) == 0) Then
                          strpmf(1) = strpmf(1) - tmp*pmf%pxx(ipmf)*pmf%pxx(ipmf)
                          strpmf(2) = strpmf(2) - tmp*pmf%pxx(ipmf)*pmf%pyy(ipmf)
                          strpmf(3) = strpmf(3) - tmp*pmf%pxx(ipmf)*pmf%pzz(ipmf)
                          strpmf(5) = strpmf(5) - tmp*pmf%pyy(ipmf)*pmf%pyy(ipmf)
                          strpmf(6) = strpmf(6) - tmp*pmf%pyy(ipmf)*pmf%pzz(ipmf)
                          strpmf(9) = strpmf(9) - tmp*pmf%pzz(ipmf)*pmf%pzz(ipmf)

                          config%parts(l)%xxx=config%parts(l)%xxx+pmf%pxx(ipmf)*gamm(jpmf)
                          config%parts(l)%yyy=config%parts(l)%yyy+pmf%pyy(ipmf)*gamm(jpmf)
                          config%parts(l)%zzz=config%parts(l)%zzz+pmf%pzz(ipmf)*gamm(jpmf)
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
           Do ipmf=1,pmf%ntpmf
             If (esig(ipmf) >= tolnce) Then
                 Write(message,'(3(a,i10),a)') &
                 'global PMF constraint number', ipmf, &
                 ' , with head particle numbers, U1:', pmf%listpmf(1,1,ipmf), &
                 ' & U2:', pmf%listpmf(1,2,ipmf)
               Call info(message)
               Write(message,'(a,f8.2,a,1p,e12.4)') &
                 'converges to a length of', Sqrt(pt2(ipmf)+pmf%dis2), &
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
     stat%passpmf(1,1,1)=Real(icyc-1,wp)
     stat%passpmf(3,1,1)=stat%passpmf(2,1,1)*stat%passpmf(3,1,1)
     stat%passpmf(2,1,1)=stat%passpmf(2,1,1)+1.0_wp
     stat%passpmf(3,1,1)=stat%passpmf(3,1,1)/stat%passpmf(2,1,1)+stat%passpmf(1,1,1)/stat%passpmf(2,1,1)
     stat%passpmf(4,1,1)=Min(stat%passpmf(1,1,1),stat%passpmf(4,1,1))
     stat%passpmf(5,1,1)=Max(stat%passpmf(1,1,1),stat%passpmf(5,1,1))

     stat%passpmf(1,2,1)=stat%passpmf(1,2,1)+stat%passpmf(1,1,1)
     stat%passpmf(1,1,1)=0.0_wp ! Reset
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
           config,stat,pmf,comm)

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
  Type( configuration_type ), Intent( InOut ) :: config
  Type( stats_type ), Intent( InOut ) :: stat
  Type(pmf_type), Intent( Inout ) :: pmf
  Type( comms_type ), Intent( InOut ) :: comm


  Logical                 :: safe
  Integer                 :: fail(1:2),ipmf,jpmf,k,l,icyc
  Real( Kind = wp )       :: dis,esig,gamma,gamm(1:2)

  Real( Kind = wp ), Dimension( : )   , Allocatable :: vxt,vyt,vzt
  Real( Kind = wp ), Dimension( :, : ), Allocatable :: xpmf,ypmf,zpmf
  Character( Len = 256 ) :: message

  fail=0
  Allocate (vxt(1:config%mxatms),vyt(1:config%mxatms),vzt(1:config%mxatms),             Stat=fail(1))
  Allocate (xpmf(1:2,1:pmf%mxpmf),ypmf(1:2,1:pmf%mxpmf),zpmf(1:2,1:pmf%mxpmf), Stat=fail(2))
  If (Any(fail > 0)) Then
    Write(message,'(a)') 'pmf_rattle allocation failure'
     Call error(0,message)
  End If

! Get PMF units' reciprocal masses

  If (pmf%newjob_rattle) Then
    pmf%newjob_rattle = .false.

     Do jpmf=1,2
        If (pmf%pmffrz(jpmf) == pmf%mxtpmf(jpmf)) Then
          pmf%amtr(jpmf)=0.0_wp
        Else
          pmf%amtr(jpmf)=pmf%pmfwg1(0,jpmf)
        End If
      End Do
    End If

! normalise PMF constraint vectors on first pass outside

    If (lfst) Then
      Do ipmf=1,pmf%ntpmf
        dis=1.0_wp/Sqrt(pmf%pxx(ipmf)**2+pmf%pyy(ipmf)**2+pmf%pzz(ipmf)**2)
        pmf%pxx(ipmf)=pmf%pxx(ipmf)*dis
        pmf%pyy(ipmf)=pmf%pyy(ipmf)*dis
        pmf%pzz(ipmf)=pmf%pzz(ipmf)*dis
      End Do
    End If

! application of PMF constraint (rattle) algorithm
! Initialise number of cycles to zero and unsafe passage of the algorithm

    safe=.false.
    icyc=0
    Do While ((.not.safe) .and. icyc < mxshak)
      icyc=icyc+1

! initialise velocity correction arrays

     Do l=1,config%natms
        vxt(l)=0.0_wp
        vyt(l)=0.0_wp
        vzt(l)=0.0_wp
      End Do

! calculate temporary COM velocity of each unit

     Call pmf_vcoms(xpmf,ypmf,zpmf,pmf,config,comm)

! calculate PMF velocity corrections

     esig=0.0_wp
     Do ipmf=1,pmf%ntpmf

! calculate constraint force parameter - gamma

        gamma = pmf%pxx(ipmf)*(xpmf(1,ipmf)-xpmf(2,ipmf)) + &
                pmf%pyy(ipmf)*(ypmf(1,ipmf)-ypmf(2,ipmf)) + &
                pmf%pzz(ipmf)*(zpmf(1,ipmf)-zpmf(2,ipmf))

        esig=Max(esig,0.5_wp*tstep*Abs(gamma))

        gamma = gamma / (pmf%amtr(1)+pmf%amtr(2))

        Do jpmf=1,2

! If this unit is present on my domain

          If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then
            gamm(jpmf) = Real(1-2*Mod(jpmf,2),wp)*gamma*pmf%amtr(jpmf)
            Do k=1,pmf%mxtpmf(jpmf)
              l=pmf%indpmf(k,jpmf,ipmf)

! improve approximate PMF particles velocity and force (if non-frozen)

                 If (l > 0 .and. l <= config%natms) Then ! l is a domain particle
                    If (config%lfrzn(l) == 0) Then
                       vxt(l)=vxt(l)+pmf%pxx(ipmf)*gamm(jpmf)
                       vyt(l)=vyt(l)+pmf%pyy(ipmf)*gamm(jpmf)
                       vzt(l)=vzt(l)+pmf%pzz(ipmf)*gamm(jpmf)
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
        Do ipmf=1,pmf%ntpmf
           Do jpmf=1,2
              If (pmf%listpmf(0,2,ipmf) == jpmf .or. pmf%listpmf(0,2,ipmf) == 3) Then
                 Do k=1,pmf%mxtpmf(jpmf)
                    l=pmf%indpmf(k,jpmf,ipmf)
                    If (l > 0 .and. l <= config%natms) Then ! l is a domain particle
                       If (config%lfrzn(l) == 0) Then
                          config%vxx(l)=config%vxx(l)+vxt(l)
                          config%vyy(l)=config%vyy(l)+vyt(l)
                          config%vzz(l)=config%vzz(l)+vzt(l)
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
     stat%passpmf(1,1,2)=Real(icyc-1,wp)
     stat%passpmf(3,1,2)=stat%passpmf(2,1,2)*stat%passpmf(3,1,2)
     stat%passpmf(2,1,2)=stat%passpmf(2,1,2)+1.0_wp
     stat%passpmf(3,1,2)=stat%passpmf(3,1,2)/stat%passpmf(2,1,2)+stat%passpmf(1,1,2)/stat%passpmf(2,1,2)
     stat%passpmf(4,1,2)=Min(stat%passpmf(1,1,2),stat%passpmf(4,1,2))
     stat%passpmf(5,1,2)=Max(stat%passpmf(1,1,2),stat%passpmf(5,1,2))

     stat%passpmf(1,2,2)=stat%passpmf(1,2,2)+stat%passpmf(1,1,2)
     If (lcol) Then ! Collect
        stat%passpmf(3,2,2)=stat%passpmf(2,2,2)*stat%passpmf(3,2,2)
        stat%passpmf(2,2,2)=stat%passpmf(2,2,2)+1.0_wp
        stat%passpmf(3,2,2)=stat%passpmf(3,2,2)/stat%passpmf(2,2,2)+stat%passpmf(1,2,2)/stat%passpmf(2,2,2)
        stat%passpmf(4,2,2)=Min(stat%passpmf(1,2,2),stat%passpmf(4,2,2))
        stat%passpmf(5,2,2)=Max(stat%passpmf(1,2,2),stat%passpmf(5,2,2))
        stat%passpmf(1,2,2)=0.0_wp ! Reset
     End If
     stat%passpmf(1,1,2)=0.0_wp ! Reset
  End If

  Deallocate (vxt,vyt,vzt,    Stat=fail(1))
  Deallocate (xpmf,ypmf,zpmf, Stat=fail(2))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'pmf_rattle deallocation failure'
     Call error(0,message)
  End If

End Subroutine pmf_rattle


End module pmf
