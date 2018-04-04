Module greenkubo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring arrays for green-kubo relations
! calculations to calculate transport properties
!
! copyright - daresbury laboratory
! author    - m.a.seaton june 2014
! contrib   - i.t.todorov july 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp
  
  Use comms,     Only : comms_type,gsum,gcheck
  Use setup,     Only : nrite,mxatyp,mxbuff,zero_plus,nvafdt,mxatms
  Use configuration,    Only : natms,ltype,lfrzn,vxx,vyy,vzz,cfgname
  Use site,      Only : unqatm,numtypnf

  Use errors_warnings, Only : error
  Implicit None

  Integer           :: isvaf   = 1 , & ! VAF sampling frequency in steps
                       nsvaf   = 0 , & ! VAF sample size
                       vaftsts =-1 , & ! VAF timestep start
                       vafsamp = 0     ! VAF simultaneously overlapping samples

  Real( Kind = wp ) :: vafcount = 0.0_wp

  Integer,           Allocatable, Save :: vafstep(:)
  Real( Kind = wp ), Allocatable, Save :: vxi(:,:),vyi(:,:),vzi(:,:)
  Real( Kind = wp ), Allocatable, Save :: vafdata(:,:),vaftime(:),vaf(:,:)

!  Real( Kind = wp ), Allocatable, Save :: stxx(:),stxy(:),stxz(:),styy(:),styz(:),stzz(:)
!  Real( Kind = wp ), Allocatable, Save :: gkpot(:),tcond(:),tctime(:)

  Public :: allocate_greenkubo_arrays

Contains

  Subroutine allocate_greenkubo_arrays()

    Integer :: i
    Integer, Dimension( 1:7 ) :: fail

    fail = 0

    Allocate (vafstep(1:vafsamp),                                                           Stat = fail(1))
    Allocate (vxi(1:mxatms,1:vafsamp),vyi(1:mxatms,1:vafsamp),vzi(1:mxatms,1:vafsamp),      Stat = fail(2))
    Allocate (vafdata(0:nsvaf,1:vafsamp*(mxatyp+1)),vaftime(0:nsvaf),vaf(0:nsvaf,1:mxatyp), Stat = fail(3))
!    Allocate (gkpot(1:mxatms),                                                              Stat = fail(4))
!    Allocate (stxx(1:mxatms),stxy(1:mxatms),stxz(1:mxatms),                                 Stat = fail(5))
!    Allocate (styy(1:mxatms),styz(1:mxatms),stzz(1:mxatms),                                 Stat = fail(6))
!    Allocate (tcond(0:mxgkstk),tctime(0:mxgkstk),                                           Stat = fail(7))

    If (Any(fail > 0)) Call error(1080)

    vxi = 0.0_wp ; vyi = 0.0_wp ; vzi = 0.0_wp
    vaf = 0.0_wp ; vaftime = 0.0_wp
    vafdata = 0.0_wp

! setup step counters for samples, starting after equilibration

    Do i=1,vafsamp
      vafstep(i)=(1-i)*isvaf-1
    End Do

!    sxx = 0.0_wp ; sxy = 0.0_wp ; sxz = 0.0_wp ; syy = 0.0_wp; syz = 0.0_wp; szz = 0.0_wp
!    tcond = 0.0_wp ; tctime = 0.0_wp

  End Subroutine allocate_greenkubo_arrays
  
  Subroutine vaf_collect(lvafav,leql,nsteql,nstep,time,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistics for velocity
! autocorrelation functions
!
! copyright - daresbury laboratory
! author    - m.a.seaton june 2014
! amended   - i.t.todorov july 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( In    ) :: lvafav,leql
  Integer,           Intent( In    ) :: nsteql,nstep
  Real( Kind = wp ), Intent( In    ) :: time
  Type( comms_type ), Intent( InOut ) :: comm
  Integer                        :: fail,i,j,k,l,nsum

  Real( Kind = wp ), Allocatable :: vafcoll(:)

  Character ( Len = 256 )  ::  message

  fail=0
  Allocate (vafcoll(1:mxatyp), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'vaf_collect allocation failure'
     Call error(0,message)
  End If

! set VAF timestep start

  If (vaftsts < 0) vaftsts=Merge(nsteql,nstep,leql)

! set reference velocities only for first sampling passes

  If (nstep <= vaftsts+(vafsamp-1)*isvaf) Then
     Do j=1,vafsamp
        If (nstep == (vaftsts+(j-1)*isvaf)) Then
           vxi(1:natms,j) = vxx(1:natms)
           vyi(1:natms,j) = vyy(1:natms)
           vzi(1:natms,j) = vzz(1:natms)
        End If
     End Do
  End If

! advance time step, calculate vaf contribution for each sample

  Do j=1,vafsamp
     vafstep(j) = vafstep(j) + 1
     k=vafstep(j)

     If (k >= 0 .and. k <= nsvaf) Then
        vafcoll=0.0_wp
        Do i=1,natms
           If (lfrzn(i) == 0) Then
              l=ltype(i)
              vafcoll(l) = vafcoll(l) + vxx(i)*vxi(i,j)+vyy(i)*vyi(i,j)+vzz(i)*vzi(i,j)
           End If
        End Do

        Do l=1,mxatyp
           vafdata(k,(j-1)*(mxatyp+1)+l) = vafcoll(l)
        End Do

        vafdata(k,j*(mxatyp+1)) = time
     End If

     If (k == nsvaf) Then
        vafcount = vafcount + 1.0_wp

        If (comm%mxnode > 1) Then
           nsum = mxbuff/(nsvaf+1)

           l=(j-1)*(mxatyp+1) ! avoid summing up timing information
           Do i=1,mxatyp,nsum
              Call gsum(comm,vafdata(:,l+i:l+Min(i+nsum-1,mxatyp)))
           End Do
        End If

        If (lvafav) Then ! if time-averaging, add vaf data to sampling array
           Do i=0,nsvaf
              vaftime(i) = vafdata(i,j*(mxatyp+1))
              Do l=1,mxatyp
                 vaf(i,l) = vaf(i,l) + vafdata(i,(j-1)*(mxatyp+1)+l)
              End Do
           End Do
        Else             ! if not time-averaging, move vaf data to sampling array
           vaftime(0:nsvaf) = vafdata(0:nsvaf,j*(mxatyp+1))
           Do l=1,mxatyp
              vaf(0:nsvaf,l) = vafdata(0:nsvaf,(j-1)*(mxatyp+1)+l)
           End Do
        End If
     End If

! reset counter and reference velocities and get vaf
! at end of time span or at start of sample

     If (k == nsvaf .or. (isvaf > nsvaf .and. k == isvaf)) Then
        vafstep(j) = 0

        vxi(1:natms,j) = vxx(1:natms)
        vyi(1:natms,j) = vyy(1:natms)
        vzi(1:natms,j) = vzz(1:natms)

        vafcoll = 0.0_wp
        Do i=1,natms
           If (lfrzn(i) == 0) Then
              l=ltype(i)
              vafcoll(l) = vafcoll(l) + vxi(i,j)*vxi(i,j)+vyi(i,j)*vyi(i,j)+vzi(i,j)*vzi(i,j)
           End If
        End Do

        vafdata(0,j*(mxatyp+1)) = time
        Do l=1,mxatyp
           vafdata(0,(j-1)*(mxatyp+1)+l) = vafcoll(l)
        End Do
     End If
  End Do

  Deallocate (vafcoll, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'vaf_collect deallocation failure'
     Call error(0,message)
  End If

End Subroutine vaf_collect

Subroutine vaf_compute(lvafav,tstep,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating velocity autocorrelation
! functions from accumulated data
!
! copyright - daresbury laboratory
! author    - m.a.seaton june 2014
! amended   - i.t.todorov july 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Logical,           Intent( In    ) :: lvafav
  Real( Kind = wp ), Intent( In    ) :: tstep
  Type( comms_type ), Intent( InOut ) :: comm

  Integer           :: i
  Real( Kind = wp ) :: factor,gvaf,ovaf,time0,timei,numt

  If (comm%idnode == 0) Then
     If (lvafav) Then
        Write(nrite,"(/,/,12x,'VELOCITY AUTOCORRELATION FUNCTIONS',/,/, &
             & 'calculated using ',i8,' samples')") Nint(vafcount)
     Else
        Write(nrite,"(/,/,12x,'VELOCITY AUTOCORRELATION FUNCTIONS',/,/, &
             & 'calculated using sample ',i8,' starting at ',f10.4,' ps')") Nint(vafcount), vaftime(0)
     End If
  End If

  time0 = vaftime(0)
  numt = Sum(numtypnf(1:mxatyp))
  factor = 1.0_wp/Sum(vaf(0,1:mxatyp))
  ovaf = Sum(vaf(0,1:mxatyp))/Real(numt,Kind=wp)
  If (lvafav) ovaf = ovaf/vafcount

  If (comm%idnode == 0) Write(nrite,"(12x,'absolute value at origin (3kT/m) = ',1p,e16.8)") ovaf

! loop over time steps

  Do i=0,nsvaf

! determine time

     If (lvafav) Then
        timei = tstep*Real(i,Kind=wp)
     Else
        timei = vaftime(i)-time0
     End If

! null it if number of particles is zero

     If (numt > zero_plus) Then
        gvaf = Sum(vaf(i,1:mxatyp))*factor
     Else ! THIS CAN NEVER HAPPEN AS TOTAL NON-FROZEN PARTICLES IS ALWAYS > 1
        gvaf = 0.0_wp
     End If

! print out information

     If (comm%idnode == 0) Write(nrite,"(12x,f10.4,1p,e14.6)") timei,gvaf

  End Do

End Subroutine vaf_compute

Subroutine vaf_write(lvafav,keyres,nstep,tstep,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing VAFDAT file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov & m.a.seaton february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( In    ) :: lvafav
  Integer,           Intent( In    ) :: keyres,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep
  Type( comms_type), Intent( InOut ) :: comm

  Logical,     Save :: newjob = .true.

  Logical           :: lexist
  Integer           :: i,j
  Real( Kind = wp ) :: factor,gvaf,ovaf,time0,timei,numt

  If (vaftsts < 0) Return
  If (.not.(nstep >= vaftsts+nsvaf .and. Mod(nstep-vaftsts-nsvaf,isvaf) == 0)) Return

  If (newjob) Then
     newjob = .false.

! If the keyres=1, is VAFDAT old (does it exist) and
! how many frames and records are in there

     Do i=1,mxatyp
        lexist=.true.
        If (keyres == 1) Then
           If (comm%idnode == 0) Inquire(File='VAFDAT_'//unqatm(i), Exist=lexist)
           Call gcheck(comm,lexist)
        Else
           lexist=.false.
        End If

! Generate file if non-existent or replace it if outputting time-averaged VAF

        If ((.not.lexist) .or. lvafav) Then
           If (comm%idnode == 0) Then
              Open(Unit=nvafdt, File='VAFDAT_'//unqatm(i), Status='replace')
              Write(nvafdt,'(a)') cfgname
              Close(Unit=nvafdt)
           End If
        End If
     End Do
  End If

  time0 = vaftime(0)

! loop over species types

  Do j=1,mxatyp
     numt = numtypnf(j)
     factor = 1.0_wp
     If (Abs(vaf(0,j)) > 0.0e-6_wp) factor = 1.0_wp/vaf(0,j)
     ovaf = vaf(0,j)/Real(numt,Kind=wp)
     If (lvafav) ovaf=ovaf/vafcount

! replace file (if outputting time-averaged VAF); prepare for data set

     If (comm%idnode == 0) Then
        If (lvafav) Then
           Open(Unit=nvafdt, File='VAFDAT_'//unqatm(j), Status='replace')
           Write(nvafdt,'(a)') cfgname
           Close(Unit=nvafdt)
        End If
        Open(Unit=nvafdt, File='VAFDAT_'//unqatm(j), Position='append')
        Write(nvafdt,'(a8,i10,1p,e16.8,f20.6)') unqatm(j),nsvaf,ovaf,time0
     End If

! Then loop over time steps

     Do i=0,nsvaf

! determine time

        If (lvafav) Then
           timei = tstep*Real(i,Kind=wp)
        Else
           timei = vaftime(i)-time0
        End If

! null it if number of particles is zero

        If (numt > zero_plus) Then
           gvaf = vaf(i,j)*factor
        Else
           gvaf = 0.0_wp
        End If

! print out information

        If (comm%idnode == 0) Write(nvafdt,"(1p,2e14.6)") timei,gvaf

     End Do

     If (comm%idnode == 0) Then
        If (.not. lvafav) Write(nvafdt,'(2/)')
        Close(Unit=nvafdt)
     End If
  End Do

End Subroutine vaf_write


End Module greenkubo
