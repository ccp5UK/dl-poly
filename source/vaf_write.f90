Subroutine vaf_write(lvafav,keyres,nstep,tstep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing VAFDAT file at selected intervals
! in simulation
!
! copyright - daresbury laboratory
! author    - i.t.todorov & m.a.seaton february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,     Only : idnode,mxnode,gcheck
  Use setup_module,     Only : mxatyp,nvafdt,zero_plus
  Use site_module,      Only : unqatm,numtypnf
  Use config_module,    Only : cfgname
  Use greenkubo_module, Only : isvaf,nsvaf,vaftsts,vafcount,vaftime,vaf

  Implicit None

  Logical,           Intent( In    ) :: lvafav
  Integer,           Intent( In    ) :: keyres,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep

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
           If (idnode == 0) Inquire(File='VAFDAT_'//unqatm(i), Exist=lexist)
           If (mxnode > 1) Call gcheck(lexist)
        Else
           lexist=.false.
        End If

! Generate file if non-existent or replace it if outputting time-averaged VAF

        If ((.not.lexist) .or. lvafav) Then
           If (idnode == 0) Then
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

     If (idnode == 0) Then
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

        If (idnode == 0) Write(nvafdt,"(1p,2e14.6)") timei,gvaf

     End Do

     If (idnode == 0) Then
        If (.not. lvafav) Write(nvafdt,'(2/)')
        Close(Unit=nvafdt)
     End If
  End Do

End Subroutine vaf_write