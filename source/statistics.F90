Module statistics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module declaring global simulation property variables and
! arrays
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi,li
  Use setup, Only : mxbfss,nstats,statis,zero_plus,&
                   prsunt,tenunt,boltz,engunit,eu_ev,eu_kcpm,&
                   eu_kjpm,mxatyp,pi,mxatms

  Use comms,   Only : comms_type,gsum,Spread_tag,wp_mpi,gtime,gmax,gsend, &
                      gwait,girecv,gcheck
  Use site, Only : site_type
  Use configuration,  Only : cfgname,imcon,cell,volm,natms,ltype, &
                             xxx,yyy,zzz,vxx,vyy,vzz,ixyz,lsa,lsi,ltg
  Use domains,    Only : nprx,npry,nprz,map,r_nprx,r_npry,r_nprz,&
                                nprx_r,npry_r,nprz_r,idx,idy,idz
  Use neighbours,  Only  : neighbours_type

  Use core_shell,  Only : passshl
  Use z_density,   Only : z_density_type,z_density_collect
  Use msd,         Only : msd_type
  Use greenkubo,   Only : greenkubo_type
  Use errors_warnings, Only : error,warning,info
  Use numerics,    Only : dcell,invert,shellsort,shellsort2,pbcshfrc,pbcshfrl
  Use thermostat, Only : thermostat_type

  Implicit None
  Type, Public :: stats_type

  Integer( Kind = wi )           :: numacc = 0 , &
                                    natms0 = 0
  Integer( Kind = wi )           :: mxnstk,mxstak,intsta
  Logical                        :: statis_file_open = .false.

  Real( Kind = wp )              :: consv = 0.0_wp,shlke = 0.0_wp,engke = 0.0_wp,&
                                    engrot = 0.0_wp,engcpe = 0.0_wp,engsrp = 0.0_wp,&
                                    engter  = 0.0_wp,engtbp = 0.0_wp,engfbp = 0.0_wp,&
                                    engshl = 0.0_wp,engtet = 0.0_wp,engbnd = 0.0_wp,&
                                    engang = 0.0_wp,engdih = 0.0_wp,enginv = 0.0_wp,&
                                    engfld = 0.0_wp,engcon = 0.0_wp,engpmf = 0.0_wp

  Real( Kind = wp )              :: stptmp = 0.0_wp,stpprs = 0.0_wp,stpvol = 0.0_wp,&
                                    stpcfg = 0.0_wp,stpeng = 0.0_wp,stpeth = 0.0_wp,&
                                    stpvir = 0.0_wp

  Real( Kind = wp )              :: virtot = 0.0_wp,vircom = 0.0_wp,vircpe = 0.0_wp,&
                                    virsrp = 0.0_wp,virshl = 0.0_wp,virter = 0.0_wp,&
                                    virtbp = 0.0_wp,virfbp = 0.0_wp,vircon = 0.0_wp,&
                                    virpmf = 0.0_wp,virtet = 0.0_wp,virbnd = 0.0_wp,&
                                    virang = 0.0_wp,virdih = 0.0_wp,virinv = 0.0_wp,&
                                    virfld = 0.0_wp,virdpd = 0.0_wp

  Real( Kind = wp )              :: strtot(1:9) = 0.0_wp,strkin(1:9) = 0.0_wp,strknf(1:9) = 0.0_wp,&
                                    strknt(1:9) = 0.0_wp,strcom(1:9) = 0.0_wp,strcon(1:9) = 0.0_wp,&
                                    strpmf(1:9) = 0.0_wp, stress(1:9) = 0.0_wp,strdpd(1:9) = 0.0_wp

  Real( Kind = wp )              :: clin(1:9) = 0.0_wp
  ! constraints accumulators
  Real( Kind = wp ),              Public :: passcnq(1:5) = (/ & ! QUENCHING per call
    0.0_wp         ,  & ! cycles counter
    0.0_wp         ,  & ! access counter
    0.0_wp         ,  & ! average cycles
    999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
    0.0_wp /)           ! maximum cycles
  Real( Kind = wp ),              Public :: passcon(1:5,1:2,1:2) = Reshape( (/ & ! dim::1-shake, dim:1:-per-call
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::1-shake, dim:2:-per-tst
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:1:-per-call
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:2:-per-tst
    0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp /) , (/5,2,2/) )

  Real( Kind = wp ),              Public :: passpmq(1:5) = (/ & ! QUENCHING per call
                                          0.0_wp         ,  & ! cycles counter
                                          0.0_wp         ,  & ! access counter
                                          0.0_wp         ,  & ! average cycles
                                          999999999.0_wp ,  & ! minimum cycles : ~Huge(1)
                                          0.0_wp /)           ! maximum cycles
  Real( Kind = wp ),              Public :: passpmf(1:5,1:2,1:2) = Reshape( (/ & ! dim::1-shake, dim:1:-per-call
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::1-shake, dim:2:-per-tst
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:1:-per-call
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp , & ! dim::2-rattle, dim:2:-per-tst
                            0.0_wp, 0.0_wp, 0.0_wp, 999999999.0_wp, 0.0_wp /) , (/5,2,2/) )

  Real( Kind = wp ), Allocatable :: xin(:),yin(:),zin(:)
  Real( Kind = wp ), Allocatable :: xto(:),yto(:),zto(:),rsd(:)

  Real( Kind = wp ), Allocatable :: stpval(:),stpvl0(:),sumval(:),ssqval(:)
  Real( Kind = wp ), Allocatable :: zumval(:),ravval(:),stkval(:,:)

  Integer,           Allocatable :: found(:), found0(:)
  Integer,           Allocatable :: lsi0(:),lsa0(:),ltg0(:)

  Real( Kind = wp ), Allocatable :: xin0(:),yin0(:),zin0(:)
  Real( Kind = wp ), Allocatable :: xto0(:),yto0(:),zto0(:)

  Real( Kind = wp ), Allocatable :: stpval0(:),stpvl00(:),sumval0(:),ssqval0(:)
  Real( Kind = wp ), Allocatable :: zumval0(:),ravval0(:),stkval0(:,:)

End Type

  Public :: allocate_statistics_arrays, allocate_statistics_connect, &
            deallocate_statistics_connect

Contains

  Subroutine allocate_statistics_arrays(mxatdm,stats)
    Integer, Intent( In ) ::           mxatdm
    Type( stats_type ), Intent( InOut ) ::  stats

    Integer, Dimension( 1:4 ) :: fail
    Integer ::  mxnstk ,mxstak

    fail = 0

    mxnstk = stats%mxnstk
    mxstak = stats%mxstak
    
    Allocate (stats%xin(1:mxatdm),stats%yin(1:mxatdm),stats%zin(1:mxatdm),                           Stat = fail(1))
    Allocate (stats%xto(1:mxatdm),stats%yto(1:mxatdm),stats%zto(1:mxatdm),stats%rsd(1:mxatdm),             Stat = fail(2))
    Allocate (stats%stpval(0:mxnstk),stats%stpvl0(0:mxnstk),stats%sumval(0:mxnstk),stats%ssqval(0:mxnstk), Stat = fail(3))
    Allocate (stats%zumval(0:mxnstk),stats%ravval(0:mxnstk),stats%stkval(1:mxstak,0:mxnstk),         Stat = fail(4))
    If (Any(fail > 0)) Call error(1016)

    stats%xin = 0.0_wp ; stats%yin = 0.0_wp ; stats%zin = 0.0_wp
    stats%xto = 0.0_wp ; stats%yto = 0.0_wp ; stats%zto = 0.0_wp ; stats%rsd = 0.0_wp

    stats%stpval = 0.0_wp ; stats%stpvl0 = 0.0_wp ; stats%sumval = 0.0_wp ; stats%ssqval = 0.0_wp
    stats%zumval = 0.0_wp ; stats%ravval = 0.0_wp ; stats%stkval = 0.0_wp

  End Subroutine allocate_statistics_arrays

  Subroutine allocate_statistics_connect(mxatdm,stats)
    Integer, Intent( InOut ) ::           mxatdm
    Type( stats_type ), Intent( InOut ) ::  stats

    Integer, Dimension( 1:6 ) :: fail
    Integer :: mxstak

    mxstak=stats%mxstak
    fail = 0

    Allocate (stats%found(1:mxatdm),stats%found0(1:mxatdm),                                                Stat = fail(1))
    Allocate (stats%lsi0(1:mxatdm),stats%lsa0(1:mxatdm),stats%ltg0(1:mxatdm),                                    Stat = fail(2))
    Allocate (stats%xin0(1:mxatdm),stats%yin0(1:mxatdm),stats%zin0(1:mxatdm),                                    Stat = fail(3))
    Allocate (stats%xto0(1:mxatdm),stats%yto0(1:mxatdm),stats%zto0(1:mxatdm),                                    Stat = fail(4))
    Allocate (stats%stpval0(1:2*mxatdm),stats%stpvl00(1:2*mxatdm),stats%sumval0(1:2*mxatdm),stats%ssqval0(1:2*mxatdm), &
      Stat = fail(5))
    Allocate (stats%zumval0(1:2*mxatdm),stats%ravval0(1:2*mxatdm),stats%stkval0(1:mxstak,1:2*mxatdm),            Stat = fail(6))

    If (Any(fail > 0)) Call error(1060)

  End Subroutine allocate_statistics_connect

  Subroutine deallocate_statistics_connect(stats)
    Type( stats_type ), Intent( InOut ) ::  stats

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Deallocate (stats%found,stats%found0,                    Stat = fail(1))
    Deallocate (stats%lsi0,stats%lsa0,stats%ltg0,                  Stat = fail(2))
    Deallocate (stats%xin0,stats%yin0,stats%zin0,                  Stat = fail(3))
    Deallocate (stats%xto0,stats%yto0,stats%zto0,                  Stat = fail(4))
    Deallocate (stats%stpval0,stats%stpvl00,stats%sumval0,stats%ssqval0, Stat = fail(5))
    Deallocate (stats%zumval0,stats%ravval0,stats%stkval0,         Stat = fail(6))

    If (Any(fail > 0)) Call error(1061)

  End Subroutine deallocate_statistics_connect

  Subroutine statistics_collect           &
           (lsim,leql,nsteql,lmsd, &
           keyres,                 &
           degfre,degshl,degrot,          &
           nstep,tstep,time,tmst,         &
           mxatdm,stats,thermo,zdensity,site,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating periodic data during the
! molecular dynamics simulation and computing the rolling averages
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov march 2016
! contrib   - a.m.elena february 2017
! contrib   - i.t.todorov february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( In    ) :: lsim,leql,lmsd
  Integer,           Intent( In    ) :: nsteql,keyres, &
                                        nstep

  Integer(Kind=li),  Intent( In    ) :: degfre,degshl,degrot

  Real( Kind = wp ), Intent( In    ) :: tstep,time

  Real( Kind = wp ), Intent( InOut ) :: tmst
  Integer( Kind = wi),Intent( In    ) :: mxatdm
  Type( stats_type ), Intent( InOut ) :: stats
  Type( thermostat_type ), Intent( In    ) :: thermo
  Type( z_density_type ), Intent( InOut ) :: zdensity
  Type( site_type ), Intent( In    ) :: site
  Type( comms_type ), Intent( InOut ) :: comm

  Logical,           Save :: newjob = .true.

  Logical                 :: l_tmp
  Integer                 :: fail,i,j,k,iadd,kstak
  Real( Kind = wp )       :: stpcns,stpshl,stprot,stpipv,celprp(1:10), &
                             zistk,sclnv1,sclnv2,h_z

  Real( Kind = wp ), Allocatable :: amsd(:)
  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)
  Character( Len = 256 ) :: message


  fail=0
  Allocate (amsd(1:mxatyp), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'statistics_collect allocation failure'
     Call error(0,message)
  End If

! open statistics file and put header

  If (newjob .and. comm%idnode == 0) Then
     newjob = .false.

! If the keyres=1 is the file old (does it exist)?

     l_tmp=.false.
     If (keyres == 1) Inquire(File=Trim(statis), Exist=l_tmp)

     If (.not.l_tmp) Then
        Open(Unit=nstats, File=Trim(statis), Status='replace')
        stats%statis_file_open = .true.

        Write(nstats,'(a)') cfgname

        If      (Abs(engunit - eu_ev)   <= zero_plus) Then
           Write(nstats,'(1x,a)') 'ENERGY UNITS = electron Volts'
        Else If (Abs(engunit - eu_kcpm) <= zero_plus) Then
           Write(nstats,'(1x,a)') 'ENERGY UNITS = kcal/mol'
        Else If (Abs(engunit - eu_kjpm) <= zero_plus) Then
           Write(nstats,'(1x,a)') 'ENERGY UNITS = kjoule/mol'
        Else If (Abs(engunit - 1.0_wp)  <= zero_plus) Then
           Write(nstats,'(1x,a)') 'ENERGY UNITS = DL_POLY Internal UNITS (10 J/mol)'
        Else If (Abs(engunit - boltz)   <= zero_plus) Then
           Write(nstats,'(1x,a)') 'ENERGY UNITS = Kelvin/Boltzmann'
        Else ! once in a blue moon
           Write(nstats,'(1x,a)') 'ENERGY UNITS = DPD (Unknown)'
        End If

     End If
  End If

! instantaneous properties of system

! configurational energy

  stats%stpcfg = stats%engcpe + stats%engsrp + stats%engter + stats%engtbp + stats%engfbp + &
           stats%engfld + stats%engshl +                            &
           stats%engtet + stats%engbnd + stats%engang + stats%engdih + stats%enginv

! system energy

  stats%stpeng = stats%stpcfg + stats%engke + stats%engrot

! energy + conserved quantity (for true ensembles)

  stpcns = stats%stpeng + stats%consv

! rotational temperature

  stprot = 2.0_wp*(stats%engrot) / (boltz*Max(1.0_wp,Real(degrot,wp)))

! core-shell units temperature

  stpshl = 2.0_wp*(stats%shlke) / (boltz*Max(1.0_wp,Real(degshl,wp)))

! system temperature

  stats%stptmp = 2.0_wp*(stats%engke+stats%engrot) / (boltz*Real(degfre,wp))

! system virial
! Note: originally, purely angle dependent interactions have zero virial!!!
! So, stats%virfbp, stats%virinv and stats%virdih are allegedly always zero!  virdih has an exception!

  stats%stpvir = stats%vircpe + stats%virsrp + stats%virter + stats%virtbp + stats%virfbp + &
           stats%virfld + stats%virshl + stats%vircon + stats%virpmf + stats%vircom + &
           stats%virtet + stats%virbnd + stats%virang + stats%virdih + stats%virinv + stats%virdpd

! system volume

  stats%stpvol = volm

! system pressure

  stats%stpprs = (2.0_wp*stats%engke-stats%stpvir) / (3.0_wp*stats%stpvol)

! system PV

  stpipv = stats%stpprs*stats%stpvol

! system enthalpy

  If (thermo%variable_cell) Then             ! P_target*V_instantaneous
     stats%stpeth = stats%stpeng + thermo%press*stats%stpvol
  Else                               ! for thermo%variable_cell=.false. V_instantaneous=V_target
     stats%stpeth = stats%stpeng + stpipv        ! and there is only P_instantaneous
  End If

  Call dcell(cell,celprp)

! store current values in statistics array

  stats%stpval(0) =stats%consv/engunit
  stats%stpval(1) =stpcns/engunit
  stats%stpval(2) =stats%stptmp
  stats%stpval(3) =stats%stpcfg/engunit
  stats%stpval(4) =(stats%engsrp+stats%engter)/engunit
  stats%stpval(5) =stats%engcpe/engunit
  stats%stpval(6) =stats%engbnd/engunit
  stats%stpval(7) =(stats%engang+stats%engtbp)/engunit
  stats%stpval(8) =(stats%engdih+stats%enginv+stats%engfbp)/engunit
  stats%stpval(9) =stats%engtet/engunit
  stats%stpval(10)=stats%stpeth/engunit
  stats%stpval(11)=stprot
  stats%stpval(12)=stats%stpvir/engunit
  stats%stpval(13)=(stats%virsrp+stats%virter)/engunit
  stats%stpval(14)=stats%vircpe/engunit
  stats%stpval(15)=stats%virbnd/engunit
  stats%stpval(16)=(stats%virtbp+stats%virang)/engunit
  stats%stpval(17)=stats%vircon/engunit
  stats%stpval(18)=stats%virtet/engunit
  stats%stpval(19)=stats%stpvol
  stats%stpval(20)=stpshl
  stats%stpval(21)=stats%engshl/engunit
  stats%stpval(22)=stats%virshl/engunit
  stats%stpval(23)=Acos(celprp(6))*180.0_wp/pi
  stats%stpval(24)=Acos(celprp(5))*180.0_wp/pi
  stats%stpval(25)=Acos(celprp(4))*180.0_wp/pi
  stats%stpval(26)=stats%virpmf/engunit
  stats%stpval(27)=stats%stpprs*prsunt

  iadd = 27

! iadd = iadd + 1 ! for the stpval(0)!!! Thus to account for in printing

! mean squared displacements per species, dependent on
! particle displacements from initial positions (at t=0)

  amsd = 0.0_wp ! initialise

  If (nstep == nsteql+1) Then ! re-initialise
     Do i=1,natms
        stats%xto(i) = 0.0_wp
        stats%yto(i) = 0.0_wp
        stats%zto(i) = 0.0_wp
     End Do
  End If

  If (nstep > 0) Then
     If (lsim)  Then ! real dynamics is happening
        Do i=1,natms
           stats%xto(i)=stats%xto(i)+vxx(i)*tstep
           stats%yto(i)=stats%yto(i)+vyy(i)*tstep
           stats%zto(i)=stats%zto(i)+vzz(i)*tstep
        End Do
     Else            ! HISTORY is replayed
        Allocate (xxt(1:mxatms),yyt(1:mxatms),zzt(1:mxatms), Stat=fail)
        If (fail > 0) Then
           Write(message,'(a)') 'statistics_collect allocation failure 1'
           Call error(0,message)
        End If
        Do i=1,natms
           xxt(i)=xxx(i)
           yyt(i)=yyy(i)
           zzt(i)=zzz(i)
        End Do
        Call pbcshfrc(imcon,cell,natms,xxt,yyt,zzt)
        Call pbcshfrc(imcon,stats%clin,natms,stats%xin,stats%yin,stats%zin)
        Do i=1,natms
           stats%xin(i)=xxt(i)-stats%xin(i)
           stats%yin(i)=yyt(i)-stats%yin(i)
           stats%zin(i)=zzt(i)-stats%zin(i)
        End Do
        Deallocate (xxt,yyt,zzt, Stat=fail)
        If (fail > 0) Then
           Write(message,'(a)') 'statistics_collect deallocation failure 1'
           Call error(0,message)
        End If
        Call pbcshfrl(imcon,cell,natms,stats%xin,stats%yin,stats%zin)
        Do i=1,natms
           stats%xto(i)=stats%xto(i)+stats%xin(i)
           stats%yto(i)=stats%yto(i)+stats%yin(i)
           stats%zto(i)=stats%zto(i)+stats%zin(i)
        End Do
     End If

     Do i=1,natms
        stats%rsd(i)=Sqrt(stats%xto(i)**2+stats%yto(i)**2+stats%zto(i)**2)

        k=ltype(i)
        amsd(k)=amsd(k)+stats%rsd(i)**2
     End Do
     Call gsum(comm,amsd(1:site%ntype_atom))
  End If

  If (lmsd) Then
     Do i=1,natms
        j=2*i
        stats%stpval(iadd+j-1)=stats%rsd(i)**2
        stats%stpval(iadd+j  )=vxx(i)**2+vyy(i)**2+vzz(i)**2
     End Do
     iadd = iadd + 2*mxatdm
  End If

  Do k=1,site%ntype_atom
     If (site%num_type_nf(k) > zero_plus) Then
        stats%stpval(iadd+k)=amsd(k)/site%num_type_nf(k)
     Else
        stats%stpval(iadd+k)=0.0_wp
     End If
  End Do
  iadd = iadd + site%ntype_atom

! pressure tensor (derived for the stress tensor)

  Do i=1,9
     stats%stpval(iadd+i)=stats%strtot(i)*prsunt/stats%stpvol
  End Do
  iadd = iadd + 9

  If (thermo%variable_cell) Then

! cell parameters

     Do i=1,9
       stats%stpval(iadd+i)=cell(i)
     End Do
     iadd = iadd + 9

! instantaneous PV

     stats%stpval(iadd+1)=stpipv/engunit
     iadd = iadd + 1

     If (thermo%iso > 0) Then
        h_z=celprp(9)

        stats%stpval(iadd+1)=h_z
        stats%stpval(iadd+2)=stats%stpvol/h_z
        iadd = iadd + 2

        If (thermo%iso > 1) Then
           stats%stpval(iadd+1)= -h_z*(stats%strtot(1)-(thermo%press+thermo%stress(1)))*tenunt
           stats%stpval(iadd+2)= -h_z*(stats%strtot(5)-(thermo%press+thermo%stress(5)))*tenunt
           iadd = iadd + 2
        End If
     End If
  End If

! write statistics file

  If (comm%idnode == 0 .and. Mod(nstep,stats%intsta) == 0) Then
     If (.not. stats%statis_file_open) Then
         Open(Unit=nstats, File=Trim(STATIS), Position='append')
         stats%statis_file_open = .true.
     End If

     If (lmsd) Then
        Write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))') &
             nstep,time,iadd+1-2*mxatdm,stats%stpval(1:  27),stats%stpval(0),stats%stpval(28+2*mxatdm:iadd)
     Else
        Write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))') &
             nstep,time,iadd+1,         stats%stpval(1:  27),stats%stpval(0),stats%stpval(28         :iadd)
     End If

  End If

! check on number of variables for stack

  If (iadd > stats%mxnstk) Call error(170)

! No totals for timestep zero

  If (nstep == 0) Go To 10

! current stack value

  kstak=Mod(nstep-1,stats%mxstak)+1

! subtract old stack value from the stack average

  If (nstep > stats%mxstak) Then
     Do i=0,stats%mxnstk
        stats%zumval(i)=stats%zumval(i)-stats%stkval(kstak,i)
     End Do
  End If

! store quantities in stack and update the stack average

  Do i=0,stats%mxnstk
     stats%stkval(kstak,i)=stats%stpval(i)
     stats%zumval(i)=stats%zumval(i)+stats%stpval(i)
  End Do

! calculate rolling averages

  zistk=Real(Min(stats%mxstak,nstep),wp)

  Do i=0,stats%mxnstk
     stats%ravval(i)=stats%zumval(i)/zistk
  End Do

! accumulate totals over steps

  If ((.not.leql) .or. nstep > nsteql) Then
     stats%numacc=stats%numacc+1
     sclnv2=1.0_wp/Real(stats%numacc,wp)
     sclnv1=Real(stats%numacc-1,wp)/Real(stats%numacc,wp)

! average squared sum and sum (keep in this order!!!)

     If (nstep == nsteql+1 .or. ((.not.leql) .and. nstep == 1)) stats%stpvl0=stats%stpval
     stats%stpval=stats%stpval-stats%stpvl0
     Do i=0,stats%mxnstk
        stats%ssqval(i)=sclnv1*(stats%ssqval(i)+sclnv2*(stats%stpval(i)-stats%sumval(i))**2)

!stats%sumval has to be shifted back tostats%sumval+stpvl0 in statistics_result
! when averaging is printed since stpval is only shifted back and forth
! which does not affect the fluctuations Sqrtstats%ssqval) only their accuracy

        stats%sumval(i)=sclnv1*stats%sumval(i)+sclnv2*stats%stpval(i)
     End Do
     stats%stpval=stats%stpval+stats%stpvl0
  End If

10 Continue

! z-density collection

  If ( zdensity%l_collect .and. ((.not.leql) .or. nstep >= nsteql) .and. &
       Mod(nstep,zdensity%frequency) == 0 ) Call z_density_collect(zdensity)

! Catch time of starting statistical averages

  If (((.not.leql) .or. nstep == nsteql) .and. tmst < tstep) tmst=time

  Deallocate (amsd, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'statistics_collect deallocation failure'
     Call error(0,message)
  End If

End Subroutine statistics_collect


Subroutine statistics_connect_frames(megatm,mxatdm,lmsd,stats,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of data between neighbouring
! domains/nodes in order to reconnect some statistical information
! between replayed frames of history
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Integer( Kind = wi ), Intent ( In    ) :: megatm
  Integer( Kind = wi ), Intent ( In    ) :: mxatdm
  Logical,            Intent( In    ) :: lmsd
  Type( stats_type ), Intent( InOut ) :: stats
  Type( comms_type ), Intent( InOut ) :: comm

  Integer :: icyc,nres
  Character( Len = 256 ) :: message

  stats%found = 0 ; icyc = 0 ; nres = 1
  Do While (icyc <= Max(nprx,npry,nprz)/2 .and. nres > 0)
     Call match_compress_spread_sort(-1,mxatdm) ! -x direction spread
     Call match_compress_spread_sort( 1,mxatdm) ! +x direction spread

     Call match_compress_spread_sort(-2,mxatdm) ! -y direction spread
     Call match_compress_spread_sort( 2,mxatdm) ! +y direction spread

     Call match_compress_spread_sort(-3,mxatdm) ! -z direction spread
     Call match_compress_spread_sort( 3,mxatdm) ! +z direction spread

     Call match_compress_spread_sort( 0,mxatdm) ! no spreading

     nres=stats%natms0
     Call gsum(comm,nres)
     If (nres > 0) Then
        nres=Merge(0,Sum(stats%found(1:natms)),natms > 0)
        Call gsum(comm,nres)
        If (nres /= megatm) icyc = icyc + 1
     End If
  End Do

  If (nres > 0) Then
    Write(message,'(a)') ' particles dynamics properties will be corrupted'
    Call warning(message,.true.)
  End If

Contains

  Subroutine match_compress_spread_sort(mdir,mxatdm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to simplify the repetition of the procedures above
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Integer( Kind = wi ), Intent( In    ) :: mdir ! +/-1,+/-2,+/-3,0 is the direction of spread
    Integer( Kind = wi ), Intent( In    ) :: mxatdm

    Integer( Kind = wi ) :: fail,i,i0,j,j0,kk

    Integer( Kind = wi ), Allocatable, Save :: lsa00(:)

! Search for matches
    
    stats%found0 = 0
    If (stats%natms0 > 0) Then
       i0 = 1
       Do i=1,natms
          Do While (i0 <= stats%natms0)
             If      (lsa(i) <  stats%lsa0(i0)) Then
                Exit
             Else If (lsa(i) == stats%lsa0(i0)) Then
                If (stats%found(lsi(i)) > 0) Then ! ghost arrival
                   stats%found0(stats%lsi0(i0)) = 1     ! erase at compression
                Else                        ! new arrival to claim
                   stats%found(lsi(i)) = 1 ; stats%found0(stats%lsi0(i0)) = 1

                   stats%xin(lsi(i)) = stats%xin0(stats%lsi0(i0))
                   stats%yin(lsi(i)) = stats%yin0(stats%lsi0(i0))
                   stats%zin(lsi(i)) = stats%zin0(stats%lsi0(i0))

                   stats%xto(lsi(i)) = stats%xto0(stats%lsi0(i0))
                   stats%yto(lsi(i)) = stats%yto0(stats%lsi0(i0))
                   stats%zto(lsi(i)) = stats%zto0(stats%lsi0(i0))

                   If (lmsd) Then
                      j =27+2*lsi(i)
                      j0=2*stats%lsi0(i0)
                      stats%stpvl0(j-1)=stats%stpvl00(j0-1)
                      stats%stpvl0(j  )=stats%stpvl00(j0  )
                      stats%stpval(j-1)=stats%stpval0(j0-1)
                      stats%stpval(j  )=stats%stpval0(j0  )
                      stats%zumval(j-1)=stats%zumval0(j0-1)
                      stats%zumval(j  )=stats%zumval0(j0  )
                      stats%ravval(j-1)=stats%ravval0(j0-1)
                      stats%ravval(j  )=stats%ravval0(j0  )
                      stats%ssqval(j-1)=stats%ssqval0(j0-1)
                      stats%ssqval(j  )=stats%ssqval0(j0  )
                      stats%sumval(j-1)=stats%sumval0(j0-1)
                      stats%sumval(j  )=stats%sumval0(j0  )
                      Do kk=1,stats%mxstak
                         stats%stkval(kk,j-1)=stats%stkval0(kk,j0-1)
                         stats%stkval(kk,j  )=stats%stkval0(kk,j0  )
                      End Do
                   End If
                End If
             End If
             i0 = i0 + 1                    ! move along
          End Do
       End Do
    End If

! Invalidate lazies and deallocate at last use

    If (mdir ==  0) Then
       i = 1
       Do i0=1,stats%natms0
          Do While (lsa00(i) /= 0 .and. i < mxatdm)
             If      (stats%lsa0(i0) <  lsa00(i)) Then
                Exit
             Else If (stats%lsa0(i0) == lsa00(i)) Then
                stats%found0(stats%lsi0(i0)) = 1        ! erase at compression
             End If
             i = i + 1                      ! move along
          End Do
       End Do

       Deallocate (lsa00, Stat = fail)
       If (fail > 0) Then
          Write(message,'(a)') 'match_compress_spread_sort deallocation failure'
          Call error(0,message)
       End If
    End If

! Compress remainder

    i0 = 1
    Do While (i0 <= stats%natms0 .and. stats%natms0 > 0)
       If (stats%found0(i0) == 0) Then          ! Not claimed
          If (ixyz(i0) == 0) Then         ! pass along
             If      (mdir == -1) Then    ! -x to do
                ixyz(i0) = 333            ! all since b4 0
             Else If (mdir ==  1) Then    !  x to do
                ixyz(i0) = 331            ! not back to -1
             Else If (mdir == -2) Then    ! -y to do
                ixyz(i0) = 332            ! not back to  1
             Else If (mdir ==  2) Then    !  y to do
                ixyz(i0) = 313            ! not back to -2
             Else If (mdir == -3) Then    ! -z to do
                ixyz(i0) = 323            ! not back to  2
             Else If (mdir ==  3) Then    !  z to do
                ixyz(i0) = 133            ! not back to -3
             Else If (mdir ==  0) Then    ! end of cycle to do
                ixyz(i0) = 233            ! not back to  3
             Else                         ! abort
                Call error(160)
             End If
          End If

          i0 = i0 + 1                     ! Increase lower bound marker
       Else                               ! claimed, to erase entry
          If (stats%found0(stats%natms0) == 0) Then   ! try to refill with the last unclaimed entry
             ixyz(i0) = ixyz(stats%natms0)
             If (ixyz(i0) == 0) Then      ! pass along
                If      (mdir == -1) Then ! -x to do
                   ixyz(i0) = 333         ! all since b4 0
                Else If (mdir ==  1) Then !  x to do
                   ixyz(i0) = 331         ! not back to -1
                Else If (mdir == -2) Then ! -y to do
                   ixyz(i0) = 332         ! not back to  1
                Else If (mdir ==  2) Then !  y to do
                   ixyz(i0) = 313         ! not back to -2
                Else If (mdir == -3) Then ! -z to do
                   ixyz(i0) = 323         ! not back to  2
                Else If (mdir ==  3) Then !  z to do
                   ixyz(i0) = 133         ! not back to -3
                Else If (mdir ==  0) Then ! end of cycle to do
                   ixyz(i0) = 233         ! not back to  3
                Else                      ! abort
                   Call error(160)
                End If
             End If
             stats%ltg0(i0) = stats%ltg0(stats%natms0)

             stats%xin0(i0) = stats%xin0(stats%natms0)
             stats%yin0(i0) = stats%yin0(stats%natms0)
             stats%zin0(i0) = stats%zin0(stats%natms0)

             stats%xto0(i0) = stats%xto0(stats%natms0)
             stats%yto0(i0) = stats%yto0(stats%natms0)
             stats%zto0(i0) = stats%zto0(stats%natms0)

             If (lmsd) Then
                j =2*i0
                j0=2*stats%natms0
                stats%stpvl00(j-1)=stats%stpvl00(j0-1)
                stats%stpvl00(j  )=stats%stpvl00(j0  )
                stats%stpval0(j-1)=stats%stpval0(j0-1)
                stats%stpval0(j  )=stats%stpval0(j0  )
                stats%zumval0(j-1)=stats%zumval0(j0-1)
                stats%zumval0(j  )=stats%zumval0(j0  )
                stats%ravval0(j-1)=stats%ravval0(j0-1)
                stats%ravval0(j  )=stats%ravval0(j0  )
                stats%ssqval0(j-1)=stats%ssqval0(j0-1)
                stats%ssqval0(j  )=stats%ssqval0(j0  )
                stats%sumval0(j-1)=stats%sumval0(j0-1)
                stats%sumval0(j  )=stats%sumval0(j0  )
                Do kk=1,stats%mxstak
                   stats%stkval0(kk,j-1)=stats%stkval0(kk,j0-1)
                   stats%stkval0(kk,j  )=stats%stkval0(kk,j0  )
                End Do
             End If

             i0 = i0 + 1                  ! increase lower bound marker if entry is refilled
          End If

!! Erase upper holdings in either case
!
!          ixyz(stats%natms0) = 0
!         stats%ltg0(stats%natms0) = 0
!
!          stats%xin0(stats%natms0) = 0
!          stats%yin0(stats%natms0) = 0
!          stats%zin0(stats%natms0) = 0
!
!         stats%xto0(stats%natms0) = 0
!         stats%yto0(stats%natms0) = 0
!         stats%zto0(stats%natms0) = 0
!
!          If (lmsd) Then
!             j0=2*stats%natms0
!            stats%stpvl00(j0-1)=0.0_wp
!            stats%stpvl00(j0  )=0.0_wp
!            stats%stpval0(j0-1)=0.0_wp
!            stats%stpval0(j0  )=0.0_wp
!            stats%zumval0(j0-1)=0.0_wp
!            stats%zumval0(j0  )=0.0_wp
!            stats%ravval0(j0-1)=0.0_wp
!            stats%ravval0(j0  )=0.0_wp
!            stats%ssqval0(j0-1)=0.0_wp
!            stats%ssqval0(j0  )=0.0_wp
!            stats%sumval0(j0-1)=0.0_wp
!            stats%sumval0(j0  )=0.0_wp
!             Do kk=1,mxstak
!               stats%stkval0(kk,j0-1)=0.0_wp
!               stats%stkval0(kk,j0  )=0.0_wp
!             End Do
!          End If
          stats%natms0 = stats%natms0 - 1             ! Decrease upper bound marker
       End If
    End Do

! Allocate and initialise at first use
! Detect unknown lazies sort them lsa like

    If (mdir == -1) Then
       fail = 0
       Allocate (lsa00(1:mxatdm), Stat = fail)
       If (fail > 0) Then
          Write(message,'(a)') 'match_compress_spread_sort allocation failure'
          Call error(0,message)
       End If
       lsa00=0

       i=0
       Do i0=1,stats%natms0
          If (ixyz(i0) == 333) Then
             i=i+1
             lsa00(i)=stats%ltg0(i0)
          End If
       End Do
       Call shellsort(i,lsa00)
    End If

! Spread atom data in the mdir direction

    If (mdir /= 0) Call statistics_connect_spread(mdir,mxatdm,lmsd,stats,comm)

! Sort past frame remainder of global atom indices

!    lsi0=0 ; lsa0=0
    Do i0=1,stats%natms0
       stats%lsi0(i0)=i0
       stats%lsa0(i0)=stats%ltg0(i0)
    End Do
    Call shellsort2(stats%natms0,stats%lsi0,stats%lsa0)

  End Subroutine match_compress_spread_sort

End Subroutine statistics_connect_frames

  Subroutine statistics_connect_set(rcut,mxatdm,lmsd,stats,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of data between neighbouring
! domains/nodes in order to reconnect some statistical information
! between replayed frames of history
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Real( Kind = wp ), Intent( In    ) :: rcut
  Integer( Kind = wi), Intent( In   ) :: mxatdm
  Logical,            Intent( In    ) :: lmsd
  Type( stats_type ), Intent( InOut ) :: stats
  Type( comms_type ), Intent( InOut ) :: comm

  Real( Kind = wp ), Save :: cut

  Integer           :: nlx,nly,nlz,i,i0,kk
  Real( Kind = wp ) :: det,celprp(1:10),rcell(1:9),x,y,z, &
                       xdc,ydc,zdc,cwx,cwy,cwz,ecwx,ecwy,ecwz

  If (comm%mxnode > 1) Then

! Define cut

     cut=rcut+1.0e-6_wp

     Call dcell(cell,celprp)
     Call invert(cell,rcell,det)

! calculate link cell dimensions per node

     nlx=Int(celprp(7)/(cut*nprx_r))
     nly=Int(celprp(8)/(cut*npry_r))
     nlz=Int(celprp(9)/(cut*nprz_r))

! Get the total number of link-cells in MD cell per direction

     xdc=Real(nlx*nprx,wp)
     ydc=Real(nly*npry,wp)
     zdc=Real(nlz*nprz,wp)

! link-cell widths in reduced space

     cwx=1.0_wp/xdc
     cwy=1.0_wp/ydc
     cwz=1.0_wp/zdc

! Distance from the - edge of this domain

     ecwx=Nearest( (-0.5_wp+cwx)+Real(idx,wp)*r_nprx , +1.0_wp)+zero_plus
     ecwy=Nearest( (-0.5_wp+cwy)+Real(idy,wp)*r_npry , +1.0_wp)+zero_plus
     ecwz=Nearest( (-0.5_wp+cwz)+Real(idz,wp)*r_nprz , +1.0_wp)+zero_plus

! Distance from the + edge of this domain with a possible
! extension strip for the one linked cell per domain scenario

     cwx=Nearest( (-0.5_wp-cwx)+Real(idx+1,wp)*r_nprx , -1.0_wp)-zero_plus-Merge( cwx*1.0e-10_wp , 0.0_wp , nlx == 1 )
     cwy=Nearest( (-0.5_wp-cwy)+Real(idy+1,wp)*r_npry , -1.0_wp)-zero_plus-Merge( cwy*1.0e-10_wp , 0.0_wp , nly == 1 )
     cwz=Nearest( (-0.5_wp-cwz)+Real(idz+1,wp)*r_nprz , -1.0_wp)-zero_plus-Merge( cwz*1.0e-10_wp , 0.0_wp , nlz == 1 )

     ixyz(1:mxatdm)=0 ! Initialise move (former halo) indicator
     Do i=1,natms
        x=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        y=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        z=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

        If (x <= ecwx) ixyz(i)=ixyz(i)+1
        If (x >=  cwx) ixyz(i)=ixyz(i)+2

        If (y <= ecwy) ixyz(i)=ixyz(i)+10
        If (y >=  cwy) ixyz(i)=ixyz(i)+20

        If (z <= ecwz) ixyz(i)=ixyz(i)+100
        If (z >=  cwz) ixyz(i)=ixyz(i)+200
     End Do

     lsi=0;lsa=0 ! This is a must, unfortunately
     Do i=1,natms
        lsi(i)=i
        lsa(i)=ltg(i)
     End Do
     Call shellsort2(natms,lsi,lsa)

     stats%natms0 = natms
     stats%ltg0(1:stats%natms0) = ltg(1:stats%natms0) !;stats%ltg0(stats%natms0+1: ) = 0
     stats%lsa0(1:stats%natms0) = lsa(1:stats%natms0) !; lsa0(stats%natms0+1: ) = 0
     stats%lsi0(1:stats%natms0) = lsi(1:stats%natms0) !; lsi0(stats%natms0+1: ) = 0

     stats%xin0(1:stats%natms0) = stats%xin(1:stats%natms0) !; stats%xin0(stats%natms0+1: ) = 0 ; stats%xin = 0.0_wp
     stats%yin0(1:stats%natms0) = stats%yin(1:stats%natms0) !; stats%yin0(stats%natms0+1: ) = 0 ; stats%yin = 0.0_wp
     stats%zin0(1:stats%natms0) = stats%zin(1:stats%natms0) !; stats%zin0(stats%natms0+1: ) = 0 ; stats%zin = 0.0_wp

    stats%xto0(1:stats%natms0) =stats%xto(1:stats%natms0) !;stats%xto0(stats%natms0+1: ) = 0
    stats%yto0(1:stats%natms0) =stats%yto(1:stats%natms0) !;stats%yto0(stats%natms0+1: ) = 0
    stats%zto0(1:stats%natms0) =stats%zto(1:stats%natms0) !;stats%zto0(stats%natms0+1: ) = 0

     If (lmsd) Then
        i0=2*stats%natms0
       stats%stpvl00(1:i0)=stats%stpvl0(28:27+i0) !;stats%stpvl00(i0+1: )=0.0_wp
       stats%stpval0(1:i0)=stats%stpval(28:27+i0) !;stats%stpval0(i0+1: )=0.0_wp
       stats%zumval0(1:i0)=stats%zumval(28:27+i0) !;stats%zumval0(i0+1: )=0.0_wp
       stats%ravval0(1:i0)=stats%ravval(28:27+i0) !;stats%ravval0(i0+1: )=0.0_wp
       stats%ssqval0(1:i0)=stats%ssqval(28:27+i0) !;stats%ssqval0(i0+1: )=0.0_wp
       stats%sumval0(1:i0)=stats%sumval(28:27+i0) !;stats%sumval0(i0+1: )=0.0_wp
        Do kk=1,stats%mxstak
          stats%stkval0(kk,1:i0)=stats%stkval(kk,1:i0) !;stats%stkval0(kk,i0+1: )=0.0_wp
        End Do
     End If
  End If

End Subroutine statistics_connect_set

Subroutine statistics_connect_spread(mdir,mxatdm,lmsd,stats,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to spread atomic and topological data of particles
! leaving this domain
!
! NOTE: When executing on one node we need not get here at all!
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  Integer( Kind = wi ), Intent( In    ) :: mdir
  Integer( Kind = wi ), Intent( In    ) :: mxatdm
  Logical,            Intent( In    ) :: lmsd
  Type( stats_type ), Intent( InOut ) :: stats
  Type( comms_type ), Intent( InOut ) :: comm
  Logical           :: safe,stay,move
  Integer           :: fail,iblock,jdnode,kdnode,   &
                       imove,jmove,kmove,keep,send, &
                       i,j,l,jj,kk,jxyz,kxyz,       &
                       ix,iy,iz,kx,ky,kz,newatm

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
  Character( Len = 256 ) :: message

  If (comm%mxnode == 1) Return

  fail=0
  Allocate (buffer(1:mxbfss), Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'statistics_connect_spread allocation failure'
     Call error(0,message)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=mxbfss/2

! DIRECTION SETTINGS INITIALISATION

! define the neighbouring domains as sending and receiving with
! respect to the direction (mdir)
! k.   - direction selection factor
! jxyz - halo reduction factor
! kxyz - corrected halo reduction factor particles haloing both +&- sides
! jdnode - destination (send to), kdnode - source (receive from)

  kx = 0 ; ky = 0 ; kz = 0
  If      (mdir == -1) Then ! Direction -x
     kx  = 1
     jxyz= 1
     kxyz= Merge(3,jxyz,nprx <= 2)

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     kxyz= Merge(3,jxyz,nprx <= 2)

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     kxyz= Merge(30,jxyz,npry <= 2)

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     kxyz= Merge(30,jxyz,npry <= 2)

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     kxyz= Merge(300,jxyz,nprz <= 2)

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     kxyz= Merge(300,jxyz,nprz <= 2)

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(160)
  End If

! Initialise counters for length of sending and receiving buffers
! buffer(1) and buffer(iblock+1) contain the actual number of
! particles to get transfered, imove and jmove are the lengths of
! the buffers

  imove=1
  jmove=1

! Initialise how many particles are to be kept and sent

  keep=0
  send=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL PREVIOUS FRAME'S PARTICLES ON THIS NODE

  Do i=1,stats%natms0

! particle designated directions

     ix=Mod(ixyz(i),10)           ! [0,1,2,3]
     iy=Mod(ixyz(i)-ix,100)       ! [0,10,20,30]
     iz=Mod(ixyz(i)-(ix+iy),1000) ! [0,100,200,300]

! Filter the move index for the selected direction

     j=ix*kx+iy*ky+iz*kz

! If the particle is scheduled to be sent in the selected
! direction then indicate it in move

     move=.false.
     If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then
        move=(comm%idnode /= jdnode) ! but don't move it if back to itself

! reduce particle move index (ixyz) using the corrected halo reduction
! factor when the particle is sent to both +&- sides

        ixyz(i)=ixyz(i)-Merge(kxyz,jxyz,j /= jxyz)
     End If
     stay = (ixyz(i) /= 0) ! decide on keeping it when to be sent elsewhere

     If (stay) Then ! keep it
        keep=keep+1

! retain config indexing and move indexing arrays

       stats%ltg0(keep)=stats%ltg0(i)
        ixyz(keep)=ixyz(i)

! retain initial positions

        stats%xin0(keep)=stats%xin0(i)
        stats%yin0(keep)=stats%yin0(i)
        stats%zin0(keep)=stats%zin0(i)

! retain final displacements

       stats%xto0(keep)=stats%xto0(i)
       stats%yto0(keep)=stats%yto0(i)
       stats%zto0(keep)=stats%zto0(i)

        If (lmsd) Then
           jj=2*i
           j =2*keep
          stats%stpvl00(j-1)=stats%stpvl00(jj-1)
          stats%stpvl00(j  )=stats%stpvl00(jj  )
          stats%stpval0(j-1)=stats%stpval0(jj-1)
          stats%stpval0(j  )=stats%stpval0(jj  )
          stats%zumval0(j-1)=stats%zumval0(jj-1)
          stats%zumval0(j  )=stats%zumval0(jj  )
          stats%ravval0(j-1)=stats%ravval0(jj-1)
          stats%ravval0(j  )=stats%ravval0(jj  )
          stats%ssqval0(j-1)=stats%ssqval0(jj-1)
          stats%ssqval0(j  )=stats%ssqval0(jj  )
          stats%sumval0(j-1)=stats%sumval0(jj-1)
          stats%sumval0(j  )=stats%sumval0(jj  )
           Do kk=1,stats%mxstak
             stats%stkval0(kk,j-1)=stats%stkval0(kk,jj-1)
             stats%stkval0(kk,j  )=stats%stkval0(kk,jj  )
           End Do
        End If
     End If

     If (move) Then ! copy it
        send=send+1
        If (imove+8 <= iblock) Then ! If safe to proceed

! pack config indexing and move indexing arrays

           buffer(imove+1)=Real(stats%ltg0(i),wp)
           buffer(imove+2)=Real(ixyz(i),wp)

! pack initial positions

           buffer(imove+3)=stats%xin0(i)
           buffer(imove+4)=stats%yin0(i)
           buffer(imove+5)=stats%zin0(i)

! pack final displacements

           buffer(imove+6)=stats%xto0(i)
           buffer(imove+7)=stats%yto0(i)
           buffer(imove+8)=stats%zto0(i)
        Else
           safe=.false.
        End If
        imove=imove+8

! pack MSD arrays

        If (lmsd) Then
           If (imove+2*(6+stats%mxstak) <= iblock) Then
              jj=2*i
              buffer(imove+ 1)=stats%stpvl00(jj-1)
              buffer(imove+ 2)=stats%stpvl00(jj  )
              buffer(imove+ 3)=stats%stpval0(jj-1)
              buffer(imove+ 4)=stats%stpval0(jj  )
              buffer(imove+ 5)=stats%zumval0(jj-1)
              buffer(imove+ 6)=stats%zumval0(jj  )
              buffer(imove+ 7)=stats%ravval0(jj-1)
              buffer(imove+ 8)=stats%ravval0(jj  )
              buffer(imove+ 9)=stats%ssqval0(jj-1)
              buffer(imove+10)=stats%ssqval0(jj  )
              buffer(imove+11)=stats%sumval0(jj-1)
              buffer(imove+12)=stats%sumval0(jj  )
              Do kk=1,stats%mxstak
                 l=12+2*kk
                 buffer(imove+l-1)=stats%stkval0(kk,jj-1)
                 buffer(imove+l  )=stats%stkval0(kk,jj  )
              End Do
           Else
              safe=.false.
           End If
           imove=imove+2*(6+stats%mxstak)
        End If
     End If

  End Do

! Check for array bound overflow (have arrays coped with outgoing data)

  Call gcheck(comm,safe)
  If (.not.safe) Call error(163)

! record of number of atoms for transfer

  buffer(1)=Real(send,wp)

! exchange information on buffer sizes

  Call girecv(comm,jmove,kdnode,Spread_tag)
  Call gsend(comm,imove,jdnode,Spread_tag)
  Call gwait(comm)

! exchange buffers between nodes (this is a MUST)

  Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,Spread_tag)
  Call gsend(comm,buffer(1:imove),jdnode,Spread_tag)
  Call gwait(comm)

! check arrays can cope with incoming atom numbers

  kmove=iblock+1
  jmove=Nint(buffer(kmove))

! Test for overloading and collect how many are to really be accepted

  imove=0
  Do i=1,jmove
     l=Nint(buffer(kmove+1))
     If (All(stats%ltg0(1:stats%natms0) /= l)) imove=imove+1
     kmove=kmove+8
     If (lmsd) kmove=kmove+2*(6+stats%mxstak)
  End Do

  stats%natms0=keep+imove

! Check for array bound overflow (can arrays cope with incoming data)

  safe=(stats%natms0 <= mxatdm)
  Call gcheck(comm,safe)
  If (.not.safe) Call error(164)

! load transferred data

  kmove=iblock+1 ! restore kmove
  newatm=keep    ! restore newatm
  Do i=1,jmove
     If (imove /= jmove) Then
        l=Nint(buffer(kmove+1))
        If (Any(stats%ltg0(1:keep) == l)) Then
           kmove=kmove+8
           If (lmsd) kmove=kmove+2*(6+stats%mxstak)
           Cycle
        End If
     End If

     newatm=newatm+1

! unpack config indexing, site and move indexing arrays

    stats%ltg0(newatm)=Nint(buffer(kmove+1))
     ixyz(newatm)=Nint(buffer(kmove+2))

! unpack initial positions arrays

     stats%xin0(newatm)=buffer(kmove+3)
     stats%yin0(newatm)=buffer(kmove+4)
     stats%zin0(newatm)=buffer(kmove+5)

! unpack initial positions arrays

    stats%xto0(newatm)=buffer(kmove+6)
    stats%yto0(newatm)=buffer(kmove+7)
    stats%zto0(newatm)=buffer(kmove+8)

     kmove=kmove+8

! unpack MSD arrays

     If (lmsd) Then
        jj=2*newatm
       stats%stpvl00(jj-1)=buffer(kmove+1 )
       stats%stpvl00(jj  )=buffer(kmove+2 )
       stats%stpval0(jj-1)=buffer(kmove+3 )
       stats%stpval0(jj  )=buffer(kmove+4 )
       stats%zumval0(jj-1)=buffer(kmove+5 )
       stats%zumval0(jj  )=buffer(kmove+6 )
       stats%ravval0(jj-1)=buffer(kmove+7 )
       stats%ravval0(jj  )=buffer(kmove+8 )
       stats%ssqval0(jj-1)=buffer(kmove+9 )
       stats%ssqval0(jj  )=buffer(kmove+10)
       stats%sumval0(jj-1)=buffer(kmove+11)
       stats%sumval0(jj  )=buffer(kmove+12)
        Do kk=1,stats%mxstak
           l=12+2*kk
          stats%stkval0(kk,jj-1)=buffer(kmove+l-1)
          stats%stkval0(kk,jj  )=buffer(kmove+l  )
        End Do

        kmove=kmove+2*(6+stats%mxstak)
     End If
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'statistics_connect_spread deallocation failure'
     Call error(0,message)
  End If

End Subroutine statistics_connect_spread

Subroutine statistics_result                                    &
           (lmin,lmsd, &
           nstrun,keyshl,megcon,megpmf,              &
           nstep,tstep,time,tmst, &
           mxatdm,stats,thermo,green,neigh,site,comm,passmin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing simulation summary
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov november 2016
! contrib   - m.a.seaton june 2014
! contrib   - a.b.g.chalk january 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( In    ) :: lmin,lmsd
  Integer( Kind = wi ),    Intent( In    ) :: nstrun,keyshl,megcon,megpmf,nstep
  Real( Kind = wp ), Intent( In    ) :: tstep,time,tmst
  Integer( Kind = wi ),    Intent( In    ) :: mxatdm
  Type( stats_type ), Intent( InOut ) :: stats
  Type( thermostat_type ), Intent( In    ) :: thermo
  Type( greenkubo_type ), Intent( In    ) :: green
  Type( neighbours_type ), Intent( InOut ) :: neigh
  Type( site_type ), Intent( In    ) :: site
  Type( comms_type ), Intent( InOut ) :: comm
  Real( Kind = wp ), Intent( In    ) ::  passmin(:)

  Logical           :: check
  Integer           :: i,iadd
  Real( Kind = wp ) :: avvol,dc,srmsd,timelp,tmp,h_z,tx,ty,temp
  Character( Len = 256 ) :: message
  Character( Len = 256 ), Dimension(5) :: messages

  Integer :: mxnstk

  mxnstk = stats%mxnstk

  Call info('',.true.)

! VNL skipping statistics

  If (neigh%unconditional_update .and. nstep > 0) Then
     If (.not.neigh%update) Then ! Include the final skip in skipping statistics
        neigh%skip(3)=neigh%skip(2)*neigh%skip(3)
        neigh%skip(2)=neigh%skip(2)+1.0_wp
        neigh%skip(3)=neigh%skip(3)/neigh%skip(2)+neigh%skip(1)/neigh%skip(2)
        neigh%skip(4)=Min(neigh%skip(1),neigh%skip(4))
        neigh%skip(5)=Max(neigh%skip(1),neigh%skip(5))
     End If

     Write(message,'(a,f7.2,2(a,i4))') &
       'VNL skipping run statistics - skips per timestep: average ',neigh%skip(3), &
       ' minimum ',Nint(Merge(neigh%skip(4),neigh%skip(5),neigh%skip(4)<neigh%skip(5))), &
       ' maximum ',Nint(neigh%skip(5))
     Call info(message,.true.)
  End If

! minimisation convergence statistics

  If (lmin) Then
     Write(message,'(a,f7.2,2(a,i4))') &
       'minimisation run statistics - cycles per call: average ',passmin(3), &
       ' minimum ',Nint(passmin(4)),' maximum ',Nint(passmin(5))
    Call info(message,.true.)
  End If

! shell relaxation convergence statistics

  If (keyshl == 2) Then
     Write(message,'(a,f7.2,2(a,i4))') &
       'hell relaxation run statistics - cycles per timestep: average ',passshl(3), &
       ' minimum ',Nint(passshl(4)),' maximum ',Nint(passshl(5))
    Call info(message,.true.)
  End If

! bond constraints iterative cycles statistics

  If (megcon > 0) Then
     Call gmax(comm,stats%passcon(3:5,1,1)) ; Call gmax(comm,stats%passcon(3:5,2,1))
     If (stats%passcon(3,1,1) > 0.0_wp) Then
       Write(message,'(2(a,f5.2),4(a,i3))') &
         'constraints shake  run statistics - cycles per call/timestep: average ', &
         stats%passcon(3,1,1),' / ',stats%passcon(3,2,1), &
         ' minimum ',Nint(stats%passcon(4,1,1)),' / ',Nint(stats%passcon(4,2,1)), &
         ' maximum ',Nint(stats%passcon(5,1,1)),' / ',Nint(stats%passcon(5,2,1))
       Call info(message,.true.)
     End If

     Call gmax(comm,stats%passcon(3:5,1,2)) ; Call gmax(comm,stats%passcon(3:5,2,2))
     If (stats%passcon(3,1,2) > 0.0_wp) Then
       Write(message,'(2(a,f5.2),4(a,i3))') &
         'constraints rattle  run statistics - cycles per call/timestep: average ', &
         stats%passcon(3,1,1),' / ',stats%passcon(3,2,1), &
         ' minimum ',Nint(stats%passcon(4,1,2)),' / ',Nint(stats%passcon(4,2,2)), &
         ' maximum ',Nint(stats%passcon(5,1,2)),' / ',Nint(stats%passcon(5,2,2))
       Call info(message,.true.)
     End If
  End If

! PMF constraints iterative cycles statistics

  If (megpmf > 0) Then
     Call gmax(comm,stats%passpmf(3:5,1,1)) ; Call gmax(comm,stats%passpmf(3:5,2,1))
     If (stats%passpmf(3,1,1) > 0.0_wp) Then
       Write(message,'(2(a,f5.2),4(a,i3))') &
         'PMFs shake  run statistics - cycles per call/timestep: average ', &
         stats%passpmf(3,1,1),' / ',stats%passpmf(3,2,1), &
         ' minimum ',Nint(stats%passpmf(4,1,1)),' / ',Nint(stats%passpmf(4,2,1)), &
         ' maximum ',Nint(stats%passpmf(5,1,1)),' / ',Nint(stats%passpmf(5,2,1))
       Call info(message,.true.)
     EndIf

     Call gmax(comm,stats%passpmf(3:5,1,2)) ; Call gmax(comm,stats%passpmf(3:5,2,2))
     If (stats%passpmf(3,1,2) > 0.0_wp) Then
       Write(message,'(2(a,f5.2),4(a,i3))') &
         'PMFs rattle  run statistics - cycles per call/timestep: average ', &
         stats%passpmf(3,1,2),' / ',stats%passpmf(3,2,2), &
         ' minimum ',Nint(stats%passpmf(4,1,2)),' / ',Nint(stats%passpmf(4,2,2)), &
         ' maximum ',Nint(stats%passpmf(5,1,2)),' / ',Nint(stats%passpmf(5,2,2))
       Call info(message,.true.)
     End If
  End If

! Get elapsed time

  Call gtime(timelp)

! Get simulation time for averages

  If (stats%numacc == 0) Then
     tmp=0.0_wp
  Else
     tmp=time-tmst
  End If

! Report termination


  If ((nstep == 0 .and. nstrun == 0) .or. stats%numacc == 0) Then
    Write(message,'(a)') 'dry run terminated'
  Else
    Write(message,'(2(a,i9,a,f10.3),a)') 'run terminated after ',nstep, &
      ' steps (',time,' ps), final averages calculated over',stats%numacc, &
      ' steps (',tmp,' ps)'
  End If
  Call info(message,.true.)

! safe average volume and cell

  avvol = volm

! If dry/static/minimisation run - NO AVERAGES
! Print pressure tensor and jump to possible RDF and Z-Density

  If (nstep == 0 .and. nstrun == 0) Then
     iadd = 27+2*Merge(mxatdm,0,lmsd)+site%ntype_atom

     If (comm%idnode == 0) Then
        Write(message,'(a)') 'pressure tensor  (katms):'
        Call info(message,.true.)

        Do i=iadd,iadd+6,3
           Write(message,'(2x,1p,3e12.4)') stats%stpval(i+1:i+3)
           Call info(message,.true.)
        End Do

        Write(message,'(2x,a,1p,e12.4)') 'trace/3  ', (stats%stpval(iadd+1)+stats%stpval(iadd+5)+stats%stpval(iadd+9))/3.0_wp
        Call info(message,.true.)
     End If

     Go To 10
  End If

  ! If still running in the pure equilibration regime - NO AVERAGES
  If (stats%numacc /= 0) Then
    ! shift back statistical averages as from statistics_collect

    Do i=0,stats%mxnstk
      stats%sumval(i)=stats%sumval(i)+stats%stpvl0(i)
    End Do

    ! calculate final fluctuations

    Do i=0,stats%mxnstk
      stats%ssqval(i)=Sqrt(stats%ssqval(i))
    End Do

    ! average volume

    avvol =stats%sumval(19)

    ! final averages and fluctuations
    Write(messages(1),'(a)') Repeat('-',130)
    Write(messages(2),'(9x,a4,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
      'step','eng_tot','temp_tot','eng_cfg','eng_src','eng_cou','eng_bnd','eng_ang','eng_dih','eng_tet'
    Write(messages(3),'(5x,a8,5x,a7,4x,a8,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7,5x,a7)') &
      'time(ps)',' eng_pv','temp_rot','vir_cfg','vir_src','vir_cou','vir_bnd','vir_ang','vir_con','vir_tet'
    Write(messages(4), '(5x,a8,5x,a6,4x,a8,5x,a7,5x,a7,7x,a5,8x,a4,7x,a5,5x,a7,7x,a5)') &
      'cpu  (s)','volume','temp_shl','eng_shl','vir_shl','alpha','beta','gamma','vir_pmf','press'
    Write(messages(5),'(a)') Repeat('-',130)
    Call info(messages,5,.true.)

    Write(messages(1),'(i13,1p,9e12.4)')stats%numacc,stats%sumval(1:9)
    Write(messages(2),'(f13.5,1p,9e12.4)')tmp,stats%sumval(10:18)
    Write(messages(3),'(0p,f13.3,1p,9e12.4)') timelp,stats%sumval(19:27)
    Write(messages(4),'(a)')''
    Call info(messages,4,.true.)

    Write(messages(1),'(6x,a8,1p,9e12.4)') ' r.m.s. ',stats%ssqval(1:9)
    Write(messages(2),'(6x,a8,1p,9e12.4)') 'fluctu- ',stats%ssqval(10:18)
    Write(messages(3),'(6x,a8,1p,9e12.4)') 'ations  ',stats%ssqval(19:27)
    Write(messages(4),'(a)') Repeat('-',130)
    Call info(messages,4,.true.)

    ! Some extra information - conserved quantity=extended ensemble energy

    Write(message,"(a,1p,e12.4,5x,a,1p,e12.4)") &
      "Extended energy:       ",stats%sumval(0),     &
      " r.m.s. fluctuations:  ",stats%ssqval(0)
    Call info(message,.true.)

    ! Some extra information - <P*V> term - only matters for NP/sT ensembles

    If (thermo%variable_cell) Then
      Write(message,"(a,1p,e12.4,5x,a,1p,e12.4)")           &
        "<P*V> term:            ",stats%sumval(37+site%ntype_atom+2*Merge(mxatdm,0,lmsd)), &
        " r.m.s. fluctuations:  ",stats%ssqval(37+site%ntype_atom+2*Merge(mxatdm,0,lmsd))
      Call info(message,.true.)
    End If

    Write(messages(1),"(130('-'))")
    Write(messages(2),'(a)')''
    Call info(messages,2,.true.)

    ! Move at the end of the default 27 quantities

    iadd = 27

    If (lmsd) iadd = iadd+2*mxatdm

    ! Write out estimated diffusion coefficients

    Write(messages(1),'(a)') 'Approximate 3D Diffusion Coefficients and square root of MSDs:'
    Write(messages(2),'(6x,a4,2x,a19,6x,a15)') 'atom','DC (10^-9 m^2 s^-1)','Sqrt[MSD] (Ang)'
    Call info(messages,2,.true.)

    Do i=1,site%ntype_atom
      If (site%num_type_nf(i) > zero_plus) Then
        dc = 10.0_wp * (stats%ravval(iadd+i)-stats%sumval(iadd+i)) / &
          (3.0_wp*Real(stats%numacc-Min(stats%mxstak,stats%numacc-1),wp)*tstep)
        If (dc < 1.0e-10_wp) dc = 0.0_wp

        srmsd = Sqrt(stats%ravval(iadd+i))
        Write(message,'(2x,a8,1p,2(8x,e13.4))') site%unique_atom(i),dc,srmsd
      Else
        Write(message,'(2x,a8,1p,2(8x,e13.4))') site%unique_atom(i),0.0_wp,0.0_wp
      End If
      Call info(message,.true.)
    End Do
    Call info('',.true.)

    iadd = iadd+site%ntype_atom

    ! print out average pressure tensor

    If (comm%idnode == 0) Then
      Write(messages(1),'(a)') 'Pressure tensor:'
      Write(messages(2),'(6x,a32,5x,17x,a19)') 'Average pressure tensor  (katms)','r.m.s. fluctuations'
      Call info(messages,2,.true.)

      Do i=iadd,iadd+6,3
        Write(message,'(2x,1p,3e12.4,5x,3e12.4)')stats%sumval(i+1:i+3),stats%ssqval(i+1:i+3)
        Call info(message,.true.)
      End Do

      Write(message,'(2x,a,1p,e12.4)') 'trace/3  ', (stats%sumval(iadd+1)+stats%sumval(iadd+5)+stats%sumval(iadd+9))/3.0_wp
      Call info(message,.true.)
      Call info('',.true.)
    End If

    iadd = iadd+9

    ! Write out mean cell vectors for npt/nst

    If (thermo%variable_cell) Then

      If (comm%idnode == 0) Then
        Write(message,'(a32,33x,a19)') 'Average cell vectors     (Angs) ','r.m.s. fluctuations'
        Call info(message,.true.)


        Do i=iadd,iadd+6,3
          Write(message,'(3f20.10,5x,1p,3e12.4)')stats%sumval(i+1:i+3),stats%ssqval(i+1:i+3)
          Call info(message,.true.)
        End Do
      End If

      iadd = iadd+9

      ! PV term used above

      iadd = iadd+1

      If (thermo%iso > 0) Then
        h_z=stats%sumval(iadd+1)

        Write(message,"('Average surface area, fluctuations & mean estimate (Angs^2)')")
        Call info(message,.true.)
        Write(message,'(1p,3e12.4)')stats%sumval(iadd+2),stats%ssqval(iadd+2),avvol/h_z
        Call info(message,.true.)

        iadd = iadd+2

        If (thermo%iso > 1) Then
          tx= -h_z * (stats%sumval(iadd-9-8-2)/prsunt - (thermo%press+thermo%stress(1)) ) * tenunt
          ty= -h_z * (stats%sumval(iadd-9-7-2)/prsunt - (thermo%press+thermo%stress(5)) ) * tenunt
          Write(message,"('Average surface tension, fluctuations & mean estimate in x (dyn/cm)')")
          Call info(message,.true.)
          Write(message,'(1p,3e12.4)')stats%sumval(iadd+1),stats%ssqval(iadd+1),tx
          Call info(message,.true.)
          Write(message,"('Average surface tension, fluctuations & mean estimate in y (dyn/cm)')")
          Call info(message,.true.)
          Write(message,'(1p,3e12.4)')stats%sumval(iadd+2),stats%ssqval(iadd+2),ty
          Call info(message,.true.)

          iadd = iadd+2
        End If
      End If

    End If

    ! Write out remaining registers

    check = .false.
    Do i=iadd+1,stats%mxnstk
      If (Abs(stats%sumval(i)) > zero_plus .or. Abs(stats%ssqval(i)) > zero_plus) check=.true.
    End Do

    If (check) Then
      Write(messages(1),"('Remaining non-zero statistics registers:')")
      Write(messages(2),"(4x,'Register',7x,'Average value',8x,'r.m.s. fluc.')")
      Call info(messages,2,.true.)
    End If

    If (comm%idnode == 0) Then
      Do i=iadd+1,mxnstk
        If (Abs(stats%sumval(i)) > zero_plus .or. Abs(stats%ssqval(i)) > zero_plus) Then
          Write(message,'(2x,i10,2f20.10)') i,stats%sumval(i),stats%ssqval(i)
          Call info(message,.true.)
        End If
      End Do
    End If

  End If
10 Continue

! print final time check

  Call gtime(timelp)

  Write(message,'("time elapsed since job start: ", f12.3, " sec")') timelp
  Call info(message,.true.)

End Subroutine statistics_result
End Module statistics
