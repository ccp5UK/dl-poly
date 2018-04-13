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

  Use kinds, Only : wp,li
  Use setup, Only : mxatdm,mxnstk,mxstak,mxbfss,zero_plus,&
                   prsunt,nstats,tenunt,boltz,engunit,eu_ev,eu_kcpm,&
                   eu_kjpm,mxatyp,statis,pi

  Use comms,   Only : comms_type,gsum,Spread_tag,wp_mpi,gtime,gmax,gsend, &
                      gwait,girecv
  Use site,    Only : ntpatm,numtypnf,unqatm,dens
  Use configuration,  Only : cfgname,imcon,cell,volm,natms,ltype, &
                             xxx,yyy,zzz,vxx,vyy,vzz,ixyz,lsa,lsi,ltg
  Use domains,    Only : nprx,npry,nprz,map,r_nprx,r_npry,r_nprz,&
                                nprx_r,npry_r,nprz_r,idx,idy,idz
  Use msd,    Only : l_msd

  Use vnl
  Use core_shell,  Only : passshl
  Use constraints, Only : passcon
  Use pmf,         Only : passpmf
  Use bonds,       Only : ncfbnd,mxgbnd1
  Use angles,      Only : ncfang,mxgang1
  Use dihedrals,   Only : ncfdih,mxgdih1
  Use inversions,  Only : ncfinv,mxginv1
  Use rdfs,         Only : ncfrdf,l_errors_jack,l_errors_block,ncfusr, &
                           calculate_errors,calculate_errors_jackknife, &
                           rdf_compute,usr_compute
  Use z_density,   Only : ncfzdn,z_density_compute,z_density_collect
  Use msd
  Use greenkubo,   Only : vafsamp,vafcount,vaf_compute
  Use bonds,       Only : bonds_compute
  Use angles,      Only : angles_compute
  Use dihedrals,   Only : dihedrals_compute
  Use inversions,  Only : inversions_compute
  Use errors_warnings, Only : error,warning,info
  Use numerics,    Only : dcell,invert,shellsort,shellsort2,pbcshfrc,pbcshfrl

  Implicit None

  Integer,                        Save :: numacc = 0 , &
                                          natms0 = 0

  Real( Kind = wp ),              Save :: clin(1:9) = 0.0_wp

  Real( Kind = wp ), Allocatable, Save :: xin(:),yin(:),zin(:)
  Real( Kind = wp ), Allocatable, Save :: xto(:),yto(:),zto(:),rsd(:)

  Real( Kind = wp ), Allocatable, Save :: stpval(:),stpvl0(:),sumval(:),ssqval(:)
  Real( Kind = wp ), Allocatable, Save :: zumval(:),ravval(:),stkval(:,:)

  Integer,           Allocatable, Save :: found(:), found0(:)
  Integer,           Allocatable, Save :: lsi0(:),lsa0(:),ltg0(:)

  Real( Kind = wp ), Allocatable, Save :: xin0(:),yin0(:),zin0(:)
  Real( Kind = wp ), Allocatable, Save :: xto0(:),yto0(:),zto0(:)

  Real( Kind = wp ), Allocatable, Save :: stpval0(:),stpvl00(:),sumval0(:),ssqval0(:)
  Real( Kind = wp ), Allocatable, Save :: zumval0(:),ravval0(:),stkval0(:,:)

  Logical,                        Save :: statis_file_open = .false.

  Public :: allocate_statistics_arrays, allocate_statistics_connect, &
            deallocate_statistics_connect

Contains

  Subroutine allocate_statistics_arrays()

    Integer, Dimension( 1:4 ) :: fail

    fail = 0

    Allocate (xin(1:mxatdm),yin(1:mxatdm),zin(1:mxatdm),                           Stat = fail(1))
    Allocate (xto(1:mxatdm),yto(1:mxatdm),zto(1:mxatdm),rsd(1:mxatdm),             Stat = fail(2))
    Allocate (stpval(0:mxnstk),stpvl0(0:mxnstk),sumval(0:mxnstk),ssqval(0:mxnstk), Stat = fail(3))
    Allocate (zumval(0:mxnstk),ravval(0:mxnstk),stkval(1:mxstak,0:mxnstk),         Stat = fail(4))

    If (Any(fail > 0)) Call error(1016)

    xin = 0.0_wp ; yin = 0.0_wp ; zin = 0.0_wp
    xto = 0.0_wp ; yto = 0.0_wp ; zto = 0.0_wp ; rsd = 0.0_wp

    stpval = 0.0_wp ; stpvl0 = 0.0_wp ; sumval = 0.0_wp ; ssqval = 0.0_wp
    zumval = 0.0_wp ; ravval = 0.0_wp ; stkval = 0.0_wp

  End Subroutine allocate_statistics_arrays

  Subroutine allocate_statistics_connect()

    

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (found(1:mxatdm),found0(1:mxatdm),                                                Stat = fail(1))
    Allocate (lsi0(1:mxatdm),lsa0(1:mxatdm),ltg0(1:mxatdm),                                    Stat = fail(2))
    Allocate (xin0(1:mxatdm),yin0(1:mxatdm),zin0(1:mxatdm),                                    Stat = fail(3))
    Allocate (xto0(1:mxatdm),yto0(1:mxatdm),zto0(1:mxatdm),                                    Stat = fail(4))
    Allocate (stpval0(1:2*mxatdm),stpvl00(1:2*mxatdm),sumval0(1:2*mxatdm),ssqval0(1:2*mxatdm), Stat = fail(5))
    Allocate (zumval0(1:2*mxatdm),ravval0(1:2*mxatdm),stkval0(1:mxstak,1:2*mxatdm),            Stat = fail(6))

    If (Any(fail > 0)) Call error(1060)

  End Subroutine allocate_statistics_connect

  Subroutine deallocate_statistics_connect()

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Deallocate (found,found0,                    Stat = fail(1))
    Deallocate (lsi0,lsa0,ltg0,                  Stat = fail(2))
    Deallocate (xin0,yin0,zin0,                  Stat = fail(3))
    Deallocate (xto0,yto0,zto0,                  Stat = fail(4))
    Deallocate (stpval0,stpvl00,sumval0,ssqval0, Stat = fail(5))
    Deallocate (zumval0,ravval0,stkval0,         Stat = fail(6))

    If (Any(fail > 0)) Call error(1061)

  End Subroutine deallocate_statistics_connect

  Subroutine statistics_collect             &
           (lsim,leql,nsteql,lzdn,nstzdn, &
           keyres,keyens,iso,intsta,      &
           degfre,degshl,degrot,          &
           nstep,tstep,time,tmst,         &
           engcpe,vircpe,engsrp,virsrp,   &
           engter,virter,                 &
           engtbp,virtbp,engfbp,virfbp,   &
           engshl,virshl,shlke,           &
           vircon,virpmf,                 &
           engtet,virtet,engfld,virfld,   &
           engbnd,virbnd,engang,virang,   &
           engdih,virdih,enginv,virinv,   &
           engke,engrot,consv,vircom,     &
           strtot,press,strext,           &
           stpeng,stpvir,stpcfg,stpeth,   &
           stptmp,stpprs,stpvol,comm,virdpd)

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

  Logical,           Intent( In    ) :: lsim,leql,lzdn
  Integer,           Intent( In    ) :: nsteql,nstzdn,keyres, &
                                        keyens,iso,intsta,nstep

  Integer(Kind=li),  Intent( In    ) :: degfre,degshl,degrot

  Real( Kind = wp ), Intent( In    ) :: tstep,time,                  &
                                        engcpe,vircpe,engsrp,virsrp, &
                                        engter,virter,               &
                                        engtbp,virtbp,engfbp,virfbp, &
                                        engshl,virshl,shlke,         &
                                        vircon,virpmf,               &
                                        engtet,virtet,engfld,virfld, &
                                        engbnd,virbnd,engang,virang, &
                                        engdih,virdih,enginv,virinv, &
                                        engke,engrot,consv,vircom,   &
                                        strtot(1:9),press,strext(1:9)

  Real( Kind = wp ), Intent( InOut ) :: tmst
  Real( Kind = wp ), Intent(   Out ) :: stpeng,stpvir,stpcfg,stpeth, &
                                        stptmp,stpprs,stpvol
  Type( comms_type ), Intent( InOut ) :: comm
  Real( Kind = wp ), Intent( In ) :: virdpd

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
        Open(Unit=nstats, File=Trim(STATIS), Status='replace')
        statis_file_open = .true.

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

  stpcfg = engcpe + engsrp + engter + engtbp + engfbp + &
           engfld + engshl +                            &
           engtet + engbnd + engang + engdih + enginv

! system energy

  stpeng = stpcfg + engke + engrot

! energy + conserved quantity (for true ensembles)

  stpcns = stpeng + consv

! rotational temperature

  stprot = 2.0_wp*(engrot) / (boltz*Max(1.0_wp,Real(degrot,wp)))

! core-shell units temperature

  stpshl = 2.0_wp*(shlke) / (boltz*Max(1.0_wp,Real(degshl,wp)))

! system temperature

  stptmp = 2.0_wp*(engke+engrot) / (boltz*Real(degfre,wp))

! system virial
! Note: originally, purely angle dependent interactions have zero virial!!!
! So, virfbp, virinv and virdih are allegedly always zero!  virdih has an exception!

  stpvir = vircpe + virsrp + virter + virtbp + virfbp + &
           virfld + virshl + vircon + virpmf + vircom + &
           virtet + virbnd + virang + virdih + virinv + virdpd

! system volume

  stpvol = volm

! system pressure

  stpprs = (2.0_wp*engke-stpvir) / (3.0_wp*stpvol)

! system PV

  stpipv = stpprs*stpvol

! system enthalpy

  If (keyens >= 20) Then             ! P_target*V_instantaneous
     stpeth = stpeng + press*stpvol
  Else                               ! for keyens < 20 V_instantaneous=V_target
     stpeth = stpeng + stpipv        ! and there is only P_instantaneous
  End If

  Call dcell(cell,celprp)

! store current values in statistics array

!  If (idnode == 0) Write(nrite,'(5(/,5e12.4))')   &
!     engke  , engrot , stpprs , stpprs , stpvol , &
!     vircpe , virsrp , virter , virtbp , virfbp , &
!     virfld , virshl , vircon , virpmf , vircom , &
!     virtet , virbnd , virang , virdih , virinv , &
!     virdpd

  stpval(0) =consv/engunit
  stpval(1) =stpcns/engunit
  stpval(2) =stptmp
  stpval(3) =stpcfg/engunit
  stpval(4) =(engsrp+engter)/engunit
  stpval(5) =engcpe/engunit
  stpval(6) =engbnd/engunit
  stpval(7) =(engang+engtbp)/engunit
  stpval(8) =(engdih+enginv+engfbp)/engunit
  stpval(9) =engtet/engunit
  stpval(10)=stpeth/engunit
  stpval(11)=stprot
  stpval(12)=stpvir/engunit
  stpval(13)=(virsrp+virter)/engunit
  stpval(14)=vircpe/engunit
  stpval(15)=virbnd/engunit
  stpval(16)=(virtbp+virang)/engunit
  stpval(17)=vircon/engunit
  stpval(18)=virtet/engunit
  stpval(19)=stpvol
  stpval(20)=stpshl
  stpval(21)=engshl/engunit
  stpval(22)=virshl/engunit
  stpval(23)=Acos(celprp(6))*180.0_wp/pi
  stpval(24)=Acos(celprp(5))*180.0_wp/pi
  stpval(25)=Acos(celprp(4))*180.0_wp/pi
  stpval(26)=virpmf/engunit
  stpval(27)=stpprs*prsunt

  iadd = 27

! iadd = iadd + 1 ! for the stpval(0)!!! Thus to account for in printing

! mean squared displacements per species, dependent on
! particle displacements from initial positions (at t=0)

  amsd = 0.0_wp ! initialise

  If (nstep == nsteql+1) Then ! re-initialise
     Do i=1,natms
        xto(i) = 0.0_wp
        yto(i) = 0.0_wp
        zto(i) = 0.0_wp
     End Do
  End If

  If (nstep > 0) Then
     If (lsim)  Then ! real dynamics is happening
        Do i=1,natms
           xto(i)=xto(i)+vxx(i)*tstep
           yto(i)=yto(i)+vyy(i)*tstep
           zto(i)=zto(i)+vzz(i)*tstep
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
        Call pbcshfrc(imcon,clin,natms,xin,yin,zin)
        Do i=1,natms
           xin(i)=xxt(i)-xin(i)
           yin(i)=yyt(i)-yin(i)
           zin(i)=zzt(i)-zin(i)
        End Do
        Deallocate (xxt,yyt,zzt, Stat=fail)
        If (fail > 0) Then
           Write(message,'(a)') 'statistics_collect deallocation failure 1'
           Call error(0,message)
        End If
        Call pbcshfrl(imcon,cell,natms,xin,yin,zin)
        Do i=1,natms
           xto(i)=xto(i)+xin(i)
           yto(i)=yto(i)+yin(i)
           zto(i)=zto(i)+zin(i)
        End Do
     End If

     Do i=1,natms
        rsd(i)=Sqrt(xto(i)**2+yto(i)**2+zto(i)**2)

        k=ltype(i)
        amsd(k)=amsd(k)+rsd(i)**2
     End Do
     Call gsum(comm,amsd(1:ntpatm))
  End If

  If (l_msd) Then
     Do i=1,natms
        j=2*i
        stpval(iadd+j-1)=rsd(i)**2
        stpval(iadd+j  )=vxx(i)**2+vyy(i)**2+vzz(i)**2
     End Do
     iadd = iadd + 2*mxatdm
  End If

  Do k=1,ntpatm
     If (numtypnf(k) > zero_plus) Then
        stpval(iadd+k)=amsd(k)/numtypnf(k)
     Else
        stpval(iadd+k)=0.0_wp
     End If
  End Do
  iadd = iadd + ntpatm

! pressure tensor (derived for the stress tensor)

  Do i=1,9
     stpval(iadd+i)=strtot(i)*prsunt/stpvol
  End Do
  iadd = iadd + 9

  If (keyens >= 20) Then

! cell parameters

     Do i=1,9
       stpval(iadd+i)=cell(i)
     End Do
     iadd = iadd + 9

! instantaneous PV

     stpval(iadd+1)=stpipv/engunit
     iadd = iadd + 1

     If (iso > 0) Then
        h_z=celprp(9)

        stpval(iadd+1)=h_z
        stpval(iadd+2)=stpvol/h_z
        iadd = iadd + 2

        If (iso > 1) Then
           stpval(iadd+1)= -h_z*(strtot(1)-(press+strext(1)))*tenunt
           stpval(iadd+2)= -h_z*(strtot(5)-(press+strext(5)))*tenunt
           iadd = iadd + 2
        End If
     End If
  End If

! write statistics file

  If (comm%idnode == 0 .and. Mod(nstep,intsta) == 0) Then
     If (.not. statis_file_open) Then
         Open(Unit=nstats, File=Trim(STATIS), Position='append')
         statis_file_open = .true.
     End If

     If (l_msd) Then
        Write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))') &
             nstep,time,iadd+1-2*mxatdm,stpval(1:  27),stpval(0),stpval(28+2*mxatdm:iadd)
     Else
        Write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))') &
             nstep,time,iadd+1,         stpval(1:  27),stpval(0),stpval(28         :iadd)
     End If

  End If

! check on number of variables for stack

  If (iadd > mxnstk) Call error(170)

! No totals for timestep zero

  If (nstep == 0) Go To 10

! current stack value

  kstak=Mod(nstep-1,mxstak)+1

! subtract old stack value from the stack average

  If (nstep > mxstak) Then
     Do i=0,mxnstk
        zumval(i)=zumval(i)-stkval(kstak,i)
     End Do
  End If

! store quantities in stack and update the stack average

  Do i=0,mxnstk
     stkval(kstak,i)=stpval(i)
     zumval(i)=zumval(i)+stpval(i)
  End Do

! calculate rolling averages

  zistk=Real(Min(mxstak,nstep),wp)

  Do i=0,mxnstk
     ravval(i)=zumval(i)/zistk
  End Do

! accumulate totals over steps

  If ((.not.leql) .or. nstep > nsteql) Then
     numacc=numacc+1
     sclnv2=1.0_wp/Real(numacc,wp)
     sclnv1=Real(numacc-1,wp)/Real(numacc,wp)

! average squared sum and sum (keep in this order!!!)

     If (nstep == nsteql+1 .or. ((.not.leql) .and. nstep == 1)) stpvl0=stpval
     stpval=stpval-stpvl0
     Do i=0,mxnstk
        ssqval(i)=sclnv1*(ssqval(i)+sclnv2*(stpval(i)-sumval(i))**2)

! sumval has to be shifted back to sumval+stpvl0 in statistics_result
! when averaging is printed since stpval is only shifted back and forth
! which does not affect the fluctuations Sqrt(ssqval) only their accuracy

        sumval(i)=sclnv1*sumval(i)+sclnv2*stpval(i)
     End Do
     stpval=stpval+stpvl0
  End If

10 Continue

! z-density collection

  If ( lzdn .and. ((.not.leql) .or. nstep >= nsteql) .and. &
       Mod(nstep,nstzdn) == 0 ) Call z_density_collect()

! Catch time of starting statistical averages

  If (((.not.leql) .or. nstep == nsteql) .and. tmst < tstep) tmst=time

  Deallocate (amsd, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'statistics_collect deallocation failure'
     Call error(0,message)
  End If

End Subroutine statistics_collect


Subroutine statistics_connect_frames(megatm,comm)

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


  Integer, Intent ( In    ) :: megatm
  Type( comms_type ), Intent( InOut ) :: comm

  Integer :: icyc,nres
  Character( Len = 256 ) :: message

  found = 0 ; icyc = 0 ; nres = 1
  Do While (icyc <= Max(nprx,npry,nprz)/2 .and. nres > 0)
     Call match_compress_spread_sort(-1) ! -x direction spread
     Call match_compress_spread_sort( 1) ! +x direction spread

     Call match_compress_spread_sort(-2) ! -y direction spread
     Call match_compress_spread_sort( 2) ! +y direction spread

     Call match_compress_spread_sort(-3) ! -z direction spread
     Call match_compress_spread_sort( 3) ! +z direction spread

     Call match_compress_spread_sort( 0) ! no spreading

     nres=natms0
     Call gsum(comm,nres)
     If (nres > 0) Then
        nres=Merge(0,Sum(found(1:natms)),natms > 0)
        Call gsum(comm,nres)
        If (nres /= megatm) icyc = icyc + 1
     End If
  End Do

  If (nres > 0) Then
    Write(message,'(a)') ' particles dynamics properties will be corrupted'
    Call warning(message,.true.)
  End If

Contains

  Subroutine match_compress_spread_sort(mdir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to simplify the repetition of the procedures above
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Integer, Intent( In    ) :: mdir ! +/-1,+/-2,+/-3,0 is the direction of spread

    Integer :: fail,i,i0,j,j0,kk

    Integer, Allocatable, Save :: lsa00(:)

! Search for matches

    found0 = 0
    If (natms0 > 0) Then
       i0 = 1
       Do i=1,natms
          Do While (i0 <= natms0)
             If      (lsa(i) <  lsa0(i0)) Then
                Exit
             Else If (lsa(i) == lsa0(i0)) Then
                If (found(lsi(i)) > 0) Then ! ghost arrival
                   found0(lsi0(i0)) = 1     ! erase at compression
                Else                        ! new arrival to claim
                   found(lsi(i)) = 1 ; found0(lsi0(i0)) = 1

                   xin(lsi(i)) = xin0(lsi0(i0))
                   yin(lsi(i)) = yin0(lsi0(i0))
                   zin(lsi(i)) = zin0(lsi0(i0))

                   xto(lsi(i)) = xto0(lsi0(i0))
                   yto(lsi(i)) = yto0(lsi0(i0))
                   zto(lsi(i)) = zto0(lsi0(i0))

                   If (l_msd) Then
                      j =27+2*lsi(i)
                      j0=2*lsi0(i0)
                      stpvl0(j-1)=stpvl00(j0-1)
                      stpvl0(j  )=stpvl00(j0  )
                      stpval(j-1)=stpval0(j0-1)
                      stpval(j  )=stpval0(j0  )
                      zumval(j-1)=zumval0(j0-1)
                      zumval(j  )=zumval0(j0  )
                      ravval(j-1)=ravval0(j0-1)
                      ravval(j  )=ravval0(j0  )
                      ssqval(j-1)=ssqval0(j0-1)
                      ssqval(j  )=ssqval0(j0  )
                      sumval(j-1)=sumval0(j0-1)
                      sumval(j  )=sumval0(j0  )
                      Do kk=1,mxstak
                         stkval(kk,j-1)=stkval0(kk,j0-1)
                         stkval(kk,j  )=stkval0(kk,j0  )
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
       Do i0=1,natms0
          Do While (lsa00(i) /= 0 .and. i < mxatdm)
             If      (lsa0(i0) <  lsa00(i)) Then
                Exit
             Else If (lsa0(i0) == lsa00(i)) Then
                found0(lsi0(i0)) = 1        ! erase at compression
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
    Do While (i0 <= natms0 .and. natms0 > 0)
       If (found0(i0) == 0) Then          ! Not claimed
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
          If (found0(natms0) == 0) Then   ! try to refill with the last unclaimed entry
             ixyz(i0) = ixyz(natms0)
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
             ltg0(i0) = ltg0(natms0)

             xin0(i0) = xin0(natms0)
             yin0(i0) = yin0(natms0)
             zin0(i0) = zin0(natms0)

             xto0(i0) = xto0(natms0)
             yto0(i0) = yto0(natms0)
             zto0(i0) = zto0(natms0)

             If (l_msd) Then
                j =2*i0
                j0=2*natms0
                stpvl00(j-1)=stpvl00(j0-1)
                stpvl00(j  )=stpvl00(j0  )
                stpval0(j-1)=stpval0(j0-1)
                stpval0(j  )=stpval0(j0  )
                zumval0(j-1)=zumval0(j0-1)
                zumval0(j  )=zumval0(j0  )
                ravval0(j-1)=ravval0(j0-1)
                ravval0(j  )=ravval0(j0  )
                ssqval0(j-1)=ssqval0(j0-1)
                ssqval0(j  )=ssqval0(j0  )
                sumval0(j-1)=sumval0(j0-1)
                sumval0(j  )=sumval0(j0  )
                Do kk=1,mxstak
                   stkval0(kk,j-1)=stkval0(kk,j0-1)
                   stkval0(kk,j  )=stkval0(kk,j0  )
                End Do
             End If

             i0 = i0 + 1                  ! increase lower bound marker if entry is refilled
          End If

!! Erase upper holdings in either case
!
!          ixyz(natms0) = 0
!          ltg0(natms0) = 0
!
!          xin0(natms0) = 0
!          yin0(natms0) = 0
!          zin0(natms0) = 0
!
!          xto0(natms0) = 0
!          yto0(natms0) = 0
!          zto0(natms0) = 0
!
!          If (l_msd) Then
!             j0=2*natms0
!             stpvl00(j0-1)=0.0_wp
!             stpvl00(j0  )=0.0_wp
!             stpval0(j0-1)=0.0_wp
!             stpval0(j0  )=0.0_wp
!             zumval0(j0-1)=0.0_wp
!             zumval0(j0  )=0.0_wp
!             ravval0(j0-1)=0.0_wp
!             ravval0(j0  )=0.0_wp
!             ssqval0(j0-1)=0.0_wp
!             ssqval0(j0  )=0.0_wp
!             sumval0(j0-1)=0.0_wp
!             sumval0(j0  )=0.0_wp
!             Do kk=1,mxstak
!                stkval0(kk,j0-1)=0.0_wp
!                stkval0(kk,j0  )=0.0_wp
!             End Do
!          End If
          natms0 = natms0 - 1             ! Decrease upper bound marker
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
       Do i0=1,natms0
          If (ixyz(i0) == 333) Then
             i=i+1
             lsa00(i)=ltg0(i0)
          End If
       End Do
       Call shellsort(i,lsa00)
    End If

! Spread atom data in the mdir direction

    If (mdir /= 0) Call statistics_connect_spread(mdir,comm)

! Sort past frame remainder of global atom indices

!    lsi0=0 ; lsa0=0
    Do i0=1,natms0
       lsi0(i0)=i0
       lsa0(i0)=ltg0(i0)
    End Do
    Call shellsort2(natms0,lsi0,lsa0)

  End Subroutine match_compress_spread_sort

End Subroutine statistics_connect_frames

  Subroutine statistics_connect_set(rcut,comm)

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

     natms0 = natms
     ltg0(1:natms0) = ltg(1:natms0) !; ltg0(natms0+1: ) = 0
     lsa0(1:natms0) = lsa(1:natms0) !; lsa0(natms0+1: ) = 0
     lsi0(1:natms0) = lsi(1:natms0) !; lsi0(natms0+1: ) = 0

     xin0(1:natms0) = xin(1:natms0) !; xin0(natms0+1: ) = 0 ; xin = 0.0_wp
     yin0(1:natms0) = yin(1:natms0) !; yin0(natms0+1: ) = 0 ; yin = 0.0_wp
     zin0(1:natms0) = zin(1:natms0) !; zin0(natms0+1: ) = 0 ; zin = 0.0_wp

     xto0(1:natms0) = xto(1:natms0) !; xto0(natms0+1: ) = 0
     yto0(1:natms0) = yto(1:natms0) !; yto0(natms0+1: ) = 0
     zto0(1:natms0) = zto(1:natms0) !; zto0(natms0+1: ) = 0

     If (l_msd) Then
        i0=2*natms0
        stpvl00(1:i0)=stpvl0(28:27+i0) !; stpvl00(i0+1: )=0.0_wp
        stpval0(1:i0)=stpval(28:27+i0) !; stpval0(i0+1: )=0.0_wp
        zumval0(1:i0)=zumval(28:27+i0) !; zumval0(i0+1: )=0.0_wp
        ravval0(1:i0)=ravval(28:27+i0) !; ravval0(i0+1: )=0.0_wp
        ssqval0(1:i0)=ssqval(28:27+i0) !; ssqval0(i0+1: )=0.0_wp
        sumval0(1:i0)=sumval(28:27+i0) !; sumval0(i0+1: )=0.0_wp
        Do kk=1,mxstak
           stkval0(kk,1:i0)=stkval(kk,1:i0) !; stkval0(kk,i0+1: )=0.0_wp
        End Do
     End If
  End If

End Subroutine statistics_connect_set

Subroutine statistics_connect_spread(mdir,comm)

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



  Integer, Intent( In    ) :: mdir
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

  Do i=1,natms0

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

        ltg0(keep)=ltg0(i)
        ixyz(keep)=ixyz(i)

! retain initial positions

        xin0(keep)=xin0(i)
        yin0(keep)=yin0(i)
        zin0(keep)=zin0(i)

! retain final displacements

        xto0(keep)=xto0(i)
        yto0(keep)=yto0(i)
        zto0(keep)=zto0(i)

        If (l_msd) Then
           jj=2*i
           j =2*keep
           stpvl00(j-1)=stpvl00(jj-1)
           stpvl00(j  )=stpvl00(jj  )
           stpval0(j-1)=stpval0(jj-1)
           stpval0(j  )=stpval0(jj  )
           zumval0(j-1)=zumval0(jj-1)
           zumval0(j  )=zumval0(jj  )
           ravval0(j-1)=ravval0(jj-1)
           ravval0(j  )=ravval0(jj  )
           ssqval0(j-1)=ssqval0(jj-1)
           ssqval0(j  )=ssqval0(jj  )
           sumval0(j-1)=sumval0(jj-1)
           sumval0(j  )=sumval0(jj  )
           Do kk=1,mxstak
              stkval0(kk,j-1)=stkval0(kk,jj-1)
              stkval0(kk,j  )=stkval0(kk,jj  )
           End Do
        End If
     End If

     If (move) Then ! copy it
        send=send+1
        If (imove+8 <= iblock) Then ! If safe to proceed

! pack config indexing and move indexing arrays

           buffer(imove+1)=Real(ltg0(i),wp)
           buffer(imove+2)=Real(ixyz(i),wp)

! pack initial positions

           buffer(imove+3)=xin0(i)
           buffer(imove+4)=yin0(i)
           buffer(imove+5)=zin0(i)

! pack final displacements

           buffer(imove+6)=xto0(i)
           buffer(imove+7)=yto0(i)
           buffer(imove+8)=zto0(i)
        Else
           safe=.false.
        End If
        imove=imove+8

! pack MSD arrays

        If (l_msd) Then
           If (imove+2*(6+mxstak) <= iblock) Then
              jj=2*i
              buffer(imove+ 1)=stpvl00(jj-1)
              buffer(imove+ 2)=stpvl00(jj  )
              buffer(imove+ 3)=stpval0(jj-1)
              buffer(imove+ 4)=stpval0(jj  )
              buffer(imove+ 5)=zumval0(jj-1)
              buffer(imove+ 6)=zumval0(jj  )
              buffer(imove+ 7)=ravval0(jj-1)
              buffer(imove+ 8)=ravval0(jj  )
              buffer(imove+ 9)=ssqval0(jj-1)
              buffer(imove+10)=ssqval0(jj  )
              buffer(imove+11)=sumval0(jj-1)
              buffer(imove+12)=sumval0(jj  )
              Do kk=1,mxstak
                 l=12+2*kk
                 buffer(imove+l-1)=stkval0(kk,jj-1)
                 buffer(imove+l  )=stkval0(kk,jj  )
              End Do
           Else
              safe=.false.
           End If
           imove=imove+2*(6+mxstak)
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
     If (All(ltg0(1:natms0) /= l)) imove=imove+1
     kmove=kmove+8
     If (l_msd) kmove=kmove+2*(6+mxstak)
  End Do

  natms0=keep+imove

! Check for array bound overflow (can arrays cope with incoming data)

  safe=(natms0 <= mxatdm)
  Call gcheck(comm,safe)
  If (.not.safe) Call error(164)

! load transferred data

  kmove=iblock+1 ! restore kmove
  newatm=keep    ! restore newatm
  Do i=1,jmove
     If (imove /= jmove) Then
        l=Nint(buffer(kmove+1))
        If (Any(ltg0(1:keep) == l)) Then
           kmove=kmove+8
           If (l_msd) kmove=kmove+2*(6+mxstak)
           Cycle
        End If
     End If

     newatm=newatm+1

! unpack config indexing, site and move indexing arrays

     ltg0(newatm)=Nint(buffer(kmove+1))
     ixyz(newatm)=Nint(buffer(kmove+2))

! unpack initial positions arrays

     xin0(newatm)=buffer(kmove+3)
     yin0(newatm)=buffer(kmove+4)
     zin0(newatm)=buffer(kmove+5)

! unpack initial positions arrays

     xto0(newatm)=buffer(kmove+6)
     yto0(newatm)=buffer(kmove+7)
     zto0(newatm)=buffer(kmove+8)

     kmove=kmove+8

! unpack MSD arrays

     If (l_msd) Then
        jj=2*newatm
        stpvl00(jj-1)=buffer(kmove+1 )
        stpvl00(jj  )=buffer(kmove+2 )
        stpval0(jj-1)=buffer(kmove+3 )
        stpval0(jj  )=buffer(kmove+4 )
        zumval0(jj-1)=buffer(kmove+5 )
        zumval0(jj  )=buffer(kmove+6 )
        ravval0(jj-1)=buffer(kmove+7 )
        ravval0(jj  )=buffer(kmove+8 )
        ssqval0(jj-1)=buffer(kmove+9 )
        ssqval0(jj  )=buffer(kmove+10)
        sumval0(jj-1)=buffer(kmove+11)
        sumval0(jj  )=buffer(kmove+12)
        Do kk=1,mxstak
           l=12+2*kk
           stkval0(kk,jj-1)=buffer(kmove+l-1)
           stkval0(kk,jj  )=buffer(kmove+l  )
        End Do

        kmove=kmove+2*(6+mxstak)
     End If
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'statistics_connect_spread deallocation failure'
     Call error(0,message)
  End If

End Subroutine statistics_connect_spread

Subroutine statistics_result                                    &
           (rcut,lmin,lpana,lrdf,lprdf,lzdn,lpzdn,lvafav,lpvaf, &
           nstrun,keyens,keyshl,megcon,megpmf,iso,              &
           press,strext,nstep,tstep,time,tmst,comm,passmin)

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

  Logical,           Intent( In    ) :: lmin,lpana,lrdf,lprdf,lzdn,lpzdn,lvafav,lpvaf
  Integer,           Intent( In    ) :: nstrun,keyens,keyshl,megcon,megpmf,iso,nstep
  Real( Kind = wp ), Intent( In    ) :: rcut,press,strext(1:9),tstep,time,tmst
  Type( comms_type ), Intent( InOut ) :: comm
  Real( Kind = wp ), Intent( In    ) ::  passmin(:)
  Logical           :: check
  Integer           :: i,iadd
  Real( Kind = wp ) :: avvol,avcel(1:9),dc,srmsd,timelp,tmp,h_z,tx,ty,temp
  Character( Len = 256 ) :: message
  Character( Len = 256 ), Dimension(5) :: messages

! VNL skipping statistics

  If (llvnl .and. nstep > 0) Then
     If (.not.l_vnl) Then ! Include the final skip in skipping statistics
        skipvnl(3)=skipvnl(2)*skipvnl(3)
        skipvnl(2)=skipvnl(2)+1.0_wp
        skipvnl(3)=skipvnl(3)/skipvnl(2)+skipvnl(1)/skipvnl(2)
        skipvnl(4)=Min(skipvnl(1),skipvnl(4))
        skipvnl(5)=Max(skipvnl(1),skipvnl(5))
     End If

     Write(message,"(                                        &
       & ' VNL skipping run statistics - skips per timestep: average ', f7.2, &
       & ' minimum ', i4, ' maximum ', i4)")                                  &
       skipvnl(3),Nint(Merge(skipvnl(4),skipvnl(5),skipvnl(4)<skipvnl(5))),Nint(skipvnl(5))
     Call info(message,.true.)
  End If

! minimisation convergence statistics

  If (lmin) Then
    Write(message,"(                          &
      & ' minimisation run statistics - cycles per call: average ', f7.2, &
      & ' minimum ', i4, ' maximum ', i4)") passmin(3),Nint(passmin(4)),Nint(passmin(5))
    Call info(message,.true.)
  End If

! shell relaxation convergence statistics

  If (keyshl == 2) Then
    Write(message,"(                           &
      & ' shell relaxation run statistics - cycles per timestep: average ', f7.2, &
      & ' minimum ', i4, ' maximum ', i4)") passshl(3),Nint(passshl(4)),Nint(passshl(5))
    Call info(message,.true.)
  End If

! bond constraints iterative cycles statistics

  If (megcon > 0) Then
     Call gmax(comm,passcon(3:5,1,1)) ; Call gmax(comm,passcon(3:5,2,1))
     If (passcon(3,1,1) > 0.0_wp) Then
       Write(message,"(                                   &
         & ' constraints shake  run statistics - cycles per call/timestep: average ', f5.2, ' / ', f5.2, &
         & ' minimum ', i3, ' / ', i3, ' maximum ', i3, ' / ', i3)")                                     &
         passcon(3,1,1),passcon(3,2,1),Nint(passcon(4,1,1)),Nint(passcon(4,2,1)),Nint(passcon(5,1,1)),Nint(passcon(5,2,1))
       Call info(message,.true.)
     End If

     Call gmax(comm,passcon(3:5,1,2)) ; Call gmax(comm,passcon(3:5,2,2))
     If (passcon(3,1,2) > 0.0_wp) Then
       Write(message,"(                                      &
         & ' constraints rattle run statistics - cycles per call/timestep: average ', f5.2, ' / ', f5.2, &
         & ' minimum ', i3, ' / ', i3, ' maximum ', i3, ' / ', i3)")                                     &
         passcon(3,1,2),passcon(3,2,2),Nint(passcon(4,1,2)),Nint(passcon(4,2,2)),Nint(passcon(5,1,2)),Nint(passcon(5,2,2))
       Call info(message,.true.)
     End If
  End If

! PMF constraints iterative cycles statistics

  If (megpmf > 0) Then
     Call gmax(comm,passpmf(3:5,1,1)) ; Call gmax(comm,passpmf(3:5,2,1))
     If (passpmf(3,1,1) > 0.0_wp) Then
       Write(message,"(                            &
         & ' PMFs shake  run statistics - cycles per call/timestep: average ', f5.2, ' / ', f5.2, &
         & ' minimum ', i3, ' / ', i3, ' maximum ', i3, ' / ', i3)")                              &
         passpmf(3,1,1),passpmf(3,2,1),Nint(passpmf(4,1,1)),Nint(passpmf(4,2,1)),Nint(passpmf(5,1,1)),Nint(passpmf(5,2,1))
       Call info(message,.true.)
     EndIf

     Call gmax(comm,passpmf(3:5,1,2)) ; Call gmax(comm,passpmf(3:5,2,2))
     If (passpmf(3,1,2) > 0.0_wp) Then
       Write(message,"(                               &
         & ' PMFs rattle run statistics - cycles per call/timestep: average ', f5.2, ' / ', f5.2, &
         & ' minimum ', i3, ' / ', i3, ' maximum ', i3, ' / ', i3)")                              &
         passpmf(3,1,2),passpmf(3,2,2),Nint(passpmf(4,1,2)),Nint(passpmf(4,2,2)),Nint(passpmf(5,1,2)),Nint(passpmf(5,2,2))
       Call info(message,.true.)
     End If
  End If

! Get elapsed time

  Call gtime(timelp)

! Get simulation time for averages

  If (numacc == 0) Then
     tmp=0.0_wp
  Else
     tmp=time-tmst
  End If

! Report termination


  If ((nstep == 0 .and. nstrun == 0) .or. numacc == 0) Then
    Write(message,"('dry run terminated')")
  Else
    Write(message,"('run terminated after',i9,' steps (',f10.3,   &
      & ' ps), final averages calculated over',i9,' steps (',f10.3, &
      & ' ps).')") nstep,time,numacc,tmp
  End If
  Call info(message,.true.)

! safe average volume and cell

  avvol = volm
  avcel = cell

! If dry/static/minimisation run - NO AVERAGES
! Print pressure tensor and jump to possible RDF and Z-Density

  If (nstep == 0 .and. nstrun == 0) Then
     iadd = 27+2*Merge(mxatdm,0,l_msd)+ntpatm

     If (comm%idnode == 0) Then
        Write(message,"(16x,'pressure tensor  (katms)')")
        Call info(message,.true.)

        Do i=iadd,iadd+6,3
           Write(message,'(9x,1p,3e12.4)') stpval(i+1:i+3)
           Call info(message,.true.)
        End Do

        Write(message,'(12x,a,1p,e12.4)') 'trace/3  ', (stpval(iadd+1)+stpval(iadd+5)+stpval(iadd+9))/3.0_wp
        Call info(message,.true.)
     End If

     Go To 10
  End If

! If still running in the pure equilibration regime - NO AVERAGES

  If (numacc == 0) Go To 20

! shift back statistical averages as from statistics_collect

  Do i=0,mxnstk
     sumval(i)=sumval(i)+stpvl0(i)
  End Do

! calculate final fluctuations

  Do i=0,mxnstk
     ssqval(i)=Sqrt(ssqval(i))
  End Do

! average volume

  avvol = sumval(19)

! final averages and fluctuations
  Write(messages(1),"(1x,130('-'))")
  Write(messages(2),"(10x,'step',5x,'eng_tot',4x,'temp_tot',5x,'eng_cfg',     &
    & 5x,'eng_src',5x,'eng_cou',5x,'eng_bnd',5x,'eng_ang',    &
    & 5x,'eng_dih',5x,'eng_tet')")
  Write(messages(3),"(6x,'time(ps)',5x,' eng_pv', &
    & 4x,'temp_rot',5x,'vir_cfg',5x,'vir_src',5x,'vir_cou',   &
    & 5x,'vir_bnd',5x,'vir_ang',5x,'vir_con',5x,'vir_tet',  &
    & 6x,'cpu  (s)',6x,'volume',4x,'temp_shl',5x,'eng_shl',   &
    & 5x,'vir_shl',7x,'alpha',8x,'beta',7x,'gamma',           &
    & 5x,'vir_pmf',7x,'press')")
  Write(messages(4),"(1x,130('-'))")
  Call info(messages,4,.true.)

  Write(messages(1),'(1x,i13,1p,9e12.4)')numacc,sumval(1:9)
  Write(messages(2),'(f14.5,1p,9e12.4)')tmp,sumval(10:18)
  Write(messages(3),'(1x,0p,f13.3,1p,9e12.4)') timelp,sumval(19:27)
  Call info(messages,3,.true.)

  Write(messages(1),"(6x,' r.m.s. ',1p,9e12.4)")ssqval(1:9)
  Write(messages(2),"(6x,'fluctu- ',1p,9e12.4)") ssqval(10:18)
  Write(messages(3),"(6x,'ations  ',1p,9e12.4)") ssqval(19:27)
  Call info(messages,3,.true.)

  Write(message,"(1x,130('-'))")
  Call info(message,.true.)

  ! Some extra information - conserved quantity=extended ensemble energy

  Write(message,"(a,1p,e12.4,5x,a,1p,e12.4)") &
    "Extended energy:       ", sumval(0),     &
    " r.m.s. fluctuations:  ", ssqval(0)
  Call info(message,.true.)

  ! Some extra information - <P*V> term - only matters for NP/sT ensembles

  If (keyens >= 20) Write(message,"(a,1p,e12.4,5x,a,1p,e12.4)")           &
    "<P*V> term:            ", sumval(37+ntpatm+2*Merge(mxatdm,0,l_msd)), &
    " r.m.s. fluctuations:  ", ssqval(37+ntpatm+2*Merge(mxatdm,0,l_msd))
  Call info(message,.true.)

  Write(message,"(130('-'))")
  Call info(message,.true.)

! Move at the end of the default 27 quantities

  iadd = 27

  If (l_msd) iadd = iadd+2*mxatdm

! Write out estimated diffusion coefficients

  Write(messages(1),"(12x,a)") 'Approximate 3D Diffusion Coefficients and square root of MSDs'
  Write(messages(2),"(12x,'atom',9x,'DC (10^-9 m^2 s^-1)',3x,'Sqrt[MSD] (Ang)')")
  Call info(messages,2,.true.)

  Do i=1,ntpatm
     If (numtypnf(i) > zero_plus) Then
        dc = 10.0_wp * (ravval(iadd+i)-sumval(iadd+i)) / &
             (3.0_wp*Real(numacc-Min(mxstak,numacc-1),wp)*tstep)
        If (dc < 1.0e-10_wp) dc = 0.0_wp

        srmsd = Sqrt(ravval(iadd+i))
        Write(message,'(12x,a8,1p,2(7x,e13.4))') unqatm(i),dc,srmsd
     Else
        Write(message,'(12x,a8,1p,2(7x,e13.4))') unqatm(i),0.0_wp,0.0_wp
     End If
     Call info(message,.true.)
  End Do

  iadd = iadd+ntpatm

! print out average pressure tensor

  If (comm%idnode == 0) Then
    Write(message,"(16x,'Average pressure tensor  (katms)',30x,'r.m.s. fluctuations')")
    Call info(message,.true.)

    Do i=iadd,iadd+6,3
      Write(message,'(9x,1p,3e12.4,24x,3e12.4)') sumval(i+1:i+3),ssqval(i+1:i+3)
      Call info(message,.true.)
    End Do

    Write(message,'(12x,a,1p,e12.4)') 'trace/3  ', (sumval(iadd+1)+sumval(iadd+5)+sumval(iadd+9))/3.0_wp
    Call info(message,.true.)
  End If

  iadd = iadd+9

! Write out mean cell vectors for npt/nst

  If (keyens >= 20) Then

! average cell (again)

     Do i=1,9
        avcel(i) = sumval(iadd+i)
     End Do

     If (comm%idnode == 0) Then
       Write(message,"(16x,'Average cell vectors     (Angs) ',30x,'r.m.s. fluctuations')")
       Call info(message,.true.)

       Do i=iadd,iadd+6,3
         Write(message,'(3f20.10,9x,1p,3e12.4)') sumval(i+1:i+3),ssqval(i+1:i+3)
         Call info(message,.true.)
       End Do
     End If

     iadd = iadd+9

! PV term used above

     iadd = iadd+1

     If (iso > 0) Then
        h_z=sumval(iadd+1)

        Write(message,"(16x,'Average surface area, fluctuations & mean estimate (Angs^2)')")
        Call info(message,.true.)
        Write(message,'(1p,3e12.4)') sumval(iadd+2),ssqval(iadd+2),avvol/h_z
        Call info(message,.true.)

        iadd = iadd+2

        If (iso > 1) Then
           tx= -h_z * ( sumval(iadd-9-8-2)/prsunt - (press+strext(1)) ) * tenunt
           ty= -h_z * ( sumval(iadd-9-7-2)/prsunt - (press+strext(5)) ) * tenunt
           Write(message,"(16x,'Average surface tension, fluctuations & mean estimate in x (dyn/cm)')")
           Call info(message,.true.)
           Write(message,'(1p,3e12.4)') sumval(iadd+1),ssqval(iadd+1),tx
           Call info(message,.true.)
           Write(message,"(16x,'Average surface tension, fluctuations & mean estimate in y (dyn/cm)')")
           Call info(message,.true.)
           Write(message,'(1p,3e12.4)') sumval(iadd+2),ssqval(iadd+2),ty
           Call info(message,.true.)

           iadd = iadd+2
        End If
     End If

  End If

! Write out remaining registers

  check = .false.
  Do i=iadd+1,mxnstk
     If (Abs(sumval(i)) > zero_plus .or. Abs(ssqval(i)) > zero_plus) check=.true.
  End Do

  If (check) Then
     Write(messages(1),"(12x,'Remaining non-zero statistics registers ')")
     Write(messages(2),"(12x,'Register',7x,'Average value',8x,'r.m.s. fluc.')")
     Call info(messages,2,.true.)
   End If

   If (comm%idnode == 0) Then
     Do i=iadd+1,mxnstk
       If (Abs(sumval(i)) > zero_plus .or. Abs(ssqval(i)) > zero_plus) Then
         Write(message,'(10x,i10,2f20.10)') i,sumval(i),ssqval(i)
         Call info(message,.true.)
       End If
     End Do
   End If

10 Continue

! scale densities for average volume and average volume and cell

  Do i=1,ntpatm
     dens(i)=dens(i)*(volm/avvol)
  End Do

! volm and cell become the averaged ones, as is the local temp

  volm = avvol
  cell = avcel
  temp = sumval(2)

! calculate and print radial distribution functions

!If block average errors, output that, else if jackknife errors output those, else just RDF.
  If (lrdf .and. lprdf .and. ncfrdf > 0 .and. l_errors_block) Then
    Call calculate_errors(temp, rcut, nstep, comm)
  End If
  If (lrdf .and. lprdf .and. ncfrdf > 0 .and. l_errors_jack .and. .not. l_errors_block) Then
    Call calculate_errors_jackknife(temp, rcut, nstep, comm)
  End If
  If (lrdf .and. lprdf .and. ncfrdf > 0 .and. .not.(l_errors_block .or. l_errors_jack)) Then
    Call rdf_compute(lpana,rcut,temp,comm)
  End IF
  If (ncfusr > 0) Call usr_compute(comm)

! calculate and print z-density profile
  If (lzdn .and. lpzdn .and. ncfzdn > 0) Then
    Call z_density_compute(comm)
  End If

! calculate and print velocity autocorrelation function
  If (vafsamp > 0 .and. lpvaf .and. vafcount > zero_plus) Then
    Call vaf_compute(lvafav,tstep,comm)
  End If

! Calculate and print PDFs
  If (lpana) Then
     If (mxgbnd1 > 0 .and. ncfbnd > 0) Call bonds_compute(temp,comm)
     If (mxgang1 > 0 .and. ncfang > 0) Call angles_compute(temp,comm)
     If (mxgdih1 > 0 .and. ncfdih > 0) Call dihedrals_compute(temp,comm)
     If (mxginv1 > 0 .and. ncfinv > 0) Call inversions_compute(temp,comm)
  End If

20 Continue

! print final time check

  Call gtime(timelp)

  Write(message,'("time elapsed since job start: ", f12.3, " sec")') timelp
  Call info(message,.true.)

End Subroutine statistics_result
End Module statistics
