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
           stptmp,stpprs,stpvol)

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

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module
  Use site_module,    Only : ntpatm,numtypnf
  Use config_module,  Only : cfgname,imcon,cell,volm,natms,ltype, &
                             xxx,yyy,zzz,vxx,vyy,vzz
  Use dpd_module,     Only : virdpd
  Use statistics_module
  Use msd_module

  Implicit None

  Logical,           Intent( In    ) :: lsim,leql,lzdn
  Integer,           Intent( In    ) :: nsteql,nstzdn,keyres, &
                                        keyens,iso,intsta,nstep

  Integer(Kind=ip),  Intent( In    ) :: degfre,degshl,degrot

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

  Logical,           Save :: newjob = .true.

  Logical                 :: l_tmp
  Integer                 :: fail,i,j,k,iadd,kstak
  Real( Kind = wp )       :: stpcns,stpshl,stprot,stpipv,celprp(1:10), &
                             zistk,sclnv1,sclnv2,h_z

  Real( Kind = wp ), Allocatable :: amsd(:)
  Real( Kind = wp ), Allocatable :: xxt(:),yyt(:),zzt(:)

  fail=0
  Allocate (amsd(1:mxatyp), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'statistics_collect allocation failure, node: ', idnode
     Call error(0)
  End If

! open statistics file and put header

  If (newjob .and. idnode == 0) Then
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
           Write(nrite,'(/,1x,a,i0)') 'statistics_collect allocation failure 1, node: ', idnode
           Call error(0)
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
           Write(nrite,'(/,1x,a,i0)') 'statistics_collect deallocation failure 1, node: ', idnode
           Call error(0)
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
     If (mxnode > 1) Call gsum(amsd(1:ntpatm))
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

  If (idnode == 0 .and. Mod(nstep,intsta) == 0) Then
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
     Write(nrite,'(/,1x,a,i0)') 'statistics_collect deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine statistics_collect
