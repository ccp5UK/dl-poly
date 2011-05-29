Subroutine statistics_collect                &
           (leql,nsteql,lzdn,nstzdn,         &
           keyres,keyens,intsta,imcon,       &
           degfre,degshl,degrot,             &
           nstep,tstep,time,tmst,            &
           engcpe,vircpe,engsrp,virsrp,      &
           engter,virter,                    &
           engtbp,virtbp,engfbp,virfbp,      &
           engshl,virshl,shlke,              &
           vircon,virpmf,                    &
           engtet,virtet,engfld,virfld,      &
           engbnd,virbnd,engang,virang,      &
           engdih,virdih,enginv,virinv,      &
           engke,engrot,consv,vircom,strtot, &
           stpeng,stpvir,stpcfg,stpeth,      &
           stptmp,stpprs,stpvol)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating periodic data during the
! molecular dynamics simulation and computing the rolling averages
!
! copyright - daresbury laboratory
! author    - w.smith august 1992
! amended   - i.t.todorov april 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gsum
  Use setup_module
  Use site_module,    Only : ntpatm,numtyp
  Use config_module,  Only : cfgname,cell,volm,natms,ltype, &
                             xxx,yyy,zzz,vxx,vyy,vzz
  Use statistics_module
  Use msd_module

  Implicit None

  Logical,           Intent( In    ) :: leql,lzdn
  Integer,           Intent( In    ) :: nsteql,nstzdn,keyres, &
                                        keyens,intsta,imcon,nstep

  Integer(Kind=ip),  Intent( In    ) :: degfre,degshl,degrot

  Real( Kind = wp ), Intent( In    ) :: tstep,time,          &
                                        engcpe,vircpe,engsrp,virsrp, &
                                        engter,virter,               &
                                        engtbp,virtbp,engfbp,virfbp, &
                                        engshl,virshl,shlke,         &
                                        vircon,virpmf,               &
                                        engtet,virtet,engfld,virfld, &
                                        engbnd,virbnd,engang,virang, &
                                        engdih,virdih,enginv,virinv, &
                                        engke,engrot,consv,vircom,strtot(1:9)

  Real( Kind = wp ), Intent( InOut ) :: tmst
  Real( Kind = wp ), Intent(   Out ) :: stpeng,stpvir,stpcfg,stpeth, &
                                        stptmp,stpprs,stpvol

  Logical,           Save :: newjob = .true.

  Logical                 :: l_tmp
  Integer                 :: fail,i,j,k,iadd,kstak
  Real( Kind = wp )       :: stpcns,stpshl,stprot,sclnv1,sclnv2,zistk,celprp(1:10)


  Real( Kind = wp ), Allocatable :: amsd(:)

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
     If (keyres == 1) Inquire(File='STATIS', Exist=l_tmp)

     If (.not.l_tmp) Then
        Open(Unit=nstats, File='STATIS', Status='replace')

        Write(nstats,'(a)') cfgname

        If (Abs(engunit-9648.530821_wp) <= zero_plus) &
           Write(nstats,'(1x,a)') 'ENERGY UNITS = electron Volts'
        If (Abs(engunit-418.4_wp)       <= zero_plus) &
           Write(nstats,'(1x,a)') 'ENERGY UNITS = kcal/mol'
        If (Abs(engunit-1.0e2_wp)       <= zero_plus) &
           Write(nstats,'(1x,a)') 'ENERGY UNITS = kjoule/mol'
        If (Abs(engunit-1.0_wp)         <= zero_plus) &
           Write(nstats,'(1x,a)') 'ENERGY UNITS = DL_POLY Internal UNITS'

        Close(Unit=nstats)
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
! Note: virfbp, virinv and virdih are allegedly always zero!!!

  stpvir = vircpe + virsrp + virter + virtbp + virfbp + &
           virfld + virshl + vircon + virpmf + vircom + &
           virtet + virbnd + virang + virdih + virinv

! system volume

  stpvol = volm

! system pressure

  stpprs = (2.0_wp*engke-stpvir) / (3.0_wp*stpvol)

! system enthalpy

  stpeth = stpeng + stpprs*stpvol

  Call dcell(cell,celprp)

! store current values in statistics array

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
     If (stptmp < 1.0_wp) Then ! HISTORY is replayed and no velocity field exists
        Do i=1,natms
           xto(i)=xxx(i)-xin(i)
           yto(i)=yyy(i)-yin(i)
           zto(i)=zzz(i)-zin(i)
        End Do
        Call images(imcon,cell,natms,xto,yto,xto)
     Else                      ! velocity field exists
        Do i=1,natms
           xto(i)=xto(i)+vxx(i)*tstep
           yto(i)=yto(i)+vyy(i)*tstep
           zto(i)=zto(i)+vzz(i)*tstep
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
     If (numtyp(k) > zero_plus) stpval(iadd+k)=amsd(k)/Max(numtyp(k),1.0_wp)
  End Do

  iadd = iadd + ntpatm

! pressure tensor (derived for the stress tensor)

  Do i=1,9
     stpval(iadd+i)=strtot(i)*prsunt/(stpvol)
  End Do

  iadd = iadd + 9

! cell parameters

  If (keyens >= 20) Then
     Do i=1,9
       stpval(iadd+i)=cell(i)
     End Do

     iadd = iadd + 9
  End If

! write statistics file

  If (idnode == 0 .and. Mod(nstep,intsta) == 0) Then
     Open(Unit=nstats, File='STATIS', Position='append')

     If (l_msd) Then
        Write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))') &
             nstep,time,iadd-2*mxatdm,(stpval(k),k=1,27),(stpval(k),k=28+2*mxatdm,iadd)
     Else
        Write(nstats,'(i10,1p,e14.6,0p,i10,/,(1p,5e14.6))') &
             nstep,time,iadd,(stpval(k),k=1,iadd)
     End If

     Close(Unit=nstats)
  End If

! check on number of variables for stack

  If (iadd > mxnstk) Call error(170)

! No totals for timestep zero

  If (nstep == 0) Go To 10

! current stack value

  kstak=Mod(nstep-1,mxstak)+1

! subtract old stack value from the stack average

  If (nstep > mxstak) Then
     Do i=1,mxnstk
        zumval(i)=zumval(i)-stkval(kstak,i)
     End Do
  End If

! store quantities in stack and update the stack average

  Do i=1,mxnstk
     stkval(kstak,i)=stpval(i)
     zumval(i)=zumval(i)+stpval(i)
  End Do

! calculate rolling averages

  zistk=Min(mxstak,nstep)

  Do i=1,mxnstk
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
     Do i=1,mxnstk
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
