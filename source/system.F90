Module system

  Use kinds, Only : wp,li
  Use comms, Only : comms_type, gbcast,SysExpand_tag,Revive_tag,wp_mpi,gsync, &
                    gsend
  Use setup
  Use site,        Only : ntpatm,numtyp,numtypnf,dens,ntpmls,numsit,&
                                 nummols
  Use configuration,      Only : volm,natms,ltg,ltype,lfrzn,xxx,yyy,zzz, &
                                 cfgname,imcon,cell,lsi,lsa,atmnam, &
                                 write_config
  Use statistics
  Use rdfs,         Only : ncfrdf,rdf,ncfusr,rusr,usr
  Use z_density,   Only : ncfzdn,zdens
  Use bonds,       Only : ldfbnd,ncfbnd,dstbnd,numbonds,lstbnd,keybnd
  Use angles,      Only : ldfang,ncfang,dstang,numang,lstang,keyang
  Use dihedrals,   Only : ldfdih,ncfdih,dstdih,numdih,lstdih
  Use inversions,  Only : ldfinv,ncfinv,dstinv,numinv,lstinv
  Use vdw,         Only : ls_vdw,ntpvdw
  Use metal,       Only : ntpmet
  Use greenkubo,   Only : nsvaf,vafsamp,vafcount,vafstep, &
                                 vxi,vyi,vzi,vafdata,vaf,vaftime

  Use development, Only : l_rin, l_rout
  Use core_shell,   Only : numshl,lstshl
  Use constraints,  Only : numcon,lstcon
  Use rigid_bodies, Only : numrgd,lstrgd
  Use parse,        Only : tabs_2_blanks, get_word, strip_blanks, &
                                  lower_case, word_2_real
  Use io,           Only : io_set_parameters,        &
                                  io_get_parameters,        &
                                  io_init, io_nc_create,    &
                                  io_open, io_write_record, &
                                  io_nc_put_var,            &
                                  io_write_sorted_file,     &
                                  io_close, io_finalize,    &
                                  io_delete,                &
                                  io_close, io_finalize,    &
                                  IO_RESTART,               &
                                  IO_BASE_COMM_NOT_SET,     &
                                  IO_ALLOCATION_ERROR,      &
                                  IO_UNKNOWN_WRITE_OPTION,  &
                                  IO_UNKNOWN_WRITE_LEVEL,   &
                                  IO_WRITE_UNSORTED_MPIIO,  &
                                  IO_WRITE_UNSORTED_DIRECT, &
                                  IO_WRITE_UNSORTED_MASTER, &
                                  IO_WRITE_SORTED_MPIIO,    &
                                  IO_WRITE_SORTED_DIRECT,   &
                                  IO_WRITE_SORTED_NETCDF,   &
                                  IO_WRITE_SORTED_MASTER
  Use vdw,             Only : vdw_lrc
  Use metal,           Only : metal_lrc
  Use errors_warnings, Only : error, info
  Implicit None
  Private
  Public :: system_revive
  Public :: system_init
  Public :: system_expand
  Contains
  
  Subroutine system_init                                             &
           (levcfg,rcut,rvdw,rbin,rmet,lrdf,lzdn,keyres,megatm,    &
           time,tmst,nstep,tstep,chit,cint,chip,eta,virtot,stress, &
           vircon,strcon,virpmf,strpmf,elrc,virlrc,elrcm,vlrcm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading the REVIVE file data and defining the
! initial thermodynamic and structural accumulators
!
! copyright - daresbury laboratory
! author    - i.t.todorov november 2016
! contrib   - m.a.seaton june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( In    ) :: lrdf,lzdn
  Integer,           Intent( InOut ) :: levcfg,keyres
  Integer,           Intent( In    ) :: megatm
  Real( Kind = wp ), Intent( In    ) :: rcut,rvdw,rbin,rmet

  Integer,           Intent(   Out ) :: nstep
  Real( Kind = wp ), Intent( InOut ) :: tstep
  Real( Kind = wp ), Intent(   Out ) :: time,tmst,chit,cint,chip,eta(1:9),     &
                                        virtot,stress(1:9),vircon,strcon(1:9), &
                                        virpmf,strpmf(1:9),elrc,virlrc,        &
                                        elrcm(0:mxatyp),vlrcm(0:mxatyp)
  Type( comms_type ), Intent( InOut ) :: comm

  Character( Len = 40 ) :: forma  = ' '

  Logical               :: l_tmp
  Integer               :: i,j,k,l,keyio,i_tmp,gidx
  Real( Kind = wp )     :: dnstep,dtstep,dnumacc,dncfrdf,dmxgusr,drusr,dncfusr,dncfzdn, &
                           dncfbnd,dncfang,dncfdih,dncfinv,r_mxnode,xyz(0:6),dvafstep(1:vafsamp)


! Define format for REVOLD reading in ASCII

  If (l_rin) Then
     i = 64/4 - 1 ! Bit_Size(0.0_wp)/4 - 1
     j = Max(mxstak*mxnstk,mxgrdf*mxrdf,mxgusr,mxgana*mxtana)

     Write(forma ,10) j/4+1,i+9,i
10   Format('(1p,',i0,'(/,4e',i0,'.',i0,'E3))')
  End If

! Initialise read failure flag

  keyio=0

50 Continue
  If (keyres /= 1 .or. comm%idnode /= 0) Then

! initialise step and time related accumulators

     nstep =0
     dtstep=0.0_wp
     time  =0.0_wp
     tmst  =0.0_wp

! initialise temperature and pressure coupling parameters
! and integral for conserved quantity

     chit = 0.0_wp
     cint = 0.0_wp
     chip = 0.0_wp
     eta  = 0.0_wp

! initialise strcon,stress,virtot and vircon

     virtot = 0.0_wp
     stress = 0.0_wp
     vircon = 0.0_wp
     strcon = 0.0_wp
     virpmf = 0.0_wp
     strpmf = 0.0_wp

! initialise accumulator arrays if reading failure occurred

     If (keyio > 0) Then

        numacc=0
        stpval=0.0_wp
        stpvl0=0.0_wp
        sumval=0.0_wp
        ssqval=0.0_wp
        zumval=0.0_wp
        ravval=0.0_wp
        stkval=0.0_wp

        If (lrdf) Then
           ncfrdf=0
           rdf   =0.0_wp
        End If

        If (mxgusr > 0) Then
           rusr  =0.0_wp
           ncfusr=0
           usr   =0
        End If

        If (lzdn) Then
           ncfzdn=0
           zdens =0.0_wp
        End If

        If (vafsamp > 0) Then
           vafcount=0.0_wp
           vafstep =0
           vafdata =0.0_wp
           vaf     =0.0_wp
           vaftime =0.0_wp
        End If

        If (mxgbnd1 > 0) Then
           ncfbnd=0
           dstbnd=0.0_wp
        End If

        If (mxgang1 > 0) Then
           ncfang=0
           dstang=0.0_wp
        End If

        If (mxgdih1 > 0) Then
           ncfdih=0
           dstdih=0.0_wp
        End If

        If (mxginv1 > 0) Then
           ncfinv=0
           dstinv=0.0_wp
        End If

     End If

  End If

! restart simulation and continue

  If (keyres == 1) Then

! If REVOLD doesn't exist then abort (mishmashed REVOLD is handled separately)

     l_tmp=.true.
     If (comm%idnode == 0) Inquire(File=Trim(revold), Exist=l_tmp)
     Call gcheck(comm,l_tmp,"enforce")
     If (.not.l_tmp) Call error(519)

! Check REVOLD restart compatibility: rcut,rbin,megatm

     xyz(1:3)=0.0_wp
     If (comm%idnode == 0) Then
        If (l_rin) Then
           Open(Unit=nrest, file=Trim(revold), form='formatted', IOStat=keyio)
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) xyz(1),xyz(2),xyz(3)
        Else
           Open(Unit=nrest, file=Trim(revold), form='unformatted', IOStat=keyio)
           Read(Unit=nrest, IOStat=keyio, End=100) xyz(1),xyz(2),xyz(3)
        End If
     End If
     Call gbcast(comm,xyz,0)
     If (Abs(xyz(1)-rcut) > 1.0e-6_wp .or. Abs(xyz(2)-rbin) > 1.0e-6_wp .or. &
         Nint(xyz(3)) /= megatm) Call error(519)

! read the rest of the accumulator data from dump file

     If (comm%idnode == 0) Then
        If (l_rin) Then
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) &
               dnstep,dtstep,time,tmst,dnumacc,chit,chip,cint
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) eta
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) stpval
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) stpvl0
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) sumval
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) ssqval
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) zumval
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) ravval
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) stkval
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) strcon
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) strpmf
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) stress

           If (lrdf) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfrdf,rdf
           If (mxgusr > 0) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dmxgusr,drusr,dncfusr,usr
           If (lzdn) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfzdn,zdens
           If (vafsamp > 0) Then
             Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) vafcount
             Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dvafstep
             Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) vafdata
             Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) vaf
             Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) vaftime
           End If

           If (mxgbnd1 > 0) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfbnd,dstbnd
           If (mxgang1 > 0) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfang,dstang
           If (mxgdih1 > 0) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfdih,dstdih
           If (mxginv1 > 0) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) dncfinv,dstinv
        Else
           Read(Unit=nrest, IOStat=keyio, End=100) &
               dnstep,dtstep,time,tmst,dnumacc,chit,chip,cint
           Read(Unit=nrest, IOStat=keyio, End=100) eta
           Read(Unit=nrest, IOStat=keyio, End=100) stpval
           Read(Unit=nrest, IOStat=keyio, End=100) stpvl0
           Read(Unit=nrest, IOStat=keyio, End=100) sumval
           Read(Unit=nrest, IOStat=keyio, End=100) ssqval
           Read(Unit=nrest, IOStat=keyio, End=100) zumval
           Read(Unit=nrest, IOStat=keyio, End=100) ravval
           Read(Unit=nrest, IOStat=keyio, End=100) stkval
           Read(Unit=nrest, IOStat=keyio, End=100) strcon
           Read(Unit=nrest, IOStat=keyio, End=100) strpmf
           Read(Unit=nrest, IOStat=keyio, End=100) stress

           If (lrdf) Read(Unit=nrest, IOStat=keyio, End=100) dncfrdf,rdf
           If (mxgusr > 0) Read(Unit=nrest, IOStat=keyio, End=100) dmxgusr,drusr,dncfusr,usr
           If (lzdn) Read(Unit=nrest, IOStat=keyio, End=100) dncfzdn,zdens
           If (vafsamp > 0) Then
             Read(Unit=nrest, IOStat=keyio, End=100) vafcount
             Read(Unit=nrest, IOStat=keyio, End=100) dvafstep
             Read(Unit=nrest, IOStat=keyio, End=100) vafdata
             Read(Unit=nrest, IOStat=keyio, End=100) vaf
             Read(Unit=nrest, IOStat=keyio, End=100) vaftime
           End If

           If (mxgbnd1 > 0) Read(Unit=nrest, IOStat=keyio, End=100) dncfbnd,dstbnd
           If (mxgang1 > 0) Read(Unit=nrest, IOStat=keyio, End=100) dncfang,dstang
           If (mxgdih1 > 0) Read(Unit=nrest, IOStat=keyio, End=100) dncfdih,dstdih
           If (mxginv1 > 0) Read(Unit=nrest, IOStat=keyio, End=100) dncfinv,dstinv
        End If

        nstep =Nint(dnstep)

        numacc=Nint(dnumacc)

        If (lrdf) ncfrdf=Nint(dncfrdf)
        If (mxgusr > 0) Then
           mxgusr=Nint(dmxgusr)
           rusr  =drusr
           ncfusr=Nint(dncfusr)
        End If
        If (lzdn) ncfzdn=Nint(dncfzdn)
        If (vafsamp > 0) vafstep=Nint(dvafstep)

        If (mxgbnd1 > 0) ncfbnd=Nint(dncfbnd)
        If (mxgang1 > 0) ncfang=Nint(dncfang)
        If (mxgdih1 > 0) ncfdih=Nint(dncfdih)
        If (mxginv1 > 0) ncfinv=Nint(dncfinv)

! calculate virtot = virtot-vircon-virpmf

        vircon = stpval(17) * engunit
        virpmf = stpval(26) * engunit
        virtot = (stpval(12)-stpval(17)-stpval(26)) * engunit
     End If

100  Continue

! If 'restart' is impossible go to 'restart noscale' and reinitialise

     Call gbcast(comm,keyio,0)
     If (keyio /= 0) Then
        If (comm%idnode == 0) Then
           Call warning(190,0.0_wp,0.0_wp,0.0_wp)
           Close(Unit=nrest)
        End If
        keyres=3
        Go To 50
     End If

! broadcast stored variables

     If (comm%mxnode > 1) Then

        Call gbcast(comm,nstep,0)
        Call gbcast(comm,dtstep,0)
        Call gbcast(comm,time,0)
        Call gbcast(comm,tmst,0)
        Call gbcast(comm,numacc,0)
        Call gbcast(comm,chit,0)
        Call gbcast(comm,chip,0)
        Call gbcast(comm,cint,0)
        Call gbcast(comm,eta,0)
        Call gbcast(comm,stpval,0)
        Call gbcast(comm,stpvl0,0)
        Call gbcast(comm,sumval,0)
        Call gbcast(comm,ssqval,0)
        Call gbcast(comm,zumval,0)
        Call gbcast(comm,ravval,0)
        Do k=0,mxnstk
           Call gbcast(comm,stkval(:,k),0)
        End Do
        Call gbcast(comm,strcon,0)
        Call gbcast(comm,strpmf,0)
        Call gbcast(comm,stress,0)
        Call gbcast(comm,vircon,0)
        Call gbcast(comm,virpmf,0)
        Call gbcast(comm,virtot,0)

! Reset timestep

        tstep = dtstep

        r_mxnode=1.0_wp/Real(comm%mxnode,wp)

! rdf table - broadcast and normalise

        If (lrdf) Then
           Call gbcast(comm,ncfrdf,0)
           Do k=1,mxrdf
              Call gbcast(comm,rdf(:,k),0)
              rdf(:,k) = rdf(:,k) * r_mxnode
           End Do
        End If

! USR RDF table - broadcast and normalise

        If (mxgusr > 0) Then
           Call gbcast(comm,mxgusr,0)
           Call gbcast(comm,rusr,0)
           Call gbcast(comm,ncfusr,0)
           Call gbcast(comm,usr,0)
           usr(:) = usr(:) * r_mxnode
        End If

! z-density table - broadcast and normalise

        If (lzdn) Then
           Call gbcast(comm,ncfzdn,0)
           Do k=1,mxatyp
              Call gbcast(comm,zdens(:,k),0)
              zdens(:,k) = zdens(:,k) * r_mxnode
           End Do
        End If

! vafdata table - broadcast and normalise

        If (vafsamp > 0) Then
           Do j=1,vafsamp
              l=(j-1)*(mxatyp+1)
              Do k=1,mxatyp+1
                 Call gbcast(comm,vafdata(:,l+k),0)


! avoid normalising timing information

                 If (k /= mxatyp+1) vafdata(:,l+k) = vafdata(:,l+k) * r_mxnode
              End Do
           End Do
        End If

! bonds table - broadcast and normalise

        If (mxgbnd1 > 0) Then
           Call gbcast(comm,ncfbnd,0)
           Do k=1,ldfbnd(0)
              Call gbcast(comm,dstbnd(:,k),0)

              dstbnd(:,k) = dstbnd(:,k) * r_mxnode
           End Do
        End If

! angles table - broadcast and normalise

        If (mxgang1 > 0) Then
           Call gbcast(comm,ncfang,0)
           Do k=1,ldfang(0)
              Call gbcast(comm,dstang(:,k),0)

              dstang(:,k) = dstang(:,k) * r_mxnode
           End Do
        End If

! dihedrals table - broadcast and normalise

        If (mxgdih1 > 0) Then
           Call gbcast(comm,ncfdih,0)
           Do k=1,ldfdih(0)
              Call gbcast(comm,mxgdih1,0)

              dstdih(:,k) = dstdih(:,k) * r_mxnode
           End Do
        End If

! inversions table - broadcast and normalise

        If (mxginv1 > 0) Then
           Call gbcast(comm,ncfinv,0)
           Do k=1,ldfinv(0)
              Call gbcast(comm,mxginv1,0)

              dstinv(:,k) = dstinv(:,k) * r_mxnode
           End Do
        End If
     End If

  End If

! initialise initial positions to current positions
! and final displacements to zero

  Do i=1,natms
     xin(i)=xxx(i)
     yin(i)=yyy(i)
     zin(i)=zzz(i)

     xto(i)=0.0_wp
     yto(i)=0.0_wp
     zto(i)=0.0_wp
  End Do

  If (keyres == 1) Then

! Error accumulator: keyio is still zero otherwise we cannot get here

     i_tmp=0

     Do k=1,megatm
        xyz=0.0_wp

        If (comm%idnode == 0) Then
           If (l_rin) Then
              Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio) &
                  xyz(0),xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6)
           Else
              Read(Unit=nrest, IOStat=keyio) &
                  xyz(0),xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6)
           End If
        End If
        If (keyio /= 0) i_tmp=1

        Call gbcast(comm,xyz,0)
        gidx=Nint(xyz(0))

! assign particle initial positions and final displacements
! to the corresponding domains

        Do i=1,natms
           If (ltg(i) == gidx) Then
              xin(i)=xyz(1)
              yin(i)=xyz(2)
              zin(i)=xyz(3)

              xto(i)=xyz(4)
              yto(i)=xyz(5)
              zto(i)=xyz(6)
           End If
        End Do
     End Do

! If 'restart' is impossible go to 'restart noscale' and reinitialise

     Call gbcast(comm,i_tmp,0)
     If (i_tmp /= 0) Then
        If (comm%idnode == 0) Then
           Call warning(190,0.0_wp,0.0_wp,0.0_wp)
           Close(Unit=nrest)
        End If
        keyres=3
        Go To 50
     End If

! Read velocities for VAF calculations if needed

     If (vafsamp > 0) Then

       i_tmp=0

       Do j=1,vafsamp
         Do k=1,megatm
            xyz=0.0_wp

            If (comm%idnode == 0) Then
               If (l_rin) Then
                  Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio) &
                      xyz(0),xyz(1),xyz(2),xyz(3)
               Else
                  Read(Unit=nrest, IOStat=keyio) &
                      xyz(0),xyz(1),xyz(2),xyz(3)
               End If
            End If
            If (keyio /= 0) i_tmp=1

            Call gbcast(comm,xyz,0)
            gidx=Nint(xyz(0))

! assign particle velocities to the corresponding domains

            Do i=1,natms
               If (ltg(i) == gidx) Then
                  vxi(i,j)=xyz(1)
                  vyi(i,j)=xyz(2)
                  vzi(i,j)=xyz(3)
               End If
            End Do
         End Do
       End Do

! If 'restart' is impossible go to 'restart noscale' and reinitialise

       Call gbcast(comm,i_tmp,0)
       If (i_tmp /= 0) Then
          If (comm%idnode == 0) Then
             Call warning(190,0.0_wp,0.0_wp,0.0_wp)
             Close(Unit=nrest)
          End If
          keyres=3
          Go To 50
       End If

     End If

     If (comm%idnode == 0) Close(Unit=nrest)

  Else

! force force and stress recalculation when 'restart' is not on

     levcfg=1

  End If


! number densities needed for long-range corrections

! evaluate species populations in system (separate totals for non-frozen atoms)

  Do i=1,natms
     k = ltype(i)
     numtyp(k) = numtyp(k)+1.0_wp
     If (lfrzn(i) == 0) numtypnf(k) = numtypnf(k)+1.0_wp
  End Do

! global number densities

  
    Call gsum(comm,numtyp(1:ntpatm))
    Call gsum(comm,numtypnf(1:ntpatm))
  

! number densities

  Do i=1,ntpatm
     If (numtyp(i) > zero_plus) dens(i) = numtyp(i)/volm
  End Do

! Get long-range corrections

! elrc & virlrc arrays are zeroed in vdw,
! no lrc when vdw interactions are force-shifted

  If (ntpvdw > 0 .and. (.not.ls_vdw)) Call vdw_lrc(rvdw,elrc,virlrc,comm)

! elrcm & vlrcm arrays are zeroed in metal_module

  If (ntpmet > 0) Call metal_lrc(rmet,elrcm,vlrcm,comm)

End Subroutine system_init

Subroutine system_expand(l_str,rcut,nx,ny,nz,megatm,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 utility to expand the MD system by a nx*ny*nz volumetric
! replication of its contents along the MD cell lattice vectors,
! creating a new matching pair of topology-interaction (FIELD) and
! crystallographic (CONFIG) files, preserving FIELD's template intact
!
! supported image conditions: 1,2,3, 6(nz==1)
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
! contrib   - w.smith, i.j.bush
! contrib   - a.m.elena february 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Logical,           Intent( In    ) :: l_str
  Integer,           Intent( In    ) :: nx,ny,megatm
  Real( Kind = wp ), Intent( In    ) :: rcut
  Integer,           Intent( InOut ) :: nz
  Type( comms_type ), Intent( InOut ) :: comm

  Integer, Parameter     :: recsz = 73 ! default record size

  Logical                :: safex,safey,safez,safer,safel,safem,safe,safeg, &
                            lmpldt=.false.

  Character( Len = 200 ) :: record,record1
  Character( Len = 40  ) :: word,fcfg,ffld,fmpl
  Integer                :: fail(1:5),nall, i,j,ix,iy,iz,m, &
                            itmols,setspc,imols,          &
                            indatm,indatm1,nattot,mxiter, &
                            sapmpt,sapmtt,iatm,jatm,      &
                            ishls,nshels,icnst,nconst,    &
                            nrigid,irgd,lrgd,             &
                            ibond,nbonds,iang,nangle,     &
                            idih,ndihed,iinv,ninver,      &
                            idm,loc_ind,index,at_scaled
  Integer(Kind=li)       :: offset,rec
  Real( Kind = wp )      :: fx,fy,fz, x,y,z, t,c1,c2,c3,c4, &
                            celprp(1:10), hwx,hwy,hwz, r,   &
                            x1(1:1),y1(1:1),z1(1:1),        &
                            cell_vecs(1:3,1:3), lengths(1:3), angles(1:3)

! Some parameters and variables needed by io interfaces

  Integer                           :: fh, io_write
  Character( Len = recsz )          :: record2, record3
  Character                         :: lf

  Integer( Kind = MPI_OFFSET_KIND ) :: top_skip

  Real( Kind = wp ), Dimension( : ),     Allocatable :: f1,f2,f3, &
                                                        f4,f5,f6, &
                                                        f7,f8,f9, &
                                                        x_scaled, &
                                                        y_scaled, &
                                                        z_scaled, &
                                                        xm,ym,zm

  Integer,           Dimension( :,:,: ), Allocatable :: i_xyz

  Integer,           Dimension( : ),     Allocatable :: ltg_scaled

  Character( Len = Len( atmnam ) ), Dimension( : ), Allocatable :: atmnam_scaled
  Integer :: ierr

  Character ( len = 256 )  :: message

  fail=0
  Allocate (f1(1:nx),f2(1:nx),f3(1:nx),                      Stat=fail(1))
  Allocate (f4(1:ny),f5(1:ny),f6(1:ny),                      Stat=fail(2))
  Allocate (f7(1:nz),f8(1:nz),f9(1:nz),                      Stat=fail(3))
  Allocate (i_xyz(1:nx,1:ny,1:nz),                           Stat=fail(4))
  Allocate (xm(1:10*mxatms),ym(1:10*mxatms),zm(1:10*mxatms), Stat=fail(5))

  If (Any(fail > 0)) Then
     Write(message,'(a)') 'system_expand allocation failure'
     Call error(0,message)
  End If

! Get write buffer size and line feed character

  Call io_get_parameters( user_method_write      = io_write )
  Call io_get_parameters( user_line_feed         = lf       )

! Print elapsed time and option header

  Call gtime(t)
  If (comm%idnode == 0) Then
     Write(nrite,'(/, "time elapsed since job start: ", f12.3, " sec",/)') t
     Write(nrite,'(4(1x,a,/))')                                                     &
     "*** Expanding the MD system by a nx*ny*nz volumetric replication        ***", &
     "*** of its contents along the MD cell lattice vectors, creating         ***", &
     "*** a new matching pair of topology-interaction (FIELD) and             ***", &
     "*** crystallographic (CONFIG) files, preserving FIELD's template intact ***"
     Write(nrite,'(1x,a,3i5,/)') '*** Replication dimensions (nx,ny,nz):', nx,ny,nz
  End If

! Holt or change execution if imcon is unsupported

  If (imcon == 0) Call error(570)
  If (imcon == 6 .and. nz > 1) Then
     nz=1
     Call warning(350,0.0_wp,0.0_wp,0.0_wp)
     Write(nrite,'(1x,a,3i5,/)') '*** Replication dimensions (nx,ny,nz):', nx,ny,nz
  End If

! Create names for the expanded CONFIG and FIELD

  record= ' ' ; Write(record,'(3(a1,i0))') '_',nx,'_',ny,'_',nz
  fcfg=' '
  fcfg=Trim(config) // record(1:Len_Trim(record))
  ffld=' '
  ffld=Trim(field) // record(1:Len_Trim(record))
  fmpl=' '
  fmpl="MPOLES" // record(1:Len_Trim(record))

! netCDF CONFIG name convention

  If (io_write == IO_WRITE_SORTED_NETCDF) fcfg=fcfg(1:Len_Trim(fcfg)) // '.nc'

  fx=Real(nx,wp)
  fy=Real(ny,wp)
  fz=Real(nz,wp)

  nall=nx*ny*nz

! Define cell vector displacement in z direction

  Do iz=1,nz
     z=Real(2*iz-nz-1,wp)
     f7(iz)=cell(7)*z/2.0_wp
     f8(iz)=cell(8)*z/2.0_wp
     f9(iz)=cell(9)*z/2.0_wp
  End Do

! Define cell vector displacement in y direction

  Do iy=1,ny
     y=Real(2*iy-ny-1,wp)
     f4(iy)=cell(4)*y/2.0_wp
     f5(iy)=cell(5)*y/2.0_wp
     f6(iy)=cell(6)*y/2.0_wp
  End Do

! Define cell vector displacement in x direction

  Do ix=1,nx
     x=Real(2*ix-nx-1,wp)
     f1(ix)=cell(1)*x/2.0_wp
     f2(ix)=cell(2)*x/2.0_wp
     f3(ix)=cell(3)*x/2.0_wp
  End Do

! Define hypercube counter

  Do iz=1,nz
     Do iy=1,ny
        Do ix=1,nx
           i_xyz(ix,iy,iz)=(ix-1)+nx*((iy-1)+ny*(iz-1))
        End Do
     End Do
  End Do

  Call dcell(cell,celprp) ! get cell properties

! define half cell widths and bond-length limit

  hwx = celprp(7)/2.0_wp
  hwy = celprp(8)/2.0_wp
  hwz = celprp(9)/2.0_wp
  c1  = Min(rcut/2.0_wp , 1.75_wp)
  c2  = c1*4.0_wp/3.0_wp
  c3  = c2*4.0_wp/3.0_wp
  c4  = c3*4.0_wp/3.0_wp

  If (comm%idnode == 0) Then

! Make sure CONFIG(new) is empty and open it

     If (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
         io_write == IO_WRITE_UNSORTED_DIRECT .or. &
         io_write == IO_WRITE_SORTED_MPIIO    .or. &
         io_write == IO_WRITE_SORTED_DIRECT   .or. &
         io_write == IO_WRITE_SORTED_NETCDF) Then

        Call io_set_parameters( user_comm = MPI_COMM_SELF )
        Call io_init( recsz )
        Call io_delete( fcfg(1:Len_Trim(fcfg) ),comm )
        If (io_write == IO_WRITE_SORTED_NETCDF) Call io_nc_create( MPI_COMM_SELF, fcfg(1:Len_Trim(fcfg)), cfgname, megatm*nall )
        Call io_open( io_write, MPI_COMM_SELF, fcfg(1:Len_Trim(fcfg)), MPI_MODE_WRONLY + MPI_MODE_CREATE, fh )

     Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
              io_write == IO_WRITE_SORTED_MASTER ) Then

        Open(Unit=nconf, File=fcfg(1:Len_Trim(fcfg)), Status='replace')
        Close(Unit=nconf)
        Open(Unit=nconf, File=fcfg(1:Len_Trim(fcfg)), Form='formatted', Access='direct', Recl=recsz)
     End If

! Write configuration file headers

     Write(nrite,'(1x,2a)') '*** Expanding CONFIG in file ',fcfg(1:Len_Trim(fcfg))
     Write(nrite,'(1x,a)') '***'

     If      (io_write == IO_WRITE_UNSORTED_MPIIO  .or. &
              io_write == IO_WRITE_UNSORTED_DIRECT .or. &
              io_write == IO_WRITE_SORTED_MPIIO    .or. &
              io_write == IO_WRITE_SORTED_DIRECT) Then

        Write(record2, Fmt='(a72,a1)') cfgname(1:72),lf
        Call io_write_record( fh, Int(0,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3i10,a42,a1)') 0,imcon,nall*megatm,Repeat(' ',42),lf
        Call io_write_record( fh, Int(1,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3f20.10,a12,a1)') fx*cell(1),fx*cell(2),fx*cell(3),Repeat(' ',12),lf
        Call io_write_record( fh, Int(2,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3f20.10,a12,a1)') fy*cell(4),fy*cell(5),fy*cell(6),Repeat(' ',12),lf
        Call io_write_record( fh, Int(3,MPI_OFFSET_KIND), record2 )

        Write(record2, Fmt='(3f20.10,a12,a1)') fz*cell(7),fz*cell(8),fz*cell(9),Repeat(' ',12),lf
        Call io_write_record( fh, Int(4,MPI_OFFSET_KIND), record2 )

     Else If (io_write == IO_WRITE_SORTED_NETCDF) Then

        i=1 ! For config there is only one frame

        Call io_nc_put_var( 'time'           , fh, 0.0_wp, i, 1 )
        Call io_nc_put_var( 'step'           , fh,      0, i, 1 )
        Call io_nc_put_var( 'datalevel'      , fh,      0, i, 1 )
        Call io_nc_put_var( 'imageconvention', fh,  imcon, i, 1 )

        cell_vecs( :, 1 ) = fx * cell( 1:3 )
        cell_vecs( :, 2 ) = fy * cell( 4:6 )
        cell_vecs( :, 3 ) = fz * cell( 7:9 )

        lengths( 1 ) = fx * celprp( 1 )
        lengths( 2 ) = fy * celprp( 2 )
        lengths( 3 ) = fz * celprp( 3 )

        angles ( 1 ) = Acos( celprp( 5 ) )
        angles ( 2 ) = Acos( celprp( 6 ) )
        angles ( 3 ) = Acos( celprp( 4 ) )
        angles = angles * 180.0_wp / ( 4.0_wp * Atan( 1.0_wp ) ) ! Convert to degrees

        Call io_nc_put_var( 'cell'        , fh, cell_vecs, (/ 1, 1, i /), (/ 3, 3, 1 /) )
        Call io_nc_put_var( 'cell_lengths', fh, lengths  , (/    1, i /), (/    3, 1 /) )
        Call io_nc_put_var( 'cell_angles' , fh, angles   , (/    1, i /), (/    3, 1 /) )

     Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
              io_write == IO_WRITE_SORTED_MASTER ) Then

        Write(Unit=nconf, Fmt='(a72,a1)',         Rec=Int(1,li)) cfgname(1:72),lf
        Write(Unit=nconf, Fmt='(3i10,a42,a1)',    Rec=Int(2,li)) 0,imcon,nall*megatm,Repeat(' ',42),lf
        Write(Unit=nconf, Fmt='(3f20.12,a12,a1)', Rec=Int(3,li)) fx*cell(1),fx*cell(2),fx*cell(3),Repeat(' ',12),lf
        Write(Unit=nconf, Fmt='(3f20.12,a12,a1)', Rec=Int(4,li)) fy*cell(4),fy*cell(5),fy*cell(6),Repeat(' ',12),lf
        Write(Unit=nconf, Fmt='(3f20.12,a12,a1)', Rec=Int(5,li)) fz*cell(7),fz*cell(8),fz*cell(9),Repeat(' ',12),lf

     End If

  End If

  If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
      io_write == IO_WRITE_SORTED_MPIIO   .or. &
      io_write == IO_WRITE_SORTED_NETCDF) Then

     If (comm%idnode == 0) Then
        Call io_close( fh )
        Call io_finalize
     End If

     Allocate ( atmnam_scaled( 1:natms * nall ), ltg_scaled( 1:natms * nall ), Stat = fail(1) )
     Allocate ( x_scaled( 1:natms * nall ), y_scaled( 1:natms * nall ), z_scaled( 1:natms * nall ), &
               Stat = fail(2) )
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'system_expand allocation failure 0 '
        Call error(0,message)
     End If

  End If

! Line counter in CONFIG(new):levcfg=0

  offset=Int(5,li)

! Global atom counter

  nattot=0

! Local atom counter

  indatm=1

! local counter in scaled system

  at_scaled = 0

! running site and topology related intra-indices

  nshels=0
  nconst=0
  nrigid=0
  nbonds=0
  nangle=0
  ndihed=0
  ninver=0

  Write(message,'(a)') 'Checking topological contiguity of molecules...'
  call info(message,.true.)

  safeg=.true. ! topology presumed safe

  sapmpt=0
  Do itmols=1,ntpmls
     setspc=nummols(itmols)*numsit(itmols)

     sapmtt=0
     Do imols=1,nummols(itmols)
        If (numsit(itmols) > 10*mxatms) Call error(0,message)

! Grab the coordinates of the atoms constituting this molecule

        indatm1=indatm
        Do m=1,numsit(itmols)
           nattot=nattot+1 ! Increase global atom counter in CONFIG(old)

           If (lsa(indatm1) == nattot) Then  ! If a local atom has a global index nattot
              loc_ind=lsi(indatm1)
              xm(m)=xxx(loc_ind)
              ym(m)=yyy(loc_ind)
              zm(m)=zzz(loc_ind)
              indatm1=indatm1+1 ! Increase local atom counter
           Else
              xm(m)=0.0_wp
              ym(m)=0.0_wp
              zm(m)=0.0_wp
           End If
        End Do
        nattot=nattot-numsit(itmols)

           Call gsum(comm,xm(1:numsit(itmols)))
           Call gsum(comm,ym(1:numsit(itmols)))
           Call gsum(comm,zm(1:numsit(itmols)))


! Start unwrapping - not safe at start for each molecule

        indatm1=nattot-sapmpt-sapmtt
        safe=.false. ; mxiter=0
        Do While ((.not.safe) .and. mxiter < 42) ! meaning of LUEE is the limit
           If (.not.safe) mxiter=mxiter+1

           If ((mxiter == 42 .and. (.not.safe)) .and. (l_str .and. comm%idnode == 0)) &
              Write(nrite,Fmt='(/,1x,2(a,i10),/)') 'MOLECULAR TYPE #: ',itmols, ' MOLECULE #: ',imols

           safe=.true.

           safel=.true.
           Do ishls=1,numshl(itmols)
              nshels=nshels+1

              iatm=lstshl(1,nshels)-indatm1
              jatm=lstshl(2,nshels)-indatm1

              safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
              safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
              safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
              safer=(safex .and. safey .and. safez)
              If (.not.safer) Then
                 x1(1)=xm(jatm)-xm(iatm)
                 y1(1)=ym(jatm)-ym(iatm)
                 z1(1)=zm(jatm)-zm(iatm)
                 Call images(imcon,cell,1,x1,y1,z1)
                 xm(jatm)=x1(1)+xm(iatm)
                 ym(jatm)=y1(1)+ym(iatm)
                 zm(jatm)=z1(1)+zm(iatm)
                 safer=.true.
              End If
              x=Abs(xm(jatm)-xm(iatm))
              y=Abs(ym(jatm)-ym(iatm))
              z=Abs(zm(jatm)-zm(iatm))
              t=c1
              safex=(x < t)
              safey=(y < t)
              safez=(z < t)
              r=Sqrt(x**2+y**2+z**2)
              If (safex .and. safey .and. safez) Then
                 safer=(r < t)
              Else
                 safer=.false.
              End If
              If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. comm%idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', r, ' > ', t, ' Angstroms'

                 t=c2
                 If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'

                 Write(nrite,Fmt='(1x,a,3i10)') 'CORE_SHELL UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',ishls,itmols,imols
                 Write(nrite,Fmt='(1x,a,3(1x,l1))') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z',safex,safey,safez
                 Write(nrite,Fmt='(1x,a,i10,3f10.1)')   'CORE  ',nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                 Write(nrite,Fmt='(1x,a,i10,3f10.1,/)') 'SHELL ',nattot+jatm,xm(jatm),ym(jatm),zm(jatm)
              End If
              safel=(safel .and. safer)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do icnst=1,numcon(itmols)
              nconst=nconst+1

              iatm=lstcon(1,nconst)-indatm1
              jatm=lstcon(2,nconst)-indatm1

              safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
              safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
              safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
              safer=(safex .and. safey .and. safez)
              If (.not.safer) Then
                 x1(1)=xm(jatm)-xm(iatm)
                 y1(1)=ym(jatm)-ym(iatm)
                 z1(1)=zm(jatm)-zm(iatm)
                 Call images(imcon,cell,1,x1,y1,z1)
                 xm(jatm)=x1(1)+xm(iatm)
                 ym(jatm)=y1(1)+ym(iatm)
                 zm(jatm)=z1(1)+zm(iatm)
                 safer=.true.
              End If
              x=Abs(xm(jatm)-xm(iatm))
              y=Abs(ym(jatm)-ym(iatm))
              z=Abs(zm(jatm)-zm(iatm))
              t=c2
              safex=(x < t)
              safey=(y < t)
              safez=(z < t)
              r=Sqrt(x**2+y**2+z**2)
              If (safex .and. safey .and. safez) Then
                 safer=(r < t)
              Else
                 safer=.false.
              End If
              If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. comm%idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', r, ' > ', t, ' Angstroms'

                 t=c3
                 If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'

                 Write(nrite,Fmt='(1x,a,3i10)') 'CONSTRAINT UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',icnst,itmols,imols
                 Write(nrite,Fmt='(1x,a,3(1x,l1))') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z',safex,safey,safez
                 Write(nrite,Fmt='(2i10,3f10.1)')   1,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                 Write(nrite,Fmt='(2i10,3f10.1,/)') 2,nattot+jatm,xm(jatm),ym(jatm),zm(jatm)
              End If
              safel=(safel .and. safer)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do irgd=1,numrgd(itmols)
              nrigid=nrigid+1

              safem=.true.
              lrgd=lstrgd(0,nrigid)
              Do i=1,lrgd-1
                 iatm=lstrgd(i,nrigid)-indatm1
                 Do j=i+1,lrgd
                    jatm=lstrgd(j,nrigid)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                          xm(jatm)=x1(1)+xm(iatm)
                          ym(jatm)=y1(1)+ym(iatm)
                          zm(jatm)=z1(1)+zm(iatm)
                          safer=.true.
                    End If
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    t=rcut
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. comm%idnode == 0)) &
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') &
                            '*** WARNING **** WARNING *** DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'
                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. comm%idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'RIGID BODY UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',irgd,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,lrgd
                    iatm=lstrgd(i,nrigid)-indatm1
                    If (i < lrgd) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do ibond=1,numbonds(itmols)
              nbonds=nbonds+1

              iatm=lstbnd(1,nbonds)-indatm1
              jatm=lstbnd(2,nbonds)-indatm1

              safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
              safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
              safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
              safer=(safex .and. safey .and. safez)
              If (.not.safer) Then
                 x1(1)=xm(jatm)-xm(iatm)
                 y1(1)=ym(jatm)-ym(iatm)
                 z1(1)=zm(jatm)-zm(iatm)
                 Call images(imcon,cell,1,x1,y1,z1)
                 xm(jatm)=x1(1)+xm(iatm)
                 ym(jatm)=y1(1)+ym(iatm)
                 zm(jatm)=z1(1)+zm(iatm)
                 safer=.true.
              End If
              x=Abs(xm(jatm)-xm(iatm))
              y=Abs(ym(jatm)-ym(iatm))
              z=Abs(zm(jatm)-zm(iatm))
              If (keybnd(nbonds) > 0) Then
                 t=c2
              Else
                 t=3.0_wp*c1
              End If
              safex=(x < t)
              safey=(y < t)
              safez=(z < t)
              r=Sqrt(x**2+y**2+z**2)
              If (safex .and. safey .and. safez) Then
                 safer=(r < t)
              Else
                 safer=.false.
              End If
              If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. comm%idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', r, ' > ', t, ' Angstroms'

                 If (keybnd(nbonds) > 0) Then
                    t=c3
                 Else
                    t=3.0_wp*c2
                 End If
                 If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'

                 Write(nrite,Fmt='(1x,a,3i10)') 'BOND UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',ibond,itmols,imols
                 Write(nrite,Fmt='(1x,a,3(1x,l1))') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z',safex,safey,safez
                 Write(nrite,Fmt='(2i10,3f10.1)')   1,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                 Write(nrite,Fmt='(2i10,3f10.1,/)') 2,nattot+jatm,xm(jatm),ym(jatm),zm(jatm)
              End If
              safel=(safel .and. safer)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do iang=1,numang(itmols)
              nangle=nangle+1

              safem=.true.
              Do i=1,2
                 iatm=lstang(i,nangle)-indatm1
                 Do j=i+1,3
                    jatm=lstang(j,nangle)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                       xm(jatm)=x1(1)+xm(iatm)
                       ym(jatm)=y1(1)+ym(iatm)
                       zm(jatm)=z1(1)+zm(iatm)
                       safer=.true.
                    End If
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    If (keyang(nangle) > 0) Then
                       t=c1*Real(j-i+1,wp)
                    Else
                       t=c3*Real(j-i+1,wp)
                    End If
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. comm%idnode == 0)) Then
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'

                       If (keyang(nangle) > 0) Then
                          t=c2*Real(j-i+1,wp)
                       Else
                          t=c4*Real(j-i+1,wp)
                       End If
                       If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'
                    End If
                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. comm%idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'ANGLE UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',iang,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,3
                    iatm=lstang(i,nangle)-indatm1
                    If (i < 3) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do idih=1,numdih(itmols)
              ndihed=ndihed+1

              safem=.true.
              Do i=1,3
                 iatm=lstdih(i,ndihed)-indatm1
                 Do j=i+1,4
                    jatm=lstdih(j,ndihed)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                       xm(jatm)=x1(1)+xm(iatm)
                       ym(jatm)=y1(1)+ym(iatm)
                       zm(jatm)=z1(1)+zm(iatm)
                       safer=.true.
                    End If
                    t=(c1+c2)*Real(j-i+1,wp)/2.0_wp
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. comm%idnode == 0)) Then
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'

                       t=(c2+c3)*Real(j-i+1,wp)/2.0_wp
                       If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'
                    End If
                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. comm%idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'DIHEDRAL UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',idih,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,4
                    iatm=lstdih(i,ndihed)-indatm1
                    If (i < 4) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           safel=.true.
           Do iinv=1,numinv(itmols)
              ninver=ninver+1

              safem=.true.
              Do i=1,3
                 iatm=lstinv(i,ninver)-indatm1
                 Do j=i+1,4
                    jatm=lstinv(j,ninver)-indatm1

                    safex=(Abs(xm(jatm)-xm(iatm)) < hwx)
                    safey=(Abs(ym(jatm)-ym(iatm)) < hwy)
                    safez=(Abs(zm(jatm)-zm(iatm)) < hwz)
                    safer=(safex .and. safey .and. safez)
                    If (.not.safer) Then
                       x1(1)=xm(jatm)-xm(iatm)
                       y1(1)=ym(jatm)-ym(iatm)
                       z1(1)=zm(jatm)-zm(iatm)
                       Call images(imcon,cell,1,x1,y1,z1)
                       xm(jatm)=x1(1)+xm(iatm)
                       ym(jatm)=y1(1)+ym(iatm)
                       zm(jatm)=z1(1)+zm(iatm)
                       safer=.true.
                    End If
                    x=Abs(xm(jatm)-xm(iatm))
                    y=Abs(ym(jatm)-ym(iatm))
                    z=Abs(zm(jatm)-zm(iatm))
                    t=c2
                    safex=(x < t)
                    safey=(y < t)
                    safez=(z < t)
                    r=Sqrt(x**2+y**2+z**2)
                    If (safex .and. safey .and. safez) Then
                       safer=(r < t)
                    Else
                       safer=.false.
                    End If
                    If ((mxiter == 42 .and. (.not.safer)) .and. (l_str .and. comm%idnode == 0)) Then
                       Write(nrite,Fmt='(1x,a,2i10,2(f7.2,a))') 'POSSIBLE DISTANCE VIOLATION: ', i,j,r, ' > ', t, ' Angstroms'

                       t=c3
                       If (r > t) Write(nrite,Fmt='(1x,a,f7.2,a)') '*** WARNING **** WARNING *** CUTOFF: ', t, ' Angstroms'
                    End If

                    safem=(safem .and. safer)
                 End Do
              End Do

              If ((mxiter == 42 .and. (.not.safem)) .and. (l_str .and. comm%idnode == 0)) Then
                 Write(nrite,Fmt='(1x,a,3i10)') 'INVERSION UNIT #(LOCAL) -> M. TYPE # -> MOLECULE #:',iinv,itmols,imols
                 Write(nrite,Fmt='(1x,a)') 'MEMBER :: GLOBAL INDEX :: X ::      Y ::      Z'
                 Do i=1,4
                    iatm=lstinv(i,ninver)-indatm1
                    If (i < 4) Then
                       Write(nrite,Fmt='(2i10,3f10.1)')   i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    Else
                       Write(nrite,Fmt='(2i10,3f10.1,/)') i,nattot+iatm,xm(iatm),ym(iatm),zm(iatm)
                    End If
                 End Do
              End If
              safel=(safel .and. safem)
           End Do
           safe=(safe .and. safel)

           If ( ((.not.safe) .and. imols <= nummols(itmols) .and. mxiter < 42) .or. &
                imols < nummols(itmols) ) Then
              nshels=nshels-numshl(itmols)
              nconst=nconst-numcon(itmols)
              nrigid=nrigid-numrgd(itmols)
              nbonds=nbonds-numbonds(itmols)
              nangle=nangle-numang(itmols)
              ndihed=ndihed-numdih(itmols)
              ninver=ninver-numinv(itmols)
           End If
        End Do

        safeg=(safeg.and.safe)

        Do m=1,numsit(itmols)
           nattot=nattot+1 ! Increase global atom counter in CONFIG(old)

           If (lsa(indatm) == nattot) Then ! If a local atom has a global index nattot

! Determine sending node for UN/SORTED MASTER

              If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                  io_write == IO_WRITE_SORTED_MASTER) Then
                 idm=comm%idnode
                 Call gsum(comm,idm)
              End If

! Get the local index of the particle

              loc_ind=lsi(indatm)

! Do particle replication by vector displacements in cyclic (z,y,x) directions

              Do iz=1,nz
                 Do iy=1,ny
                    Do ix=1,nx

                       x=xm(m)+f1(ix)+f4(iy)+f7(iz)
                       y=ym(m)+f2(ix)+f5(iy)+f8(iz)
                       z=zm(m)+f3(ix)+f6(iy)+f9(iz)

! Write 2 records @ line 'rec+1' and 'rec+2' for particle 'index' in CONFIG(new)

                       index = i_xyz(ix,iy,iz)*setspc + m
                       rec   = offset + Int(2,li)*Int(index,li) - Int(2,li)
                       index = index + Int((offset - Int(5,li))/Int(2,li))

                       If (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
                           io_write == IO_WRITE_SORTED_MPIIO   .or. &
                           io_write == IO_WRITE_SORTED_NETCDF) Then

                          at_scaled = at_scaled + 1

                          atmnam_scaled( at_scaled ) = atmnam( loc_ind )

                          ltg_scaled( at_scaled ) = index

                          x_scaled( at_scaled ) = x
                          y_scaled( at_scaled ) = y
                          z_scaled( at_scaled ) = z

                       Else

                          Write(record2, Fmt='(a8,i10,a54,a1)') atmnam(loc_ind),index,Repeat(' ',54),lf
                          Write(record3, Fmt='(3g20.12,a12,a1)') x,y,z,Repeat(' ',12),lf

                          If (io_write == IO_WRITE_UNSORTED_DIRECT .or. &
                              io_write == IO_WRITE_SORTED_DIRECT) Then

                             Call io_write_record( fh, Int(rec,MPI_OFFSET_KIND), record2 )
                             rec=rec+Int(1,li)
                             Call io_write_record( fh, Int(rec,MPI_OFFSET_KIND), record3 )

                          Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                                   io_write == IO_WRITE_SORTED_MASTER) Then

                             If (comm%idnode == 0) Then
                                rec=rec+Int(1,li)
                                Write(Unit=nconf, Fmt='(a73)', Rec=rec) record2
                                rec=rec+Int(1,li)
                                Write(Unit=nconf, Fmt='(a73)', Rec=rec) record3
                             Else
                                Call gsend(comm,record2,0,SysExpand_tag)
                                Call gsend(comm,record3,0,SysExpand_tag)
                             End If

                          End If

                       End If

                    End Do
                 End Do
              End Do

! Increase local atom counter

              indatm=indatm+1

           Else

! Determine sending node for UN/SORTED MASTER

              If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
                  io_write == IO_WRITE_SORTED_MASTER) Then
                 idm=0 ! Initialise node number
                 Call gsum(comm,idm)

                 Do iz=1,nz
                    Do iy=1,ny
                       Do ix=1,nx
                          rec   = offset + Int(2,li)*(Int(i_xyz(ix,iy,iz),li)*Int(setspc,li) + Int(m,li)) - Int(2,li)

                          If (comm%idnode == 0) Then
                             Call MPI_RECV(record2,recsz,MPI_CHARACTER,idm,SysExpand_tag,comm%comm,comm%status,comm%ierr)
                             Call MPI_RECV(record3,recsz,MPI_CHARACTER,idm,SysExpand_tag,comm%comm,comm%status,comm%ierr)

                             rec=rec+Int(1,li)
                             Write(Unit=nconf, Fmt='(a73)', Rec=rec) record2
                             rec=rec+Int(1,li)
                             Write(Unit=nconf, Fmt='(a73)', Rec=rec) record3
                          End If
                       End Do
                    End Do
                 End Do
              End If

           End If
        End Do

        sapmtt = sapmtt + numsit(itmols)
        offset = offset + Int(2,li)*Int(numsit(itmols),li)
     End Do

     sapmpt = sapmpt + sapmtt
     offset = offset + Int(2,li)*Int(nall-1,li)*Int(setspc,li)
  End Do

  If ((.not.safeg) .and. comm%idnode == 0) Write(nrite,Fmt='(/,1x,a)') &
     '*** warning - possible topological contiguity failures occurred !!! ***'

  If      (io_write == IO_WRITE_UNSORTED_MPIIO .or. &
           io_write == IO_WRITE_SORTED_MPIIO   .or. &
           io_write == IO_WRITE_SORTED_NETCDF) Then

     Call io_set_parameters( user_comm = comm%comm )
     Call io_init( recsz )
     Call io_open( io_write, comm%comm, fcfg(1:Len_Trim(fcfg)), MPI_MODE_WRONLY, fh )

     If (io_write /= IO_WRITE_SORTED_NETCDF) Then
        top_skip = Int(5,MPI_OFFSET_KIND)
     Else
        top_skip = Int(1,MPI_OFFSET_KIND) ! netCDF frame
     End If

     Call io_write_sorted_file( fh, 0, IO_RESTART, top_skip, at_scaled, &
          ltg_scaled, atmnam_scaled,                                    &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /),                     &
          x_scaled, y_scaled, z_scaled,                                 &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /),                     &
          (/ 0.0_wp /), (/ 0.0_wp /), (/ 0.0_wp /), ierr )

     If ( ierr /= 0 ) Then
        Select Case( ierr )
        Case( IO_BASE_COMM_NOT_SET )
           Call error( 1050 )
        Case( IO_ALLOCATION_ERROR )
           Call error( 1053 )
        Case( IO_UNKNOWN_WRITE_OPTION )
           Call error( 1056 )
        Case( IO_UNKNOWN_WRITE_LEVEL )
           Call error( 1059 )
        End Select
     End If

     Call io_close( fh )
     Call io_finalize

     Deallocate ( atmnam_scaled, ltg_scaled,    Stat = fail(1) )
     Deallocate ( x_scaled, y_scaled, z_scaled, Stat = fail(2) )
     If (Any(fail > 0)) Then
        Write(message,'(a)') 'system_expand deallocation failure 0 '
        Call error(0,message)
     End If

  Else If (io_write == IO_WRITE_UNSORTED_DIRECT .or. &
           io_write == IO_WRITE_SORTED_DIRECT) Then

     Call io_close( fh )
     Call io_finalize

  Else If (io_write == IO_WRITE_UNSORTED_MASTER .or. &
           io_write == IO_WRITE_SORTED_MASTER) Then

     Close(Unit=nconf)

  End If

  Call gtime(t)

! Write summary data and proceed with FIELD

  x=0.5_wp*Min(fx*celprp(7),fy*celprp(8),fz*celprp(9))

  If (comm%idnode == 0) Then
     Write(nrite,'(/,1x,3a)') '*** ', fcfg(1:Len_Trim(fcfg)), ' expansion completed !'
     Write(nrite,'(1x,a,i10,a)') '*** Size: ', nall*megatm, ' particles'
     Write(nrite,'(1x,a,f10.2,a)') '*** Maximum radius of cutoff: ', x, ' Angstroms'
     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec")') t
     Open(Unit=nfield, File='FIELD', Status='old')
     Open(Unit=nconf, File=ffld(1:Len_Trim(ffld)), Status='replace')
     Write(nrite,'(/,1x,2a)')'*** Expanding FIELD in file ', ffld(1:Len_Trim(ffld))
     Write(nrite,'(1x,a)') '***'

! omit first line

     record=' '
     Read(Unit=nfield, Fmt='(a)', End=10) record
     Call tabs_2_blanks(record) ; Call strip_blanks(record)
     Write(nconf,'(a)') record(1:Len_Trim(record))

! read and process directives from field file

     Do
        record=' '
        Read(Unit=nfield, Fmt='(a)', End=10) record
        Call tabs_2_blanks(record) ; Call strip_blanks(record)
        record1=record
        Call get_word(record,word)
        Call lower_case(word)

        If      (word(1:5) == 'multi') Then

! MPOLES should exist

           lmpldt=.true.

! number of molecules of this type

        Else If (word(1:6) == 'nummol') Then

           Call get_word(record,word)
           index=Nint(word_2_real(word))
           Call get_word(record1,word)
           Write(Unit=nconf,Fmt='(a,i10)') word(1:Len_Trim(word)), nall*index

! close force field file

        Else If (word(1:5) == 'close') Then

           Call gtime(t)

           Write(Unit=nconf,Fmt='(a)') record1(1:Len_Trim(record1))
           Close(Unit=nfield)
           Close(Unit=nconf)
           Write(nrite,'(1x,3a)') '*** ', ffld(1:Len_Trim(ffld)), ' expansion done !'
           Exit

! just paste the copy

        Else

           Write(nconf,'(a)') record1(1:Len_Trim(record1))

        End If
     End Do

10   Continue

     If (lmpldt) Inquire(File='MPOLES', Exist=lmpldt)
     If (lmpldt) Then
        Open(Unit=nmpldt, File='MPOLES', Status='old')
        Open(Unit=nconf, File=fmpl(1:Len_Trim(fmpl)), Status='replace')
        Write(nrite,'(/,1x,2a)')'*** Expanding MPOLES in file ', fmpl(1:Len_Trim(fmpl))
        Write(nrite,'(1x,a)') '***'

! omit first line

        record=' '
        Read(Unit=nmpldt, Fmt='(a)', End=10) record
        Call tabs_2_blanks(record) ; Call strip_blanks(record)
        Write(nconf,'(a)') record(1:Len_Trim(record))

! read and process directives from mpoles file

        Do
           record=' '
           Read(Unit=nmpldt, Fmt='(a)', End=20) record
           Call tabs_2_blanks(record) ; Call strip_blanks(record)
           record1=record
           Call get_word(record,word)
           Call lower_case(word)

           If      (word(1:6) == 'nummol') Then

              Call get_word(record,word)
              index=Nint(word_2_real(word))
              Call get_word(record1,word)
              Write(Unit=nconf,Fmt='(a,i10)') word(1:Len_Trim(word)), nall*index

! close mpoles file

           Else If (word(1:5) == 'close') Then

              Call gtime(t)

              Write(Unit=nconf,Fmt='(a)') record1(1:Len_Trim(record1))
              Close(Unit=nmpldt)
              Close(Unit=nconf)
              Write(nrite,'(1x,3a)') '*** ', fmpl(1:Len_Trim(fmpl)), ' expansion done !'
              Exit

! just paste the copy

           Else

              Write(nconf,'(a)') record1(1:Len_Trim(record1))

           End If
        End Do
     End If

20   Continue

     Write(nrite,'(/,1x, "time elapsed since job start: ", f12.3, " sec",/)') t
     Write(nrite,'(1x,a)') '*** Simulation continues as scheduled...'
  End If
  Call gsync(comm)

  Deallocate (f1,f2,f3, Stat=fail(1))
  Deallocate (f4,f5,f6, Stat=fail(2))
  Deallocate (f7,f8,f9, Stat=fail(3))
  Deallocate (i_xyz,    Stat=fail(4))
  Deallocate (xm,ym,zm, Stat=fail(5))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a)') 'system_expand dellocation failure '
     Call error(0,message)
  End If

End Subroutine system_expand

Subroutine system_revive                                      &
           (rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing restart files at job termination or
! selected intervals in simulation
!
! copyright - daresbury laboratory
! author    - w.smith december 1992
! amended   - i.t.todorov november 2016
! contrib   - m.a.seaton june 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer,           Intent( In    ) :: megatm,nstep
  Logical,           Intent( In    ) :: lrdf,lzdn
  Real( Kind = wp ), Intent( In    ) :: rcut,rbin,tstep,time,tmst, &
                                        chit,cint,chip,eta(1:9),   &
                                        strcon(1:9),strpmf(1:9),stress(1:9)
  Type( comms_type ), Intent( InOut ) :: comm

  Logical               :: ready
  Character( Len = 42 )  :: name
  Character( Len = 40 ) :: forma  = ' '

  Integer               :: fail(1:3),i,j,l,levcfg,jdnode,jatms,nsum
  Real( Kind = wp )     :: r_mxnode

  Integer,           Dimension( : ), Allocatable :: iwrk
  Real( Kind = wp ), Dimension( : ), Allocatable :: axx,ayy,azz
  Real( Kind = wp ), Dimension( : ), Allocatable :: bxx,byy,bzz
  Character ( Len = 256 )  :: message 

  fail=0
  Allocate (iwrk(1:mxatms),                            Stat=fail(1))
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
  Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'system_revive allocation failure '
     Call error(0,message)
  End If


! Define format for REVIVE printing in ASCII

  If (l_rout) Then
     i = 64/4 - 1 ! Bit_Size(0.0_wp)/4 - 1
     j = Max(mxstak*mxnstk+1,mxgrdf*mxrdf,mxgusr,mxgana*mxtana)

     Write(forma ,10) j/4+1,i+9,i
10   Format('(1p,',i0,'(/,4e',i0,'.',i0,'E3))')
  End If

  If (comm%mxnode > 1) Then

! globally sum RDF information before saving

     If (lrdf) Then

! maximum rdfs that can be summed in each step

        nsum = mxbuff/mxgrdf
        If (nsum == 0) Call error(200)

        Do i=1,mxrdf,nsum
           Call gsum(comm,rdf(:,i:Min(i+nsum-1,mxrdf)))
        End Do

     End If

! globally sum USR RDF information before saving

     If (mxgusr > 0) Call gsum(comm,usr(1:mxgusr))

! globally sum z-density information before saving

     If (lzdn) Then

! maximum zdens that can be summed in each step

        nsum = mxbuff/mxgrdf
        If (nsum == 0) Call error(200)

        Do i=1,mxatyp,nsum
           Call gsum(comm,zdens(:,i:Min(i+nsum-1,mxatyp)))
        End Do

     End If

! globally sum vafdata information before saving

     If (vafsamp > 0) Then

! maximum vafdata that can be summed in each step

        nsum = mxbuff/(nsvaf+1)
        If (nsum == 0) Call error(200)

        Do j=1,vafsamp
           l=(j-1)*(mxatyp+1) ! avoid summing up timing information
           Do i=1,mxatyp,nsum
              Call gsum(comm,vafdata(:,l+i:l+Min(i+nsum-1,mxatyp)))
           End Do
        End Do

     End If

! globally sum bonds' distributions information before saving

     If (mxgbnd1 > 0) Then

! maximum dstbnd that can be summed in each step

        nsum = mxbuff/(mxgbnd1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfbnd(0),nsum
           Call gsum(comm,dstbnd(:,i:Min(i+nsum-1,ldfbnd(0))))
        End Do

     End If

! globally sum angles' distributions information before saving

     If (mxgang1 > 0) Then

! maximum dstang that can be summed in each step

        nsum = mxbuff/(mxgang1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfang(0),nsum
           Call gsum(comm,dstang(:,i:Min(i+nsum-1,ldfang(0))))
        End Do

     End If

! globally sum dihedrals' distributions information before saving

     If (mxgdih1 > 0) Then

! maximum dstdih that can be summed in each step

        nsum = mxbuff/(mxgdih1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfdih(0),nsum
           Call gsum(comm,dstdih(:,i:Min(i+nsum-1,ldfdih(0))))
        End Do

     End If

! globally sum inversions' distributions information before saving

     If (mxginv1 > 0) Then

! maximum dstinv that can be summed in each step

        nsum = mxbuff/(mxginv1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfinv(0),nsum
           Call gsum(comm,dstinv(:,i:Min(i+nsum-1,ldfinv(0))))
        End Do

     End If

  End If

! Write REVCON

  name = Trim(revcon) ! file name
  levcfg = 2      ! define level of information in REVCON

  Call write_config(name,levcfg,megatm,nstep,tstep,time,comm)

! node 0 handles I/O

  If (comm%idnode == 0) Then

! Write accumulator data to dump file

     If (l_rout) Then
        Open(Unit=nrest, File=Trim(revive), Form='formatted', Status='replace')

        Write(Unit=nrest, Fmt=forma, Advance='No') rcut,rbin,Real(megatm,wp)
        Write(Unit=nrest, Fmt=forma, Advance='No') &
             Real(nstep,wp),tstep,time,tmst,Real(numacc,wp),chit,chip,cint
        Write(Unit=nrest, Fmt=forma, Advance='No') eta
        Write(Unit=nrest, Fmt=forma, Advance='No') stpval
        Write(Unit=nrest, Fmt=forma, Advance='No') stpvl0
        Write(Unit=nrest, Fmt=forma, Advance='No') sumval
        Write(Unit=nrest, Fmt=forma, Advance='No') ssqval
        Write(Unit=nrest, Fmt=forma, Advance='No') zumval
        Write(Unit=nrest, Fmt=forma, Advance='No') ravval
        Write(Unit=nrest, Fmt=forma, Advance='No') stkval
        Write(Unit=nrest, Fmt=forma, Advance='No') strcon
        Write(Unit=nrest, Fmt=forma, Advance='No') strpmf
        Write(Unit=nrest, Fmt=forma, Advance='No') stress

        If (lrdf) Write(Unit=nrest, Fmt=forma, Advance='No') Real(ncfrdf,wp),rdf
        If (mxgusr > 0) Write(Unit=nrest, Fmt=forma, Advance='No') Real(mxgusr),rusr,Real(ncfusr,wp),usr
        If (lzdn) Write(Unit=nrest, Fmt=forma, Advance='No') Real(ncfzdn,wp),zdens
        If (vafsamp > 0) Then
          Write(Unit=nrest, Fmt=forma, Advance='No') vafcount
          Write(Unit=nrest, Fmt=forma, Advance='No') Real(vafstep,wp)
          Write(Unit=nrest, Fmt=forma, Advance='No') vafdata
          Write(Unit=nrest, Fmt=forma, Advance='No') vaf
          Write(Unit=nrest, Fmt=forma, Advance='No') vaftime
        End If

        If (mxgbnd1 > 0) Write(Unit=nrest, Fmt=forma, Advance='No') Real(ncfbnd,wp),dstbnd
        If (mxgang1 > 0) Write(Unit=nrest, Fmt=forma, Advance='No') Real(ncfang,wp),dstang
        If (mxgdih1 > 0) Write(Unit=nrest, Fmt=forma, Advance='No') Real(ncfdih,wp),dstdih
        If (mxginv1 > 0) Write(Unit=nrest, Fmt=forma, Advance='No') Real(ncfinv,wp),dstinv
     Else
        Open(Unit=nrest, File=Trim(revive), Form='unformatted', Status='replace')

        Write(Unit=nrest) rcut,rbin,Real(megatm,wp)
        Write(Unit=nrest) &
             Real(nstep,wp),tstep,time,tmst,Real(numacc,wp),chit,chip,cint
        Write(Unit=nrest) eta
        Write(Unit=nrest) stpval
        Write(Unit=nrest) stpvl0
        Write(Unit=nrest) sumval
        Write(Unit=nrest) ssqval
        Write(Unit=nrest) zumval
        Write(Unit=nrest) ravval
        Write(Unit=nrest) stkval
        Write(Unit=nrest) strcon
        Write(Unit=nrest) strpmf
        Write(Unit=nrest) stress

        If (lrdf) Write(Unit=nrest) Real(ncfrdf,wp),rdf
        If (mxgusr > 0) Write(Unit=nrest) Real(mxgusr),rusr,Real(ncfusr,wp),usr
        If (lzdn) Write(Unit=nrest) Real(ncfzdn,wp),zdens
        If (vafsamp > 0) Then
          Write(Unit=nrest) vafcount
          Write(Unit=nrest) Real(vafstep,wp)
          Write(Unit=nrest) vafdata
          Write(Unit=nrest) vaf
          Write(Unit=nrest) vaftime
        End If

        If (mxgbnd1 > 0) Write(Unit=nrest) Real(ncfbnd,wp),dstbnd
        If (mxgang1 > 0) Write(Unit=nrest) Real(ncfang,wp),dstang
        If (mxgdih1 > 0) Write(Unit=nrest) Real(ncfdih,wp),dstdih
        If (mxginv1 > 0) Write(Unit=nrest) Real(ncfinv,wp),dstinv
     End If

! Write initial position and final displacement data to REVIVE

     jatms=natms

     Do i=1,natms
        iwrk(i)=ltg(i)

        axx(i)=xin(i)
        ayy(i)=yin(i)
        azz(i)=zin(i)

        bxx(i)=xto(i)
        byy(i)=yto(i)
        bzz(i)=zto(i)
     End Do

     ready=.true.
     Do jdnode=0,comm%mxnode-1
        If (jdnode > 0) Then
           Call gsend(comm,ready,jdnode,Revive_tag)

           Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)

           Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)

           Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
           Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
           Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)

           Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
           Call MPI_RECV(byy,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
           Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
        End If

        If (l_rout) Then
           Do i=1,jatms
              Write(Unit=nrest, Fmt=forma, Advance='No') &
                   Real(iwrk(i),wp),axx(i),ayy(i),azz(i),bxx(i),byy(i),bzz(i)
           End Do
        Else
           Do i=1,jatms
              Write(Unit=nrest) Real(iwrk(i),wp),axx(i),ayy(i),azz(i),bxx(i),byy(i),bzz(i)
           End Do
        End If
     End Do

  Else

     Call MPI_RECV(ready,1,MPI_LOGICAL,0,Revive_tag,comm%comm,comm%status,comm%ierr)

     Call gsend(comm,natms,0,Revive_tag)

     Call gsend(comm,ltg(1:natms),0,Revive_tag)

     Call gsend(comm,xin(1:natms),0,Revive_tag)
     Call gsend(comm,yin(1:natms),0,Revive_tag)
     Call gsend(comm,zin(1:natms),0,Revive_tag)

     Call gsend(comm,xto(1:natms),0,Revive_tag)
     Call gsend(comm,yto(1:natms),0,Revive_tag)
     Call gsend(comm,zto(1:natms),0,Revive_tag)

  End If

! Write initial velocities for VAF calculations if needed

  If (vafsamp > 0) Then

    If (comm%idnode == 0) Then

       jatms=natms
       Do j=1,vafsamp
         Do i=1,natms
            iwrk(i)=ltg(i)
            axx(i)=vxi(i,j)
            ayy(i)=vyi(i,j)
            azz(i)=vzi(i,j)
         End Do

         ready=.true.
         Do jdnode=0,comm%mxnode-1
            If (jdnode > 0) Then
               Call gsend(comm,ready,jdnode,Revive_tag)

               Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)

               Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)

               Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
               Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
               Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Revive_tag,comm%comm,comm%status,comm%ierr)
            End If


            If (l_rout) Then
               Do i=1,jatms
                  Write(Unit=nrest, Fmt=forma, Advance='No') &
                       Real(iwrk(i),wp),axx(i),ayy(i),azz(i)
               End Do
            Else
               Do i=1,jatms
                  Write(Unit=nrest) Real(iwrk(i),wp),axx(i),ayy(i),azz(i)
               End Do
            End If
         End Do
       End Do

    Else

       Do j=1,vafsamp
         Call MPI_RECV(ready,1,MPI_LOGICAL,0,Revive_tag,comm%comm,comm%status,comm%ierr)

         Call gsend(comm,natms,0,Revive_tag)

         Call gsend(comm,ltg(1:natms),0,Revive_tag)

         Call gsend(comm,vxi(1,j),0,Revive_tag)
         Call gsend(comm,vyi(1,j),0,Revive_tag)
         Call gsend(comm,vzi(1,j),0,Revive_tag)
       End Do

    End If

    Call gsync(comm)

  End If

  Call gsync(comm)
  If (comm%idnode == 0) Close(Unit=nrest)

       r_mxnode=1.0_wp/Real(comm%mxnode,wp)

! globally divide rdf data between nodes

     If (lrdf) rdf = rdf * r_mxnode

! globally divide USR rdf data between nodes

     If (mxgusr > 0) usr = usr * r_mxnode

! globally divide z-density data between nodes

     If (lzdn) zdens = zdens * r_mxnode

! globally divide bonds' distributions data between nodes

     If (mxgbnd1 > 0) dstbnd = dstbnd * r_mxnode

! globally divide angles' distributions data between nodes

     If (mxgang1 > 0) dstang = dstang * r_mxnode

! globally divide dihedrals' distributions data between nodes

     If (mxgdih1 > 0) dstdih = dstdih * r_mxnode

! globally divide inversions' distributions data between nodes

     If (mxginv1 > 0) dstinv = dstinv * r_mxnode

  Deallocate (iwrk,        Stat=fail(1))
  Deallocate (axx,ayy,azz, Stat=fail(2))
  Deallocate (bxx,byy,bzz, Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'system_revive deallocation failure'
     Call error(0,message)
  End If

End Subroutine system_revive



End Module system
