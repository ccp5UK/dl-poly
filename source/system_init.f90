Subroutine system_init                                       &
           (levcfg,imcon,rcut,rvdw,rbin,rmet,                &
           lrdf,lzdn,keyres,megatm,                          &
           time,tmst,nstep,chit,cint,chip,eta,virtot,stress, &
           vircon,strcon,virpmf,strpmf,elrc,virlrc,elrcm,vlrcm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for reading the REVIVE file data and defining the
! initial thermodynamic and structural accumulators
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use site_module,   Only : ntpatm,numtyp,dens
  Use config_module, Only : volm,natms,ltg,ltype,xxx,yyy,zzz
  Use langevin_module
  Use vdw_module,    Only : ls_vdw,ntpvdw
  Use metal_module,  Only : ntpmet

  Use statistics_module
  Use development_module

  Implicit None

  Logical,           Intent( In    ) :: lrdf,lzdn
  Integer,           Intent( InOut ) :: levcfg,keyres
  Integer,           Intent( In    ) :: imcon,megatm
  Real( Kind = wp ), Intent( In    ) :: rcut,rvdw,rbin,rmet

  Integer,           Intent(   Out ) :: nstep
  Real( Kind = wp ), Intent(   Out ) :: time,tmst,chit,cint,chip,eta(1:9),     &
                                        virtot,stress(1:9),vircon,strcon(1:9), &
                                        virpmf,strpmf(1:9),elrc,virlrc,        &
                                        elrcm(0:mxatyp),vlrcm(0:mxatyp)

  Logical              , Save :: newjob = .true.
  Character( Len = 40 ), Save :: forma  = ' '

  Logical           :: l_tmp
  Integer           :: i,j,k,keyio,i_tmp,gidx
  Real( Kind = wp ) :: dnstep,dnumac,dnumrd,dnumzd,r_mxnode,xyz(0:6)


  If (newjob .and. l_rin) Then
     newjob = .false.

! Define format for REVOLD reading in ASCII

     i = 64/4 - 1 ! Bit_Size(0.0_wp)/4 - 1

     Write(forma ,10) Max(mxgrdf*mxrdf,mxstak*mxnstk)/4+1,i+9,i
10   Format('(1p,',i0,'(/,4e',i0,'.',i0,'E3))')
  End If

! Initialise read failure flag

  keyio=0

50 Continue
  If (keyres /= 1 .or. idnode /= 0) Then

! initialise step counters

     numacc=0
     numrdf=0
     numzdn=0

     nstep =0
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

! initialise accumulator arrays if reading failure occured

     If (keyio > 0) Then

        stpval=0.0_wp
        stpvl0=0.0_wp
        sumval=0.0_wp
        ssqval=0.0_wp
        zumval=0.0_wp
        ravval=0.0_wp

        stkval=0.0_wp

        If (lrdf) rdf=0.0_wp
        If (lzdn) zdens=0.0_wp

     End If

  End If

! restart simulation and continue

  If (keyres == 1) Then

! If REVOLD doesn't exist then abort (mishmashed REVOLD is handled separately)

     l_tmp=.true.
     If (idnode == 0) Inquire(File='REVOLD', Exist=l_tmp)
     If (mxnode > 1) Call gcheck(l_tmp)
     If (.not.l_tmp) Call error(519)

! Check REVOLD restart compatibility: rcut,rbin,megatm

     xyz(1:3)=0.0_wp
     If (idnode == 0) Then
        If (l_rin) Then
           Open(Unit=nrest, file='REVOLD', form='formatted', IOStat=keyio)
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) xyz(1),xyz(2),xyz(3)
        Else
           Open(Unit=nrest, file='REVOLD', form='unformatted', IOStat=keyio)
           Read(Unit=nrest, IOStat=keyio, End=100) xyz(1),xyz(2),xyz(3)
        End If
     End If
     If (mxnode > 1) Call MPI_BCAST(xyz(1:3), 3, wp_mpi, 0, dlp_comm_world, ierr)
     If (Abs(xyz(1)-rcut) > 1.0e-6_wp .or. Abs(xyz(2)-rbin) > 1.0e-6_wp .or. &
         Nint(xyz(3)) /= megatm) Call error(519)

! read the rest of the accumulator data from dump file

     If (idnode == 0) Then
        If (l_rin) Then
           Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) &
               dnstep,dnumac,dnumrd,dnumzd,time,tmst,chit,chip,cint
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

           If (lrdf) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) rdf
           If (lzdn) Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio, End=100) zdens
        Else
           Read(Unit=nrest, IOStat=keyio, End=100) &
               dnstep,dnumac,dnumrd,dnumzd,time,tmst,chit,chip,cint
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

           If (lrdf) Read(Unit=nrest, IOStat=keyio, End=100) rdf
           If (lzdn) Read(Unit=nrest, IOStat=keyio, End=100) zdens
        End If

        nstep =Nint(dnstep)
        numacc=Nint(dnumac)
        numrdf=Nint(dnumrd)
        numzdn=Nint(dnumzd)

! calculate virtot = virtot-vircon-virpmf

        vircon = stpval(17) * engunit
        virpmf = stpval(26) * engunit
        virtot = (stpval(12)-stpval(17)-stpval(26)) * engunit
     End If

100  Continue

! If 'restart' is impossible go to 'restart noscale' and reinitialise

     If (mxnode > 1) Call MPI_BCAST(keyio, 1, MPI_INTEGER, 0, dlp_comm_world, ierr)
     If (keyio /= 0) Then
        If (idnode == 0) Then
           Call warning(190,0.0_wp,0.0_wp,0.0_wp)
           Close(Unit=nrest)
        End If
        keyres=3
        Go To 50
     End If

! broadcast stored variables

     If (mxnode > 1) Then
        Call MPI_BCAST(nstep,                 1, MPI_INTEGER, 0, dlp_comm_world, ierr)
        Call MPI_BCAST(numacc,                1, MPI_INTEGER, 0, dlp_comm_world, ierr)
        Call MPI_BCAST(numrdf,                1, MPI_INTEGER, 0, dlp_comm_world, ierr)
        Call MPI_BCAST(numzdn,                1, MPI_INTEGER, 0, dlp_comm_world, ierr)
        Call MPI_BCAST(time,                  1, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(tmst,                  1, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(chit,                  1, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(chip,                  1, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(cint,                  1, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(eta(1:9),              9, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(stpval(1:mxnstk), mxnstk, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(stpvl0(1:mxnstk), mxnstk, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(sumval(1:mxnstk), mxnstk, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(ssqval(1:mxnstk), mxnstk, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(zumval(1:mxnstk), mxnstk, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(ravval(1:mxnstk), mxnstk, wp_mpi,      0, dlp_comm_world, ierr)
        Do k=1,mxnstk
           Call MPI_BCAST(stkval(1:mxstak,k), mxstak, wp_mpi,      0, dlp_comm_world, ierr)
        End Do
        Call MPI_BCAST(strcon(1:9),           9, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(strpmf(1:9),           9, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(stress(1:9),           9, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(vircon,                1, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(virpmf,                1, wp_mpi,      0, dlp_comm_world, ierr)
        Call MPI_BCAST(virtot,                1, wp_mpi,      0, dlp_comm_world, ierr)

! rdf table - broadcast and normalise

        r_mxnode=1.0_wp/Real(mxnode,wp)
        If (lrdf) Then
           Do k=1,mxrdf
              Call MPI_BCAST(rdf(1:mxgrdf,k), mxgrdf, wp_mpi, 0, dlp_comm_world, ierr)

              Do j=1,mxgrdf
                 rdf(j,k) = rdf(j,k)*r_mxnode
              End Do
           End Do
        End If

! z-density table - broadcast and normalise

        If (lzdn) Then
           Do k=1,mxatyp
              Call MPI_BCAST(zdens(1:mxgrdf,k), mxgrdf, wp_mpi, 0, dlp_comm_world, ierr)

              Do j=1,mxgrdf
                 zdens(j,k) = zdens(j,k)*r_mxnode
              End Do
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

        If (idnode == 0) Then
           If (l_rin) Then
              Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio) &
                  xyz(0),xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6)
           Else
              Read(Unit=nrest, IOStat=keyio) &
                  xyz(0),xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6)
           End If
        End If
        If (keyio /= 0) i_tmp=1

        If (mxnode > 1) Call MPI_BCAST(xyz(0:6), 7, wp_mpi, 0, dlp_comm_world, ierr)
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

     If (mxnode > 1) Call MPI_BCAST(i_tmp, 1, MPI_INTEGER, 0, dlp_comm_world, ierr)
     If (i_tmp /= 0) Then
        If (idnode == 0) Then
           Call warning(190,0.0_wp,0.0_wp,0.0_wp)
           Close(Unit=nrest)
        End If
        keyres=3
        Go To 50
     End If

! Read Langevin arrays if needed

     If (l_lan) Then

        i_tmp=0 ! Error accumulator: keyio is still zero otherwise we cannot get here

        Do k=1,megatm
           xyz(0:3)=0.0_wp

           If (idnode == 0) Then
              If (l_rin) Then
                 Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio) &
                     xyz(0),xyz(1),xyz(2),xyz(3)
              Else
                 Read(Unit=nrest, IOStat=keyio) &
                     xyz(0),xyz(1),xyz(2),xyz(3)
              End If
           End If
           If (keyio /= 0) i_tmp=1

           If (mxnode > 1) Call MPI_BCAST(xyz(0:3), 4, wp_mpi, 0, dlp_comm_world, ierr)
           gidx=Nint(xyz(0))

! assign random particle forces to the corresponding domains

           Do i=1,natms
              If (ltg(i) == gidx) Then
                 fxl(i)=xyz(1)
                 fyl(i)=xyz(2)
                 fzl(i)=xyz(3)
              End If
           End Do
        End Do
        If (idnode == 0) Then
           If (l_rin) Then
              Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio) fpl(1:9)
           Else
              Read(Unit=nrest, IOStat=keyio) fpl(1:9)
           End If
        End If
        If (keyio /= 0) i_tmp=1
        If (mxnode > 1) Call MPI_BCAST(fpl(1:9), 9, wp_mpi, 0, dlp_comm_world, ierr)

! If anything is wrong

        If (mxnode > 1) Call MPI_BCAST(i_tmp, 1, MPI_INTEGER, 0, dlp_comm_world, ierr)
        If (i_tmp /= 0) Then
           If (idnode == 0) Then
              Call warning(190,0.0_wp,0.0_wp,0.0_wp)
              Close(Unit=nrest)
           End If
           l_lan_s=.true.
           Go To 50
        Else
           l_lan_s=.false.
        End If

     End If

! Read Langevin process gaussian variable if needed

     If (l_gst) Then

        i_tmp=0 ! Error accumulator: keyio is still zero otherwise we cannot get here

        If (idnode == 0) Then
           If (l_rin) Then
              Read(Unit=nrest, Fmt=forma, Advance='No', IOStat=keyio) r_0
           Else
              Read(Unit=nrest, IOStat=keyio) r_0
           End If
        End If
        If (keyio /= 0) i_tmp=1

        If (mxnode > 1) Call MPI_BCAST(r_0, 1, wp_mpi, 0, dlp_comm_world, ierr)

! If anything is wrong

        If (mxnode > 1) Call MPI_BCAST(i_tmp, 1, MPI_INTEGER, 0, dlp_comm_world, ierr)
        If (i_tmp /= 0) Then
           If (idnode == 0) Then
              Call warning(190,0.0_wp,0.0_wp,0.0_wp)
              Close(Unit=nrest)
           End If
           l_gst_s=.true.
           Go To 50
        Else
           l_gst_s=.false.
        End If

     End If

     If (idnode == 0) Close(Unit=nrest)

  Else

! force force and stress recalculation when 'restart' is not on

     levcfg=1

  End If


! number densities needed for long-range corrections

! evaluate species populations in system

  Do i=1,natms
     k = ltype(i)
     numtyp(k) = numtyp(k)+1.0_wp
  End Do

! global number densities

  If (mxnode > 1) Call gsum(numtyp(1:ntpatm))

! number densities

  Do i=1,ntpatm
     If (numtyp(i) > zero_plus) dens(i) = numtyp(i)/volm
  End Do

! Get long-range corrections

! elrc & virlrc arrays are zeroed in vdw_module,
! no lrc when vdw interactions are force-shifted

  If (ntpvdw > 0 .and. (.not.ls_vdw)) Call vdw_lrc(imcon,rvdw,elrc,virlrc)

! elrcm & vlrcm arrays are zeroed in metal_module

  If (ntpmet > 0) Call metal_lrc(imcon,rmet,elrcm,vlrcm)

End Subroutine system_init
