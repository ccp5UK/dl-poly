Subroutine system_revive                                            &
           (imcon,rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing restart files at job termination or
! selected intervals in simulation
!
! copyright - daresbury laboratory
! author    - w.smith december 1992
! amended   - i.t.todorov january 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module
  Use setup_module
  Use config_module, Only : natms,ltg
  Use langevin_module
  Use statistics_module
  Use development_module

  Implicit None

  Integer,           Intent( In    ) :: imcon,megatm,nstep
  Logical,           Intent( In    ) :: lrdf,lzdn
  Real( Kind = wp ), Intent( In    ) :: rcut,rbin,tstep,time,tmst, &
                                        chit,cint,chip,eta(1:9),   &
                                        strcon(1:9),strpmf(1:9),stress(1:9)

  Logical              , Save :: newjob = .true.
  Character( Len = 40 ), Save :: forma  = ' '

  Logical              :: ready
  Character( Len = 6 ) :: name

  Integer              :: fail(1:3),i,levcfg,jdnode,jatms,nsum

  Integer,           Dimension( : ), Allocatable :: iwrk
  Real( Kind = wp ), Dimension( : ), Allocatable :: axx,ayy,azz
  Real( Kind = wp ), Dimension( : ), Allocatable :: bxx,byy,bzz

  fail=0
  Allocate (iwrk(1:mxatms),                            Stat=fail(1))
  Allocate (axx(1:mxatms),ayy(1:mxatms),azz(1:mxatms), Stat=fail(2))
  Allocate (bxx(1:mxatms),byy(1:mxatms),bzz(1:mxatms), Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'system_revive allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob .and. l_rout) Then
     newjob = .false.

! Define format for REVIVE printing in ASCII

     i = 64/4 - 1 ! Bit_Size(0.0_wp)/4 - 1

     Write(forma ,10) Max(mxgrdf*mxrdf,mxstak*mxnstk)/4+1,i+9,i
10   Format('(1p,',i0,'(/,4e',i0,'.',i0,'E3))')
  End If

  If (mxnode > 1) Then

! globally sum rdf information before saving

     If (lrdf) Then

! maximum rdfs that can be summed in each step

        nsum = mxbuff/mxgrdf
        If (nsum == 0) Call error(200)

        Do i=1,mxrdf,nsum
           Call gsum(rdf(:,i:Min(i+nsum-1,mxrdf)))
        End Do
     End If

! globally sum zden information before saving

     If (lzdn) Then

! maximum rdfs that can be summed in each step

        nsum = mxbuff/mxgrdf
        If (nsum == 0) Call error(200)

        Do i=1,mxatyp,nsum
           Call gsum(zdens(:,i:Min(i+nsum-1,mxatyp)))
        End Do
     End If
  End If

! Write REVCON

  name = 'REVCON' ! file name
  levcfg = 2      ! define level of information in REVCON

  Call write_config(name,imcon,levcfg,megatm,nstep,tstep,time)

! node 0 handles I/O

  If (idnode == 0) Then

! Write accumulator data to dump file

     If (l_rout) Then
        Open(Unit=nrest, File='REVIVE', Form='formatted', Status='replace')

        Write(Unit=nrest, Fmt=forma, Advance='No') rcut,rbin,Real(megatm,wp)
        Write(Unit=nrest, Fmt=forma, Advance='No')           &
             Real(nstep,wp),Real(numacc,wp),Real(numrdf,wp), &
             Real(numzdn,wp),time,tmst,chit,chip,cint
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

        If (lrdf) Write(Unit=nrest, Fmt=forma, Advance='No') rdf
        If (lzdn) Write(Unit=nrest, Fmt=forma, Advance='No') zdens
     Else
        Open(Unit=nrest, File='REVIVE', Form='unformatted', Status='replace')

        Write(Unit=nrest) rcut,rbin,Real(megatm,wp)
        Write(Unit=nrest) Real(nstep,wp),Real(numacc,wp),Real(numrdf,wp), &
                          Real(numzdn,wp),time,tmst,chit,chip,cint
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

        If (lrdf) Write(Unit=nrest) rdf
        If (lzdn) Write(Unit=nrest) zdens
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
     Do jdnode=0,mxnode-1
        If (jdnode > 0) Then
           Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Revive_tag,dlp_comm_world,ierr)

           Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Revive_tag,dlp_comm_world,status,ierr)

           Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Revive_tag,dlp_comm_world,status,ierr)

           Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
           Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
           Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)

           Call MPI_RECV(bxx,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
           Call MPI_RECV(byy,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
           Call MPI_RECV(bzz,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
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

     Call MPI_RECV(ready,1,MPI_LOGICAL,0,Revive_tag,dlp_comm_world,status,ierr)

     Call MPI_SEND(natms,1,MPI_INTEGER,0,Revive_tag,dlp_comm_world,ierr)

     Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Revive_tag,dlp_comm_world,ierr)

     Call MPI_SEND(xin,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
     Call MPI_SEND(yin,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
     Call MPI_SEND(zin,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)

     Call MPI_SEND(xto,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
     Call MPI_SEND(yto,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
     Call MPI_SEND(zto,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)

  End If

  If (mxnode > 1) Call gsync()

! Write Langevin arrays if needed

  If (l_lan) Then
     If (idnode == 0) Then
        jatms=natms

        Do i=1,natms
           iwrk(i)=ltg(i)

           axx(i)=fxl(i)
           ayy(i)=fyl(i)
           azz(i)=fzl(i)
        End Do

        ready=.true.
        Do jdnode=0,mxnode-1
           If (jdnode > 0) Then
              Call MPI_SEND(ready,1,MPI_LOGICAL,jdnode,Revive_tag,dlp_comm_world,ierr)

              Call MPI_RECV(jatms,1,MPI_INTEGER,jdnode,Revive_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(iwrk,jatms,MPI_INTEGER,jdnode,Revive_tag,dlp_comm_world,status,ierr)

              Call MPI_RECV(axx,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(ayy,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
              Call MPI_RECV(azz,jatms,wp_mpi,jdnode,Revive_tag,dlp_comm_world,status,ierr)
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

     Else

        Call MPI_RECV(ready,1,MPI_LOGICAL,0,Revive_tag,dlp_comm_world,status,ierr)

        Call MPI_SEND(natms,1,MPI_INTEGER,0,Revive_tag,dlp_comm_world,ierr)

        Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Revive_tag,dlp_comm_world,ierr)

        Call MPI_SEND(fxl,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
        Call MPI_SEND(fyl,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
        Call MPI_SEND(fzl,natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)

     End If

     If (mxnode > 1) Call gsync()
     If (idnode == 0) Then
        If (l_rout) Then
           Write(Unit=nrest, Fmt=forma, Advance='No') fpl(1:9)
        Else
           Write(Unit=nrest) fpl(1:9)
        End If
     End If
  End If

  If (mxnode > 1) Call gsync()

! Write Langevin process gaussian variable if needed

  If (l_gst) Then
     If (idnode == 0) Then
        If (l_rout) Then
           Write(Unit=nrest, Fmt=forma, Advance='No') r_0
        Else
           Write(Unit=nrest) r_0
        End If
     End If
  End If

  If (idnode == 0) Close(Unit=nrest)
  If (mxnode > 1) Call gsync()

! divide rdf data between nodes

  If (lrdf) rdf = rdf / Real(mxnode,wp)

! divide zdensity data between nodes

  If (lzdn) zdens = zdens / Real(mxnode,wp)

  Deallocate (iwrk,        Stat=fail(1))
  Deallocate (axx,ayy,azz, Stat=fail(2))
  Deallocate (bxx,byy,bzz, Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'system_revive deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine system_revive
