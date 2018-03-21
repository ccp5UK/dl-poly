Subroutine system_revive                                      &
           (rcut,rbin,lrdf,lzdn,megatm,nstep,tstep,time,tmst, &
           chit,cint,chip,eta,strcon,strpmf,stress)

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

  Use kinds, only : wp
  Use comms_module
  Use setup_module
  Use configuration,      Only : natms,ltg
  Use statistics_module
  Use rdf_module,         Only : ncfrdf,rdf,ncfusr,rusr,usr
  Use z_density_module,   Only : ncfzdn,zdens
  Use bonds_module,       Only : ldfbnd,ncfbnd,dstbnd
  Use angles_module,      Only : ldfang,ncfang,dstang
  Use dihedrals,   Only : ldfdih,ncfdih,dstdih
  Use inversions,  Only : ldfinv,ncfinv,dstinv
  Use greenkubo_module,   Only : nsvaf,vafsamp,vafcount,vafstep, &
                                 vxi,vyi,vzi,vafdata,vaf,vaftime

  Use development_module, Only : l_rout

  Implicit None

  Integer,           Intent( In    ) :: megatm,nstep
  Logical,           Intent( In    ) :: lrdf,lzdn
  Real( Kind = wp ), Intent( In    ) :: rcut,rbin,tstep,time,tmst, &
                                        chit,cint,chip,eta(1:9),   &
                                        strcon(1:9),strpmf(1:9),stress(1:9)

  Logical               :: ready
  Character( Len = 42 )  :: name
  Character( Len = 40 ) :: forma  = ' '

  Integer               :: fail(1:3),i,j,l,levcfg,jdnode,jatms,nsum
  Real( Kind = wp )     :: r_mxnode

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


! Define format for REVIVE printing in ASCII

  If (l_rout) Then
     i = 64/4 - 1 ! Bit_Size(0.0_wp)/4 - 1
     j = Max(mxstak*mxnstk+1,mxgrdf*mxrdf,mxgusr,mxgana*mxtana)

     Write(forma ,10) j/4+1,i+9,i
10   Format('(1p,',i0,'(/,4e',i0,'.',i0,'E3))')
  End If

  If (mxnode > 1) Then

! globally sum RDF information before saving

     If (lrdf) Then

! maximum rdfs that can be summed in each step

        nsum = mxbuff/mxgrdf
        If (nsum == 0) Call error(200)

        Do i=1,mxrdf,nsum
           Call gsum(rdf(:,i:Min(i+nsum-1,mxrdf)))
        End Do

     End If

! globally sum USR RDF information before saving

     If (mxgusr > 0) Call gsum(usr(1:mxgusr))

! globally sum z-density information before saving

     If (lzdn) Then

! maximum zdens that can be summed in each step

        nsum = mxbuff/mxgrdf
        If (nsum == 0) Call error(200)

        Do i=1,mxatyp,nsum
           Call gsum(zdens(:,i:Min(i+nsum-1,mxatyp)))
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
              Call gsum(vafdata(:,l+i:l+Min(i+nsum-1,mxatyp)))
           End Do
        End Do

     End If

! globally sum bonds' distributions information before saving

     If (mxgbnd1 > 0) Then

! maximum dstbnd that can be summed in each step

        nsum = mxbuff/(mxgbnd1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfbnd(0),nsum
           Call gsum(dstbnd(:,i:Min(i+nsum-1,ldfbnd(0))))
        End Do

     End If

! globally sum angles' distributions information before saving

     If (mxgang1 > 0) Then

! maximum dstang that can be summed in each step

        nsum = mxbuff/(mxgang1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfang(0),nsum
           Call gsum(dstang(:,i:Min(i+nsum-1,ldfang(0))))
        End Do

     End If

! globally sum dihedrals' distributions information before saving

     If (mxgdih1 > 0) Then

! maximum dstdih that can be summed in each step

        nsum = mxbuff/(mxgdih1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfdih(0),nsum
           Call gsum(dstdih(:,i:Min(i+nsum-1,ldfdih(0))))
        End Do

     End If

! globally sum inversions' distributions information before saving

     If (mxginv1 > 0) Then

! maximum dstinv that can be summed in each step

        nsum = mxbuff/(mxginv1+1)
        If (nsum == 0) Call error(200)

        Do i=1,ldfinv(0),nsum
           Call gsum(dstinv(:,i:Min(i+nsum-1,ldfinv(0))))
        End Do

     End If

  End If

! Write REVCON

  name = Trim(revcon) ! file name
  levcfg = 2      ! define level of information in REVCON

  Call write_config(name,levcfg,megatm,nstep,tstep,time)

! node 0 handles I/O

  If (idnode == 0) Then

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

! Write initial velocities for VAF calculations if needed

  If (vafsamp > 0) Then

    If (idnode == 0) Then

       jatms=natms
       Do j=1,vafsamp
         Do i=1,natms
            iwrk(i)=ltg(i)
            axx(i)=vxi(i,j)
            ayy(i)=vyi(i,j)
            azz(i)=vzi(i,j)
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
       End Do

    Else

       Do j=1,vafsamp
         Call MPI_RECV(ready,1,MPI_LOGICAL,0,Revive_tag,dlp_comm_world,status,ierr)

         Call MPI_SEND(natms,1,MPI_INTEGER,0,Revive_tag,dlp_comm_world,ierr)

         Call MPI_SEND(ltg,natms,MPI_INTEGER,0,Revive_tag,dlp_comm_world,ierr)

         Call MPI_SEND(vxi(1,j),natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
         Call MPI_SEND(vyi(1,j),natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
         Call MPI_SEND(vzi(1,j),natms,wp_mpi,0,Revive_tag,dlp_comm_world,ierr)
       End Do

    End If

    If (mxnode > 1) Call gsync()

  End If

  If (mxnode > 1) Call gsync()
  If (idnode == 0) Close(Unit=nrest)

  If (mxnode > 1) Then
     r_mxnode=1.0_wp/Real(mxnode,wp)

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
  End If

  Deallocate (iwrk,        Stat=fail(1))
  Deallocate (axx,ayy,azz, Stat=fail(2))
  Deallocate (bxx,byy,bzz, Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'system_revive deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine system_revive
