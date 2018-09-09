Subroutine ttm_system_revive    &
           (dumpfile,nstep,time,freq,nstrun)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for writing electronic temperature restart files
! at job termination or selected intervals in simulation
!
! copyright - daresbury laboratory
! authors   - s.l.daraszewicz & m.a.seaton september 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use setup_module
  Use ttm_module
  Use comms_module

  Implicit None

  Character (Len = *), Intent ( In ) :: dumpfile
  Integer, Intent ( In ) :: nstep,freq,nstrun
  Real(Kind=wp), Intent ( In ) :: time
  Integer :: iounit = 117
  Integer :: id,ii,jj,kk,imin,jmin,kmin,imax,jmax,kmax,i,j,k,ijk,ix,iy,iz
  Logical :: lrange

  If (freq /=0) Then
    If (Mod(nstep,freq)==0 .or. nstep==nstrun) Then

      If (idnode==0) Then
        Open(Unit=iounit, File=dumpfile, Status='replace')
        Write(iounit,'(3i8)') eltsys(1),eltsys(2),eltsys(3)
        Write(iounit,'(i12,3(2x,es24.15))') nstep,time,depostart,depoend
        Close(iounit)
      End If
      If (mxnode > 1) Call gsync()

      Do id=0,mxnode-1
        If (idnode==id) Then
          Open(Unit=iounit, File=dumpfile, Status='old', Position='append')
          Do kk = -eltcell(3), eltcell(3)
            If (eltcell(3)>0 .and. kk == -eltcell(3) .and. ttmbcmap(5)>=0) Then
              kmin = ttmbc(5)
            Else
              kmin = 1
            End If
            If (eltcell(3)>0 .and. kk == eltcell(3) .and. ttmbcmap(6)>=0) Then
              kmax = ttmbc(6)
            Else
              kmax = ntcell(3)
            End If
            Do jj = -eltcell(2), eltcell(2)
              If (eltcell(2)>0 .and. jj == -eltcell(2) .and. ttmbcmap(3)>=0) Then
                jmin = ttmbc(3)
              Else
                jmin = 1
              End If
              If (eltcell(2)>0 .and. jj == eltcell(2) .and. ttmbcmap(4)>=0) Then
                jmax = ttmbc(4)
              Else
                jmax = ntcell(2)
              End If
              Do ii = -eltcell(1), eltcell(1)
                If (eltcell(1)>0 .and. ii == -eltcell(1) .and. ttmbcmap(1)>=0) Then
                  imin = ttmbc(1)
                Else
                  imin = 1
                End If
                If (eltcell(1)>0 .and. ii == eltcell(1) .and. ttmbcmap(2)>=0) Then
                  imax = ttmbc(2)
                Else
                  imax = ntcell(1)
                End If
                Do k = kmin, kmax
                  iz = k + ntcelloff(3) + (kk + eltcell(3)) * ntsys(3) - zeroE(3)
                  Do j = jmin, jmax
                    iy = j + ntcelloff(2) + (jj + eltcell(2)) * ntsys(2) - zeroE(2)
                    Do i = imin, imax
                      ix = i + ntcelloff(1) + (ii + eltcell(1)) * ntsys(1) - zeroE(1)
                      lrange = (ix>0 .and. ix<=eltsys(1) .and. iy>0 .and. iy<=eltsys(2) .and. iz>0 .and. iz<=eltsys(3))
                      ijk = 1 + i + (ntcell(1)+2) * (j + (ntcell(2)+2) * k)
                      If (lrange) Write(iounit,'(3i8,2x,es24.15)') ix-midE(1),iy-midE(2),iz-midE(3),eltemp(ijk,ii,jj,kk)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
          Close(iounit)
        End If
        IF (mxnode>1) Call gsync()
      End Do

    End If
  End If

End Subroutine ttm_system_revive
