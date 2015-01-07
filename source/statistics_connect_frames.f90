Subroutine statistics_connect_frames()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to arrange exchange of data between neighbouring
! domains/nodes in order to reconnect some statistical information
! between replayed frames of history
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use comms_module,  Only : mxnode,idnode,gsum
  Use setup_module,  Only : nrite,mxstak
  Use config_module, Only : natms,ixyz,lsa

  Use statistics_module

  Use msd_module,    Only : l_msd

  Implicit None

  Integer :: icyc,nres,i,i0,j,j0,kk

  icyc = 0 ; nres = 1
  Do While (mxnode > 1 .and. (icyc < 4 .and. nres > 0))
     found = 0

! Search for matches

     found0 = 0
     i0=1
     Do i=1,natms
        If      (lsa(i) <  lsa0(i0)) Then
           Cycle
        Else If (lsa(i) == lsa0(i0)) Then
           found(i) = 1 ; found0(i0) = 1
           xin(i) = xin0(i0) ; yin(i) = yin0(i0) ; zin(i) = zin0(i0)
           xto(i) = xto0(i0) ; yto(i) = yto0(i0) ; zto(i) = zto0(i0)

           If (l_msd) Then
              j =27+2*i
              j0=2*i0
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

           If (i0 < natms0) Then
              i0 = i0 + 1
           Else
              Exit
           End If
        Else
           If (i0 < natms0) Then
              i0 = i0 + 1
           Else
              Exit
           End If
        End If
     End Do

! Compress remainder

     i0=0
     Do While (i0 < natms0)
        i0 = i0 + 1
        If (found0(i0) == 0) Then
10         Continue
           If (found0(natms0) > 0) Then
              ixyz(i0) = ixyz(natms0)
              If (ixyz(i0) == 0) ixyz(i0) = 333 ! Bloody send all directions
              ltg0(i0) = ltg0(natms0)

              xin0(i0) = xin0(natms0) ; yin0(i0) = yin0(natms0) ; zin0(i0) = zin0(natms0)
              xto0(i0) = xto0(natms0) ; yto0(i0) = yto0(natms0) ; zto0(i0) = zto0(natms0)

              If (l_msd) Then
                 j =27+2*i0
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

              natms0 = natms0 - 1
           Else
              natms0 = natms0 - 1
              If (i0 < natms0) Go To 10
           End If
        End If
     End Do

! spread atom data in -/+ x directions

     Call statistics_connect_spread(-1)

! Record global atom indices for halo sorting

     Do i0=1,natms0
        lsi0(i0)=i0
        lsa0(i0)=ltg0(i0)
     End Do
     Call shellsort2(natms0,lsi0,lsa0)

! Search for matches

     found0 = 0
     i0=1
     Do i=1,natms
        If (found(i) == 0) Then
           If      (lsa(i) <  lsa0(i0)) Then
              Cycle
           Else If (lsa(i) == lsa0(i0)) Then
              found(i) = 1 ; found0(i0) = 1
              xin(i) = xin0(i0) ; yin(i) = yin0(i0) ; zin(i) = zin0(i0)
              xto(i) = xto0(i0) ; yto(i) = yto0(i0) ; zto(i) = zto0(i0)

              If (l_msd) Then
                 j =27+2*i
                 j0=2*i0
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

              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           Else
              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           End If
        End If
     End Do

! Compress remainder

     i0=0
     Do While (i0 < natms0)
        i0 = i0 + 1
        If (found0(i0) == 0) Then
20         Continue
           If (found0(natms0) > 0) Then
              ixyz(i0) = ixyz(natms0)
              If (ixyz(i0) == 0) ixyz(i0) = 333 ! Bloody send all directions
              ltg0(i0) = ltg0(natms0)

              xin0(i0) = xin0(natms0) ; yin0(i0) = yin0(natms0) ; zin0(i0) = zin0(natms0)
              xto0(i0) = xto0(natms0) ; yto0(i0) = yto0(natms0) ; zto0(i0) = zto0(natms0)

              If (l_msd) Then
                 j =27+2*i0
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

              natms0 = natms0 - 1
           Else
              natms0 = natms0 - 1
              If (i0 < natms0) Go To 20
           End If
        End If
     End Do

     Call statistics_connect_spread( 1)

! Record global atom indices for halo sorting

     Do i0=1,natms0
        lsi0(i0)=i0
        lsa0(i0)=ltg0(i0)
     End Do
     Call shellsort2(natms0,lsi0,lsa0)

! Search for matches

     found0 = 0
     i0=1
     Do i=1,natms
        If (found(i) == 0) Then
           If      (lsa(i) <  lsa0(i0)) Then
              Cycle
           Else If (lsa(i) == lsa0(i0)) Then
              found(i) = 1 ; found0(i0) = 1
              xin(i) = xin0(i0) ; yin(i) = yin0(i0) ; zin(i) = zin0(i0)
              xto(i) = xto0(i0) ; yto(i) = yto0(i0) ; zto(i) = zto0(i0)

              If (l_msd) Then
                 j =27+2*i
                 j0=2*i0
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

              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           Else
              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           End If
        End If
     End Do

! Compress remainder

     i0=0
     Do While (i0 < natms0)
        i0 = i0 + 1
        If (found0(i0) == 0) Then
30         Continue
           If (found0(natms0) > 0) Then
              ixyz(i0) = ixyz(natms0)
              If (ixyz(i0) == 0) ixyz(i0) = 333 ! Bloody send all directions
              ltg0(i0) = ltg0(natms0)

              xin0(i0) = xin0(natms0) ; yin0(i0) = yin0(natms0) ; zin0(i0) = zin0(natms0)
              xto0(i0) = xto0(natms0) ; yto0(i0) = yto0(natms0) ; zto0(i0) = zto0(natms0)

              If (l_msd) Then
                 j =27+2*i0
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

              natms0 = natms0 - 1
           Else
              natms0 = natms0 - 1
              If (i0 < natms0) Go To 30
           End If
        End If
     End Do

! spread atom data in -/+ y directions

     Call statistics_connect_spread(-2)

! Record global atom indices for halo sorting

     Do i0=1,natms0
        lsi0(i0)=i0
        lsa0(i0)=ltg0(i0)
     End Do
     Call shellsort2(natms0,lsi0,lsa0)

! Search for matches

     found0 = 0
     i0=1
     Do i=1,natms
        If (found(i) == 0) Then
           If      (lsa(i) <  lsa0(i0)) Then
              Cycle
           Else If (lsa(i) == lsa0(i0)) Then
              found(i) = 1 ; found0(i0) = 1
              xin(i) = xin0(i0) ; yin(i) = yin0(i0) ; zin(i) = zin0(i0)
              xto(i) = xto0(i0) ; yto(i) = yto0(i0) ; zto(i) = zto0(i0)

              If (l_msd) Then
                 j =27+2*i
                 j0=2*i0
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

              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           Else
              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           End If
        End If
     End Do

! Compress remainder

     i0=0
     Do While (i0 < natms0)
        i0 = i0 + 1
        If (found0(i0) == 0) Then
40         Continue
           If (found0(natms0) > 0) Then
              ixyz(i0) = ixyz(natms0)
              If (ixyz(i0) == 0) ixyz(i0) = 333 ! Bloody send all directions
              ltg0(i0) = ltg0(natms0)

              xin0(i0) = xin0(natms0) ; yin0(i0) = yin0(natms0) ; zin0(i0) = zin0(natms0)
              xto0(i0) = xto0(natms0) ; yto0(i0) = yto0(natms0) ; zto0(i0) = zto0(natms0)

              If (l_msd) Then
                 j =27+2*i0
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

              natms0 = natms0 - 1
           Else
              natms0 = natms0 - 1
              If (i0 < natms0) Go To 40
           End If
        End If
     End Do

     Call statistics_connect_spread( 2)

! Record global atom indices for halo sorting

     Do i0=1,natms0
        lsi0(i0)=i0
        lsa0(i0)=ltg0(i0)
     End Do
     Call shellsort2(natms0,lsi0,lsa0)

! Search for matches

     found0 = 0
     i0=1
     Do i=1,natms
        If (found(i) == 0) Then
           If      (lsa(i) <  lsa0(i0)) Then
              Cycle
           Else If (lsa(i) == lsa0(i0)) Then
              found(i) = 1 ; found0(i0) = 1
              xin(i) = xin0(i0) ; yin(i) = yin0(i0) ; zin(i) = zin0(i0)
              xto(i) = xto0(i0) ; yto(i) = yto0(i0) ; zto(i) = zto0(i0)

              If (l_msd) Then
                 j =27+2*i
                 j0=2*i0
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

              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           Else
              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           End If
        End If
     End Do

! Compress remainder

     i0=0
     Do While (i0 < natms0)
        i0 = i0 + 1
        If (found0(i0) == 0) Then
50         Continue
           If (found0(natms0) > 0) Then
              ixyz(i0) = ixyz(natms0)
              If (ixyz(i0) == 0) ixyz(i0) = 333 ! Bloody send all directions
              ltg0(i0) = ltg0(natms0)

              xin0(i0) = xin0(natms0) ; yin0(i0) = yin0(natms0) ; zin0(i0) = zin0(natms0)
              xto0(i0) = xto0(natms0) ; yto0(i0) = yto0(natms0) ; zto0(i0) = zto0(natms0)

              If (l_msd) Then
                 j =27+2*i0
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

              natms0 = natms0 - 1
           Else
              natms0 = natms0 - 1
              If (i0 < natms0) Go To 50
           End If
        End If
     End Do

! spread atom data in -/+ z directions

     Call statistics_connect_spread(-3)

! Record global atom indices for halo sorting

     Do i0=1,natms0
        lsi0(i0)=i0
        lsa0(i0)=ltg0(i0)
     End Do
     Call shellsort2(natms0,lsi0,lsa0)

! Search for matches

     found0 = 0
     i0=1
     Do i=1,natms
        If (found(i) == 0) Then
           If      (lsa(i) <  lsa0(i0)) Then
              Cycle
           Else If (lsa(i) == lsa0(i0)) Then
              found(i) = 1 ; found0(i0) = 1
              xin(i) = xin0(i0) ; yin(i) = yin0(i0) ; zin(i) = zin0(i0)
              xto(i) = xto0(i0) ; yto(i) = yto0(i0) ; zto(i) = zto0(i0)

              If (l_msd) Then
                 j =27+2*i
                 j0=2*i0
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

              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           Else
              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           End If
        End If
     End Do

! Compress remainder

     i0=0
     Do While (i0 < natms0)
        i0 = i0 + 1
        If (found0(i0) == 0) Then
60         Continue
           If (found0(natms0) > 0) Then
              ixyz(i0) = ixyz(natms0)
              If (ixyz(i0) == 0) ixyz(i0) = 333 ! Bloody send all directions
              ltg0(i0) = ltg0(natms0)

              xin0(i0) = xin0(natms0) ; yin0(i0) = yin0(natms0) ; zin0(i0) = zin0(natms0)
              xto0(i0) = xto0(natms0) ; yto0(i0) = yto0(natms0) ; zto0(i0) = zto0(natms0)

              If (l_msd) Then
                 j =27+2*i0
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

              natms0 = natms0 - 1
           Else
              natms0 = natms0 - 1
              If (i0 < natms0) Go To 60
           End If
        End If
     End Do

     Call statistics_connect_spread( 3)

! Record global atom indices for halo sorting

     Do i0=1,natms0
        lsi0(i0)=i0
        lsa0(i0)=ltg0(i0)
     End Do
     Call shellsort2(natms0,lsi0,lsa0)

! Search for matches

     found0 = 0
     i0=1
     Do i=1,natms
        If (found(i) == 0) Then
           If      (lsa(i) <  lsa0(i0)) Then
              Cycle
           Else If (lsa(i) == lsa0(i0)) Then
              found(i) = 1 ; found0(i0) = 1
              xin(i) = xin0(i0) ; yin(i) = yin0(i0) ; zin(i) = zin0(i0)
              xto(i) = xto0(i0) ; yto(i) = yto0(i0) ; zto(i) = zto0(i0)

              If (l_msd) Then
                 j =27+2*i
                 j0=2*i0
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

              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           Else
              If (i0 < natms0) Then
                 i0 = i0 + 1
              Else
                 Exit
              End If
           End If
        End If
     End Do

! Compress remainder

     i0=0
     Do While (i0 < natms0)
        i0 = i0 + 1
        If (found0(i0) == 0) Then
70         Continue
           If (found0(natms0) > 0) Then
              ixyz(i0) = ixyz(natms0)
              If (ixyz(i0) == 0) ixyz(i0) = 333 ! Bloody send all directions
              ltg0(i0) = ltg0(natms0)

              xin0(i0) = xin0(natms0) ; yin0(i0) = yin0(natms0) ; zin0(i0) = zin0(natms0)
              xto0(i0) = xto0(natms0) ; yto0(i0) = yto0(natms0) ; zto0(i0) = zto0(natms0)

              If (l_msd) Then
                 j =27+2*i0
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

              natms0 = natms0 - 1
           Else
              natms0 = natms0 - 1
              If (i0 < natms0) Go To 70
           End If
        End If
     End Do

     nres=natms0
     Call gsum(nres)
     If (nres > 0) icyc = icyc + 1
  End Do

  If (icyc == 3 .and. nres > 0 .and. idnode == 0) &
        Write(nrite,'(/,1x,a)') '*** warning - particles dynamics properties will be corrupted!!! ***'

End Subroutine statistics_connect_frames
