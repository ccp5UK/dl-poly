Subroutine statistics_connect_frames(megatm)

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

  Use comms_module,      Only : mxnode,idnode,gsum
  Use setup_module,      Only : nrite
  Use domains_module,    Only : nprx,npry,nprz
  Use config_module,     Only : natms
  Use statistics_module, Only : found,natms0

  Implicit None

  Integer, Intent ( In    ) :: megatm

  Integer :: icyc,nres

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
     Call gsum(nres)
     If (nres > 0) Then
        nres=Merge(0,Sum(found(1:natms)),natms > 0)
        Call gsum(nres)
        If (nres /= megatm) icyc = icyc + 1
     End If
  End Do

  If (nres > 0 .and. idnode == 0) &
     Write(nrite,'(/,1x,a)') '*** warning - particles dynamics properties will be corrupted!!! ***'

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

    Use setup_module,  Only : mxatdm,mxstak
    Use config_module, Only : ixyz,lsa,lsi

    Use statistics_module

    Use msd_module,    Only : l_msd

    Implicit None

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
          Write(nrite,'(/,1x,a,i0)') 'match_compress_spread_sort deallocation failure, node: ', idnode
          Call error(0)
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
          Write(nrite,'(/,1x,a,i0)') 'match_compress_spread_sort allocation failure, node: ', idnode
          Call error(0)
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

    If (mdir /= 0) Call statistics_connect_spread(mdir)

! Sort past frame remainder of global atom indices

!    lsi0=0 ; lsa0=0
    Do i0=1,natms0
       lsi0(i0)=i0
       lsa0(i0)=ltg0(i0)
    End Do
    Call shellsort2(natms0,lsi0,lsa0)

  End Subroutine match_compress_spread_sort

End Subroutine statistics_connect_frames
