Module angular_distribution
!> Module to calculate the angular distribution function
!>
!> Copyright - Daresbury Laboratory
!>
!> Author:  A.J.Diver and O.A.Dicks,  December 2019
!> Contrib: a.m.elena, January 2020 - cleanup and error handling

  Use comms,           Only: comms_type,&
                             grecv,&
                             gsend
  Use configuration,   Only: configuration_type
  Use constants,       Only: nchadf,&
                             pi
  Use coord,           Only: coord_type
  Use errors_warnings, Only: error
  Use flow_control,    Only: flow_type
  Use kinds,           Only: wp
  Use site,            Only: site_type

  Implicit None
  Type, Public :: adf_type

    Real(Kind=wp)        :: rij(1:100), rik, rjk, prec
    Integer, Allocatable :: astat(:, :), coordlist(:, :)
    Integer              :: interval
    Logical              :: adfon
!        real(wp), allocatable :: rij(1:10)

  End Type adf_type

Contains

  Subroutine adf_calculate(config, sites, flow, crd, adf, comm)
    Type(configuration_type), Intent(In   ) :: config
    Type(site_type),          Intent(In   ) :: sites
    Type(flow_type),          Intent(In   ) :: flow
    Type(coord_type),         Intent(In   ) :: crd
    Type(adf_type),           Intent(InOUT) :: adf
    Type(comms_type),         Intent(InOut) :: comm

    Integer              :: i, ii, iii, j, jj, k, kk, nab, numbins
    Integer, Allocatable :: adfbuff(:)
    Logical              :: itsopen
    Real                 :: costheta, temptheta

    If (.not. adf%adfon) Return
    If (.not. crd%coordon) Return
    If (Mod(flow%step, adf%interval) /= 0) Return
    Open (Unit=nchadf, File='ADFDAT', Form='formatted')
    If (flow%step .eq. 0) Then
      numbins = 180.0 / adf%prec
      Allocate (adf%astat(-1:numbins, 1:2 * crd%ncoordpairs))
    Endif
    adf%astat(:, :) = 0
    Do i = 1, crd%ncoordpairs
      adf%astat(-1, (2 * i) - 1) = crd%ltype(i, 1)
      adf%astat(0, (2 * i) - 1) = crd%ltype(i, 2)
      adf%astat(-1, (2 * i)) = crd%ltype(i, 2)
      adf%astat(0, (2 * i)) = crd%ltype(i, 1)
    Enddo
    numbins = 180.0 / adf%prec
    Do i = 1, config%natms
      adf%rij(:) = 0_wp
      adf%rjk = 0_wp

      Do j = 1, crd%adfcoordlist(0, i)
        adf%rij(j) = (config%parts(i)%xxx - config%parts(crd%adfcoordlist(j, i))%xxx)**2 &
                     + (config%parts(i)%yyy - config%parts(crd%adfcoordlist(j, i))%yyy)**2 &
                     + (config%parts(i)%zzz - config%parts(crd%adfcoordlist(j, i))%zzz)**2
      End Do

      Do j = 1, crd%adfcoordlist(0, i) - 1
        Do jj = 1 + j, crd%adfcoordlist(0, i)
          If (config%ltype(j) .eq. config%ltype(jj)) Then
            Do ii = 1, 2 * crd%ncoordpairs
              If (adf%astat(-1, ii) == config%ltype(crd%adfcoordlist(j, i)) .and. adf%astat(0, ii) == config%ltype(i)) Then
                adf%rjk = (config%parts(crd%adfcoordlist(j, i))%xxx - config%parts(crd%adfcoordlist(jj, i))%xxx)**2 &
                          + (config%parts(crd%adfcoordlist(j, i))%yyy - config%parts(crd%adfcoordlist(jj, i))%yyy)**2 &
                          + (config%parts(crd%adfcoordlist(j, i))%zzz - config%parts(crd%adfcoordlist(jj, i))%zzz)**2
                costheta = ((adf%rij(j) + adf%rij(jj) - adf%rjk) / (2 * Sqrt(adf%rij(j)) * Sqrt(adf%rij(jj))))
                temptheta = Acos(costheta) * (180 / pi)

                Do iii = 1, numbins

                  If (temptheta .ge. (iii - 1) * adf%prec .and. temptheta .lt. (iii) * adf%prec) Then

                    adf%astat(iii, ii) = adf%astat(iii, ii) + 1
                  End If
                End Do
              End If
            End Do
          End If
        End Do
      End Do

    End Do
    nab = 0
    Do i = 1, 2 * crd%ncoordpairs
      nab = nab + 2 + numbins
    Enddo

    If (comm%idnode == 0) Then
      Inquire (unit=nchadf, opened=itsopen)
      If (.not. itsopen) Then
        Open (Unit=nchadf, File='ADFDAT', Form='formatted')
      Endif

      Do j = 1, comm%mxnode - 1
        Call grecv(comm, nab, j, j)
        If (nab > 0) Then
          Allocate (adfbuff(nab))
          Call grecv(comm, adfbuff, j, j)
          Do ii = 1, 2 * crd%ncoordpairs
            If (adf%astat(-1, ii) /= adfbuff(1 + (numbins + 2) * (ii - 1)) .or. &
                adf%astat(0, ii) /= adfbuff(2 + (numbins + 2) * (ii - 1))) Then
              Call error(0, 'ERROR: adf pairs do not match in MPI')
            Endif
            Do kk = 1, numbins
              adf%astat(kk, ii) = adf%astat(kk, ii) + adfbuff(2 + kk + (numbins + 2) * (ii - 1))
            Enddo
          Enddo
          Deallocate (adfbuff)
        End If
      Enddo

!        do ii=1,2*crd%ncoordpairs
!        do kk=1,180
!        adf%astatavg(kk,ii)=adf%astatavg(kk,ii)+adf%astat(kk,ii)
!        enddo
!        enddo

      Write (nchadf, '(A29,I10,F20.6)') "Angular distribution function", flow%step, flow%time
      Do i = 1, 2 * crd%ncoordpairs
        Write (nchadf, *) Trim(sites%unique_atom(adf%astat(-1, i))), '-', Trim(sites%unique_atom(adf%astat(0, i))) &
          , '-', Trim(sites%unique_atom(adf%astat(-1, i)))
        Do ii = 1, numbins
          Write (nchadf, *) (adf%prec * ii) - adf%prec / 2, adf%astat(ii, i)
        End Do
      End Do

    Else
      Allocate (adfbuff(nab))
      k = 0
      Do i = 1, 2 * crd%ncoordpairs
        k = k + 1
        adfbuff(k) = adf%astat(-1, i)
        k = k + 1
        adfbuff(k) = adf%astat(0, i)
        Do ii = 1, numbins
          k = k + 1
          adfbuff(K) = adf%astat(ii, i)
        Enddo
      Enddo
      Call gsend(comm, nab, 0, comm%idnode)
      If (nab > 0) Then
        Call gsend(comm, adfbuff, 0, comm%idnode)
      Endif
      Deallocate (adfbuff)
    Endif
  End Subroutine adf_calculate

End Module angular_distribution

