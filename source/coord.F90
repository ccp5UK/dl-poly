Module coord

!> Module to calculate the coordination number distributions during a DLPOLY run
!>
!> Author  Oliver Dicks and Aaron Diver, February 2018
!> Contrib a.m.elena January 2020 - cleanup and errors

  Use comms,           Only: comms_type,&
                             grecv,&
                             gsend
  Use configuration,   Only: configuration_type
  Use constants,       Only: nccrdt,&
                             nicrdt
  Use errors_warnings, Only: error
  Use flow_control,    Only: flow_type
  Use kinds,           Only: wp
  Use neighbours,      Only: neighbours_type
  Use site,            Only: site_type
  Use statistics,      Only: stats_type

  Implicit None

  Type, Public :: coord_type

    Integer                       :: ncoordpairs, coordinterval, coordstart, coordops, &
                                     ncoorddis, ncoordab, maxlist = 1000
    Real(wp), Allocatable         :: arraycuts(:), discuts(:)
    Character(Len=8), Allocatable :: arraypairs(:, :), disatms(:)
    Integer, Allocatable          :: coordlist(:, :), icoordlist(:, :), defectlist(:), adfcoordlist(:, :)
    Integer, Allocatable          :: ltype(:, :), ltypeA(:, :), ltypeB(:, :), cstat(:, :), disltype(:)
    Logical                       :: coordon
    Real(wp)                      :: coordis

  Contains

    Procedure :: init => init_coord
    Procedure :: init_coordlist
    final     :: clean_coord

  End Type coord_type
  Private
  Public :: init_coord_list, checkcoord

Contains

  Subroutine init_coord(T)
    class(coord_type) :: T

    Allocate (T%arraycuts(1:T%ncoordpairs))
    Allocate (T%discuts(1:T%ncoorddis))
    Allocate (T%arraypairs(1:T%ncoordpairs, 1:2))
    Allocate (T%ltype(1:T%ncoordpairs, 1:2))
    Allocate (T%ltypeA(1:2 * T%ncoordab, 0:2 * T%ncoordpairs))
    Allocate (T%ltypeB(1:2 * T%ncoordab, 0:2 * T%ncoordpairs))
    Allocate (T%cstat(-3:T%maxlist, 1:(2 * (T%ncoordpairs + T%ncoordab))))
    Allocate (T%disltype(1:T%ncoorddis))
!   allocate(T%coordatoms(0:2*T%ncoordpairs))
  End Subroutine init_coord

  Subroutine init_coordlist(T, n, m)
    class(coord_type)      :: T
    Integer, Intent(In   ) :: n, m

    Allocate (T%coordlist(0:n, 1:m))
    Allocate (T%icoordlist(0:n, 1:m))
    Allocate (T%defectlist(0:m))
    Allocate (T%adfcoordlist(0:n, 1:m))
  End Subroutine init_coordlist

  Subroutine clean_coord(T)
    Type(coord_type) :: T

    If (Allocated(T%arraycuts)) Then
      Deallocate (T%arraycuts)
    End If
    If (Allocated(T%discuts)) Then
      Deallocate (T%discuts)
    End If
    If (Allocated(T%coordlist)) Then
      Deallocate (T%coordlist)
    End If
    If (Allocated(T%icoordlist)) Then
      Deallocate (T%icoordlist)
    End If
    If (Allocated(T%adfcoordlist)) Then
      Deallocate (T%adfcoordlist)
    End If
    If (Allocated(T%defectlist)) Then
      Deallocate (T%defectlist)
    End If
  End Subroutine clean_coord

  Subroutine init_coord_list(config, neigh, crd, sites, flow, comm)
    Type(configuration_type), Intent(In   ) :: config
    Type(neighbours_type),    Intent(In   ) :: neigh
    Type(coord_type),         Intent(InOut) :: crd
    Type(site_type),          Intent(In   ) :: sites
    Type(flow_type),          Intent(In   ) :: flow
    Type(comms_type),         Intent(InOut) :: comm

    Character(len=1000), Allocatable :: cbuff(:)
    Character(len=8)                 :: aux
    Integer                          :: en, i, ii, j, jj, k, kk, m, mcoord, ncb, ncoord, nmax
    Integer, Allocatable             :: buff(:), coordbuff(:), dbuff(:), ltypeAB(:, :)
    Logical                          :: itsopen
    Real(kind=wp)                    :: rab, rcut

!ncb size coordbuff

    !Check whether option is called
    If (.not. crd%coordon) Return
    If (crd%ncoordpairs == 0) Return
    If (crd%coordstart > flow%step) Return
    If (Mod(flow%step, crd%coordinterval) /= 0) Return
    crd%coordlist(0, :) = 0
    Do i = 1, crd%ncoordab
!    print*,crd%ltypeA(i,:)
!    print*,crd%ltypeB(i,:)
    End Do
    ncb = 0
    Do j = 1, config%natms
      ncoord = 0
      k = neigh%list(0, j)
      Do i = 1, k
        kk = neigh%list(i, j)

        rab = (config%parts(j)%xxx - config%parts(kk)%xxx)**2 + (config%parts(j)%yyy - config%parts(kk)%yyy)**2 &
              + (config%parts(j)%zzz - config%parts(kk)%zzz)**2
!        if(j.eq.2)then
!        print*,j,kk,rab

!        endif
        rcut = 0.00

        Do ii = 1, crd%ncoordpairs
          If (((config%ltype(j) == crd%ltype(ii, 1)) .and. (config%ltype(kk) == crd%ltype(ii, 2))) &
              .or. &
              ((config%ltype(j) == crd%ltype(ii, 2)) .and. (config%ltype(kk) == crd%ltype(ii, 1)))) Then
            rcut = crd%arraycuts(ii) * crd%arraycuts(ii)
          Endif
        End Do

        If (rab <= rcut) Then
          crd%coordlist(0, j) = crd%coordlist(0, j) + 1
          crd%coordlist(crd%coordlist(0, j), j) = kk
          If (kk <= config%natms) Then
            crd%coordlist(0, kk) = crd%coordlist(0, kk) + 1
            crd%coordlist(crd%coordlist(0, kk), kk) = j
          Endif
        End If
      End Do
    End Do

    !Create coordination statistics array
    crd%cstat(-3:, :) = 0
    Do i = 1, crd%ncoordpairs
      crd%cstat(-3, (2 * i) - 1) = crd%ltype(i, 1)
      crd%cstat(-2, (2 * i) - 1) = crd%ltype(i, 2)
      crd%cstat(-3, (2 * i)) = crd%ltype(i, 2)
      crd%cstat(-2, (2 * i)) = crd%ltype(i, 1)
    End Do

    !Collect coordination statistics
    Do i = 1, config%natms
      Do j = 1, 2 * crd%ncoordpairs
        mcoord = 0
        If ((config%ltype(i)) == crd%cstat(-3, j)) Then
          Do k = 1, crd%coordlist(0, i)
            If (config%ltype(crd%coordlist(k, i)) == crd%cstat(-2, j)) Then
              mcoord = mcoord + 1
            End If
          End Do
          crd%cstat(mcoord, j) = crd%cstat(mcoord, j) + 1
          If (mcoord > crd%cstat(-1, j)) Then
            crd%cstat(-1, j) = mcoord
          End If
        End If
      End Do
    End Do

    !Create coordination AB statistics array
    Allocate (ltypeAB(1:2 * crd%ncoordab, 0:2 * crd%ncoordpairs))
    Do i = 1, crd%ncoordab
!     Do ii=1,crd%ncoordpairs
      ltypeAB(2 * i - 1, :) = crd%ltypeA(i, :)
      ltypeAB(2 * i, :) = crd%ltypeB(i, :)
      !   End Do
    End Do
    If (crd%ncoordab > 0) Then
      Do i = 1, crd%ncoordab
        crd%cstat(-3, 2 * crd%ncoordpairs + (2 * i - 1)) = 2 * i - 1
        crd%cstat(-2, 2 * crd%ncoordpairs + (2 * i - 1)) = 2 * i
        crd%cstat(-3, 2 * crd%ncoordpairs + (2 * i)) = 2 * i
        crd%cstat(-2, 2 * crd%ncoordpairs + (2 * i)) = 2 * i - 1
      End Do
    End If

!    Do i=1, 2*(crd%ncoordpairs+crd%ncoordab)
!      print*,crd%cstat(-3:0,i)
!    End do
    !   Do i = 1, 2 * crd%ncoordab
    !     Print *, ltypeAB(i, :)
    !   End Do
    !Collect AB coordination statistics
    Do i = 1, config%natms
      Do j = 1, 2 * crd%ncoordab
        mcoord = 0
!        print*,ltypeAB(crd%cstat(-3,j+(2*crd%ncoordpairs)),1:ltypeAB(crd%cstat(-3,j+(2*crd%ncoordpairs)),0))
!        print*,config%ltype(i)
!        if ((config%ltype(i)) == Any(ltypeAB(crd%cstat(-3,j+(2*crd%ncoordpairs)),1:))) then
        If (Any(ltypeAB(crd%cstat(-3, j + (2 * crd%ncoordpairs)), 1:ltypeAB(crd%cstat(-3, j + (2 * crd%ncoordpairs)), 0)) &
                == config%ltype(i))) Then
          Do k = 1, crd%coordlist(0, i)
!            If (config%ltype(crd%coordlist(k,i)) == Any(ltypeAB(crd%cstat(-2,j+(2*crd%ncoordpairs)),:))) then
            If (Any(ltypeAB(crd%cstat(-2, j + (2 * crd%ncoordpairs)), 1:ltypeAB(crd%cstat(-2, j + (2 * crd%ncoordpairs)), 0)) &
                    == config%ltype(crd%coordlist(k, i)))) Then
              mcoord = mcoord + 1
            End If
          End Do
          crd%cstat(mcoord, j + (2 * crd%ncoordpairs)) = crd%cstat(mcoord, (j + (2 * crd%ncoordpairs))) + 1
          If (mcoord > crd%cstat(-1, j + (2 * crd%ncoordpairs))) Then
            crd%cstat(-1, j + (2 * crd%ncoordpairs)) = mcoord
          End If
        End If
      End Do
    End Do

!     Do i=1, 2*(crd%ncoordpairs+crd%ncoordab)
!      print*,crd%cstat(-3:10,i)
!    End do

!Set coordbuff size
    Do i = 1, 2 * (crd%ncoordpairs + crd%ncoordab)
      ncb = ncb + (crd%cstat(-1, i) + 4)
    End Do

    Allocate (buff(2 * config%mxatms))
    Allocate (cbuff(config%mxatms))
    Allocate (dbuff(config%mxatms))
    If (comm%idnode == 0) Then
      Inquire (Unit=nicrdt, opened=itsopen)
      If (.not. itsopen) Then
        Open (Unit=nicrdt, File='ICOORD', Form='formatted')
      End If
      If (flow%step == crd%coordstart) Then
        Write (Unit=nicrdt, Fmt='(a72)') config%cfgname(1:72)
        Write (Unit=nicrdt, Fmt='(a40)') 'Initial coordination between atoms'
      Endif

      If ((crd%coordops == 2) .or. crd%coordstart == flow%step) Then
      Do i = 1, config%natms
        m = crd%coordlist(0, i)
        Write (nicrdt, Fmt='(i12,1x,a8,i12,1x)', advance="no") &
          config%ltg(i), Trim(sites%unique_atom(config%ltype(i))), crd%coordlist(0, i)
        Do ii = 1, m
          Write (nicrdt, '(i0,1x)', advance="no") config%ltg(crd%coordlist(ii, i))
        Enddo
        Do ii = 1, m
          Write (nicrdt, '(a,1x)', advance="no") Trim(sites%unique_atom(config%ltype(crd%coordlist(ii, i))))
        Enddo
        Write (nicrdt, *)
      Enddo

      Do j = 1, comm%mxnode - 1
        Call grecv(comm, en, j, j)
        If (en > 0) Then
          Call grecv(comm, buff, j, j)
          Call grecv(comm, cbuff, j, j)
          Call grecv(comm, dbuff, j, j)
          Do i = 1, en / 2
            Write (nicrdt, Fmt='(i12,1x,a8,i12,a)') &
              !              buff(2*i-1),trim(sites%unique_atom(config%ltype(config%ltg(buff(2*i-1))))),buff(2*i),trim(cbuff(i))
              buff(2 * i - 1), sites%unique_atom(dbuff(i)), buff(2 * i), Trim(cbuff(i))
          Enddo
        Endif
      Enddo
      Endif
    Else
      If ((crd%coordops == 2) .or. crd%coordstart == flow%step) Then
        en = 2 * config%natms
        Do i = 1, config%natms
          buff(2 * i - 1) = config%ltg(i)
          buff(2 * i) = crd%coordlist(0, i)
          cbuff(i) = ''
          dbuff(i) = config%ltype(i)
          Do ii = 1, crd%coordlist(0, i)
            Write (aux, '(i0)') config%ltg(crd%coordlist(ii, i))
            cbuff(i) = Trim(cbuff(i))//" "//Trim(aux)
          Enddo
          Do ii = 1, crd%coordlist(0, i)
            Write (aux, '(a)') Trim(sites%unique_atom(config%ltype(crd%coordlist(ii, i))))
            cbuff(i) = Trim(cbuff(i))//" "//Trim(aux)

          Enddo
        Enddo
        Call gsend(comm, en, 0, comm%idnode)
        If (en > 0) Then
          Call gsend(comm, buff, 0, comm%idnode)
          Call gsend(comm, cbuff, 0, comm%idnode)
          Call gsend(comm, dbuff, 0, comm%idnode)
        Endif
      Endif
    Endif
    Deallocate (buff)
    Deallocate (cbuff)
    Deallocate (dbuff)

    If (comm%idnode == 0) Then
      Open (Unit=nicrdt, File='ICOORD', Form='formatted')

      Do j = 1, comm%mxnode - 1
        Call grecv(comm, ncb, j, j)
        If (ncb > 0) Then
          Allocate (coordbuff(ncb))
          Call grecv(comm, coordbuff, j, j)
          jj = 1
          Do ii = 1, 2 * (crd%ncoordpairs + crd%ncoordab)
            nmax = coordbuff(jj)
            If (nmax > crd%cstat(-1, ii)) Then
              crd%cstat(-1, ii) = nmax
            End If
            If (crd%cstat(-3, ii) /= coordbuff(jj + 1) .or. crd%cstat(-2, ii) /= coordbuff(jj + 2)) Then
              Call error(0, "ERROR: coord pairs do not match in MPI")
            End If
            Do kk = 0, nmax
              crd%cstat(kk, ii) = crd%cstat(kk, ii) + coordbuff(jj + 3 + kk)
            End Do
            jj = jj + nmax + 4
          End Do
          Deallocate (coordbuff)
        End If
      Enddo

      Write (nicrdt, '(A36,I10,F20.6)') "Coordination distribution statistics", flow%step, flow%time

      Do i = 1, 2 * crd%ncoordpairs
        Do j = 0, crd%cstat(-1, i)
          Write (nicrdt, *) sites%unique_atom(crd%cstat(-3, i)), sites%unique_atom(crd%cstat(-2, i)), j, crd%cstat(j, i)
        End Do
      End Do

      Do i = 2 * crd%ncoordpairs + 1, 2 * (crd%ncoordpairs + crd%ncoordab)
        If (Mod(i, 2) == 1) Then
        Do j = 0, crd%cstat(-1, i)
          Write (nicrdt, '(A6,I0,1X,A6,I0,2X,I12,I12)') 'ListA', i - 2 * crd%ncoordpairs, &
            'ListB', i - 2 * crd%ncoordpairs, j, crd%cstat(j, i)
        End Do
        End If
        If (Mod(i, 2) == 0) Then
        Do j = 0, crd%cstat(-1, i)
          Write (nicrdt, '(A6,I0,1X,A6,I0,2X,I12,I12)') 'ListB', i - 2 * crd%ncoordpairs - 1, &
            'ListA', i - 2 * crd%ncoordpairs - 1, j, crd%cstat(j, i)
        End Do
        End If
      End Do

    Else
      Allocate (coordbuff(ncb))
      k = 1
      Do i = 1, 2 * (crd%ncoordpairs + crd%ncoordab)
        coordbuff(k) = crd%cstat(-1, i)
        k = k + 1
        coordbuff(k:k + 1) = crd%cstat(-3:-2, i)
        k = k + 2
        Do j = 0, crd%cstat(-1, i)
          coordbuff(k) = crd%cstat(j, i)
          k = k + 1
        End Do
      End Do
      Call gsend(comm, ncb, 0, comm%idnode)
      If (ncb > 0) Then
        Call gsend(comm, coordbuff, 0, comm%idnode)
      Endif
      Deallocate (coordbuff)
    Endif

    Do i = 1, config%natms
      crd%adfcoordlist(0, i) = crd%coordlist(0, i)
      Do j = 1, crd%coordlist(0, i)
        crd%adfcoordlist(j, i) = crd%coordlist(j, i)
        crd%coordlist(j, i) = config%ltg(crd%coordlist(j, i))
      End Do
    End Do

    !Create icoordlist
    If (flow%step == crd%coordstart) Then
      crd%icoordlist = crd%coordlist
    End If

  End Subroutine init_coord_list

  Subroutine checkcoord(config, crd, sites, flow, stats, comm)
    Type(configuration_type), Intent(In   ) :: config
    Type(coord_type),         Intent(InOut) :: crd
    Type(site_type),          Intent(In   ) :: sites
    Type(flow_type),          Intent(In   ) :: flow
    Type(stats_type),         Intent(InOUT) :: stats
    Type(comms_type),         Intent(InOut) :: comm

    Character(len=80)              :: aux
    Character(len=80), Allocatable :: rbuff(:)
    Integer                        :: defectcnt, defn, i, ii, j, k, totdefectcnt
    Integer, Allocatable           :: buff(:)
    Logical                        :: coordchange, coordfound, thisopen
    Real(kind=wp)                  :: rdis

    If (.not. crd%coordon) Return
    If (crd%ncoordpairs == 0) Return
    If (crd%coordops == 0) Return
    If (crd%coordstart > flow%step) Return
    If (Mod(flow%step, crd%coordinterval) /= 0) Return

    defectcnt = 0
    crd%defectlist(:) = 0

    Do i = 1, config%natms
      Do ii = 1, crd%ncoorddis
        If ((config%ltype(i) == crd%disltype(ii))) Then
          rdis = crd%discuts(ii)

        Endif
      End Do

      If (stats%rsd(i) > rdis) Then
        coordchange = .false.
        coordfound = .false.
        If (crd%coordlist(0, i) /= crd%icoordlist(0, i)) Then
          coordchange = .true.
        Else
          Do j = 1, crd%coordlist(0, i)
            coordfound = .false.
            Do k = 1, crd%icoordlist(0, i)
              If (crd%coordlist(j, i) == crd%icoordlist(k, i)) Then
                coordfound = .true.
              Endif
            Enddo
            If (.not. coordfound) Then
              coordchange = .true.
            Endif
          Enddo
        End If
        If (coordchange) Then
          defectcnt = defectcnt + 1
          crd%defectlist(0) = defectcnt
          crd%defectlist(defectcnt) = i
        Endif
      Endif
    Enddo

    If (comm%idnode == 0) Then
      totdefectcnt = 0
      Do j = 1, comm%mxnode - 1
        Call grecv(comm, defn, j, j)
        totdefectcnt = totdefectcnt + defn
      End Do
      totdefectcnt = totdefectcnt + crd%defectlist(0)

      Inquire (Unit=nccrdt, opened=thisopen)
      If (.not. thisopen) Then
        Open (Unit=nccrdt, File='CCOORD', Form='formatted')
      End If
      If (flow%step <= crd%coordstart) Then
        Write (Unit=nccrdt, Fmt='(a60)') config%cfgname(1:60)
        Write (Unit=nccrdt, Fmt='(a20,I10)') 'Number of frames', (flow%run_steps - crd%coordstart) / crd%coordinterval + 1
      Endif
      Write (Unit=nccrdt, Fmt='(A30,I10,I10,f20.6)') 'Number of coordination changes', totdefectcnt, flow%step, flow%time

      Do i = 0, 2
        Write (Unit=nccrdt, fmt='( 3f20.10 )') &
          config%cell(1 + i * 3), config%cell(2 + i * 3), config%cell(3 + i * 3)
      Enddo

      Do i = 1, crd%defectlist(0)
        Write (Unit=nccrdt, Fmt='(a2,I10,3f20.10,f10.5)') &
          Trim(sites%unique_atom(config%ltype(crd%defectlist(i)))), config%ltg(crd%defectlist(i)), &
          config%parts(crd%defectlist(i))%xxx, config%parts(crd%defectlist(i))%yyy, config%parts(crd%defectlist(i))%zzz, &
          stats%rsd(crd%defectlist(i))
      Enddo

      Do j = 1, comm%mxnode - 1
        Call grecv(comm, defectcnt, j, j)
        If (defectcnt > 0) Then
          Allocate (buff(2 * defectcnt))
          Allocate (rbuff(defectcnt))
          Call grecv(comm, buff, j, j)
          Call grecv(comm, rbuff, j, j)
          Do i = 1, defectcnt
            Write (nccrdt, Fmt='(a2,I10,a)') Trim(sites%unique_atom(buff(2 * i - 1))), buff(2 * i), &
              rbuff(i)
          End Do
          Deallocate (buff)
          Deallocate (rbuff)
        End If
      End Do

    Else
      defectcnt = crd%defectlist(0)
      defn = crd%defectlist(0)
      Allocate (buff(2 * crd%defectlist(0)))
      Allocate (rbuff(crd%defectlist(0)))
      Call gsend(comm, defn, 0, comm%idnode)
      Do i = 1, defectcnt
        buff(2 * i - 1) = config%ltype(crd%defectlist(i))
        buff(2 * i) = config%ltg(crd%defectlist(i))
        Write (aux, '(3f20.10,f10.5)') config%parts(crd%defectlist(i))%xxx, config%parts(crd%defectlist(i))%yyy, &
          config%parts(crd%defectlist(i))%zzz, stats%rsd(crd%defectlist(i))
        rbuff(i) = aux
      End Do
      Call gsend(comm, defectcnt, 0, comm%idnode)
      If (defectcnt > 0) Then
        Call gsend(comm, buff, 0, comm%idnode)
        Call gsend(comm, rbuff, 0, comm%idnode)
      End If
      Deallocate (buff)
      Deallocate (rbuff)
    End If

  End Subroutine checkcoord

End Module coord
