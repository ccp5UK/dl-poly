Module currents
  Use comms,         Only: comms_type,&
                           gsum
  Use configuration, Only: configuration_type
  Use constants,     Only: czero
  Use filename,      Only: file_type
  Use kinds,         Only: wp
  Use site,          Only: site_type

  Implicit None

  Type, Public :: current_type

    Complex(Kind=wp), Allocatable :: jlk(:, :, :)
    Integer                       :: nkpoints, lag
    Integer                       :: file_handle = -2
    Logical                       :: on = .false., io_yaml = .false.

  Contains
    Private
    Procedure, Public :: init
    Procedure, Public :: compute
    Final             :: cleanup
  End Type

Contains

  Subroutine init(T, nk, lag, fcurrent, comm, config, types, io_yaml)
    Class(current_type)                               :: T
    Integer,                            Intent(In   ) :: nk, lag, types
    Type(file_type),                    Intent(InOut) :: fcurrent
    Type(comms_type),                   Intent(In   ) :: comm
    Type(configuration_type),           Intent(In   ) :: config
    Logical,                  Optional, Intent(In   ) :: io_yaml


    Allocate (T%jlk(nk, 3, types))
    T%nkpoints = nk
    T%lag = lag
    If (Present(io_yaml)) Then 
      T%io_yaml = io_yaml
    Else 
      T%io_yaml = .false.
    End IF
    If (comm%idnode == 0) Then
      Open (Newunit=fcurrent%unit_no, File=fcurrent%filename, Status='unknown', Action="Write")
      If (T%io_yaml) Then
        Write (fcurrent%unit_no,'(a)') "%YAML 1.2"
        Write (fcurrent%unit_no,'(a)') "---"
        Write (fcurrent%unit_no, '(a,a)') "title: ", Trim(config%cfgname)
        Write (fcurrent%unit_no, '(a)') "timesteps: "
      End If
      T%file_handle = fcurrent%unit_no
    End If
  End Subroutine init

  Subroutine compute(T, config, time, comm, sites)
    Class(current_type),      Intent(InOut) :: T
    Type(configuration_type), Intent(In   ) :: config
    Real(Kind=wp),            Intent(In   ) :: time
    Type(comms_type),         Intent(InOut) :: comm
    Type(site_type),          Intent(In   ) :: sites

    Complex(Kind=wp) :: h(3, sites%mxatyp)
    Integer          :: i, k, atype
    Real(Kind=wp)    :: tmp

    Do k = 1, config%k%n
      T%jlk(k, :, :) = czero
      h = czero
      Do i = 1, config%natms
        atype = config%ltype(i)
        tmp = Dot_product(config%k%r(:, k), [config%parts(i)%xxx, config%parts(i)%yyy, config%parts(i)%zzz])
        h(:, atype) = h(:, atype) + config%weight(i) * [config%vxx(i), config%vyy(i), config%vzz(i)] * Exp(Cmplx(0.0_wp, tmp, wp))
      End Do
      !current%jlk(:,k,j)=kp%u(:,k)*Dot_product(kp%u(:,k),h)
      Call gsum(comm, h)
      T%jlk(k, :, :) = h
    End Do

    If (comm%idnode == 0) Then
      If (T%io_yaml) Then 
        Write (T%file_handle, '(a, g16.8, a)') "  - { time: ", time, ","
        Write (T%file_handle, '(a)')           "      atoms: {"
      End If
      Do atype = 1, sites%mxatyp
        If (T%io_yaml) Then 
          Write (T%file_handle, '(3a, *(g16.8, ","))', advance="no") &
            "      ", Trim(sites%unique_atom(atype)), ": [", T%jlk(:, :, atype)
          If (atype < sites%mxatyp) Then
            Write (T%file_handle, '(a)') "],"
          Else 
            Write (T%file_handle, '(a)') "]"
            Write (T%file_handle, '(a)') "      }"
            Write (T%file_handle, '(a)') "    }"
          End If
        Else
          Write (T%file_handle, '(g16.8,3a,*(g16.8, ","))') time, ", ", Trim(sites%unique_atom(atype)), ", ", T%jlk(:, :, atype)
        End If
      End Do
    End If

  End Subroutine compute

  Subroutine cleanup(T)
    Type(current_type) :: T

    If (Allocated(T%jlk)) Deallocate (T%jlk)
    !Close(t%file_handle)
  End Subroutine cleanup

End Module currents

