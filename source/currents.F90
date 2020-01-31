Module currents
  Use comms,         Only: comms_type,&
                           gsum
  Use configuration, Only: configuration_type
  Use constants,     Only: czero
  Use filename,      Only: file_type
  Use kinds,         Only: wp

  Implicit None

  Type, Public :: current_type

    Complex(Kind=wp), Allocatable :: jlk(:, :)
    Complex(Kind=wp), Allocatable :: c(:, :, :) ! lags,kpoints,xyz
    Complex(Kind=wp), Allocatable :: fc(:, :, :) ! lags,kpoints,xyz
    Integer                       :: nkpoints, lag
    Integer                       :: file_handle = -2
    Logical                       :: on = .false.

  Contains
    Private
    Procedure, Public :: init
    Procedure, Public :: compute
    Final             :: cleanup
  End Type

Contains

  Subroutine init(T, nk, lag, fcurrent, comm)
    Class(current_type)             :: T
    Integer,          Intent(In   ) :: nk, lag
    Type(file_type),  Intent(InOut) :: fcurrent
    Type(comms_type), Intent(In   ) :: comm

    Allocate (T%jlk(nk, 3))
    T%nkpoints = nk
    T%lag = lag
    If (comm%idnode == 0) Then
      Open (Newunit=fcurrent%unit_no, File=fcurrent%filename, Status='unknown', Action="Write")
      T%file_handle = fcurrent%unit_no
    End If
  End Subroutine init

  Subroutine compute(T, config, time, comm)
    Class(current_type),      Intent(InOut) :: T
    Type(configuration_type), Intent(In   ) :: config
    Real(Kind=wp),            Intent(In   ) :: time
    Type(comms_type),         Intent(InOut) :: comm

    Complex(Kind=wp) :: h(3)
    Integer          :: i, k
    Real(Kind=wp)    :: tmp

    Do k = 1, config%k%n
      T%jlk(k, :) = czero
      h = czero
      Do i = 1, config%natms
        tmp = Dot_product(config%k%r(:, k), [config%parts(i)%xxx, config%parts(i)%yyy, config%parts(i)%zzz])
        h = h + [config%vxx(i), config%vyy(i), config%vzz(i)] * Exp(Cmplx(0.0_wp, tmp, wp))
      End Do
      !current%jlk(:,k,j)=kp%u(:,k)*Dot_product(kp%u(:,k),h)
      Call gsum(comm, h)
      T%jlk(k, :) = h
    End Do

    If (comm%idnode == 0) Then
      Write (T%file_handle, '(g0.8,1x,*(g0.8,1x))') time, T%jlk(:, :)
    End If

  End Subroutine compute

  Subroutine cleanup(T)
    Type(current_type) :: T

    If (Allocated(T%jlk)) Deallocate (T%jlk)
    !Close(t%file_handle)
  End Subroutine cleanup

End Module currents

