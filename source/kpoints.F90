Module kpoints

  Use comms,     Only: comms_type,&
                       gbcast
  Use constants, Only: epsilon_wp
  Use kinds,     Only: wp

  Implicit None
  Private

  Type, Public :: kpoints_type

    Integer                    :: n
    Real(Kind=wp), Allocatable :: r(:, :)
    Real(Kind=wp), Allocatable :: u(:, :)
    Integer                    :: uf

  Contains
    Private

    Procedure, Public :: init
    Final             :: cleanup
  End Type

Contains

  Subroutine init(T, filename, comm)
    Class(kpoints_type)             :: T
    Character(Len=*), Intent(In   ) :: filename
    Type(comms_type), Intent(InOut) :: comm

    Integer       :: i
    Real(Kind=wp) :: h

    If (comm%idnode == 0) Then
      Open (File=Trim(filename), Newunit=T%uf, Action="Read", Status="old")
      Read (t%uf, *) T%n
    End If

    Call gbcast(comm, T%n, 0)

    If (T%n > 0) Allocate (T%r(3, t%n))
    If (T%n > 0) Allocate (T%u(3, t%n))

    If (comm%idnode == 0) Then
      Do i = 1, t%n
        Read (t%uf, *) t%r(:, i)
      End Do
    End If
    Do i = 1, T%n
      Call gbcast(comm, T%r(:,i), 0)
      h = Norm2(T%r(:, i))
      If (Abs(h) > epsilon_wp) &
        T%u(:, i) = T%r(:, i) / h
    End Do
    If (comm%idnode == 0) Then
      Close (t%uf)
    End If

  End Subroutine init

  Subroutine cleanup(T)
    Type(kpoints_type) :: T

    If (Allocated(t%r)) Deallocate (t%r)
    If (Allocated(t%u)) Deallocate (t%u)

  End Subroutine cleanup

End Module kpoints
