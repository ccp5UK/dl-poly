!> Core particle type module
!>
!> Copyright - Daresbury Laboratory
!>
!> Author A.B.G Chalk July 2018
Module particle
  Use kinds, Only: wi,&
                   wp

  Implicit None
  Private

  Type, Public :: corePart

    Real(Kind=wp)    :: xxx, yyy, zzz
    Real(Kind=wp)    :: fxx, fyy, fzz
    Real(Kind=wp)    :: chge
    Integer(Kind=wi) :: pad1, pad2

  End Type
End Module particle
