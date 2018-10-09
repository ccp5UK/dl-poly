Module particle

  Use kinds, Only : wp

  Type, Public :: corePart
    Real( Kind = wp ) :: xxx,yyy,zzz
    Real( Kind = wp ) :: fxx,fyy,fzz
    Real( Kind = wp ) :: chge
    Integer           :: pad1, pad2
  End Type


End Module particle
