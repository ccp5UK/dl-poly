Subroutine pmf_pseudo_bonds(indpmf,pxx,pyy,pzz,gxx,gyy,gzz,engpmf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for treating PMF constraints as stiff harmonic
! springs for use with the conjugate gradient method (minimise_relax.f90)
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov august 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,    Only : mxnode,gsum
  Use setup_module
  Use config_module,   Only : natms,lfrzn
  Use pmf_module

  Implicit None

  Integer,           Intent( In    ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent( In    ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)
  Real( Kind = wp ), Intent( InOut ) :: gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: engpmf

  Real( Kind = wp ), Parameter :: rigid=1.0e6_wp

  Integer           :: ipmf,jpmf,k,l
  Real( Kind = wp ) :: r,r0,ebond,gamma,tmp

  r0=prmpmf

  engpmf=0.0_wp
  Do ipmf=1,ntpmf

     r=Sqrt(pxx(ipmf)**2+pyy(ipmf)**2+pzz(ipmf)**2)

     gamma=rigid*(r-r0)/Real(mxtpmf(1)+mxtpmf(2),wp)
     ebond=gamma*0.5_wp*(r-r0)
     gamma=gamma/r

     Do jpmf=1,2
        tmp=Real(1-2*Mod(jpmf,2),wp)*gamma

! If this unit is present on my domain

        If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then

           Do k=1,mxtpmf(jpmf)
              l=indpmf(k,jpmf,ipmf)

! For domain particles

              If (l > 0 .and. l <= natms) Then

! Accumulate energy

                 engpmf=engpmf+ebond

! Add forces

                 If (lfrzn(l) == 0) Then
                    gxx(l)=gxx(l)+pxx(ipmf)*tmp
                    gyy(l)=gyy(l)+pyy(ipmf)*tmp
                    gzz(l)=gzz(l)+pzz(ipmf)*tmp
                 End If

              End If

           End Do

        End If

     End Do

  End Do

! Global sum of energy

  If (mxnode > 1) Call gsum(engpmf)

End Subroutine pmf_pseudo_bonds
