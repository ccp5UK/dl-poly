Subroutine constraints_pseudo_bonds(lstopt,dxx,dyy,dzz,gxx,gyy,gzz,engcon)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for treating constraint bonds as stiff harmonic
! springs for use with the conjugate gradient method (minimise_relax.f90)
!
! copyright - daresbury laboratory
! author    - w.smith & i.t.todorov november 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,    Only : mxnode,gsum
  Use setup_module
  Use config_module,   Only : natms,lfrzn
  Use constraints_module

  Implicit None

  Integer,           Intent( In    ) :: lstopt(0:2,1:mxcons)
  Real( Kind = wp ), Intent( In    ) :: dxx(1:mxcons),dyy(1:mxcons),dzz(1:mxcons)
  Real( Kind = wp ), Intent( InOut ) :: gxx(1:mxatms),gyy(1:mxatms),gzz(1:mxatms)
  Real( Kind = wp ), Intent(   Out ) :: engcon

  Real( Kind = wp ), Parameter :: rigid=1.0e6_wp

  Integer           :: k,i,j
  Real( Kind = wp ) :: r,r0,ebond,gamma

  engcon=0.0_wp
  Do k=1,ntcons
     If (lstopt(0,k) == 0) Then
        i=lstopt(1,k)
        j=lstopt(2,k)

! if a pair is frozen and constraint bonded, it is more frozen
! than constrained (users!!!)

        r=Sqrt(dxx(k)**2+dyy(k)**2+dzz(k)**2)
        r0=prmcon(listcon(0,k))

        gamma=rigid*(r-r0)
        ebond=gamma*0.5_wp*(r-r0)
        gamma=gamma/r

! Accumulate energy and add forces

        If (i <= natms) Then
           engcon=engcon+ebond

           If (lfrzn(i) == 0) Then
              gxx(i)=gxx(i)-dxx(k)*gamma
              gyy(i)=gyy(i)-dyy(k)*gamma
              gzz(i)=gzz(i)-dzz(k)*gamma
           End If
        End If

        If (j <= natms .and. lfrzn(j) == 0) Then
           gxx(j)=gxx(j)+dxx(k)*gamma
           gyy(j)=gyy(j)+dyy(k)*gamma
           gzz(j)=gzz(j)+dzz(k)*gamma
        End If
     End If
  End Do

! global sum of energy

  If (mxnode > 1) Call gsum(engcon)

End Subroutine constraints_pseudo_bonds
