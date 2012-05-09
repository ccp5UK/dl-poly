Subroutine dihedrals_14_vdw(rvdw,ai,aj,rad,rad2,eng,gamma)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating dihedrals 1-4 vdw interaction:
! adjust by weighting factor
!
! copyright - daresbury laboratory
! amended   - i.t.todorov march 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : mxgrid,zero_plus
  Use vdw_module

  Implicit None

  Integer,           Intent( In    ) :: ai,aj
  Real( Kind = wp ), Intent( In    ) :: rvdw,rad,rad2
  Real( Kind = wp ), Intent(   Out ) :: eng,gamma

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: dlrpot,rdr

  Integer           :: key,k,l,ityp
  Real( Kind = wp ) :: rsq,rrr,ppp,            &
                       r0,r0rn,r0rm,r_6,sor6,  &
                       rho,a,b,c,d,e0,kk,      &
                       n,m,sig,eps,alpha,beta, &
                       gk,gk1,gk2,vk,vk1,vk2,t1,t2,t3

  If (newjob) Then
     newjob = .false.

! define grid resolution for potential arrays and interpolation spacing

     dlrpot = rvdw/Real(mxgrid-4,wp)
     rdr    = 1.0_wp/dlrpot
  End If

! Zero energy and force components

  eng   = 0.0_wp
  gamma = 0.0_wp

! potential function indices

  If (ai > aj) Then
     key=ai*(ai-1)/2 + aj
  Else
     key=aj*(aj-1)/2 + ai
  End If

  k=lstvdw(key)

! validity of potential

  ityp=ltpvdw(k)
  If (ityp >= 0) Then

! Get separation distance

     rrr = rad
     rsq = rad2

     If (ld_vdw) Then ! direct calculation

        If      (ityp == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

           a=prmvdw(1,k)
           b=prmvdw(2,k)

           r_6=rrr**(-6)

           eng   = r_6*(a*r_6-b)
           gamma = 6.0_wp*r_6*(2.0_wp*a*r_6-b)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma + afs(k)/rrr
           End If

        Else If (ityp == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

           eps=prmvdw(1,k)
           sig=prmvdw(2,k)

           sor6=(sig/rrr)**6

           eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)
           gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma + afs(k)/rrr
           End If

        Else If (ityp == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

           e0=prmvdw(1,k)
           n =prmvdw(2,k)
           m =prmvdw(3,k)
           r0=prmvdw(4,k)

           a=r0/rrr
           b=1.0_wp/(n-m)
           r0rn=a**n
           r0rm=a**m

           eng   = e0*(m*r0rn-n*r0rm)*b
           gamma = e0*m*n*(r0rn-r0rm)*b/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma + afs(k)/rrr
           End If

        Else If (ityp == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

           a  =prmvdw(1,k)
           rho=prmvdw(2,k)
           c  =prmvdw(3,k)

           If (Abs(rho) <= zero_plus) Then
              If (Abs(a) <= zero_plus) Then
                 rho=1.0_wp
              Else
                 Call error(467)
              End If
           End If

           b=rrr/rho
           t1=a*Exp(-b)
           t2=-c/rrr**6

           eng   = t1+t2
           gamma = (t1*b+6.0_wp*t2)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma + afs(k)/rrr
           End If

        Else If (ityp == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

           a  =prmvdw(1,k)
           b  =prmvdw(2,k)
           sig=prmvdw(3,k)
           c  =prmvdw(4,k)
           d  =prmvdw(5,k)

           t1=a*Exp(b*(sig-rrr))
           t2=-c/rrr**6
           t3=-d/rrr**8

           eng   = t1+t2+t3
           gamma = (t1*rrr*b+6.0_wp*t2+8.0_wp*t3)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma + afs(k)/rrr
           End If

        Else If (ityp == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

           a=prmvdw(1,k)
           b=prmvdw(2,k)

           t1=a/rrr**12
           t2=-b/rrr**10

           eng   = t1+t2
           gamma = (12.0_wp*t1+10.0_wp*t2)/rsq

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma + afs(k)/rrr
           End If

        Else If (ityp == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

           e0=prmvdw(1,k)
           n =prmvdw(2,k)
           m =prmvdw(3,k)
           r0=prmvdw(4,k)

           If (n <= m) Call error(470)

           a=r0/rrr
           b=1.0_wp/(n-m)
           c=rvdw/r0 ; If (c < 1.0_wp) Call error(468)

           beta = c*( (c**(m+1.0_wp)-1.0_wp) / (c**(n+1.0_wp)-1.0_wp) )**b
           alpha= -(n-m) / ( m*(beta**n)*(1.0_wp+(n/c-n-1.0_wp)/c**n) &
                             -n*(beta**m)*(1.0_wp+(m/c-m-1.0_wp)/c**m) )
           e0 = e0*alpha

           eng   = e0*( m*(beta**n)*(a**n-(1.0_wp/c)**n)  &
                        -n*(beta**m)*(a**m-(1.0_wp/c)**m) &
                        +n*m*((rrr/rvdw-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
           gamma = e0*m*n*( (beta**n)*a**n-(beta**m)*a**m &
                            -rrr/rvdw*((beta/c)**n-(beta/c)**m) )*b/rsq

        Else If (ityp == 8) Then

! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}

           e0=prmvdw(1,k)
           r0=prmvdw(2,k)
           kk=prmvdw(3,k)

           t1=Exp(-kk*(rrr-r0))

           eng   = e0*t1*(t1-2.0_wp)
           gamma = -2.0_wp*e0*kk*t1*(1.0_wp-t1)/rrr

           If (ls_vdw) Then ! force-shifting
              eng   = eng + afs(k)*rrr + bfs(k)
              gamma = gamma + afs(k)/rrr
           End If

        Else If (ityp == 9) Then

! Weeks-Chandler-Anderson (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

           eps=prmvdw(1,k)
           sig=prmvdw(2,k)
           d  =prmvdw(3,k)

           If (rrr < prmvdw(4,k) .or. Abs(rrr-d) < 1.0e-10_wp) Then ! Else leave them zeros
              sor6=(sig/(rrr-d))**6

              eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
              gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(rrr*(rrr-d))

              If (ls_vdw) Then ! force-shifting
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma + afs(k)/rrr
              End If
           End If

        Else If (Abs(vvdw(1,k)) > zero_plus) Then ! potential read from TABLE - (ityp == 0)

           l   = Int(rrr*rdr)
           ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

           vk  = vvdw(l,k)
           vk1 = vvdw(l+1,k)
           vk2 = vvdw(l+2,k)

           t1 = vk  + (vk1 - vk )*ppp
           t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

           eng = t1 + (t2-t1)*ppp*0.5_wp
           If (ls_vdw) eng = eng - gvdw(mxgrid,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgrid,k) ! force-shifting

! calculate forces using 3-point interpolation

           gk  = gvdw(l,k)
           gk1 = gvdw(l+1,k)
           gk2 = gvdw(l+2,k)

           t1 = gk  + (gk1 - gk )*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           gamma = (t1 + (t2-t1)*ppp*0.5_wp)/rsq
           If (ls_vdw) gamma = gamma - gvdw(mxgrid,k)/(rrr*rvdw) ! force-shifting

        End If

     Else If (Abs(vvdw(1,k)) > zero_plus) Then ! no direct = fully tabulated calculation

        l   = Int(rrr*rdr)
        ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

        vk  = vvdw(l,k)
        vk1 = vvdw(l+1,k)
        vk2 = vvdw(l+2,k)

        t1 = vk  + (vk1 - vk )*ppp
        t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

        eng = t1 + (t2-t1)*ppp*0.5_wp
        If (ls_vdw) eng = eng - gvdw(mxgrid,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgrid,k) ! force-shifting

! calculate forces using 3-point interpolation

        gk  = gvdw(l,k)
        gk1 = gvdw(l+1,k)
        gk2 = gvdw(l+2,k)

        t1 = gk  + (gk1 - gk )*ppp
        t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

        gamma = (t1 + (t2-t1)*ppp*0.5_wp)/rsq
        If (ls_vdw) gamma = gamma - gvdw(mxgrid,k)/(rrr*rvdw) ! force-shifting

     End If

  End If

End Subroutine dihedrals_14_vdw