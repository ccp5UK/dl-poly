Subroutine vdw_direct_fs_generate(rvdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating force-shifted constant arrays for
! direct vdw evaluation
!
! copyright - daresbury laboratory
! amended   - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : zero_plus
  Use vdw_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rvdw

  Integer           :: ivdw,keypot
  Real( Kind = wp ) :: r0,r0rn,r0rm,r_6,sor6, &
                       rho,a,b,c,d,e0,kk,n,m, &
                       sig,eps,t1,t2,t3

! allocate arrays for force-shifted corrections

  Call allocate_vdw_direct_fs_arrays()

! construct arrays for all types of vdw potential

  Do ivdw=1,ntpvdw

     keypot=ltpvdw(ivdw)
     If      (keypot == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        r_6=rvdw**(-6)

        afs(ivdw) = 6.0_wp*r_6*(2.0_wp*a*r_6-b)
        bfs(ivdw) =-r_6*(a*r_6-b) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        sor6=(sig/rvdw)**6

        afs(ivdw) = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)
        bfs(ivdw) =-4.0_wp*eps*sor6*(sor6-1.0_wp) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0=prmvdw(1,ivdw)
        n =prmvdw(2,ivdw)
        m =prmvdw(3,ivdw)
        r0=prmvdw(4,ivdw)

        a=r0/rvdw
        b=1.0_wp/(n-m)
        r0rn=a**n
        r0rm=a**m

        afs(ivdw) = e0*m*n*(r0rn-r0rm)*b
        bfs(ivdw) =-e0*(m*r0rn-n*r0rm)*b - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

        a  =prmvdw(1,ivdw)
        rho=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)

        If (Abs(rho) <= zero_plus) Then
           If (Abs(a) <= zero_plus) Then
              rho=1.0_wp
           Else
              Call error(467)
           End If
        End If

        b=rvdw/rho
        t1=a*Exp(-b)
        t2=-c/rvdw**6

        afs(ivdw) = (t1*b+6.0_wp*t2)
        bfs(ivdw) =-(t1+t2) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a  =prmvdw(1,ivdw)
        b  =prmvdw(2,ivdw)
        sig=prmvdw(3,ivdw)
        c  =prmvdw(4,ivdw)
        d  =prmvdw(5,ivdw)

        t1=a*Exp(b*(sig-rvdw))
        t2=-c/rvdw**6
        t3=-d/rvdw**8

        afs(ivdw) = (t1*rvdw*b+6.0_wp*t2+8.0_wp*t3)
        bfs(ivdw) =-(t1+t2+t3) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        t1=a/rvdw**12
        t2=-b/rvdw**10

        afs(ivdw) = (12.0_wp*t1+10.0_wp*t2)
        bfs(ivdw) =-(t1+t2) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

     Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0=prmvdw(1,ivdw)
        r0=prmvdw(2,ivdw)
        kk=prmvdw(3,ivdw)

        t1=Exp(-kk*(rvdw-r0))

        afs(ivdw) =-2.0_wp*e0*kk*t1*(1.0_wp-t1)*rvdw
        bfs(ivdw) =-e0*t1*(t1-2.0_wp) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 9) Then

! Weeks-Chandler-Anderson (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        d  =prmvdw(3,ivdw)

        sor6=(sig/(rvdw-d))**6

        afs(ivdw) = (24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(rvdw-d))*rvdw
        bfs(ivdw) =-(4.0_wp*eps*sor6*(sor6-1.0_wp)+eps) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 10) Then ! all zeroed in vdw_module

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

!       afs(ivdw) = 0.0_wp !initialised in vdw_module
!       bfs(ivdw) = 0.0_wp !initialised in vdw_module

     Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        rho=sig/rvdw
        t1=1.0_wp/(0.07_wp+rho)
        t2=1.0_wp/(0.12_wp+rho**7)
        t3=eps*(1.07_wp/t1**7)

        afs(ivdw) =-7.0_wp*t3*rho*(((1.12_wp/t2)-2.0_wp)/t1 + (1.12_wp/t2**2)*rho**6)
        bfs(ivdw) =-t3*((1.12_wp/t2)-2.0_wp) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else

        Call error(150)

     End If

  End Do

End Subroutine vdw_direct_fs_generate
