Subroutine vdw_direct_fs_generate(rvdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating force-shifted constant arrays for
! direct vdw evaluation
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena september 2017 (ryd)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (zblb)
! contrib   - a.m.elena april 2018 (mlj/mbuc)
! contrib   - a.m.elena may 2018 (m126)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : zero_plus, r4pie0
  Use vdw_module
  Use m_zbl, Only : zbl,ab,zbls,zblb,mlj,mbuck,mlj126

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rvdw

  Integer           :: ivdw,keypot,n,m
  Real( Kind = wp ) :: r0,r0rn,r0rm,r_6,sor6,   &
                       rho,a,b,c,d,e0,kk,nr,mr, &
                       sig,eps,t1,t2,t3,z,dz,   &
                       z1,z2,rm,ic,k,ri

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
        n =Nint(prmvdw(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(prmvdw(3,ivdw)) ; mr=Real(m,wp)
        r0=prmvdw(4,ivdw)

        a=r0/rvdw
        b=1.0_wp/Real(n-m,wp)
        r0rn=a**n
        r0rm=a**m

        afs(ivdw) = e0*mr*nr*(r0rn-r0rm)*b
        bfs(ivdw) =-e0*(mr*r0rn-nr*r0rm)*b - afs(ivdw)
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

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
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

      Else If (keypot == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)

        sor6=(sig/rvdw)**6

        afs(ivdw) = 24.0_wp*eps*sor6*(2.0_wp*sor6-c)
        bfs(ivdw) =-4.0_wp*eps*sor6*(sor6-c) - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 13) Then

! Morse potential with twelve term:: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0=prmvdw(1,ivdw)
        r0=prmvdw(2,ivdw)
        kk=prmvdw(3,ivdw)
        c=prmvdw(4,ivdw)

        t1=Exp(-kk*(rvdw-r0))
        sor6 = c/rvdw**12

        afs(ivdw) =-2.0_wp*e0*kk*t1*(1.0_wp-t1)*rvdw + 12.0_wp*sor6
        bfs(ivdw) =-e0*t1*(t1-2.0_wp) + sor6 - afs(ivdw)
        afs(ivdw) = afs(ivdw)/rvdw

     Else If (keypot == 14) Then

! Morse potential with twelve term:: u=(a+b*r)Exp(-r/c)

        a = prmvdw(1,ivdw)
        b = prmvdw(2,ivdw)
        c = prmvdw(3,ivdw)

        kk=1.0_wp/c
        t1=Exp(-rvdw*kk)
        afs(ivdw) = (a+b*rvdw)*kk*t1-b*t1
        bfs(ivdw) = -(a*c+a*rvdw+b*rvdw*rvdw)*kk*t1

     Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        call zbl(rvdw,kk,a,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        ic = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0
        Call zbls(rvdw,kk,a,rm,ic,e0,k,r0,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        ic = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0
        Call zblb(rvdw,kk,a,rm,ic,e0,r0,k,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else If (keypot == 18) Then

! LJ tappered with MDF:: u=f(r)LJ(r)

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        ri=prmvdw(3,ivdw)

        Call mlj(rvdw,eps,sig,ri,rvdw,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else If (keypot == 19) Then

! Buckingham tappered with MDF:: u=f(r)Buck(r)
        a  =prmvdw(1,ivdw)
        rho=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)
        ri  =prmvdw(4,ivdw)
        Call mbuck(rvdw,A,rho,c,ri,rvdw,z,dz)

        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else If (keypot == 20) Then

! LJ tappered with MDF:: u=f(r)LJ(r)
         a = prmvdw(1,ivdw)
         b = prmvdw(2,ivdw)
        ri = prmvdw(3,ivdw)
        Call mlj126(rvdw,A,B,ri,rvdw,z,dz)
        afs(ivdw) = dz/rvdw
        bfs(ivdw) = -z-dz

     Else

        Call error(150)

     End If

  End Do

End Subroutine vdw_direct_fs_generate
