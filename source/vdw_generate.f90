Subroutine vdw_generate(rvdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating potential energy and force arrays
! for van der waals forces only
!
! copyright - daresbury laboratory
! author    - w.smith may 1992
! amended   - i.t.todorov march 2016
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena september 2017 (rydberg)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (zblb)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : mxgvdw,zero_plus,r4pie0
  Use vdw_module
  Use m_zbl, Only : ab,zbl,zbls,zblb

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rvdw

  Integer           :: i,ivdw,keypot,n,m
  Real( Kind = wp ) :: dlrpot,r,r0,r0rn,r0rm,r_6,sor6,  &
                       rho,a,b,c,d,e0,kk,nr,mr,rc,sig,eps, &
                       alpha,beta,t1,t2,t3,t,z1,z2,dphi,phi, &
                       rm,k

! allocate arrays for tabulating

  Call allocate_vdw_table_arrays()

! define grid resolution for potential arrays

  dlrpot=rvdw/Real(mxgvdw-4,wp)

! construct arrays for all types of vdw potential

  Do ivdw=1,ntpvdw

     keypot=ltpvdw(ivdw)
     If      (keypot == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           r_6=r**(-6)

           vvdw(i,ivdw)=r_6*(a*r_6-b)
           gvdw(i,ivdw)=6.0_wp*r_6*(2.0_wp*a*r_6-b)
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           If (a*b > zero_plus) Then
              sigeps(1,ivdw)=(a/b)**(1.0_wp/6.0_wp)
              sigeps(2,ivdw)=b**2/(4.0_wp*a)
           End If ! else leave undetermined
        End If

     Else If (keypot == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           sor6=(sig/r)**6

           vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)
           gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=sig
           sigeps(2,ivdw)=eps
        End If

     Else If (keypot == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0=prmvdw(1,ivdw)
        n =Nint(prmvdw(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(prmvdw(3,ivdw)) ; mr=Real(m,wp)
        r0=prmvdw(4,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           a=r0/r
           b=1.0_wp/(nr-mr)
           r0rn=(a)**n
           r0rm=(a)**m

           vvdw(i,ivdw)=e0*(mr*r0rn-nr*r0rm)*b
           gvdw(i,ivdw)=e0*mr*nr*(r0rn-r0rm)*b
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=r0*(nr/mr)**(1.0_wp/(mr-nr))
           sigeps(2,ivdw)=e0
        End If

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

! Sigma-epsilon initialisation

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           b=r/rho
           t1=a*Exp(-b)
           t2=-c/r**6

           vvdw(i,ivdw)=t1+t2
           gvdw(i,ivdw)=t1*b+6.0_wp*t2

! Sigma-epsilon search

           If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                       vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                      (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                       vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                       vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                    sigeps(2,ivdw)=-vvdw(i-1,ivdw)
              End If
           End If
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a  =prmvdw(1,ivdw)
        b  =prmvdw(2,ivdw)
        sig=prmvdw(3,ivdw)
        c  =prmvdw(4,ivdw)
        d  =prmvdw(5,ivdw)

! Sigma-epsilon initialisation

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=a*Exp(b*(sig-r))
           t2=-c/r**6
           t3=-d/r**8

           vvdw(i,ivdw)=t1+t2+t3
           gvdw(i,ivdw)=t1*r*b+6.0_wp*t2+8.0_wp*t3

! Sigma-epsilon search

           If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                       vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                      (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                       vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                       vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                    sigeps(2,ivdw)=-vvdw(i-1,ivdw)
              End If
           End If
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=a/r**12
           t2=-b/r**10

           vvdw(i,ivdw)=t1+t2
           gvdw(i,ivdw)=12.0_wp*t1+10.0_wp*t2
        End Do

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=Sqrt(a/b)
           sigeps(2,ivdw)=((b/6.0_wp)**6)*((5.0_wp/a)**5)
        End If
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

        e0=prmvdw(1,ivdw)
        n =Nint(prmvdw(2,ivdw)) ; nr=Real(n,wp)
        m =Nint(prmvdw(3,ivdw)) ; mr=Real(m,wp)
        r0=prmvdw(4,ivdw)
        rc=prmvdw(5,ivdw) ; If (rc < 1.0e-6_wp) rc=rvdw

        If (n <= m) Call error(470)

! Sigma-epsilon initialisation

        sigeps(1,ivdw)=-1.0_wp
        sigeps(2,ivdw)= 0.0_wp

        t=Real(n-m,wp)

        b=1.0_wp/t
        c = rc/r0 ; If (c < 1.0_wp) Call error(468)

        beta = c*( (c**(m+1)-1.0_wp) / (c**(n+1)-1.0_wp) )**b
        alpha= -t / (  mr*(beta**n)*(1.0_wp+(nr/c-nr-1.0_wp)/c**n) &
                      -nr*(beta**m)*(1.0_wp+(mr/c-mr-1.0_wp)/c**m) )
        e0 = e0*alpha

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot
           If (r <= rc) Then
              a=r0/r

              vvdw(i,ivdw)=e0*(  mr*(beta**n)*(a**n-(1.0_wp/c)**n) &
                                -nr*(beta**m)*(a**m-(1.0_wp/c)**m) &
                                +nr*mr*((r/rc-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
              gvdw(i,ivdw)=e0*mr*nr*( (beta**n)*a**n-(beta**m)*a**m &
                                    -r/rc*((beta/c)**n-(beta/c)**m) )*b

! Sigma-epsilon search

              If (i > 20) Then ! Assumes some safety against numeric black holes!!!
                 If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                       sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                          vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                         (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                          vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                          vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                       sigeps(2,ivdw)=-vvdw(i-1,ivdw)
                 End If
              End If
           End If ! The else condition is satisfied by the vdw_module initialisation
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0=prmvdw(1,ivdw)
        r0=prmvdw(2,ivdw)
        kk=prmvdw(3,ivdw)

        Do i=0,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=Exp(-kk*(r-r0))

           vvdw(i,ivdw)=e0*((1.0_wp-t1)**2-1.0_wp)
           gvdw(i,ivdw)=-2.0_wp*r*e0*kk*(1.0_wp-t1)*t1
        End Do
        t1=Exp(+kk*r0)
        gvdw(0,ivdw)=-2.0_wp*e0*kk*(1.0_wp-t1)*t1

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=r0-log(2.0_wp)/kk
           sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 9) Then

! Weeks-Chandler-Andersen (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        d  =prmvdw(3,ivdw)

! Sigma-epsilon initialisation

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=-1.0_wp
           sigeps(2,ivdw)= 0.0_wp
        End If

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           If (r < prmvdw(4,ivdw) .or. Abs(r-d) < 1.0e-10_wp) Then ! Else leave them zeros
              sor6=(sig/(r-d))**6

              vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
              gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*r/(r-d)

! Sigma-epsilon search

              If ((.not.ls_vdw) .and. i > 20) Then ! Assumes some safety against numeric black holes!!!
                 If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                    If (Nint(Sign(1.0_wp,vvdw(i-1,ivdw))) == -Nint(Sign(1.0_wp,vvdw(i,ivdw)))) &
                       sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
                 Else                                           ! find epsilon
                    If ( (vvdw(i-2,ivdw) >= vvdw(i-1,ivdw) .and.  &
                          vvdw(i-1,ivdw) <= vvdw(i  ,ivdw)) .and. &
                         (vvdw(i-2,ivdw) /= vvdw(i-1,ivdw) .or.   &
                          vvdw(i-2,ivdw) /= vvdw(i  ,ivdw) .or.   &
                          vvdw(i-1,ivdw) /= vvdw(i  ,ivdw)) )     &
                       sigeps(2,ivdw)=-vvdw(i-1,ivdw)
                 End If
              End If
           End If
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

     Else If (keypot == 10) Then

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.rc.(1-r/rc)^2

        a =prmvdw(1,ivdw)
        rc=prmvdw(2,ivdw)

        Do i=0,mxgvdw
           r=Real(i,wp)*dlrpot

           If (r < rc) Then
              t1=0.5_wp*a*rc
              t2=1.0_wp-r/rc

              vvdw(i,ivdw)=t1*t2**2
              gvdw(i,ivdw)=a*t2*r
           End If
        End Do
        gvdw(0,ivdw)=a

        sigeps(1,ivdw)=rc
        sigeps(2,ivdw)=a

     Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((r/sig)+0.07)]^7 * [(1.12/((r/sig)^7+0.12))-2]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           rho=r/sig
           t1=1.0_wp/(0.07_wp+rho)
           t2=1.0_wp/(0.12_wp+rho**7)
           t3=eps*(1.07_wp*t1)**7

           t=t3*((1.12_wp*t2) - 2.0_wp)

           vvdw(i,ivdw)=t
           gvdw(i,ivdw)=7.0_wp*(t1*t + 1.12_wp*t3*t2**2*rho**6)*rho
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=sig
           sigeps(2,ivdw)=eps
        End If

      Else If (keypot == 12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        c  =prmvdw(3,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           sor6=(sig/r)**6

           vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-c)
           gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-c)
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=sig
           sigeps(2,ivdw)=eps
        End If

     Else If (keypot == 13) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

        e0 = prmvdw(1,ivdw)
        r0 = prmvdw(2,ivdw)
        kk = prmvdw(3,ivdw)
        c  = prmvdw(4,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           t1=Exp(-kk*(r - r0))
           sor6=c/r**12

           vvdw(i,ivdw)=e0*((1.0_wp-t1)**2-1.0_wp)+sor6
           gvdw(i,ivdw)=-2.0_wp*r*e0*kk*(1.0_wp-t1)*t1+12.0_wp*sor6
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then !???
           sigeps(1,ivdw)=r0-log(2.0_wp)/kk
           sigeps(2,ivdw)=e0
        End If

     Else If (keypot == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)
        
        a = prmvdw(1,ivdw)
        b = prmvdw(2,ivdw)
        c = prmvdw(3,ivdw)

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           kk = r/c
           t1=Exp(-kk)           

           vvdw(i,ivdw) = (a+b*r)*t1
           gvdw(i,ivdw) = t1*kk*(a-b*c+b*r)
        End Do
        vvdw(0,ivdw)= a
        gvdw(0,ivdw)= 0

        If (.not.ls_vdw) Then !???
           sigeps(1,ivdw)=1.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        
        ! this is in fact inverse a
        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           call zbl(r,kk,a,phi,dphi)

           vvdw(i,ivdw) = phi
           gvdw(i,ivdw) = dphi

        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=0.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        c = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           Call zbls(r,kk,a,rm,c,e0,k,r0,phi,dphi)
           vvdw(i,ivdw) = phi
           gvdw(i,ivdw) = dphi
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=0.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

        z1 = prmvdw(1,ivdw)
        z2 = prmvdw(2,ivdw)
        rm = prmvdw(3,ivdw)
        c = 1.0_wp/prmvdw(4,ivdw)
        e0 = prmvdw(5,ivdw)
        r0 = prmvdw(6,ivdw)
        k = prmvdw(7,ivdw)

        a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
        kk = z1*z2*r4pie0

        Do i=1,mxgvdw
           r=Real(i,wp)*dlrpot

           Call zblb(r,kk,a,rm,c,e0,k,r0,phi,dphi)
           vvdw(i,ivdw) = phi
           gvdw(i,ivdw) = dphi
        End Do
        vvdw(0,ivdw)=Huge(vvdw(1,ivdw))
        gvdw(0,ivdw)=Huge(gvdw(1,ivdw))

        If (.not.ls_vdw) Then
           sigeps(1,ivdw)=0.0_wp
           sigeps(2,ivdw)=0.0_wp
        End If

     Else

        If (.not.lt_vdw) Call error(150)

     End If

     If (ls_vdw .and. (keypot /= 7 .and. keypot /= 10)) Then ! no shifting to shifted n-m and DPD

        sigeps(1,ivdw)=-1.0_wp
        sigeps(2,ivdw)= 0.0_wp

        Do i=1,mxgvdw-4
           t  = vvdw(i  ,ivdw) + gvdw(mxgvdw-4,ivdw)*(Real(i  ,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
           t1 = vvdw(i-1,ivdw) + gvdw(mxgvdw-4,ivdw)*(Real(i-1,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)

! Sigma-epsilon search

           If (i > 20) Then ! Assumes some safety against numeric black holes!!!
              If (Sign(1.0_wp,sigeps(1,ivdw)) < 0.0_wp) Then ! find sigma
                 If (Nint(Sign(1.0_wp,t1)) == -Nint(Sign(1.0_wp,t))) &
                    sigeps(1,ivdw)=(Real(i,wp)-0.5_wp)*dlrpot
              Else                                           ! find epsilon
                 t2 = vvdw(i-2,ivdw) + gvdw(mxgvdw-4,ivdw)*(Real(i-2,wp)*dlrpot/rvdw-1.0_wp) - vvdw(mxgvdw-4,ivdw)
                 If ( (t2 >= t1 .and. t1 <= t) .and.         &
                      (t2 /= t1 .or. t2 /= t .or. t1 /= t) ) &
                    sigeps(2,ivdw)=-t1
              End If
           End If
        End Do
     End If

! Needed to distinguish that something has been defined

     If (Abs(vvdw(0,ivdw)) <= zero_plus) vvdw(0,ivdw) = Sign(Tiny(vvdw(0,ivdw)),vvdw(0,ivdw))

  End Do

End Subroutine vdw_generate
