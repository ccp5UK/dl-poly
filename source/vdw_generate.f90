Subroutine vdw_generate(rvdw)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating potential energy and force arrays
! for van der waals forces only
!
! copyright - daresbury laboratory
! author    - w.smith may 1992
! amended   - i.t.todorov november 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : mxgrid,zero_plus
  Use vdw_module

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rvdw

  Integer           :: i,ivdw,keypot
  Real( Kind = wp ) :: dlrpot,r,r0,r0rn,r0rm,r_6,sor6, &
                       rho,a,b,c,d,e0,k,n,m,sig,eps,   &
                       alpha,beta,t1,t2,t3

! allocate arrays for tabulating

  Call allocate_vdw_table_arrays()

! define grid resolution for potential arrays

  dlrpot=rvdw/Real(mxgrid-4,wp)

! construct arrays for all types of vdw potential

  Do ivdw=1,ntpvdw

     keypot=ltpvdw(ivdw)
     If      (keypot == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           r_6=r**(-6)

           vvdw(i,ivdw)=r_6*(a*r_6-b)
           gvdw(i,ivdw)=6.0_wp*r_6*(2.0_wp*a*r_6-b)
        End Do

     Else If (keypot == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           sor6=(sig/r)**6

           vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)
           gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)
        End Do

     Else If (keypot == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

        e0=prmvdw(1,ivdw)
        n =prmvdw(2,ivdw)
        m =prmvdw(3,ivdw)
        r0=prmvdw(4,ivdw)

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           a=r0/r
           b=1.0_wp/(n-m)
           r0rn=(a)**n
           r0rm=(a)**m

           vvdw(i,ivdw)=e0*(m*r0rn-n*r0rm)*b
           gvdw(i,ivdw)=e0*m*n*(r0rn-r0rm)*b
        End Do

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

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           b=r/rho
           t1=a*Exp(-b)
           t2=-c/r**6

           vvdw(i,ivdw)=t1+t2
           gvdw(i,ivdw)=t1*b+6.0_wp*t2
        End Do

     Else If (keypot == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

        a  =prmvdw(1,ivdw)
        b  =prmvdw(2,ivdw)
        sig=prmvdw(3,ivdw)
        c  =prmvdw(4,ivdw)
        d  =prmvdw(5,ivdw)

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           t1=a*Exp(b*(sig-r))
           t2=-c/r**6
           t3=-d/r**8

           vvdw(i,ivdw)=t1+t2+t3
           gvdw(i,ivdw)=t1*r*b+6.0_wp*t2+8.0_wp*t3
        End Do

     Else If (keypot == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

        a=prmvdw(1,ivdw)
        b=prmvdw(2,ivdw)

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           t1=a/r**12
           t2=-b/r**10

           vvdw(i,ivdw)=t1+t2
           gvdw(i,ivdw)=12.0_wp*t1+10.0_wp*t2
        End Do

     Else If (keypot == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

        e0=prmvdw(1,ivdw)
        n =prmvdw(2,ivdw)
        m =prmvdw(3,ivdw)
        r0=prmvdw(4,ivdw)

        If (n <= m) Call error(470)

        b=1.0_wp/(n-m)
        c = rvdw/r0 ; If (c < 1.0_wp) Call error(468)

        beta = c*( (c**(m+1.0_wp)-1.0_wp) / (c**(n+1.0_wp)-1.0_wp))**b
        alpha= -(n-m) / ( m*(beta**n)*(1.0_wp+(n/c-n-1.0_wp)/c**n) &
                         -n*(beta**m)*(1.0_wp+(m/c-m-1.0_wp)/c**m) )
        e0 = e0*alpha

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           a=r0/r

           vvdw(i,ivdw)=e0*( m*(beta**n)*(a**n-(1.0_wp/c)**n) &
                            -n*(beta**m)*(a**m-(1.0_wp/c)**m) &
                            +n*m*((r/rvdw-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
           gvdw(i,ivdw)=e0*m*n*( (beta**n)*a**n-(beta**m)*a**m &
                                -r/rvdw*((beta/c)**n-(beta/c)**m) )*b
        End Do

     Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

        e0=prmvdw(1,ivdw)
        r0=prmvdw(2,ivdw)
        k=prmvdw(3,ivdw)

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           t1=Exp(-k*(r-r0))

           vvdw(i,ivdw)=e0*((1.0_wp-t1)**2-1.0_wp)
           gvdw(i,ivdw)=-2.0_wp*r*e0*k*(1.0_wp-t1)*t1
        End Do

     Else If (keypot == 9) Then

! Weeks-Chandler-Anderson (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

        eps=prmvdw(1,ivdw)
        sig=prmvdw(2,ivdw)
        d  =prmvdw(3,ivdw)

        Do i=1,mxgrid
           r=Real(i,wp)*dlrpot

           If (r < prmvdw(4,ivdw) .or. Abs(r-d) < 1.0e-10_wp) Then ! Else leave them zeros
              sor6=(sig/(r-d))**6

              vvdw(i,ivdw)=4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
              gvdw(i,ivdw)=24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)*r/(r-d)
           End If
        End Do

     Else

        If (.not.lt_vdw) Call error(150)

     End If

  End Do

End Subroutine vdw_generate
