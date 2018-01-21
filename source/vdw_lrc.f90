Subroutine vdw_lrc(rvdw,elrc,virlrc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to evaluate vdw long-range corrections to
! pressure and energy in a 3D periodic system
!
! copyright - daresbury laboratory
! author    - t.forester may 1993
! amended   - i.t.todorov september 2016
! contrib   - a.m.elena september 2016 (ljc)
! contrib   - a.m.elena september 2017 (rydberg)
! contrib   - a.m.elena october 2017 (zbl/zbls)
! contrib   - a.m.elena december 2017 (zblb)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module
  Use site_module,   Only : ntpatm,numtyp
  Use config_module, Only : imcon,volm,natms,ltype,lfrzn
  Use vdw_module,    Only : ls_vdw,lstvdw,ltpvdw,prmvdw
  Use mm3_module
  Use m_zbl,         Only : ab, intRadZBL, intdRadZBL

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rvdw
  Real( Kind = wp ), Intent(   Out ) :: elrc,virlrc

  Integer           :: fail,i,j,k,ivdw,keypot,n,m
  Real( Kind = wp ) :: a,b,c,d,e0,nr,mr,r0,r,eps,sig, &
                       eadd,padd,denprd,plrc,t,kk,s9, &
                       z1,z2,rm,al

  Real( Kind = wp ), Dimension( : ), Allocatable :: numfrz

  fail=0
  Allocate (numfrz(mxatyp), Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_lrc allocation failure, node: ', idnode
     Call error(0)
  End If

! initialise long-range corrections to energy and pressure

  plrc = 0.0_wp
  elrc = 0.0_wp

  If (ls_vdw) Go To 10 ! force-shifting

! initialise counter arrays and evaluate number density in system

  numfrz = 0.0_wp
  Do i=1,natms
     k = ltype(i)
     If (lfrzn(i) /= 0) numfrz(k)=numfrz(k)+1.0_wp
  End Do
  If (mxnode > 1) Call gsum(numfrz(1:ntpatm))

! Evaluate only for 3D periodic systems

  If (imcon /= 0 .and. imcon /= 6) Then
     ivdw = 0

     Do i=1,ntpatm
        Do j=1,i

           eadd = 0.0_wp
           padd = 0.0_wp

           ivdw = ivdw + 1
           k = lstvdw(ivdw)

           keypot=ltpvdw(k)
           If      (keypot ==  0) Then

! tabulated energy and pressure lrc

              eadd = prmvdw(1,k)
              padd =-prmvdw(2,k)

           Else If (keypot ==  1) Then

! 12-6 potential :: u=a/r^12-b/r^6

              a=prmvdw(1,k)
              b=prmvdw(2,k)
              r=rvdw

              eadd = a/(9.0_wp*r**9) - b/(3.0_wp*r**3)
              padd = 12.0_wp*a/(9.0_wp*r**9) - 6.0_wp*b/(3.0_wp*r**3)

           Else If (keypot ==  2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              r  =rvdw

              eadd = 4.0_wp*eps*(sig**12/(9.0_wp*r**9) - sig**6/(3.0_wp*r**3))
              padd = 8.0_wp*eps*(6.0_wp*sig**12/(9.0_wp*r**9) - sig**6/(r**3))

           Else If (keypot ==  3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

              e0=prmvdw(1,k)
              n =Nint(prmvdw(2,k)) ; nr=Real(n,wp)
              m =Nint(prmvdw(3,k)) ; mr=Real(m,wp)
              r0=prmvdw(4,k)
              r =rvdw

              eadd = e0/(nr-mr)*( mr*r0**n/((nr-3.0_wp)*r**(n-3)) - nr*r0**m/((mr-3.0_wp)*r**(m-3)) )
              padd = e0/(nr-mr)*nr*mr*( r0**n/((nr-3.0_wp)*r**(n-3)) - r0**m/((mr-3.0_wp)*r**(m-3)) )

           Else If (keypot ==  4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

              c=prmvdw(3,k)
              r=rvdw

              eadd = -c/(3.0_wp*r**3)
              padd = -2.0_wp*c/(r**3)

           Else If (keypot ==  5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

              c=prmvdw(4,k)
              d=prmvdw(5,k)
              r=rvdw

              eadd = -c/(3.0_wp*r**3) - d/(5.0_wp*r**5)
              padd = -2.0_wp*c/(r**3) - 8.0_wp*d/(5.0_wp*r**5)

           Else If (keypot ==  6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a=prmvdw(1,k)
              b=prmvdw(2,k)
              r=rvdw

              eadd = a/(9.0_wp*r**9) - b/(7.0_wp*r**7)
              padd = 12.0_wp*a/(9.0_wp*r**9) - 10.0_wp*b/(7.0_wp*r**7)

           Else If (keypot == 8) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}

              e0=prmvdw(1,k)
              r0=prmvdw(2,k)
              kk=prmvdw(3,k)
              If (kk > Tiny(kk)) Then
                 t = Exp(-kk*(rvdw - r0))

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*rvdw+1)**2 + 1) + &
                    e0*t*t/(4.0_wp*kk*kk*kk)*((kk*rvdw+1)**2 + kk*kk*rvdw*rvdw)
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*rvdw**3 + &
                      3*kk**2*rvdw**2 +6*kk*rvdw + 6) + &
                      e0*t*t/(4.0_wp*kk*kk*kk)* & 
                      (4.0_wp*kk**3*rvdw**3 + 6*kk**2*rvdw**2 + 6*kk*rvdw + 3)
              End If

           Else If (keypot == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)

              a =0.07_wp
              b =0.12_wp
              e0=1.0e-12_wp

              eadd = intRadMM3(sig,a,b,eps,rvdw,e0)
              padd = -intRaddMM3(sig,a,b,eps,rvdw,e0)

           Else If (keypot ==  12) Then

! Lennard-Jones cohesive potential :: u=4*eps*[(sig/r)^12-c*(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              c  =prmvdw(3,k)
              r  =rvdw

              eadd = 4.0_wp*eps*(sig**12/(9.0_wp*r**9) - c*sig**6/(3.0_wp*r**3))
              padd = 8.0_wp*eps*(6.0_wp*sig**12/(9.0_wp*r**9) - c*sig**6/(r**3))

           Else If (keypot == 13) Then

! Morse potential :: u=e0*{[1-Exp(-k(r-r0))]^2-1}+c/r^12

              e0 = prmvdw(1,k)
              r0 = prmvdw(2,k)
              kk = prmvdw(3,k)
               c = prmvdw(4,k)

              If (kk>Tiny(kk)) Then

                 t = Exp(-kk*(rvdw - r0))
                 s9 = c/(9.0_wp*rvdw**9)

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*rvdw+1)**2 + 1) + &
                     e0*t*t/(4.0_wp*kk*kk*kk)*((kk*rvdw+1)**2 + & 
                     kk*kk*rvdw*rvdw) + s9
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*rvdw**3 + & 
                       3*kk**2*rvdw**2 + 6*kk*rvdw + 6) + & 
                       e0*t*t/(4.0_wp*kk*kk*kk)* (4.0_wp*kk**3*rvdw**3 + & 
                       6*kk**2*rvdw**2 + 6*kk*rvdw + 3) + 12.0_wp*s9
              End If

           Else If (keypot == 14) Then

! Rydberg potential:: u=(a+b*r)Exp(-r/c)

              a = prmvdw(1,k)
              b = prmvdw(2,k)
              c = prmvdw(3,k)
              t = exp(-rvdw/c)

              eadd = (b*c*rvdw**3+(3*b*c**2+a*c)*rvdw**2+(6*b*c**3+2*a*c**2)*rvdw&
                +6*b*c**4+2*a*c**3)*t
              padd = (b*rvdw**4+(3*b*c+a)*rvdw**3+(9*b*c**2+3*a*c)*rvdw**2+& 
                (18*b*c**3+6*a*c**2)*rvdw+18*b*c**4+6*a*c**3)*t

           Else If (keypot == 15) Then

! ZBL potential:: u=Z1Z2/(4πε0r)∑_{i=1}^4b_ie^{-c_i*r/a}

              z1 = prmvdw(1,k)
              z2 = prmvdw(2,k)

        ! this is in fact inverse a
              a = (z1**0.23_wp+z2**0.23_wp)/(ab*0.88534_wp)
              kk = z1*z2*r4pie0
              eadd = intRadZBL(kk,a,rvdw,1e-12_wp)
              padd = intdRadZBL(kk,a,rvdw,1e-12_wp)

           Else If (keypot == 16) Then

! ZBL swithched with Morse:: u=f(r)zbl(r)+(1-f(r))*morse(r)

              e0 = prmvdw(5,k)
              r0 = prmvdw(6,k)
              kk = prmvdw(7,k)

              If (kk > Tiny(kk)) Then
                 t = Exp(-kk*(rvdw - r0))

                 eadd = -2.0_wp*e0*t/(kk*kk*kk)*((kk*rvdw+1)**2 + 1) + &
                    e0*t*t/(4.0_wp*kk*kk*kk)*((kk*rvdw+1)**2 + kk*kk*rvdw*rvdw)
                 padd = -2.0_wp*e0*t/(kk*kk*kk)*(kk**3*rvdw**3 + &
                      3*kk**2*rvdw**2 +6*kk*rvdw + 6) + &
                      e0*t*t/(4.0_wp*kk*kk*kk)* & 
                      (4.0_wp*kk**3*rvdw**3 + 6*kk**2*rvdw**2 + 6*kk*rvdw + 3)
              End If

           Else If (keypot == 17) Then

! ZBL swithched with Buckingham:: u=f(r)zbl(r)+(1-f(r))*buckingham(r)

              A = prmvdw(5,k)
              r0 = prmvdw(6,k)
              c = prmvdw(7,k)

              t=A*Exp(-rvdw/r0)

              eadd = (rvdw**2+2*r0*rvdw+2*r0**2)*t*r0-c/(3.0_wp*rvdw**3)
              padd = (rvdw**3+3*r0*rvdw**2+6*r0**2*rvdw+6*r0**3)*t -2.0_wp*c/(rvdw**3)

           End If

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

           If (i /= j) Then
              eadd = eadd*2.0_wp
              padd = padd*2.0_wp
           End If

           denprd=twopi * (numtyp(i)*numtyp(j) - numfrz(i)*numfrz(j)) / volm**2

           elrc = elrc + volm*denprd*eadd
           plrc = plrc + denprd*padd/3.0_wp

        End Do
     End Do

  End If

10 Continue

  If (idnode == 0) Write(nrite,"(/,/,1x, &
     & 'long-range correction for: vdw energy  ',e15.6,/,26x, &
     & ': vdw pressure',e15.6)") elrc/engunit,plrc*prsunt

! convert plrc to a viral term

  virlrc = plrc*(-3.0_wp*volm)

  Deallocate (numfrz, Stat=fail)
  If (fail > 0) Then
     Write(nrite,'(/,1x,a,i0)') 'vdw_lrc deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine vdw_lrc
