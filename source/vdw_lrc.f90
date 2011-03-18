Subroutine vdw_lrc(imcon,rvdw,elrc,virlrc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to evaluate vdw long-range corrections to
! pressure and energy in a 3D periodic system
!
! copyright - daresbury laboratory
! author    - t.forester may 1993
! amended   - i.t.todorov february 2011
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,  Only : idnode,mxnode,gsum
  Use setup_module
  Use site_module,   Only : ntpatm,numtyp
  Use config_module, Only : volm,natms,ltype,lfrzn
  Use vdw_module,    Only : ls_vdw,lstvdw,ltpvdw,prmvdw

  Implicit None

  Integer,           Intent( In    ) :: imcon
  Real( Kind = wp ), Intent( In    ) :: rvdw
  Real( Kind = wp ), Intent(   Out ) :: elrc,virlrc

  Integer           :: fail,i,j,k,ivdw
  Real( Kind = wp ) :: twopi,a,b,c,d,e0,n,m,r0,r,eps,sig, &
                       eadd,padd,denprd,plrc

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

  numfrz = 0
  Do i=1,natms
     k = ltype(i)
     If (lfrzn(i) /= 0) numfrz(k)=numfrz(k)+1.0_wp
  End Do
  If (mxnode > 1) Call gsum(numfrz(1:ntpatm))

! twopi factor

  twopi = 2.0_wp*pi

! Evaluate only for 3D periodic systems

  If (imcon /= 0 .and. imcon /= 6) Then
     ivdw = 0

     Do i=1,ntpatm
        Do j=1,i

           eadd = 0.0_wp
           padd = 0.0_wp

           ivdw = ivdw + 1
           k = lstvdw(ivdw)

           If (ltpvdw(k) == 0) Then

! tabulated energy and pressure lrc

              eadd = prmvdw(1,k)
              padd =-prmvdw(2,k)

           Else If (ltpvdw(k) == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

              a=prmvdw(1,k)
              b=prmvdw(2,k)
              r=rvdw

              eadd = a/(9.0_wp*r**9) - b/(3.0_wp*r**3)
              padd = 12.0_wp*a/(9.0_wp*r**9)- 6.0_wp*b/(3.0_wp*r**3)

           Else If (ltpvdw(k) == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              r  =rvdw

              eadd = 4.0_wp*eps*(sig**12/(9.0_wp*r**9) - sig**6/(3.0_wp*r**3))
              padd = 4.0_wp*eps*(12.0_wp*sig**12/(9.0_wp*r**9) - 2.0_wp*sig**6/(r**3))

           Else If (ltpvdw(k) == 3) Then

! n-m potential :: u={e0/(n-m)}*[m*(r0/r)^n-n*(d/r)^c]

              e0=prmvdw(1,k)
              n =prmvdw(2,k)
              m =prmvdw(3,k)
              r0=prmvdw(4,k)
              r =rvdw

              eadd = e0/(n-m)*( m*r0**n/((n-3.0_wp)*r**(n-3.0_wp)) - n*r0**m/((m-3.0_wp)*r**(m-3.0_wp)) )
              padd = e0/(n-m)*n*m*( r0**n/((n-3.0_wp)*r**(n-3.0_wp)) - r0**m/((m-3.0_wp)*r**(m-3.0_wp)) )

           Else If (ltpvdw(k) == 4) Then

! Buckingham exp-6 potential :: u=a*Exp(-r/rho)-c/r^6

              c=prmvdw(3,k)
              r=rvdw

              eadd = -c/(3.0_wp*r**3)
              padd = -2.0_wp*c/(r**3)

           Else If (ltpvdw(k) == 5) Then

! Born-Huggins-Meyer exp-6-8 potential :: u=a*Exp(b*(sig-r))-c/r^6-d/r^8

              c=prmvdw(4,k)
              d=prmvdw(5,k)
              r=rvdw

              eadd = -c/(3.0_wp*r**3) - d/(5.0_wp*r**5)
              padd = -2.0_wp*c/(r**3) - 8.0_wp*d/(5.0_wp*r**5)

           Else If (ltpvdw(k) == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a=prmvdw(1,k)
              b=prmvdw(2,k)
              r=rvdw

              eadd = a/(9.0_wp*r**9) - b/(7.0_wp*r**7)
              padd = 12.0_wp*a/(9.0_wp*r**9) - 10.0_wp*b/(7.0_wp*r**7)

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
