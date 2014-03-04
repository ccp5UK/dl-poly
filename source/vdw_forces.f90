Subroutine vdw_forces &
           (iatm,rvdw,xdf,ydf,zdf,rsqdf,engvdw,virvdw,stress)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating vdw energy and force terms using
! verlet neighbour list
!
! copyright - daresbury laboratory
! author    - w.smith august 1998
! amended   - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module
  Use config_module, Only : natms,ltg,ltype,list,fxx,fyy,fzz
  Use vdw_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rvdw
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: xdf,ydf,zdf,rsqdf
  Real( Kind = wp ),                        Intent(   Out ) :: engvdw,virvdw
  Real( Kind = wp ), Dimension( 1:9 ),      Intent( InOut ) :: stress

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: dlrpot,rdr,rcsq

  Integer           :: mm,idi,ai,aj,jatm,key,k,l,ityp
  Real( Kind = wp ) :: rsq,rrr,ppp,gamma,eng,          &
                       r0,r0rn,r0rm,r_6,sor6,          &
                       rho,a,b,c,d,e0,kk,              &
                       n,m,rc,sig,eps,alpha,beta,      &
                       fix,fiy,fiz,fx,fy,fz,           &
                       gk,gk1,gk2,vk,vk1,vk2,t1,t2,t3, &
                       strs1,strs2,strs3,strs5,strs6,strs9

  If (newjob) Then
     newjob = .false.

! define grid resolution for potential arrays and interpolation spacing

     dlrpot = rvdw/Real(mxgrid-4,wp)
     rdr    = 1.0_wp/dlrpot

! set cutoff condition

     rcsq   = rvdw**2
  End If

! initialise potential energy and virial

  engvdw=0.0_wp
  virvdw=0.0_wp

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! global identity of iatm

  idi=ltg(iatm)

! start of primary loop for forces evaluation

  ai=ltype(iatm)

! load forces

  fix=fxx(iatm)
  fiy=fyy(iatm)
  fiz=fzz(iatm)

  Do mm=1,list(0,iatm)

! atomic and potential function indices

     jatm=list(mm,iatm)
     aj=ltype(jatm)

     If (ai > aj) Then
        key=ai*(ai-1)/2 + aj
     Else
        key=aj*(aj-1)/2 + ai
     End If

     k=lstvdw(key)

! interatomic distance

     rsq = rsqdf(mm)

! validity and truncation of potential

     ityp=ltpvdw(k)
     If (ityp >= 0 .and. rsq < rcsq) Then

! Zero energy and force components

        eng   = 0.0_wp
        gamma = 0.0_wp

! Get separation distance

        rrr = Sqrt(rsq)

        If (ld_vdw) Then ! direct calculation

           If      (ityp == 1) Then

! 12-6 potential :: u=a/r^12-b/r^6

              a=prmvdw(1,k)
              b=prmvdw(2,k)

              r_6=rrr**(-6)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = r_6*(a*r_6-b)
              gamma = 6.0_wp*r_6*(2.0_wp*a*r_6-b)/rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
              End If

           Else If (ityp == 2) Then

! Lennard-Jones potential :: u=4*eps*[(sig/r)^12-(sig/r)^6]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)

              sor6=(sig/rrr)**6

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)
              gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
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

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*(m*r0rn-n*r0rm)*b
              gamma = e0*m*n*(r0rn-r0rm)*b/rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
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

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2
              gamma = (t1*b+6.0_wp*t2)/rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
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

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2+t3
              gamma = (t1*rrr*b+6.0_wp*t2+8.0_wp*t3)/rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
              End If

           Else If (ityp == 6) Then

! Hydrogen-bond 12-10 potential :: u=a/r^12-b/r^10

              a=prmvdw(1,k)
              b=prmvdw(2,k)

              t1=a/rrr**12
              t2=-b/rrr**10

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t1+t2
              gamma = (12.0_wp*t1+10.0_wp*t2)/rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
              End If

           Else If (ityp == 7) Then

! shifted and force corrected n-m potential (w.smith) ::

              e0=prmvdw(1,k)
              n =prmvdw(2,k)
              m =prmvdw(3,k)
              r0=prmvdw(4,k)
              rc=prmvdw(5,k) ; If (rc < 1.0e-6_wp) rc=rvdw

              If (n <= m) Call error(470)

              b=1.0_wp/(n-m)
              c=rc/r0 ; If (c < 1.0_wp) Call error(468)

              beta = c*( (c**(m+1.0_wp)-1.0_wp) / (c**(n+1.0_wp)-1.0_wp) )**b
              alpha= -(n-m) / (  m*(beta**n)*(1.0_wp+(n/c-n-1.0_wp)/c**n) &
                                -n*(beta**m)*(1.0_wp+(m/c-m-1.0_wp)/c**m) )
              e0 = e0*alpha

              If (rrr <= rc) Then
                 a=r0/rrr

                 If (jatm <= natms .or. idi < ltg(jatm))           &
                    eng   = e0*(  m*(beta**n)*(a**n-(1.0_wp/c)**n) &
                                 -n*(beta**m)*(a**m-(1.0_wp/c)**m) &
                                 +n*m*((rrr/rc-1.0_wp)*((beta/c)**n-(beta/c)**m)) )*b
                 gamma = e0*m*n*( (beta**n)*a**n-(beta**m)*a**m &
                                  -rrr/rc*((beta/c)**n-(beta/c)**m) )*b/rsq
              End If

           Else If (ityp == 8) Then

! Morse potential :: u=e0*{[1-Exp(-kk(r-r0))]^2-1}

              e0=prmvdw(1,k)
              r0=prmvdw(2,k)
              kk=prmvdw(3,k)

              t1=Exp(-kk*(rrr-r0))

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = e0*t1*(t1-2.0_wp)
              gamma = -2.0_wp*e0*kk*t1*(1.0_wp-t1)/rrr

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
              End If

           Else If (ityp == 9) Then

! Weeks-Chandler-Anderson (shifted & truncated Lenard-Jones) (i.t.todorov)
! :: u=4*eps*[{sig/(r-d)}^12-{sig/(r-d)}^6]-eps

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)
              d  =prmvdw(3,k)

              If (rrr < prmvdw(4,k) .or. Abs(rrr-d) < 1.0e-10_wp) Then ! Else leave them zeros
                 sor6=(sig/(rrr-d))**6

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = 4.0_wp*eps*sor6*(sor6-1.0_wp)+eps
                 gamma = 24.0_wp*eps*sor6*(2.0_wp*sor6-1.0_wp)/(rrr*(rrr-d))

                 If (ls_vdw) Then ! force-shifting
                    If (jatm <= natms .or. idi < ltg(jatm)) &
                    eng   = eng + afs(k)*rrr + bfs(k)
                    gamma = gamma - afs(k)/rrr
                 End If
              End If

           Else If (ityp == 10) Then

! DPD potential - Groot-Warren (standard) :: u=(1/2).a.r.(1-r/rc)^2

              a =prmvdw(1,k)
              rc=prmvdw(2,k)

              If (rrr < rc) Then ! Else leave them zeros
                 t2=rrr/rc
                 t1=0.5_wp*a*rrr*(1.0_wp-t2)

                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = t1*(1.0_wp-t2)
                 gamma = t1*(3.0_wp*t2-1.0_wp)/rsq
              End If

           Else If (ityp == 11) Then

! AMOEBA 14-7 :: u=eps * [1.07/((sig/r)+0.07)]^7 * [(1.12/((sig/r)^7+0.12))-2]

              eps=prmvdw(1,k)
              sig=prmvdw(2,k)

              rho=sig/rrr
              t1=1.0_wp/(0.07_wp+rho)
              t2=1.0_wp/(0.12_wp+rho**7)
              t3=eps*(1.07_wp/t1**7)

              If (jatm <= natms .or. idi < ltg(jatm)) &
              eng   = t3*((1.12_wp/t2)-2.0_wp)
              gamma =-7.0_wp*t3*rho*(((1.12_wp/t2)-2.0_wp)/t1 + (1.12_wp/t2**2)*rho**6)/rsq

              If (ls_vdw) Then ! force-shifting
                 If (jatm <= natms .or. idi < ltg(jatm)) &
                 eng   = eng + afs(k)*rrr + bfs(k)
                 gamma = gamma - afs(k)/rrr
              End If

           Else If (Abs(vvdw(1,k)) > zero_plus) Then ! potential read from TABLE - (ityp == 0)

              l   = Int(rrr*rdr)
              ppp = rrr*rdr - Real(l,wp)

! calculate interaction energy using 3-point interpolation

              If (jatm <= natms .or. idi < ltg(jatm)) Then
                 vk  = vvdw(l,k)
                 vk1 = vvdw(l+1,k)
                 vk2 = vvdw(l+2,k)

                 t1 = vk  + (vk1 - vk )*ppp
                 t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

                 eng = t1 + (t2-t1)*ppp*0.5_wp
                 If (ls_vdw) eng = eng + gvdw(mxgrid,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgrid,k) ! force-shifting
              End If

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

           If (jatm <= natms .or. idi < ltg(jatm)) Then
              vk  = vvdw(l,k)
              vk1 = vvdw(l+1,k)
              vk2 = vvdw(l+2,k)

              t1 = vk  + (vk1 - vk )*ppp
              t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

              eng = t1 + (t2-t1)*ppp*0.5_wp
              If (ls_vdw) eng = eng + gvdw(mxgrid,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgrid,k) ! force-shifting
           End If

! calculate forces using 3-point interpolation

           gk  = gvdw(l,k)
           gk1 = gvdw(l+1,k)
           gk2 = gvdw(l+2,k)

           t1 = gk  + (gk1 - gk )*ppp
           t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

           gamma = (t1 + (t2-t1)*ppp*0.5_wp)/rsq
           If (ls_vdw) gamma = gamma - gvdw(mxgrid,k)/(rrr*rvdw) ! force-shifting

        End If

! calculate forces

        fx = gamma*xdf(mm)
        fy = gamma*ydf(mm)
        fz = gamma*zdf(mm)

        fix=fix+fx
        fiy=fiy+fy
        fiz=fiz+fz

        If (jatm <= natms) Then

           fxx(jatm)=fxx(jatm)-fx
           fyy(jatm)=fyy(jatm)-fy
           fzz(jatm)=fzz(jatm)-fz

        End If

        If (jatm <= natms .or. idi < ltg(jatm)) Then

! add interaction energy

           engvdw = engvdw + eng

! add virial

           virvdw = virvdw - gamma*rsq

! add stress tensor

           strs1 = strs1 + xdf(mm)*fx
           strs2 = strs2 + xdf(mm)*fy
           strs3 = strs3 + xdf(mm)*fz
           strs5 = strs5 + ydf(mm)*fy
           strs6 = strs6 + ydf(mm)*fz
           strs9 = strs9 + zdf(mm)*fz

        End If

     End If

  End Do

! load back forces

  fxx(iatm)=fix
  fyy(iatm)=fiy
  fzz(iatm)=fiz

! complete stress tensor

  stress(1) = stress(1) + strs1
  stress(2) = stress(2) + strs2
  stress(3) = stress(3) + strs3
  stress(4) = stress(4) + strs2
  stress(5) = stress(5) + strs5
  stress(6) = stress(6) + strs6
  stress(7) = stress(7) + strs3
  stress(8) = stress(8) + strs6
  stress(9) = stress(9) + strs9

End Subroutine vdw_forces
