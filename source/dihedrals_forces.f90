Subroutine dihedrals_forces(imcon,engdih,virdih,stress, &
           rcut,rvdw,keyfce,alpha,epsq,engcpe,vircpe,engsrp,virsrp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating dihedral energy and force terms
!
! Note: scale factors for reduces electrostatic and vdw 1-4 interactions
!       assumes 1-4 interactions are in the exclude list
!
! copyright - daresbury laboratory
! author    - w.smith march 1992
! amended   - i.t.todorov february 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,     Only : idnode,mxnode,gsync,gsum,gcheck
  Use setup_module,     Only : nrite,mxdihd,mxgrid, &
                               pi,sqrpi,r4pie0,zero_plus
  Use config_module,    Only : cell,natms,nlast,lsi,lsa,ltg,lfrzn,ltype, &
                               chge,xxx,yyy,zzz,fxx,fyy,fzz
  Use dihedrals_module, Only : ntdihd,keydih,listdih,prmdih
  Use vdw_module,       Only : ntpvdw,lstvdw,vvdw,gvdw

  Implicit None

  Integer,                             Intent( In    ) :: imcon
  Real( Kind = wp ),                   Intent(   Out ) :: engdih,virdih
  Real( Kind = wp ), Dimension( 1:9 ), Intent( InOut ) :: stress
  Real( Kind = wp ),                   Intent( In    ) :: rcut,rvdw, &
                                                          alpha,epsq
  Integer,                             Intent( In    ) :: keyfce
  Real( Kind = wp ),                   Intent( InOut ) :: engcpe,vircpe, &
                                                          engsrp,virsrp

  Real( Kind = wp ), Parameter :: aa1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: aa2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: aa3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: aa4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: aa5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp  =  0.3275911_wp

  Logical,           Save :: newjob = .true. , damp
  Real( Kind = wp ), Save :: dlrpot,twopi,rtwopi, aa,bb, rfld0,rfld1

  Logical                 :: safe(1:3)
  Integer                 :: fail(1:4),i,j,ia,ib,ic,id,kk,k,ka,kb,l,local_index
  Real( Kind = wp )       :: xab,yab,zab, xac,yac,zac,                 &
                             xad,yad,zad,rad,rad2,rrad,rrad2,          &
                             xbc,ybc,zbc,rrbc,                         &
                             xcd,ycd,zcd,                              &
                             fax,fay,faz, fb1x,fb1y,fb1z,              &
                             fcx,fcy,fcz, fd1x,fd1y,fd1z,              &
                             fx,fy,fz,                                 &
                             pbx,pby,pbz,pb2,rpb1,rpb2,                &
                             pcx,pcy,pcz,pc2,rpc1,rpc2,                &
                             pbpc,cost,sint,rsint,theta,theta0,dtheta, &
                             a,a0,a1,a2,a3,d,m,term,pterm,gamma,       &
                             scale,chgprd,coul,fcoul,                  &
                             exp1,tt,erc,fer,b0,                       &
                             vk,vk1,vk2,gk,gk1,gk2,t1,t2,ppp,          &
                             engc14,virc14,engs14,virs14,              &
                             strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:5)

  Logical,           Allocatable :: lunsafe(:)
  Integer,           Allocatable :: lstopt(:,:)
  Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
  Real( Kind = wp ), Allocatable :: xdbc(:),ydbc(:),zdbc(:)
  Real( Kind = wp ), Allocatable :: xdcd(:),ydcd(:),zdcd(:)

  fail=0
  Allocate (lunsafe(1:mxdihd),lstopt(0:4,1:mxdihd),       Stat=fail(1))
  Allocate (xdab(1:mxdihd),ydab(1:mxdihd),zdab(1:mxdihd), Stat=fail(2))
  Allocate (xdbc(1:mxdihd),ydbc(1:mxdihd),zdbc(1:mxdihd), Stat=fail(3))
  Allocate (xdcd(1:mxdihd),ydcd(1:mxdihd),zdcd(1:mxdihd), Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'dihedrals_forces allocation failure, node: ', idnode
     Call error(0)
  End If


  If (newjob) Then
     newjob = .false.

! Define constants

     dlrpot= rvdw/Real(mxgrid-4,wp)
     twopi = 2.0_wp*pi
     rtwopi= 1.0_wp/twopi

! Check for damped force-shifted coulombic and reaction field interactions
! and set force and potential shifting parameters dependingly

     damp=.false.
     If (alpha > zero_plus) Then
        damp=.true.

        exp1= Exp(-(alpha*rcut)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rcut)

        erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1/rcut
        fer = (erc + 2.0_wp*(alpha/sqrpi)*exp1)/rcut**2

        aa  = fer*rcut
        bb  = -(erc + aa*rcut)
     Else If (keyfce == 8) Then
        aa =  1.0_wp/rcut**2
        bb = -2.0_wp/rcut ! = -(1.0_wp/rcut+aa*rcut)
     End If

! set reaction field terms for RFC

     If (keyfce == 10) Then
        b0    = 2.0_wp*(epsq - 1.0_wp)/(2.0_wp*epsq + 1.0_wp)
        rfld0 = b0/rcut**3
        rfld1 = (1.0_wp + 0.5_wp*b0)/rcut
     End If
  End If

! calculate atom separation vectors

  Do i=1,ntdihd
     lunsafe(i)=.false.

! indices of dihedral atoms

     ia=local_index(listdih(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
     ib=local_index(listdih(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib
     ic=local_index(listdih(3,i),nlast,lsi,lsa) ; lstopt(3,i)=ic
     id=local_index(listdih(4,i),nlast,lsi,lsa) ; lstopt(4,i)=id

     lstopt(0,i)=0
     If (ia > 0 .and. ib > 0 .and. ic > 0 .and. id > 0) Then !Tag
        If (lfrzn(ia)*lfrzn(ib)*lfrzn(ic)*lfrzn(id) == 0) Then
           If (ia <= natms .or. ib <= natms .or. ic <= natms .or. id <= natms) Then
              lstopt(0,i)=1
            End If
        End If
     Else                                                    ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= natms) .or.   &
              (ib > 0 .and. ib <= natms) .or.   &
              (ic > 0 .and. ic <= natms) .or.   &
              (id > 0 .and. id <= natms)) .and. &
             (ia == 0 .or. ib == 0 .or. ic == 0 .or. id == 0) ) lunsafe(i)=.true.
     End If

! define components of bond vectors

     If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ia)-xxx(ib)
        ydab(i)=yyy(ia)-yyy(ib)
        zdab(i)=zzz(ia)-zzz(ib)

        xdbc(i)=xxx(ib)-xxx(ic)
        ydbc(i)=yyy(ib)-yyy(ic)
        zdbc(i)=zzz(ib)-zzz(ic)

        xdcd(i)=xxx(ic)-xxx(id)
        ydcd(i)=yyy(ic)-yyy(id)
        zdcd(i)=zzz(ic)-zzz(id)
!     Else ! (DEBUG)
!        xdab(i)=0.0_wp
!        ydab(i)=0.0_wp
!        zdab(i)=0.0_wp
!
!        xdbc(i)=0.0_wp
!        ydbc(i)=0.0_wp
!        zdbc(i)=0.0_wp
!
!        xdcd(i)=0.0_wp
!        ydcd(i)=0.0_wp
!        zdcd(i)=0.0_wp
     End If
  End Do

! Check for uncompressed units

  safe(1) = .not. Any(lunsafe(1:ntdihd))
  If (mxnode > 1) Call gcheck(safe(1))
  If (.not.safe(1)) Then
     Do j=0,mxnode-1
        If (idnode == j) Then
           Do i=1,ntdihd
              If (lunsafe(i)) Write(nrite,'(/,1x,a,2(i10,a))')     &
                 '*** warning - global unit number', listdih(0,i), &
                 ' , with a head particle number', listdih(1,i),   &
                 ' contributes towards next error !!! ***'
           End Do
        End If
        Call gsync()
     End Do
     Call error(132)
  End If

! Initialise safety flags

  safe=.true.

! periodic boundary condition

  Call images(imcon,cell,ntdihd,xdab,ydab,zdab)
  Call images(imcon,cell,ntdihd,xdbc,ydbc,zdbc)
  Call images(imcon,cell,ntdihd,xdcd,ydcd,zdcd)

! initialise stress tensor accumulators

  strs1=0.0_wp
  strs2=0.0_wp
  strs3=0.0_wp
  strs5=0.0_wp
  strs6=0.0_wp
  strs9=0.0_wp

! zero dihedral energy accumulator

  engdih=0.0_wp
  virdih=0.0_wp

! zero scaled 1-4 electrostatic and short-range potential accumulators

  engc14=0.0_wp
  virc14=0.0_wp
  engs14=0.0_wp
  virs14=0.0_wp

! loop over all specified dihedrals

  Do i=1,ntdihd
     If (lstopt(0,i) > 0) Then

! indices of dihedral atoms

        ia=lstopt(1,i)
        ib=lstopt(2,i)
        ic=lstopt(3,i)
        id=lstopt(4,i)

! define components of bond vectors

        xab=xdab(i)
        yab=ydab(i)
        zab=zdab(i)

        xbc=xdbc(i)
        ybc=ydbc(i)
        zbc=zdbc(i)
        rrbc=1.0_wp/Sqrt(xbc*xbc+ybc*ybc+zbc*zbc)

        xcd=xdcd(i)
        ycd=ydcd(i)
        zcd=zdcd(i)

        xac=xab+xbc
        yac=yab+ybc
        zac=zab+zbc

! construct first dihedral vector

        pbx=yab*zbc-zab*ybc
        pby=zab*xbc-xab*zbc
        pbz=xab*ybc-yab*xbc

        pb2=pbx*pbx+pby*pby+pbz*pbz

        rpb1=1.0_wp/Sqrt(pb2)
        rpb2=rpb1*rpb1

! construct second dihedral vector

        pcx=ybc*zcd-zbc*ycd
        pcy=zbc*xcd-xbc*zcd
        pcz=xbc*ycd-ybc*xcd

        pc2=pcx*pcx+pcy*pcy+pcz*pcz

        rpc1=1.0_wp/Sqrt(pc2)
        rpc2=rpc1*rpc1

! determine dihedral angle

        pbpc=pbx*pcx+pby*pcy+pbz*pcz
        cost=pbpc*rpb1*rpc1
        If (Abs(cost) > 1.0_wp) cost=Sign(1.0_wp,cost)
        sint=( xbc*(pcy*pbz-pcz*pby) + ybc*(pbx*pcz-pbz*pcx) + &
               zbc*(pcx*pby-pcy*pbx) ) * (rpb1*rpc1*rrbc)

        theta=Atan2(sint,cost)

! avoid singularity in sint

        sint=Sign( Max(1.0e-10_wp,Abs(sint)) , sint )
        rsint=1.0_wp/sint

! selection of potential energy function type

        kk=listdih(0,i)

! calculate potential energy and scalar force term

        If      (keydih(kk) == 1) Then

! torsion dihedral potential

           a=prmdih(1,kk)
           d=prmdih(2,kk)
           m=prmdih(3,kk)

           term=m*theta-d

           pterm=a*(1.0_wp+Cos(term))
           gamma=-a*m*Sin(term)*rsint * rpb1*rpc1

        Else If (keydih(kk) == 2) Then

! harmonic improper dihedral

           a     =prmdih(1,kk)
           theta0=prmdih(2,kk)
           dtheta=theta-theta0
           dtheta=dtheta-Nint(dtheta*rtwopi)*twopi

           term  =a*dtheta

           pterm=0.5_wp*term*dtheta
           gamma=term*rsint * rpb1*rpc1

        Else If (keydih(kk) == 3) Then

! harmonic cosine dihedral (note sint is cancelled)

           a     =prmdih(1,kk)
           theta0=prmdih(2,kk)
           dtheta=Cos(theta)-Cos(theta0)

           term  =a*dtheta

           pterm=0.5_wp*term*dtheta
           gamma=-term * rpb1*rpc1

        Else If (keydih(kk) == 4) Then

! 3-term cosine dihedral

           a1=prmdih(1,kk)
           a2=prmdih(2,kk)
           a3=prmdih(3,kk)

           pterm=0.5_wp*(a1*(1.0_wp+Cos(theta)) +        &
                         a2*(1.0_wp-Cos(2.0_wp*theta)) + &
                         a3*(1.0_wp+Cos(3.0_wp*theta)) )
           gamma=-0.5_wp*(a1*Sin(theta) -               &
                          2.0_wp*a2*Sin(2.0_wp*theta) + &
                          3.0_wp*a3*Sin(3.0_wp*theta) )*rsint * rpb1*rpc1

        Else If (keydih(kk) == 5) Then

! ryckaert-bellemans potential
!
! reference: chem. phys. lett., vol. 30, p. 123 (1975)
! ATTENTION: Modified to have the transition configuration correspond
!            to theta=180 rather than theta=0 as in original form

           a=prmdih(1,kk)
           m=Cos(theta)

           pterm=a*( 1.116_wp      - 1.462_wp*m    - 1.578_wp*m**2 + &
                     0.368_wp*m**3 + 3.156_wp*m**4 + 3.788_wp*m**5 )
           gamma=a*( 1.462_wp       + 3.156_wp*m   - 1.104_wp*m**2 - &
                     12.624_wp*m**3 - 18.94_wp*m**4 ) * rpb1*rpc1

        Else If (keydih(kk) == 6) Then

! fluorinated ryckaert-bellemans potential
! reference: Rice at al., JCP 104, p. 2101 (1996)

           a=prmdih(1,kk)
           m=Cos(theta)
           d=Exp(-56.0_wp*(theta-pi)**2)
           term=-1083.04_wp*(theta-pi)*d

           pterm=a*( 3.55_wp      - 2.78_wp*m    - 3.56_wp*m**2  - &
                     1.64_wp*m**3 + 7.13_wp*m**4 + 12.84_wp*m**5 + &
                     9.67_wp*d )
           gamma=( a*(2.78_wp       + 7.12_wp*m   + 4.92_wp*m**2 - &
                      28.52_wp*m**3 - 64.2_wp*m**4) + term*rsint ) * rpb1*rpc1

        Else If (keydih(kk) == 7) Then

! opls cosine dihedral

           a0=prmdih(1,kk)
           a1=prmdih(2,kk)
           a2=prmdih(3,kk)
           a3=prmdih(6,kk)
           theta0=prmdih(7,kk)
           dtheta=theta-theta0

           pterm=a0 + 0.5_wp*( a1*(1.0_wp+Cos(dtheta)) +        &
                               a2*(1.0_wp-Cos(2.0_wp*dtheta)) + &
                               a3*(1.0_wp+Cos(3.0_wp*dtheta)) )
           gamma=-0.5_wp*(       a1*Sin(dtheta) -               &
                          2.0_wp*a2*Sin(2.0_wp*dtheta)  + &
                          3.0_wp*a3*Sin(3.0_wp*dtheta))*rsint * rpb1*rpc1

        Else

! flag undefined potential

           safe(1)=.false.
           pterm=0.0_wp
           gamma=0.0_wp

        End If

! calculate atomic forces

        fax = gamma*((-pcy*zbc+pcz*ybc)-pbpc*rpb2*(-pby*zbc+pbz*ybc))
        fay = gamma*(( pcx*zbc-pcz*xbc)-pbpc*rpb2*( pbx*zbc-pbz*xbc))
        faz = gamma*((-pcx*ybc+pcy*xbc)-pbpc*rpb2*(-pbx*ybc+pby*xbc))

        fcx = gamma*((-pcy*zab+pcz*yab)-pbpc*rpb2*(-pby*zab+pbz*yab))
        fcy = gamma*(( pcx*zab-pcz*xab)-pbpc*rpb2*( pbx*zab-pbz*xab))
        fcz = gamma*((-pcx*yab+pcy*xab)-pbpc*rpb2*(-pbx*yab+pby*xab))

        fb1x= gamma*((-pby*zcd+pbz*ycd)-pbpc*rpc2*(-pcy*zcd+pcz*ycd))
        fb1y= gamma*(( pbx*zcd-pbz*xcd)-pbpc*rpc2*( pcx*zcd-pcz*xcd))
        fb1z= gamma*((-pbx*ycd+pby*xcd)-pbpc*rpc2*(-pcx*ycd+pcy*xcd))

        fd1x= gamma*((-pby*zbc+pbz*ybc)-pbpc*rpc2*(-pcy*zbc+pcz*ybc))
        fd1y= gamma*(( pbx*zbc-pbz*xbc)-pbpc*rpc2*( pcx*zbc-pcz*xbc))
        fd1z= gamma*((-pbx*ybc+pby*xbc)-pbpc*rpc2*(-pcx*ybc+pcy*xbc))

        If (ia <= natms) Then

! sum of dihedral energy (dihedral virial is zero!!!)

           engdih=engdih+pterm

! stress tensor calculation for dihedral terms

           strs1 = strs1 + xab*fax + xbc*(fb1x-fcx) - xcd*fd1x
           strs2 = strs2 + yab*fax + ybc*(fb1x-fcx) - ycd*fd1x
           strs3 = strs3 + zab*fax + zbc*(fb1x-fcx) - zcd*fd1x
           strs5 = strs5 + yab*fay + ybc*(fb1y-fcy) - ycd*fd1y
           strs6 = strs6 + yab*faz + ybc*(fb1z-fcz) - ycd*fd1z
           strs9 = strs9 + zab*faz + zbc*(fb1z-fcz) - zcd*fd1z

           fxx(ia)=fxx(ia)+fax
           fyy(ia)=fyy(ia)+fay
           fzz(ia)=fzz(ia)+faz

        End If

        If (ib <= natms) Then

           fxx(ib)=fxx(ib)-fax-fcx+fb1x
           fyy(ib)=fyy(ib)-fay-fcy+fb1y
           fzz(ib)=fzz(ib)-faz-fcz+fb1z

        End If

        If (ic <= natms) Then

           fxx(ic)=fxx(ic)+fcx-fb1x-fd1x
           fyy(ic)=fyy(ic)+fcy-fb1y-fd1y
           fzz(ic)=fzz(ic)+fcz-fb1z-fd1z

        End If

        If (id <= natms) Then

           fxx(id)=fxx(id)+fd1x
           fyy(id)=fyy(id)+fd1y
           fzz(id)=fzz(id)+fd1z

        End If

        xad=xac+xcd
        yad=yac+ycd
        zad=zac+zcd

        rad2  = xad**2+yad**2+zad**2
        rrad2 = 1.0_wp/rad2
        rad   = Sqrt(rad2)
        rrad  = 1.0_wp/rad

! flag error if rad > cutoff

        If (rad > rcut) Then
           Write(*,*) i,ltg(ia),ltg(ib),ltg(ic),ltg(id),rcut,rad
           safe(2) = .false.
        End If

! initialise defaults for coulombic energy and force contributions

        coul =0.0_wp
        fcoul=0.0_wp

! 1-4 electrostatics: adjust by weighting factor
! assumes 1-4 interactions are in the exclude list and Rad < rcut

        scale=prmdih(4,kk)

! scaled charge product times dielectric constants

        chgprd=scale*chge(ia)*chge(id)*r4pie0/epsq
        If (Abs(chgprd) > zero_plus .and. keyfce > 0) Then

! Electrostatics by ewald sum

           If      (keyfce ==  2) Then

              coul = chgprd*rrad
              fcoul= coul*rrad2

! distance dependent dielectric

           Else If (keyfce ==  4) Then

              coul = chgprd*rrad2
              fcoul= 2.0_wp*coul*rrad2

! direct coulombic

           Else If (keyfce ==  6) Then

              coul = chgprd*rrad
              fcoul= coul*rrad2

! force shifted coulombic and reaction field

           Else If (keyfce ==  8 .or. keyfce == 10) Then

              If (damp) Then ! calculate damping contributions
                 exp1= Exp(-(alpha*rad)**2)
                 tt  = 1.0_wp/(1.0_wp+pp*alpha*rad)

                 erc = tt*(aa1+tt*(aa2+tt*(aa3+tt*(aa4+tt*aa5))))*exp1*rrad
                 fer = (erc + 2.0_wp*(alpha/sqrpi)*exp1)*rrad2

                 coul = chgprd*(erc + aa*rad + bb)
                 fcoul= chgprd*(fer - aa*rrad)
              End If

              If      (keyfce ==  8) Then ! force shifted coulombic
                 If (.not.damp) Then ! pure
                    coul = chgprd*(rrad + aa*rad + bb)
                    fcoul= chgprd*(rrad2 - aa)*rrad
                 Else                ! damped
                    coul = coul
                    fcoul= fcoul
                 End If
              Else If (keyfce == 10) Then ! reaction field
                 If (.not.damp) Then ! pure
                    coul = chgprd*(rrad + 0.5_wp*rfld0*rad2 - rfld1)
                    fcoul= chgprd*(rrad*rrad2 - rfld0)
                 Else                ! damped
                    coul = coul  + chgprd*(0.5_wp*rfld0*rad2 - rfld1)
                    fcoul= fcoul + chgprd*(- rfld0)
                 End If
              End If

           Else

              safe(3) = .false.

           End If

           fx = fcoul*xad
           fy = fcoul*yad
           fz = fcoul*zad

           If (ia <= natms) Then

! correction to electrostatic energy and virial

              engc14 = engc14 + coul
              virc14 = virc14 - fcoul*rad2

! calculate stress tensor

              strs1 = strs1 + xad*fx
              strs2 = strs2 + xad*fy
              strs3 = strs3 + xad*fz
              strs5 = strs5 + yad*fy
              strs6 = strs6 + yad*fz
              strs9 = strs9 + zad*fz

              fxx(ia) = fxx(ia) + fx
              fyy(ia) = fyy(ia) + fy
              fzz(ia) = fzz(ia) + fz

           End If

           If (id <= natms) Then

              fxx(id) = fxx(id) - fx
              fyy(id) = fyy(id) - fy
              fzz(id) = fzz(id) - fz

           End If

        End If

! 1-4 short-range (vdw) interactions: adjust by weighting factor
! assumes 1-4 interactions are in the exclude list and Rad < rvdw

        scale=prmdih(5,kk)

! initialise default force contribution

        gamma=0.0_wp

        If (ntpvdw > 0) Then

! atomic and potential function indices

           ka=Max(ltype(ia),ltype(id))
           kb=Min(ltype(ia),ltype(id))
           k=lstvdw((ka*(ka-1))/2+kb)

           If (Abs(scale*vvdw(1,k)) > zero_plus) Then

! apply truncation of potential (rad < rvdw)

              If (rad < rvdw) Then

! determine interpolation panel for force arrays

                 l=Int(rad/dlrpot)
                 ppp=rad/dlrpot-Real(l,wp)

! calculate interaction energy using 3-point interpolation

                 vk  = vvdw(l,k)
                 vk1 = vvdw(l+1,k)
                 vk2 = vvdw(l+2,k)

                 t1 = vk  + (vk1 - vk)*ppp
                 t2 = vk1 + (vk2 - vk1)*(ppp - 1.0_wp)

! calculate scaled 1-4 short-range potential energy contribution

                 If (ia <= natms) engs14=engs14 + (t1+(t2-t1)*ppp*0.5_wp)*scale

! calculate forces using 3-point interpolation

                 gk  = gvdw(l,k)
                 gk1 = gvdw(l+1,k)
                 gk2 = gvdw(l+2,k)

                 t1 = gk  + (gk1 - gk)*ppp
                 t2 = gk1 + (gk2 - gk1)*(ppp - 1.0_wp)

                 gamma = scale*(t1 +(t2-t1)*ppp*0.5_wp)*rrad2

                 fx = gamma*xad
                 fy = gamma*yad
                 fz = gamma*zad

                 If (ia <= natms) Then

! calculate scaled 1-4 short-range potential virial

                    virs14=virs14-gamma*rad2

! calculate stress tensor

                    strs1 = strs1 + xad*fx
                    strs2 = strs2 + xad*fy
                    strs3 = strs3 + xad*fz
                    strs5 = strs5 + yad*fy
                    strs6 = strs6 + yad*fz
                    strs9 = strs9 + zad*fz

                    fxx(ia)=fxx(ia)+fx
                    fyy(ia)=fyy(ia)+fy
                    fzz(ia)=fzz(ia)+fz

                 End If

                 If (id <= natms) Then

                    fxx(id)=fxx(id)-fx
                    fyy(id)=fyy(id)-fy
                    fzz(id)=fzz(id)-fz

                 End If

              End If

           End If

! End If of (ntpvdw > 0)

        End If

     End If
  End Do

! sum contributions to potentials

  If (mxnode > 1) Then

     buffer(1) = engdih
     buffer(2) = engc14
     buffer(3) = virc14
     buffer(4) = engs14
     buffer(5) = virs14

     Call gsum(buffer(1:5))

     engdih = buffer(1)
     engc14 = buffer(2)
     virc14 = buffer(3)
     engs14 = buffer(4)
     virs14 = buffer(5)

  End If

  engcpe = engcpe + engc14
  vircpe = vircpe + virc14
  engsrp = engsrp + engs14
  virsrp = virsrp + virs14

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

! check safety to continue

  If (mxnode > 1) Call gcheck(safe)
  If (.not.safe(1)) Call error(448)
  If (.not.safe(2)) Call error(445)
  If (.not.safe(3)) Call error(446)

  Deallocate (lunsafe,lstopt, Stat=fail(1))
  Deallocate (xdab,ydab,zdab, Stat=fail(2))
  Deallocate (xdbc,ydbc,zdbc, Stat=fail(3))
  Deallocate (xdcd,ydcd,zdcd, Stat=fail(4))
  If (Any(fail > 0)) Then
     Write(nrite,'(/,1x,a,i0)') 'dihedrals_forces deallocation failure, node: ', idnode
     Call error(0)
  End If

End Subroutine dihedrals_forces
