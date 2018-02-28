Subroutine intra_mcoul(keyfce,rcut,alpha,epsq,iatm,jatm,scale, &
                      rrr,xdf,ydf,zdf,coul,virele,fx,fy,fz,safe)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for calculating bond's or 1-4 dihedral
! electrostatics: adjusted by a weighting factor
!
! copyright - daresbury laboratory
! amended   - i.t.todorov & h.a.boateng february 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup_module, Only : sqrpi,r4pie0,zero_plus,mxompl,mximpl
  Use config_module,Only : natms
  Use mpoles_module

  Implicit None

  Integer,           Intent( In    ) :: keyfce,iatm,jatm
  Real( Kind = wp ), Intent( In    ) :: rcut,alpha,epsq,scale
  Real( Kind = wp ), Intent( In    ) :: xdf,ydf,zdf,rrr
  Real( Kind = wp ), Intent(   Out ) :: coul,virele,fx,fy,fz
  Logical,           Intent( InOut ) :: safe

  Integer                 :: k1,k2,k3,s1,s2,s3,n
  Integer                 :: ks1,ks2,ks3,ks11,ks21,ks31,ii,jj

  Logical,           Save :: newjob = .true. , damp
  Real( Kind = wp ), Save :: aa     = 0.0_wp , &
                             bb     = 0.0_wp , &
                             rfld0  = 0.0_wp , &
                             rfld1  = 0.0_wp , &
                             rfld2  = 0.0_wp

  Real( Kind = wp ) :: exp1,tt,erc,fer,b0,      &
                       tix,tiy,tiz,tjx,tjy,tjz, &
                       talpha,alphan,tmp,tmpi,tmpj,t1,t2,sx,sy,sz,kx,ky,kz,txyz

  Real( Kind = wp ) :: d1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: b1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: a1(-2:2*mxompl+1,-2:2*mxompl+1,-2:2*mxompl+1)
  Real( Kind = wp ) :: imp(1:mximpl),jmp(1:mximpl)
  Real( Kind = wp ) :: impx(1:mximpl),impy(1:mximpl),impz(1:mximpl)
  Real( Kind = wp ) :: jmpx(1:mximpl),jmpy(1:mximpl),jmpz(1:mximpl)

  Real( Kind = wp ), Parameter :: aa1 =  0.254829592_wp
  Real( Kind = wp ), Parameter :: aa2 = -0.284496736_wp
  Real( Kind = wp ), Parameter :: aa3 =  1.421413741_wp
  Real( Kind = wp ), Parameter :: aa4 = -1.453152027_wp
  Real( Kind = wp ), Parameter :: aa5 =  1.061405429_wp
  Real( Kind = wp ), Parameter :: pp  =  0.3275911_wp

  If (newjob) Then
     newjob = .false.

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
        rfld2 = 0.5_wp*rfld0
     End If
  End If

! initialise defaults for coulombic energy and force contributions

   coul =0.0_wp ; virele=0.0_wp
   fx =0.0_wp ; fy =0.0_wp ; fz =0.0_wp
   tix=0.0_wp ; tiy=0.0_wp ; tiz=0.0_wp
   tjx=0.0_wp ; tjy=0.0_wp ; tjz=0.0_wp

! get the multipoles for sites i and j

  imp=mplgfr(:,iatm)
  jmp=mplgfr(:,jatm)

! ignore interaction if the charge is zero

  If (Maxval(Abs(imp)) <= zero_plus .or. Maxval(Abs(jmp)) <= zero_plus) Return

  If(mxompl > 0 .and. induce) Then

     imp(2)=imp(2)+indipx(iatm)
     imp(3)=imp(3)+indipy(iatm)
     imp(4)=imp(4)+indipz(iatm)

     jmp(2)=jmp(2)+indipx(jatm)
     jmp(3)=jmp(3)+indipy(jatm)
     jmp(4)=jmp(4)+indipz(jatm)

  End If

! get the components for site i and j infinitesimal rotations

  impx=mprotx(:,iatm)
  impy=mproty(:,iatm)
  impz=mprotz(:,iatm)

  jmpx=mprotx(:,jatm)
  jmpy=mproty(:,jatm)
  jmpz=mprotz(:,jatm)

! default convergence factor and derivative of 'r'

  talpha = 1.0_wp
  a1     = 0.0_wp

! scale imp multipoles

  imp = scale*imp*r4pie0/epsq

! Electrostatics by ewald sum = direct coulombic

  If      (keyfce ==  2 .or. keyfce ==  6) Then

! compute derivatives of 1/r kernel

     Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! distance dependent dielectric

  Else If (keyfce ==  4) Then

! Compute derivatives of 1/r^2 kernel

     Call coul_deriv(2,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! force shifted coulombic and reaction field

  Else If (keyfce ==  8 .or. keyfce == 10) Then

     If (damp) Then ! calculate damping contributions

! compute derivatives of 'r'

        Call coul_deriv(-1,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

! scale the derivatives of 'r'

        a1 = aa*a1

        exp1= Exp(-(alpha*rrr)**2)
        tt  = 1.0_wp/(1.0_wp+pp*alpha*rrr)

        fer = erc/alpha

! compute derivatives of the ewald real space kernel

        Call ewald_deriv(-2,2*mxompl+1,1,fer,alpha*xdf,alpha*ydf,alpha*zdf,alpha*rrr,d1)

! scale the derivatives into the right form

        d1 = 2.0_wp*alpha*d1/sqrpi

     End If

     If      (keyfce ==  8) Then ! force shifted coulombic
        If (.not.damp) Then ! pure

! compute derivatives of '1/r'

           Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! compute derivatives of 'r'

           Call coul_deriv(-1,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

! scale the derivatives of 'r' and add to d1

           d1 = d1 + aa*a1
           a1 = 0.0_wp

        Else                ! damped

           talpha = alpha

        End If
     Else If (keyfce == 10) Then ! reaction field
        If (.not.damp) Then ! pure

! compute derivatives of '1/r'

           Call coul_deriv(1,2*mxompl+1,xdf,ydf,zdf,rrr,d1)

! compute derivatives of 'r^2'

           Call coul_deriv(-2,2*mxompl+1,xdf,ydf,zdf,rrr,a1)

! scale the derivatives of 'r' and add to d1

           d1 = d1 + rfld2*a1
           a1 = 0.0_wp

        Else

! compute derivatives of 'r^2'

           Call coul_deriv(-2,2*mxompl+1,xdf,ydf,zdf,rrr,b1)

! scale the derivatives of 'r^2' and add to scaled derivatives of 'r'

           a1 = a1 + rfld2*b1

           talpha = alpha

        End If
     End If

  Else

     safe = .false.

  End If

  If (safe) Then

     kz = 1.0_wp
     Do k3=0,mxompl

        ky = kz
        Do k2=0,mxompl-k3

           kx = ky
           Do k1=0,mxompl-k3-k2

              jj = mplmap(k1,k2,k3)

              If (Abs(jmp(jj)) > zero_plus) Then

                 txyz=kx*jmp(jj)

                 sz = 1.0_wp
                 Do s3=0,mxompl
                    ks3=k3+s3; ks31=ks3+1

                    sy = sz
                    Do s2=0,mxompl-s3
                       ks2=k2+s2; ks21=ks2+1

                       sx = sy
                       Do s1=0,mxompl-s3-s2
                          ks1=k1+s1; ks11=ks1+1

                          n      = ks1+ks2+ks3
                          alphan = talpha**n

                          ii     = mplmap(s1,s2,s3)

                          tmp    = alphan * d1(ks1,ks2,ks3) + a1(ks1,ks2,ks3)

                          tmpi   = txyz       * tmp
                          tmpj   = sx*imp(ii) * tmp

                          t2     = txyz*imp(ii)
                          t1     = alphan*t2

!  energy

                          coul    = coul + t1*d1(ks1,ks2,ks3)  + t2*a1(ks1,ks2,ks3)

!  force
                          t1      = t1*talpha

                          fx      = fx   - t1*d1(ks11,ks2,ks3) + t2*a1(ks11,ks2,ks3)
                          fy      = fy   - t1*d1(ks1,ks21,ks3) + t2*a1(ks1,ks21,ks3)
                          fz      = fz   - t1*d1(ks1,ks2,ks31) + t2*a1(ks1,ks2,ks31)

!  torque on iatm

                          tix     = tix  + impx(ii)*tmpi
                          tiy     = tiy  + impy(ii)*tmpi
                          tiz     = tiz  + impz(ii)*tmpi

!  torque on jatm

                          tjx     = tjx  + jmpx(jj)*tmpj
                          tjy     = tjy  + jmpy(jj)*tmpj
                          tjz     = tjz  + jmpz(jj)*tmpj

                          sx = -sx
                       End Do

                       sy = -sy
                    End Do

                    sz = -sz
                 End Do

              End If

              kx = -kx

           End Do

           ky = -ky

        End Do

        kz = -kz

     End Do

     If (keyfce == 2 .or. keyfce == 4 .or. keyfce == 6) Then

        virele = -coul

     Else If (keyfce ==  8) Then ! force shifted coulombic

        virele = -(fx*xdf + fy*ydf + fz*zdf)

! shift potential

        tmp    = aa*rrr + bb
        coul   = coul   + tmp*imp(1)*jmp(1)

! shift torque

        tmpi   = tmp*jmp(1)
        tix    = tix    + impx(1)*tmpi
        tiy    = tiy    + impy(1)*tmpi
        tiz    = tiz    + impz(1)*tmpi

        tmpj   = tmp*imp(1)
        tjx    = tjx    + jmpx(1)*tmpj
        tjy    = tjy    + jmpy(1)*tmpj
        tjz    = tjz    + jmpz(1)*tmpj

     Else If (keyfce == 10) Then ! reaction field

        virele = -(fx*xdf + fy*ydf + fz*zdf)

! shift potential

        coul   = coul   - rfld1*imp(1)*jmp(1)

! shift torque

        tmpi   = -rfld1*jmp(1)
        tix    = tix    + impx(1)*tmpi
        tiy    = tiy    + impy(1)*tmpi
        tiz    = tiz    + impz(1)*tmpi

        tmpj   = -rfld1*imp(1)
        tjx    = tjx    + jmpx(1)*tmpj
        tjy    = tjy    + jmpy(1)*tmpj
        tjz    = tjz    + jmpz(1)*tmpj

     End If

     tix = tix*r4pie0/epsq
     tiy = tiy*r4pie0/epsq
     tiz = tiz*r4pie0/epsq

     If (iatm <= natms) Then

        mptrqx(iatm)=mptrqx(iatm)+tix
        mptrqy(iatm)=mptrqy(iatm)+tiy
        mptrqz(iatm)=mptrqz(iatm)+tiz

     End If

     If (jatm <= natms) Then

        mptrqx(jatm)=mptrqx(jatm)+tjx
        mptrqy(jatm)=mptrqy(jatm)+tjy
        mptrqz(jatm)=mptrqz(jatm)+tjz

     End If

  End If

End Subroutine intra_mcoul
