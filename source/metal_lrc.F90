Subroutine metal_lrc(rmet,elrcm,vlrcm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine to evaluate metal long-range corrections to
! pressure and energy in a 3D periodic system
!
! copyright - daresbury laboratory
! author    - w.smith june 1995
! amended   - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use comms_module,  Only : idnode
  Use setup_module
  Use site_module,   Only : ntpatm,unqatm,dens
  Use config_module, Only : imcon,volm
  Use metal_module,  Only : lstmet,ltpmet,prmmet

  Implicit None

  Real( Kind = wp ),                        Intent( In    ) :: rmet
  Real( Kind = wp ), Dimension( 0:mxatyp ), Intent(   Out ) :: elrcm,vlrcm

  Logical, Save     :: newjob = .true.
  Integer           :: i,j,k0,k1,k2,kmet,keypot,nnn,mmm

  Real( Kind = wp ) :: elrc0,elrc1,elrc2,elrcsum,vlrc0,vlrc1,vlrc2, tmp, &
                       eps,sig,nnnr,mmmr,ccc, aaa,rr0,ppp,zet,qqq,eee


! long-range corrections to energy, pressure and density

  elrcm   = 0.0_wp
  vlrcm   = 0.0_wp
  elrcsum = 0.0_wp

  If (imcon /= 0 .and. imcon /= 6) Then
     kmet = 0

     Do i=1,ntpatm
        Do j=1,i

           elrc0=0.0_wp
           elrc1=0.0_wp
           elrc2=0.0_wp

           vlrc0=0.0_wp
           vlrc1=0.0_wp
           vlrc2=0.0_wp

           kmet = kmet + 1
           k0 = lstmet(kmet)

           keypot=ltpmet(k0)
           If      (keypot == 3) Then

! sutton-chen potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              nnn=Nint(prmmet(3,k0)) ; nnnr=Real(nnn,wp)
              mmm=Nint(prmmet(4,k0)) ; mmmr=Real(mmm,wp)
              ccc=prmmet(5,k0)

              elrc0=eps*sig**3*(sig/rmet)**(nnn-3)/(nnnr-3.0_wp)
              vlrc0=nnnr*elrc0

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

              If (i /= j) Then
                 elrc0 = elrc0*2.0_wp
                 vlrc0 = vlrc0*2.0_wp
              End If

              elrcm(0) = elrcm(0) + twopi*volm*dens(i)*dens(j)*elrc0
              vlrcm(0) = vlrcm(0) - twopi*volm*dens(i)*dens(j)*vlrc0

              tmp=sig**3*(sig/rmet)**(mmm-3)/(mmmr-3.0_wp)
              If (i == j) Then
                 elrc1=tmp*(eps*ccc)**2
                 elrcm(i)=elrcm(i)+fourpi*dens(i)*elrc1
                 elrcsum=elrcsum+twopi*volm*dens(i)**2*elrc1

                 vlrc1=mmmr*elrc1
                 vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrc1
              Else
                 k1=lstmet((i*(i+1))/2)
                 k2=lstmet((j*(j+1))/2)

                 elrc1=tmp*(prmmet(1,k1)*prmmet(5,k1))**2
                 elrc2=tmp*(prmmet(1,k2)*prmmet(5,k2))**2
                 elrcm(i)=elrcm(i)+fourpi*dens(j)*elrc1
                 elrcm(j)=elrcm(j)+fourpi*dens(i)*elrc2
                 elrcsum=elrcsum+twopi*volm*dens(i)*dens(j)*(elrc1+elrc2)

                 vlrc1=mmmr*elrc1
                 vlrc2=mmmr*elrc2
                 vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrc1
                 vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrc2
              End If

           Else If (keypot == 4) Then

! gupta potentials

              aaa=prmmet(1,k0)
              rr0=prmmet(2,k0)
              ppp=prmmet(3,k0)
              zet=prmmet(4,k0)
              qqq=prmmet(5,k0)
              eee=Exp(-ppp*(rmet-rr0)/rr0)

              elrc0=2.0_wp*aaa*(rr0/ppp)*(rmet**2+2.0_wp*rmet*(rr0/ppp)+2.0_wp*(rr0/ppp)**2)*eee
              vlrc0=2.0_wp*aaa*rmet**3*eee+3.0_wp*elrc0

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

              If (i /= j) Then
                 elrc0=elrc0*2.0_wp
                 vlrc0=vlrc0*2.0_wp
              End If

              elrcm(0)=elrcm(0)+twopi*volm*dens(i)*dens(j)*elrc0
              vlrcm(0)=vlrcm(0)-twopi*volm*dens(i)*dens(j)*vlrc0

              eee=Exp(-2.0_wp*qqq*(rmet-rr0)/rr0)

              If (i == j) Then
                 elrc1=(rmet**2+2.0_wp*rmet*(0.5_wp*rr0/qqq)+2.0_wp*(0.5_wp*rr0/qqq)**2)*(0.5_wp*rr0/qqq)*eee*zet**2
                 elrcm(i)=elrcm(i)+fourpi*dens(i)*elrc1
                 elrcsum=elrcsum+twopi*volm*dens(i)**2*elrc1

                 vlrc1=(rmet**3+3.0_wp*rmet**2*(0.5_wp*rr0/qqq)+6.0_wp*rmet*(0.5_wp*rr0/qqq)**2+(0.5_wp*rr0/qqq)**3)*eee*zet**2
                 vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrc1
              Else
                 elrc1=(rmet**2+2.0_wp*rmet*(0.5_wp*rr0/qqq)+2.0_wp*(0.5_wp*rr0/qqq)**2)*(0.5_wp*rr0/qqq)*eee*zet**2
                 elrc2=elrc2
                 elrcm(i)=elrcm(i)+fourpi*dens(j)*elrc1
                 elrcm(j)=elrcm(j)+fourpi*dens(i)*elrc2
                 elrcsum=elrcsum+twopi*volm*dens(i)*dens(j)*(elrc1+elrc2)

                 vlrc1=(rmet**3+3.0_wp*rmet**2*(0.5_wp*rr0/qqq)+6.0_wp*rmet*(0.5_wp*rr0/qqq)**2+(0.5_wp*rr0/qqq)**3)*eee*zet**2
                 vlrc2=vlrc1
                 vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrc1
                 vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrc2
              End If

           Else If (keypot == 5) Then

! many-body perturbation component only potentials

              eps=prmmet(1,k0)
              sig=prmmet(2,k0)
              mmm=Nint(prmmet(3,k0)) ; mmmr=Real(mmm,wp)

! No pairwise contributions for mbpc potentials!!!

!              elrc0=0.0_wp
!              vlrc0=0.0_wp

! Self-interaction accounted once, interaction between different species
! MUST be accounted twice!!

!              If (i /= j) Then
!                 elrc0 = elrc0*2.0_wp
!                 vlrc0 = vlrc0*2.0_wp
!              End If

!              elrcm(0) = elrcm(0) + twopi*volm*dens(i)*dens(j)*elrc0
!              vlrcm(0) = vlrcm(0) - twopi*volm*dens(i)*dens(j)*vlrc0

              tmp=sig/((mmmr-3.0_wp)*rmet**(mmm-3))
              If (i == j) Then
                 elrc1=tmp*eps**2
                 elrcm(i)=elrcm(i)+fourpi*dens(i)*elrc1
                 elrcsum=elrcsum+twopi*volm*dens(i)**2*elrc1

                 vlrc1=mmmr*elrc1
                 vlrcm(i)=vlrcm(i)+twopi*dens(i)*vlrc1
              Else
                 k1=lstmet((i*(i+1))/2)
                 k2=lstmet((j*(j+1))/2)

                 elrc1=tmp*prmmet(1,k1)**2
                 elrc2=tmp*prmmet(1,k2)**2
                 elrcm(i)=elrcm(i)+fourpi*dens(j)*elrc1
                 elrcm(j)=elrcm(j)+fourpi*dens(i)*elrc2
                 elrcsum=elrcsum+twopi*volm*dens(i)*dens(j)*(elrc1+elrc2)

                 vlrc1=mmmr*elrc1
                 vlrc2=mmmr*elrc2
                 vlrcm(i)=vlrcm(i)+twopi*dens(j)*vlrc1
                 vlrcm(j)=vlrcm(j)+twopi*dens(i)*vlrc2
              End If

           End If

        End Do
     End Do
  End If

  If (newjob) Then
     newjob =.false.

     If (idnode == 0) Then
        Write(nrite,"(/,/,1x,1p,                                  &
        & 'long-range correction to metal energy    ',e15.6,/,1x, &
        & 'lr correction for metal atom density     ',e15.6,/,1x, &
        & '1st partial lr correction to metal virial',e15.6,/)")  &
           elrcm(0)/engunit,elrcsum/engunit**2,vlrcm(0)/engunit

        Write(nrite,"(1x,'density dependent energy and virial corrections',/)")
        Do i=1,ntpatm
           kmet=lstmet((i*(i+1))/2)
           If (lstmet(kmet) > 0) Write(nrite,"(25x,a8,1p,2e15.6)") &
              unqatm(i),elrcm(i)/engunit,vlrcm(i)/engunit
        End Do
     End If
  End If

End Subroutine metal_lrc
