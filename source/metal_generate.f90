Subroutine metal_generate(rmet)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating potential energy and force arrays
! for metal potentials
!
! copyright - daresbury laboratory
! author    - w.smith june 2006
! amended   - i.t.todorov december 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module, Only : mxgmet
  Use site_module,  Only : ntpatm
  Use metal_module, Only : lstmet,ltpmet,prmmet,dmet,vmet,merf,mfer

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rmet

  Integer           :: i,imet,kmet,keypot,katom1,katom2, &
                       pmet,qmet
  Real( Kind = wp ) :: dlrpot,rrr,              &
                       eps,sig,nnn,mmm,         &
                       cc0,cc1,cc2,cc3,cc4,     &
                       aaa,bbb,ccc,ddd,ppp,qqq, &
                       bet,cut1,cut2,rr0

! define grid resolution for potential arrays

  dlrpot=rmet/Real(mxgmet-1,wp)

! construct arrays for metal potentials

  kmet=0
  Do katom1=1,ntpatm
     Do katom2=1,katom1
        kmet=kmet+1

! calculate potentials for defined interactions

        imet=lstmet(kmet)
        keypot=ltpmet(imet)
        If (keypot > 0) Then

! store array specification parameters

           vmet(1,imet,1)=Real(mxgmet,wp)
           vmet(2,imet,1)=0.0_wp          ! l_int(min) >= 1
           vmet(3,imet,1)=rmet            ! rmet=rcut
           vmet(4,imet,1)=dlrpot

           Do i=1,4
              vmet(i,imet,2)=vmet(i,imet,1)
              dmet(i,imet,1)=vmet(i,imet,1)
              dmet(i,imet,2)=0.0_wp
           End Do

           If      (keypot == 1) Then

! finnis-sinclair potentials

              cc0=prmmet(1,imet)
              cc1=prmmet(2,imet)
              cc2=prmmet(3,imet)
              ccc=prmmet(4,imet)
              ddd=prmmet(6,imet)
              bet=prmmet(7,imet)
              cut1=ccc+4.0_wp*dlrpot
              cut2=ddd+4.0_wp*dlrpot

              vmet(3,imet,1:2)=cut1
              dmet(3,imet,1)=cut2

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot

                 If (rrr <= cut1) Then
                    vmet(i,imet,1)=(cc0+cc1*rrr+cc2*rrr**2)*(rrr-ccc)**2
                    vmet(i,imet,2)=-rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2) * &
                                    (rrr-ccc)+(cc1+2.0_wp*cc2*rrr)*(rrr-ccc)**2)
                 End If

                 If (rrr <= cut2) Then
                    dmet(i,imet,1)=(rrr-ddd)**2+bet*(rrr-ddd)**3/ddd
                    dmet(i,imet,2)=-rrr*(2.0_wp*(rrr-ddd)+3.0_wp*bet*(rrr-ddd)**2/ddd)
                 End If
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=prmmet(5,imet)**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=prmmet(5,pmet)**2
                 dmet(2,imet,2)=prmmet(5,qmet)**2
              End If

           Else If (keypot == 2) Then

! extended finnis-sinclair potentials

              cc0=prmmet(1,imet)
              cc1=prmmet(2,imet)
              cc2=prmmet(3,imet)
              cc3=prmmet(4,imet)
              cc4=prmmet(5,imet)
              ccc=prmmet(6,imet)
              ddd=prmmet(8,imet)
              bbb=prmmet(9,imet)
              cut1=ccc+4.0_wp*dlrpot
              cut2=ddd+4.0_wp*dlrpot

              vmet(3,imet,1:2)=cut1
              dmet(3,imet,1)=cut2

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot

                 If (rrr <= cut1) Then
                    vmet(i,imet,1)=(cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc)**2
                    vmet(i,imet,2)=-rrr*(2.0_wp*(cc0+cc1*rrr+cc2*rrr**2+cc3*rrr**3+cc4*rrr**4)*(rrr-ccc) + &
                                         (cc1+2.0_wp*cc2*rrr+3.0_wp*cc3*rrr**2+4.0_wp*cc4*rrr**3)*(rrr-ccc)**2)
                 End If

                 If (rrr <= cut2) Then
                    dmet(i,imet,1)=(rrr-ddd)**2+bbb**2*(rrr-ddd)**4
                    dmet(i,imet,2)=-rrr*(2.0_wp*(rrr-ddd)+4.0_wp*bbb*2*(rrr-ddd)**3)
                 End If
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=prmmet(7,imet)**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=prmmet(7,pmet)**2
                 dmet(2,imet,2)=prmmet(7,qmet)**2
              End If

           Else If (keypot == 3) Then

! sutton-chen potentials

              eps=prmmet(1,imet)
              sig=prmmet(2,imet)
              nnn=prmmet(3,imet)
              mmm=prmmet(4,imet)

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot
                 vmet(i,imet,1)=eps*(sig/rrr)**nnn
                 vmet(i,imet,2)=nnn*eps*(sig/rrr)**nnn
                 dmet(i,imet,1)=(sig/rrr)**mmm
                 dmet(i,imet,2)=mmm*(sig/rrr)**mmm
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=(prmmet(1,imet)*prmmet(5,imet))**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=(prmmet(1,pmet)*prmmet(5,pmet))**2
                 dmet(2,imet,2)=(prmmet(1,qmet)*prmmet(5,qmet))**2
              End If

           Else If (keypot == 4) Then

! gupta potentials

              aaa=prmmet(1,imet)
              rr0=prmmet(2,imet)
              ppp=prmmet(3,imet)
              qqq=prmmet(5,imet)

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot

                 cut1=(rrr-rr0)/rr0
                 cut2=cut1+1.0_wp

                 vmet(i,imet,1)=2.0_wp*aaa*Exp(-ppp*cut1)
                 vmet(i,imet,2)=vmet(i,imet,1)*ppp*cut2
                 dmet(i,imet,1)=Exp(-2.0_wp*qqq*cut1)
                 dmet(i,imet,2)=2.0_wp*dmet(i,imet,1)*qqq*cut2
              End Do

              dmet(1,imet,2)=prmmet(4,imet)**2
              dmet(2,imet,2)=dmet(1,imet,2)

           Else If (keypot == 5) Then

! many-body perturbation component only potentials

              eps=prmmet(1,imet)
              sig=prmmet(2,imet)
              mmm=prmmet(3,imet)

              Do i=5,mxgmet
                 rrr=Real(i,wp)*dlrpot
!                 vmet(i,imet,1)=0.0_wp
!                 vmet(i,imet,2)=0.0_wp
                 nnn=sig/rrr**mmm
                 dmet(i,imet,1)=nnn*merf(i)
                 dmet(i,imet,2)=mmm*dmet(i,imet,1)-rrr*nnn*mfer(i)
              End Do

              If (katom1 == katom2) Then
                 dmet(1,imet,2)=prmmet(1,imet)**2
                 dmet(2,imet,2)=dmet(1,imet,2)
              Else
                 pmet=lstmet((katom1*(katom1+1))/2)
                 qmet=lstmet((katom2*(katom2+1))/2)
                 dmet(1,imet,2)=prmmet(1,pmet)**2
                 dmet(2,imet,2)=prmmet(1,qmet)**2
              End If

           Else

              Call error(151)

           End If

        End If
     End Do
  End Do

End Subroutine metal_generate
