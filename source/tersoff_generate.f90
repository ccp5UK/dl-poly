Subroutine tersoff_generate(rcter)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for generating potential energy and force arrays
! for tersoff forces only, based on potential form as defined by
! J. Tersoff, Phys. Rev. B 39 (1989) 5566
!
! copyright - daresbury laboratory
! author    - w.smith  october 2004
! amended   - i.t.todorov october 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,   Only : pi,mxgrid
  Use tersoff_module, Only : ntpter,lstter,ltpter,prmter,gmbp,vmbp

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rcter

  Integer           :: i,katom1,katom2,ipt,jpt,kpt
  Real( Kind = wp ) :: dlrpot,baij,saij,bbij,sbij,rij,sij,rrr,arg,rep,att

! define grid resolution for potential arrays

  dlrpot=rcter/Real(mxgrid-4,wp)

! construct arrays for all types of tersoff potential

  Do katom1=1,ntpter
     Do katom2=1,katom1

        If (ltpter(katom1) == 1 .and. ltpter(katom2) == 1) Then

           ipt=lstter(katom1)
           jpt=lstter(katom2)
           kpt=(Max(ipt,jpt)*(Max(ipt,jpt)-1))/2+Min(ipt,jpt)

! define tersoff parameters

           baij =    Sqrt(prmter(1,ipt)*prmter(1,jpt))
           saij = 0.5_wp*(prmter(2,ipt)+prmter(2,jpt))
           bbij =    Sqrt(prmter(3,ipt)*prmter(3,jpt))
           sbij = 0.5_wp*(prmter(4,ipt)+prmter(4,jpt))
           rij  =    Sqrt(prmter(5,ipt)*prmter(5,jpt))
           sij  =    Sqrt(prmter(6,ipt)*prmter(6,jpt))

! store potential cutoff

           vmbp(1,kpt,1)=sij

! calculate screening generic function

           Do i=2,mxgrid
              rrr=Real(i,wp)*dlrpot

              If      (rrr <= rij) Then
                 vmbp(i,kpt,1)=1.0_wp
                 gmbp(i,kpt,1)=0.0_wp
              Else
                 If (rrr <= sij) Then
                    arg=pi*(rrr-rij)/(sij-rij)

                    vmbp(i,kpt,1)=0.5_wp*(1.0_wp+Cos(arg))
                    gmbp(i,kpt,1)=0.5_wp*pi*rrr*Sin(arg)/(sij-rij)

!                Else the rest is anyway initialised to zero

                 End If
              End If
           End Do

! calculate screening repulsion & attraction functions

           Do i=2,mxgrid
              rrr=Real(i,wp)*dlrpot

! repulsion

              rep=baij*Exp(-saij*rrr)

              vmbp(i,kpt,2)=rep*  vmbp(i,kpt,1)
              gmbp(i,kpt,2)=rep*( gmbp(i,kpt,1) + saij*rrr*vmbp(i,kpt,1) )

! attraction

              att=bbij*Exp(-sbij*rrr)

              vmbp(i,kpt,3)=att*  vmbp(i,kpt,1)
              gmbp(i,kpt,3)=att*( gmbp(i,kpt,1) + sbij*rrr*vmbp(i,kpt,1) )
           End Do

        End If

     End Do
  End Do

End Subroutine tersoff_generate
