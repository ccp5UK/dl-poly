Subroutine external_field_apply(keyshl,time,engfld,virfld)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for application of an external field
!
! Note: Only one field at a time is allowed
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module,   Only : idnode,mxnode,gcheck,gsum
  Use setup_module,   Only : twopi,nrite,mxshl,mxatms
  Use config_module,  Only : imcon,cell,natms,nfree,nlast,lsi,lsa,ltg, &
                             lfrzn,lstfre,weight,chge,           &
                             xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
  Use kinetic_module, Only : getcom_mol
  Use rigid_bodies_module
  Use core_shell_module
  Use external_field_module

  Implicit None

  Integer,           Intent( In    ) :: keyshl
  Real( Kind = wp ), Intent( In    ) :: time ! for oscillating fields
  Real( Kind = wp ), Intent(   Out ) :: engfld,virfld

  Logical, Save     :: newjob = .true.

  Logical           :: safe
  Integer           :: i,j,ia,ib,fail(1:2),local_index, &
                       irgd,jrgd,lrgd,rgdtyp,megrgd
  Real( Kind = wp ) :: gamma,rrr,rz,zdif,vxt,vyt,vzt,tmp,rtmp(1:2), &
                       x(1:1),y(1:1),z(1:1),cmm(0:3)

  Integer,           Allocatable :: lstopt(:,:),list(:)
  Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)

! Recover megrgd

  megrgd=rgdmeg

! energy and virial accumulators

  engfld=0.0_wp
  virfld=0.0_wp

  If (keyfld == 1) Then

! electric field: prmfld(1-3) are field components

     Do i=1,natms
        If (lfrzn(i) == 0) Then
           fxx(i)=fxx(i) + chge(i)*prmfld(1)
           fyy(i)=fyy(i) + chge(i)*prmfld(2)
           fzz(i)=fzz(i) + chge(i)*prmfld(3)
        End If
     End Do

  Else If (keyfld == 2) Then

! oscillating shear: orthorhombic box:  Fx=a*Cos(b.2.pi.z/L)

     If (imcon /= 1 .and. imcon /= 2) Return

     rz=twopi/cell(9)

     Do i=1,natms
        If (lfrzn(i) == 0) fxx(i)=fxx(i) + prmfld(1)*Cos(prmfld(2)*zzz(i)*rz)
     End Do

  Else If (keyfld == 3) Then

! continuous shear of walls : 2D periodic box (imcon=6)

     If (imcon /= 6) Return

! shear rate=prmfld(1) angstrom per ps for non-frozen
! and non-weightless atoms at Abs(z) > prmfld(2)

     If (megrgd > 0) Then

! FPs
        Do j=1,nfree
           i=lstfre(j)

           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. Abs(zzz(i)) > prmfld(2)) &
           vxx(i)=0.5_wp*Sign(prmfld(1),zzz(i))
        End Do

! RBs

        Do irgd=1,ntrgd
           rgdtyp=listrgd(0,irgd)

! For all good RBs

           lrgd=listrgd(-1,irgd)
           If (rgdfrz(0,rgdtyp) == 0) Then
              x=rgdxxx(irgd) ; y=rgdyyy(irgd) ; z=rgdzzz(irgd)
              Call images(imcon,cell,1,x,y,z)

              rz=z(1)
              If (Abs(rz) > prmfld(2)) Then
                 tmp=0.5_wp*Sign(prmfld(1),rz)
                 vxt=tmp-rgdvxx(irgd)

                 rgdvxx(irgd)=tmp
                 Do jrgd=1,lrgd
                    i=indrgd(jrgd,irgd) ! local index of particle/site

                    If (i <= natms) vxx(i)=vxx(i)+vxt
                 End Do
              End If
           End If
        End Do

     Else

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp .and. Abs(zzz(i)) > prmfld(2)) &
           vxx(i)=0.5_wp*Sign(prmfld(1),zzz(i))
        End Do

     End If

  Else If (keyfld == 4) Then

! gravitational field: field components given by prmfld(1-3)

     If (keyshl == 1) Then

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + prmfld(1)*weight(i)
              fyy(i)=fyy(i) + prmfld(2)*weight(i)
              fzz(i)=fzz(i) + prmfld(3)*weight(i)
           End If
        End Do

     Else

        fail=0
        Allocate (lstopt(1:2,1:mxshl),                       Stat=fail(1))
        Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms), Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply allocation failure, node: ', idnode
           Call error(0)
        End If

        Do i=1,ntshl
           lstopt(1,i)=local_index(listshl(1,i),nlast,lsi,lsa)
           lstopt(2,i)=local_index(listshl(2,i),nlast,lsi,lsa)
        End Do

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              oxt(i)=prmfld(1)*weight(i)
              oyt(i)=prmfld(2)*weight(i)
              ozt(i)=prmfld(3)*weight(i)
           Else ! for the sake of massless sites of RBs
              oxt(i)=0.0_wp
              oyt(i)=0.0_wp
              ozt(i)=0.0_wp
           End If
        End Do

        If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,oxt,oyt,ozt)

! Transfer cores' forces to shells

        Do i=1,ntshl
           ia=lstopt(1,i)
           ib=lstopt(2,i)
           If (ia > 0 .and. (ib > 0 .and. ib <= natms)) Then
              oxt(ib)=oxt(ia)
              oyt(ib)=oyt(ia)
              ozt(ib)=ozt(ia)
           End If
        End Do

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + oxt(i)
              fyy(i)=fyy(i) + oyt(i)
              fzz(i)=fzz(i) + ozt(i)
           End If
        End Do

        Deallocate (lstopt,      Stat=fail(1))
        Deallocate (oxt,oyt,ozt, Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply deallocation failure, node: ', idnode
           Call error(0)
        End If

     End If

  Else If (keyfld == 5) Then

! magnetic field: field components given by prmfld(1-3)

     If (keyshl == 1) Then

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + (vyy(i)*prmfld(3)-vzz(i)*prmfld(2))*chge(i)
              fyy(i)=fyy(i) + (vzz(i)*prmfld(1)-vxx(i)*prmfld(3))*chge(i)
              fzz(i)=fzz(i) + (vxx(i)*prmfld(2)-vyy(i)*prmfld(1))*chge(i)
           End If
        End Do

     Else

        fail=0
        Allocate (lstopt(1:2,1:mxshl),                       Stat=fail(1))
        Allocate (oxt(1:mxatms),oyt(1:mxatms),ozt(1:mxatms), Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply allocation failure, node: ', idnode
           Call error(0)
        End If

        Do i=1,ntshl
           lstopt(1,i)=local_index(listshl(1,i),nlast,lsi,lsa)
           lstopt(2,i)=local_index(listshl(2,i),nlast,lsi,lsa)
        End Do

! cores' velocities

        Do i=1,natms
           If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
              oxt(i)=vxx(i)
              oyt(i)=vyy(i)
              ozt(i)=vzz(i)
           Else
              oxt(i)=0.0_wp
              oyt(i)=0.0_wp
              ozt(i)=0.0_wp
           End If
        End Do

        If (lshmv_shl) Call update_shared_units(natms,nlast,lsi,lsa,lishp_shl,lashp_shl,oxt,oyt,ozt)

! Transfer cores' velocities to shells

        Do i=1,ntshl
           ia=lstopt(1,i)
           ib=lstopt(2,i)
           If (ia > 0 .and. (ib > 0 .and. ib <= natms)) Then
              oxt(ib)=oxt(ia)
              oyt(ib)=oyt(ia)
              ozt(ib)=ozt(ia)
           End If
        End Do

        Do i=1,natms
           If (lfrzn(i) == 0) Then
              fxx(i)=fxx(i) + (oyt(i)*prmfld(3)-ozt(i)*prmfld(2))*chge(i)
              fyy(i)=fyy(i) + (ozt(i)*prmfld(1)-oxt(i)*prmfld(3))*chge(i)
              fzz(i)=fzz(i) + (oxt(i)*prmfld(2)-oyt(i)*prmfld(1))*chge(i)
           End If
        End Do

        Deallocate (lstopt,      Stat=fail(1))
        Deallocate (oxt,oyt,ozt, Stat=fail(2))
        If (Any(fail > 0)) Then
           Write(nrite,'(/,1x,a,i0)') 'external_field_apply deallocation failure, node: ', idnode
           Call error(0)
        End If

     End If

  Else If (keyfld == 6) Then

! containing sphere : r^(-n) potential

     Do i=1,natms
        If (lfrzn(i) == 0) Then
           rrr=Sqrt(xxx(i)**2+yyy(i)**2+zzz(i)**2)
           If (rrr > prmfld(4)) Then
              rrr=prmfld(2)-rrr
              If (rrr < 0.0_wp) rrr=0.1_wp

              gamma=prmfld(1)*rrr**(-prmfld(3))
              engfld=engfld + gamma
              gamma=-prmfld(3)*gamma/(rrr*rrr)

              fxx(i)=fxx(i) + gamma*xxx(i)
              fyy(i)=fyy(i) + gamma*yyy(i)
              fzz(i)=fzz(i) + gamma*zzz(i)
           End If
        End If
     End Do

     virfld=-9.0_wp*engfld

  Else If (keyfld == 7) Then

! repulsive wall (harmonic) starting at z0

     Do i=1,natms
        If (lfrzn(i) == 0 .and. prmfld(3)*zzz(i) > prmfld(3)*prmfld(2)) Then
           zdif=zzz(i)-prmfld(2)
           gamma=-prmfld(1)*zdif

           fzz(i)=fzz(i) + gamma
           engfld=engfld - gamma*zdif/2.0_wp
        End If
     End Do

  Else If (keyfld == 8) Then

! xpist - piston wall pushing down along the X=bxc direction
! prmfld(1) is the first atom of the layer of molecules (membrane) to be pushed
! prmfld(2) is the last atom of the layer of molecules (membrane) to be pushed
! prmfld(3) is the pressure applied to the layer of molecules (membrane) in the
! +X=bxc direction - i.e. left to right.  The layer plane is defined as _|_ bxc

     If (imcon /= 1 .and. imcon /= 2) Return

     ia = Nint(prmfld(1))
     ib = Nint(prmfld(2))

     If (newjob) Then
        newjob=.false.

!        mass=0.0_wp ! defined and initialise in external_field_module
        safe=.true.
        Do i=1,natms
           If ((ltg(i) >= ia .and. ltg(i) <= ib) .and. lfrzn(i) > 0) Then
              safe=.false.
           Else
              mass=mass+weight(i)
           End If
        End Do

        If (mxnode > 1) Then
           Call gsum(mass)
           Call gcheck(safe)
        End If
        If (.not.safe) Call error(456)
     End If

     rtmp=0.0_wp ! average velocity and force per atom in x direction of the piston
     Do i=1,natms
        If (ltg(i) >= ia .and. ltg(i) <= ib) Then
           rtmp(1)=rtmp(1)+vxx(i) ; rtmp(2)=rtmp(2)+fxx(i)
           vyy(i) = 0.0_wp        ; fyy(i) = 0.0_wp
           vzz(i) = 0.0_wp        ; fzz(i) = 0.0_wp
        End If
     End Do
     If (mxnode > 1) Call gsum(rtmp) ! net velocity and force to ensure solid wall behaviour

     rtmp(1)=rtmp(1)/Real(ib-ia+1)     ! averaged velocity per particle
     rtmp(2)=(rtmp(2)+prmfld(3))/mass  ! averaged acceleration of the slab

     Do i=1,natms
        If (ltg(i) >= ia .and. ltg(i) <= ib) Then
           vxx(i)=rtmp(1)
           fxx(i)=rtmp(2)*weight(i) ! force per particle
        End If
     End Do

   Else If (keyfld == 9) Then

! zres external field: restrain molecule z-position (pull in)
! prmfld(1) is the index of first atom of restrained molecule
! prmfld(2) is the index of last atom of restrained molecule
! prmfld(3) is the restraining constant
! prmfld(4) is z-min (min limit in z-direction)
! prmfld(5) is z-max (max limit in z-direction)
! where prmfld(4) < prmfld(5)

      ia = Nint(prmfld(1))
      ib = Nint(prmfld(2))

! Get molecule's weight and CoM

      Call getcom_mol(ia,ib,cmm)

! Apply force corrections

      Do i=1,natms
         If (ltg(i) >= ia .and. ltg(i) <= ib .and. lfrzn(i) == 0) Then
            If (cmm(3) < prmfld(4) .or. cmm(3) > prmfld(5)) Then
               If (cmm(3) < prmfld(4)) zdif = cmm(3) - prmfld(4)
               If (cmm(3) > prmfld(5)) zdif = prmfld(5) - cmm(3)

               gamma=-prmfld(3)*zdif*weight(i)/cmm(0)
               fzz(i)=fzz(i) + gamma
               engfld=engfld - 0.5_wp*gamma*zdif
            End If
         End If
      End Do

   Else If (keyfld == 10) Then

! extension to exzn- external field (pull out)
! prmfld(1) is the index of first atom of the water molecules to be restrained
! prmfld(2) is the index of last atom of the water molecules to be restrained
! prmfld(3) is the restraining constant
! prmfld(4) is z-min (min limit in z-direction)
! prmfld(5) is z-max (max limit in z-direction)
! where prmfld(4) < prmfld(5)
! This will keep water away from the membrane region, and allow
! the DMPC to slowly fill in the gaps before water molecules are let loose.

      ia = Nint(prmfld(1))
      ib = Nint(prmfld(2))
      Do i=1,natms
         If (ltg(i) >= ia .and. ltg(i) <= ib .and. lfrzn(i) == 0) Then
            If (zzz(i) > prmfld(4) .and. zzz(i) < prmfld(5)) Then
               tmp = prmfld(5) + prmfld(4)

               If (zzz(i) <  tmp) zdif = zzz(i) - prmfld(4)
               If (zzz(i) >= tmp) zdif = prmfld(5) - zzz(i)

               gamma=-prmfld(3)*zdif
               fzz(i)=fzz(i) + gamma
               engfld=engfld - 0.5_wp*gamma*zdif
            End If
         End If
      End Do

   Else If (keyfld == 11) Then

! extension to exzn+ external field (pull in)
! prmfld(1) is the index of first atom of the water molecules to be restrained
! prmfld(2) is the index of last atom of the water molecules to be restrained
! prmfld(3) is the restraining constant
! prmfld(4) is z-min (min limit in z-direction)
! prmfld(5) is z-max (max limit in z-direction)
! where prmfld(4) < prmfld(5)
! This will keep water within this region.

      ia = Nint(prmfld(1))
      ib = Nint(prmfld(2))
      Do i=1,natms
         If (ltg(i) >= ia .and. ltg(i) <= ib .and. lfrzn(i) == 0) Then
            If (zzz(i) < prmfld(4) .and. zzz(i) > prmfld(5)) Then
               If (zzz(i) < prmfld(4)) zdif = zzz(i) - prmfld(4)
               If (zzz(i) > prmfld(5)) zdif = prmfld(5) - zzz(i)

               gamma=-prmfld(3)*zdif
               fzz(i)=fzz(i) + gamma
               engfld=engfld - 0.5_wp*gamma*zdif
            End If
         End If
      End Do

  Else If (keyfld == 12) Then

! extension to oscillating electric field: prmfld(1-3) are field components
! prmfld(4) is the oscillating frequency defined in ps^-1!

     tmp=Sin(time*prmfld(4)*twopi)
     Do i=1,natms
        If (lfrzn(i) == 0) Then
           fxx(i)=fxx(i) + chge(i)*prmfld(1)*tmp
           fyy(i)=fyy(i) + chge(i)*prmfld(2)*tmp
           fzz(i)=fzz(i) + chge(i)*prmfld(3)*tmp
        End If
     End Do

  Else

! unidentified field potential error exit

     Call error(454)

  End If

End Subroutine external_field_apply
