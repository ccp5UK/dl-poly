Module kinetic_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 module for kinetic routines and functions
!
! Function getkin - calculates the kinetic energy
! Function getknf - calculates the kinetic energy of free particles
! Function getknt - calculates the translational kinetic energy of RBs
! Function getknr - calculates the rotational kinetic energy of RBs
! Subroutine kinstress - calculates the kinetic stress
! Subroutine kinstresf - calculates the kinetic stress of free particles
! Subroutine kinstrest - calculates the kinetic stress of RBs (t, no r)
! Subroutine getcom - calculates the centre of mass position
! Subroutine chvom - changes the behaviour of getvom
! Subroutine getvom - calculates the centre of mass momentum
! Subroutine freeze_atoms - quenches forces and velocities on 'frozen'
!                           atoms
! Subroutine cap_forces - limits the absolute magnitude of forces
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use comms_module, Only : mxnode,gsum

  Implicit None

! Remove COM motion defaults

  Logical,          Save :: l_vom = .true.
  Logical, Private, Save :: lvom  = .true.

  Public :: getkin,getknf,getknt,getknr,    &
            kinstress,kinstresf,kinstrest,  &
            getcom,getcom_mol,getvom,chvom, &
            freeze_atoms,cap_forces

  Interface getvom
     Module Procedure getvom
     Module Procedure getvom_rgd
  End Interface !getvom

Contains

  Function getkin(vxx,vyy,vzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the system kinetic energy
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use config_module, Only : natms,lfrzn,weight

    Implicit None

    Real( Kind = wp )                                    :: getkin

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz

    Integer           :: i
    Real( Kind = wp ) :: engke

    engke = 0.0_wp

    Do i=1,natms
       If (lfrzn(i) == 0) &
          engke = engke + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
    End Do

    If (mxnode > 1) Call gsum(engke)

    getkin = 0.5_wp * engke

  End Function getkin

  Function getknf(vxx,vyy,vzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the kinetic energy of free atoms
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use config_module, Only : nfree,lfrzn,lstfre,weight

    Implicit None

    Real( Kind = wp )                                    :: getknf

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz

    Integer           :: i,j
    Real( Kind = wp ) :: engke

    engke = 0.0_wp

    Do j=1,nfree
       i=lstfre(j)

       If (lfrzn(i) == 0) &
          engke = engke + weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
    End Do

    If (mxnode > 1) Call gsum(engke)

    getknf = 0.5_wp * engke

  End Function getknf

  Function getknt(rgdvxx,rgdvyy,rgdvzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the translational kinetic energy of
! rigid bodies
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use rigid_bodies_module, Only : ntrgd,rgdfrz,rgdwgt,listrgd,indrgd

    Implicit None

    Real( Kind = wp )                                    :: getknt

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: rgdvxx,rgdvyy,rgdvzz

    Integer           :: irgd,lrgd,rgdtyp
    Real( Kind = wp ) :: engtra,tmp

    engtra = 0.0_wp

    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       If (rgdfrz(0,rgdtyp) == 0) Then
          lrgd=listrgd(-1,irgd)

          tmp=rgdwgt(0,rgdtyp)*Real(indrgd(0,irgd),wp)/Real(lrgd,wp)

          engtra = engtra + tmp*(rgdvxx(irgd)**2+rgdvyy(irgd)**2+rgdvzz(irgd)**2)
       End If
    End Do

    If (mxnode > 1) Call gsum(engtra)

    getknt = 0.5_wp * engtra

  End Function getknt

  Function getknr(rgdoxx,rgdoyy,rgdozz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the rotational kinetic energy of
! rigid bodies
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use rigid_bodies_module, Only : ntrgd,rgdfrz,rgdrix,rgdriy,rgdriz,listrgd,indrgd

    Implicit None

    Real( Kind = wp )                                    :: getknr

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: rgdoxx,rgdoyy,rgdozz

    Integer           :: irgd,lrgd,rgdtyp
    Real( Kind = wp ) :: engrot,tmp

    engrot = 0.0_wp

    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) < lrgd) Then
          tmp=Real(indrgd(0,irgd),wp)/Real(lrgd,wp)

          engrot = engrot + tmp *                     &
                   (rgdrix(1,rgdtyp)*rgdoxx(irgd)**2+ &
                    rgdriy(1,rgdtyp)*rgdoyy(irgd)**2+ &
                    rgdriz(1,rgdtyp)*rgdozz(irgd)**2)
       End If
    End Do

    If (mxnode > 1) Call gsum(engrot)

    getknr = 0.5_wp * engrot

  End Function getknr

  Subroutine kinstress(vxx,vyy,vzz,strkin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate kinetic contribution to the stress
! tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use config_module, Only : natms,lfrzn,weight

    Implicit None

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz
    Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: strkin

    Integer :: i

    strkin = 0.0_wp

    Do i=1,natms
       If (lfrzn(i) == 0) Then
          strkin(1) = strkin(1) + weight(i)*vxx(i)*vxx(i)
          strkin(2) = strkin(2) + weight(i)*vxx(i)*vyy(i)
          strkin(3) = strkin(3) + weight(i)*vxx(i)*vzz(i)
          strkin(5) = strkin(5) + weight(i)*vyy(i)*vyy(i)
          strkin(6) = strkin(6) + weight(i)*vyy(i)*vzz(i)
          strkin(9) = strkin(9) + weight(i)*vzz(i)*vzz(i)
       End If
    End Do

    If (mxnode > 1) Call gsum(strkin)

! Symmetrise

    strkin(4) = strkin(2)
    strkin(7) = strkin(3)
    strkin(8) = strkin(6)

  End Subroutine kinstress

  Subroutine kinstresf(vxx,vyy,vzz,strknf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate kinetic contribution of free atoms to
! the stress tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use config_module, Only : nfree,lfrzn,lstfre,weight

    Implicit None

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz
    Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: strknf

    Integer :: i,j

    strknf = 0.0_wp

    Do j=1,nfree
       i=lstfre(j)

       If (lfrzn(i) == 0) Then
          strknf(1) = strknf(1) + weight(i)*vxx(i)*vxx(i)
          strknf(2) = strknf(2) + weight(i)*vxx(i)*vyy(i)
          strknf(3) = strknf(3) + weight(i)*vxx(i)*vzz(i)
          strknf(5) = strknf(5) + weight(i)*vyy(i)*vyy(i)
          strknf(6) = strknf(6) + weight(i)*vyy(i)*vzz(i)
          strknf(9) = strknf(9) + weight(i)*vzz(i)*vzz(i)
       End If
    End Do

    If (mxnode > 1) Call gsum(strknf)

! Symmetrise

    strknf(4) = strknf(2)
    strknf(7) = strknf(3)
    strknf(8) = strknf(6)

  End Subroutine kinstresf

  Subroutine kinstrest(rgdvxx,rgdvyy,rgdvzz,strknt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate kinetic contribution of RB
! translational motion to the stress tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use rigid_bodies_module, Only : ntrgd,rgdfrz,rgdwgt,listrgd,indrgd

    Implicit None

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: rgdvxx,rgdvyy,rgdvzz
    Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: strknt

    Integer           :: irgd,lrgd,rgdtyp
    Real( Kind = wp ) :: tmp

    strknt = 0.0_wp

    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

       If (rgdfrz(0,rgdtyp) == 0) Then
          lrgd=listrgd(-1,irgd)

          tmp=rgdwgt(0,rgdtyp)*Real(indrgd(0,irgd),wp)/Real(lrgd,wp)

          strknt(1) = strknt(1) + tmp*rgdvxx(irgd)*rgdvxx(irgd)
          strknt(2) = strknt(2) + tmp*rgdvxx(irgd)*rgdvyy(irgd)
          strknt(3) = strknt(3) + tmp*rgdvxx(irgd)*rgdvzz(irgd)
          strknt(5) = strknt(5) + tmp*rgdvyy(irgd)*rgdvyy(irgd)
          strknt(6) = strknt(6) + tmp*rgdvyy(irgd)*rgdvzz(irgd)
          strknt(9) = strknt(9) + tmp*rgdvzz(irgd)*rgdvzz(irgd)
       End If
    End Do

    If (mxnode > 1) Call gsum(strknt)

! Symmetrise

    strknt(4) = strknt(2)
    strknt(7) = strknt(3)
    strknt(8) = strknt(6)

  End Subroutine kinstrest

  Subroutine getcom(xxx,yyy,zzz,com)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate system centre of mass position
!
! copyright - daresbury laboratory
! author    - i.t.todorov october 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use setup_module,  Only : zero_plus
    Use config_module, Only : natms,lfrzn,weight

    Implicit None

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: xxx,yyy,zzz
    Real( Kind = wp ), Dimension( 1:3 ), Intent(   Out ) :: com

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: totmas

    Integer                 :: i

! total system mass

    If (newjob) Then
       newjob = .false.

       totmas = 0.0_wp
       Do i=1,natms
          If (lfrzn(i) == 0) totmas = totmas + weight(i)
       End Do

       If (mxnode > 1) Call gsum(totmas)
    End If

    com = 0.0_wp

    Do i=1,natms
       If (lfrzn(i) == 0) Then
          com(1) = com(1) + weight(i)*xxx(i)
          com(2) = com(2) + weight(i)*yyy(i)
          com(3) = com(3) + weight(i)*zzz(i)
       End If
    End Do

    If (mxnode > 1) Call gsum(com)
    If (totmas >= zero_plus) com = com/totmas

  End Subroutine getcom

  Subroutine getcom_mol(imcon,istart,ifinish,cmm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate a molecule's mass and CoM
!
! istart  - the global index of the first atom of the molecule
! ifinish - the global index of the last atom of the molecule
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use comms_module,  Only : idnode
    Use setup_module,  Only : nrite,zero_plus
    Use config_module, Only : cell,natms,ltg,lfrzn,xxx,yyy,zzz,weight

    Implicit None
    Integer,           Intent( In    ) :: imcon,istart,ifinish

    Real( Kind = wp ), Intent(   Out ) :: cmm(0:3)

    Integer           :: fail,i,j,k

    Real( Kind = wp ), Allocatable :: mol(:,:)

    fail = 0
    Allocate (mol(1:(ifinish-istart+1),0:3), Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'getcom_mol allocation failure, node: ', idnode
       Call error(0)
    End If

! Initialise

    cmm  = 0.0_wp

    mol = 0.0_wp
    Do i=1,natms
       j=ltg(i)
       If (j >= istart .and. j <= ifinish) Then
          k=j-istart+1

          mol(k,0) = weight(i)*Real(lfrzn(i)-1,wp)
          mol(k,1) = xxx(i)
          mol(k,2) = yyy(i)
          mol(k,3) = zzz(i)
       End If
    End Do

    If (mxnode > 1) Call gsum(mol)

    mol(:,1) = xxx(:)-xxx(1)
    mol(:,2) = yyy(:)-yyy(1)
    mol(:,3) = zzz(:)-zzz(1)

    k=ifinish-istart+1
    Call images(imcon,cell,k,mol(:,1),mol(:,2),mol(:,3))
    Do i=1,ifinish-istart+1
       cmm(0) = cmm(0) + mol(i,0)
       cmm(1) = cmm(1) + mol(i,0)*mol(i,1)
       cmm(2) = cmm(2) + mol(i,0)*mol(i,2)
       cmm(3) = cmm(3) + mol(i,0)*mol(i,3)
    End Do

    If (cmm(0) >= zero_plus) cmm(1:3) = cmm(1:3) / cmm(0)

    fail = 0
    Deallocate (mol, Stat = fail)
    If (fail > 0) Then
       Write(nrite,'(/,1x,a,i0)') 'getcom_mol deallocation failure, node: ', idnode
       Call error(0)
    End If

  End Subroutine getcom_mol

  Subroutine chvom(flag)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to change behaviour for COM momentum removal
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Implicit None

    Logical :: flag

    If (flag) Then
       lvom=.true.  ! Remove COM momentum
    Else
       lvom=.false. ! Don't
    End If

  End Subroutine chvom

  Subroutine getvom(vom,vxx,vyy,vzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate system centre of mass momentum
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use setup_module,  Only : mxatms,zero_plus
    Use config_module, Only : natms,lfrzn,weight

    Implicit None

    Real( Kind = wp ), Intent(   Out ) :: vom(1:3)
    Real( Kind = wp ), Intent( In    ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: totmas

    Integer                 :: i

! total system mass not including frozen mass

    If (newjob) Then
       newjob = .false.

! For all unfrozen, free particles

       totmas = 0.0_wp
       Do i=1,natms
          If (lfrzn(i) == 0) totmas = totmas + weight(i)
       End Do

       If (mxnode > 1) Call gsum(totmas)
    End If

    vom = 0.0_wp

    If (.not.lvom) Return

! For all unfrozen, free particles

    Do i=1,natms
       If (lfrzn(i) == 0) Then
          vom(1) = vom(1) + weight(i)*vxx(i)
          vom(2) = vom(2) + weight(i)*vyy(i)
          vom(3) = vom(3) + weight(i)*vzz(i)
       End If
    End Do

    If (mxnode > 1) Call gsum(vom)
    If (totmas >= zero_plus) vom = vom/totmas

  End Subroutine getvom

  Subroutine getvom_rgd(vom,vxx,vyy,vzz,rgdvxx,rgdvyy,rgdvzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate system centre of mass momentum
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use setup_module,        Only : mxatms,mxrgd,zero_plus
    Use config_module,       Only : nfree,lfrzn,lstfre,weight
    Use rigid_bodies_module, Only : ntrgd,rgdfrz,rgdwgt,listrgd,indrgd

    Implicit None

    Real( Kind = wp ), Intent(   Out ) :: vom(1:3)
    Real( Kind = wp ), Intent( In    ) :: vxx(1:mxatms),vyy(1:mxatms),vzz(1:mxatms)
    Real( Kind = wp ), Intent( In    ) :: rgdvxx(1:mxrgd),rgdvyy(1:mxrgd),rgdvzz(1:mxrgd)

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: totmas

    Integer                 :: i,j,irgd,lrgd,rgdtyp
    Real( Kind = wp )       :: tmp

! total system mass not including frozen mass

    If (newjob) Then
       newjob = .false.

       totmas = 0.0_wp

! For all unfrozen, free particles

       Do j=1,nfree
          i=lstfre(j)

          If (lfrzn(i) == 0) totmas = totmas + weight(i)
       End Do

! For all RBs without any frozen particles
! These with some frozen only turn and
! have zero net momentum enforced!

       Do irgd=1,ntrgd
          rgdtyp=listrgd(0,irgd)

          lrgd=listrgd(-1,irgd)
          If (rgdfrz(0,rgdtyp) == 0) &
             totmas = totmas + rgdwgt(0,rgdtyp)*Real(indrgd(0,irgd),wp)/Real(lrgd,wp)
       End Do

       If (mxnode > 1) Call gsum(totmas)
    End If

    vom = 0.0_wp

    If (.not.lvom) Return

! For all unfrozen, free particles

    Do j=1,nfree
       i=lstfre(j)

       If (lfrzn(i) == 0) Then
          tmp=weight(i)
          vom(1) = vom(1) + tmp*vxx(i)
          vom(2) = vom(2) + tmp*vyy(i)
          vom(3) = vom(3) + tmp*vzz(i)
       End If
    End Do

    Do irgd=1,ntrgd
       rgdtyp=listrgd(0,irgd)

! For all RBs without any frozen particles

       lrgd=listrgd(-1,irgd)
       If (rgdfrz(0,rgdtyp) == 0) Then
          tmp=rgdwgt(0,rgdtyp)*Real(indrgd(0,irgd),wp)/Real(lrgd,wp)

          vom(1) = vom(1) + tmp*rgdvxx(irgd)
          vom(2) = vom(2) + tmp*rgdvyy(irgd)
          vom(3) = vom(3) + tmp*rgdvzz(irgd)
       End If
    End Do

    If (mxnode > 1) Call gsum(vom)
    If (totmas >= zero_plus) vom = vom/totmas

  End Subroutine getvom_rgd

  Subroutine freeze_atoms()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to quench forces and velocities on 'frozen' atoms
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use config_module, Only : natms,lfrzn,vxx,vyy,vzz,fxx,fyy,fzz

    Implicit None

    Integer :: i

    Do i=1,natms
       If (lfrzn(i) /= 0) Then
           vxx(i) = 0.0_wp ; vyy(i) = 0.0_wp ; vzz(i) = 0.0_wp
           fxx(i) = 0.0_wp ; fyy(i) = 0.0_wp ; fzz(i) = 0.0_wp
       End If
    End Do

  End Subroutine freeze_atoms

  Subroutine cap_forces(fmax,temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to limit the absolute magnitude of forces.
! To be used in equilibration period only!!!
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Use setup_module,  Only : boltz
    Use config_module, Only : natms,lfrzn,weight,fxx,fyy,fzz

    Implicit None

    Real( Kind = wp ), Intent( In    ) :: fmax,temp

    Logical,           Save :: newjob = .true.
    Real( Kind = wp ), Save :: meg

    Integer           :: i
    Real( Kind = wp ) :: fmax2,fmod,scale,fcom(1:3)

    If (newjob) Then
       newjob = .false.

! Get the total number of non-frozen&non-massless particles

       meg=0.0_wp
       Do i=1,natms
          If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) &
             meg=meg+1.0_wp
       End Do
       If (mxnode > 1) Call gsum(meg)
    End If

! maximum force permitted

    fmax2 = (boltz*temp*fmax)**2

! cap forces and conserve linear momentum
! for non-frozen&non-massless particles

    fcom = 0.0_wp
    Do i=1,natms
       If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
          fmod = fxx(i)**2 + fyy(i)**2 + fzz(i)**2

          If (fmod > fmax2) Then
             scale  = Sqrt(fmax2/fmod)
             fxx(i) = fxx(i)*scale
             fyy(i) = fyy(i)*scale
             fzz(i) = fzz(i)*scale
          End If

! accummulate forces - to check on momentum conservation

          fcom(1) = fcom(1) + fxx(i)
          fcom(2) = fcom(2) + fyy(i)
          fcom(3) = fcom(3) + fzz(i)
       End If
    End Do

! ensure net forces sum to zero

    If (mxnode > 1) Call gsum(fcom)
    fcom = fcom/meg

! conserve momentum

    Do i=1,natms
       If (lfrzn(i) == 0 .and. weight(i) > 1.0e-6_wp) Then
          fxx(i) = fxx(i) - fcom(1)
          fyy(i) = fyy(i) - fcom(2)
          fzz(i) = fzz(i) - fcom(3)
       End If
    End Do

  End Subroutine cap_forces

End Module kinetic_module
