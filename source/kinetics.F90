Module kinetics

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

  Use kinds, only : wp
  Use comms, Only : comms_type, gsum
  Use setup,  Only : zero_plus,boltz
  Use configuration, Only : configuration_type,getcom_mol
  Use particle,     Only : corePart
  Use rigid_bodies, Only : rigid_bodies_type
  Implicit None

! Remove COM motion defaults


  Public :: getkin,getknf,getknt,getknr,    &
            kinstress,kinstresf,kinstrest,getvom

  Interface getvom
     Module Procedure getvom
     Module Procedure getvom_rgd
  End Interface !getvom

Contains

  Function getkin(config,vxx,vyy,vzz,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the system kinetic energy
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Real( Kind = wp )                                    :: getkin
    Type( configuration_type ), Intent( In    )                 :: config
    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz
    Type(comms_type), Intent ( InOut )                   :: comm

    Integer           :: i
    Real( Kind = wp ) :: engke

    engke = 0.0_wp

    Do i=1,config%natms
       If (config%lfrzn(i) == 0) &
          engke = engke + config%weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
    End Do

    Call gsum(comm,engke)

    getkin = 0.5_wp * engke

  End Function getkin

  Function getknf(vxx,vyy,vzz,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the kinetic energy of free atoms
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Real( Kind = wp )                                    :: getknf

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz
    Type( configuration_type ),          Intent( InOut ) :: config
    Type(comms_type),                    Intent( InOut ) :: comm

    Integer           :: i,j
    Real( Kind = wp ) :: engke

    engke = 0.0_wp

    Do j=1,config%nfree
       i=config%lstfre(j)

       If (config%lfrzn(i) == 0) &
          engke = engke + config%weight(i)*(vxx(i)**2+vyy(i)**2+vzz(i)**2)
    End Do

    Call gsum(comm,engke)

    getknf = 0.5_wp * engke

  End Function getknf

  Function getknt(rigid,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the translational kinetic energy of
! rigid bodies
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp )                                    :: getknt

    Type( rigid_bodies_type ), Intent( In    ) :: rigid
    Type(comms_type), Intent ( InOut )                     :: comm

    Integer           :: irgd,lrgd,rgdtyp
    Real( Kind = wp ) :: engtra,tmp

    engtra = 0.0_wp

    Do irgd=1,rigid%n_types
       rgdtyp=rigid%list(0,irgd)

       If (rigid%frozen(0,rgdtyp) == 0) Then
          lrgd=rigid%list(-1,irgd)

          tmp=rigid%weight(0,rgdtyp)*Real(rigid%index_local(0,irgd),wp)/Real(lrgd,wp)

          engtra = engtra + tmp*(rigid%vxx(irgd)**2+rigid%vyy(irgd)**2+rigid%vzz(irgd)**2)
       End If
    End Do

    Call gsum(comm,engtra)

    getknt = 0.5_wp * engtra

  End Function getknt

  Function getknr(rigid,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate the rotational kinetic energy of
! rigid bodies
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp )                                    :: getknr

    Type( rigid_bodies_type ), Intent( In    ) :: rigid
    Type(comms_type), Intent ( InOut )                   :: comm

    Integer           :: irgd,lrgd,rgdtyp
    Real( Kind = wp ) :: engrot,tmp

    engrot = 0.0_wp

    Do irgd=1,rigid%n_types
       rgdtyp=rigid%list(0,irgd)

       lrgd=rigid%list(-1,irgd)
       If (rigid%frozen(0,rgdtyp) < lrgd) Then
          tmp=Real(rigid%index_local(0,irgd),wp)/Real(lrgd,wp)

          engrot = engrot + tmp *                     &
                   (rigid%rix(1,rgdtyp)*rigid%oxx(irgd)**2+ &
                    rigid%riy(1,rgdtyp)*rigid%oyy(irgd)**2+ &
                    rigid%riz(1,rgdtyp)*rigid%ozz(irgd)**2)
       End If
    End Do

    Call gsum(comm,engrot)

    getknr = 0.5_wp * engrot

  End Function getknr

  Subroutine kinstress(vxx,vyy,vzz,strkin,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate kinetic contribution to the stress
! tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov august 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz
    Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: strkin
    Type( configuration_type ),          Intent( InOut ) :: config
    Type(comms_type), Intent ( InOut )                   :: comm

    Integer :: i

    strkin = 0.0_wp

    Do i=1,config%natms
       If (config%lfrzn(i) == 0) Then
          strkin(1) = strkin(1) + config%weight(i)*vxx(i)*vxx(i)
          strkin(2) = strkin(2) + config%weight(i)*vxx(i)*vyy(i)
          strkin(3) = strkin(3) + config%weight(i)*vxx(i)*vzz(i)
          strkin(5) = strkin(5) + config%weight(i)*vyy(i)*vyy(i)
          strkin(6) = strkin(6) + config%weight(i)*vyy(i)*vzz(i)
          strkin(9) = strkin(9) + config%weight(i)*vzz(i)*vzz(i)
       End If
    End Do

    Call gsum(comm,strkin)

! Symmetrise

    strkin(4) = strkin(2)
    strkin(7) = strkin(3)
    strkin(8) = strkin(6)

  End Subroutine kinstress

  Subroutine kinstresf(vxx,vyy,vzz,strknf,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate kinetic contribution of free atoms to
! the stress tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Dimension( 1:* ), Intent( In    ) :: vxx,vyy,vzz
    Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: strknf
    Type( configuration_type ),          Intent( InOut ) :: config
    Type(comms_type), Intent ( InOut )                   :: comm

    Integer :: i,j

    strknf = 0.0_wp

    Do j=1,config%nfree
       i=config%lstfre(j)

       If (config%lfrzn(i) == 0) Then
          strknf(1) = strknf(1) + config%weight(i)*vxx(i)*vxx(i)
          strknf(2) = strknf(2) + config%weight(i)*vxx(i)*vyy(i)
          strknf(3) = strknf(3) + config%weight(i)*vxx(i)*vzz(i)
          strknf(5) = strknf(5) + config%weight(i)*vyy(i)*vyy(i)
          strknf(6) = strknf(6) + config%weight(i)*vyy(i)*vzz(i)
          strknf(9) = strknf(9) + config%weight(i)*vzz(i)*vzz(i)
       End If
    End Do

    Call gsum(comm,strknf)

! Symmetrise

    strknf(4) = strknf(2)
    strknf(7) = strknf(3)
    strknf(8) = strknf(6)

  End Subroutine kinstresf

  Subroutine kinstrest(rigid,strknt,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate kinetic contribution of RB
! translational motion to the stress tensor
!
! copyright - daresbury laboratory
! author    - i.t.todorov september 2008
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( rigid_bodies_type ), Intent( In    ) :: rigid
    Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: strknt
    Type(comms_type), Intent ( InOut )                   :: comm

    Integer           :: irgd,lrgd,rgdtyp
    Real( Kind = wp ) :: tmp

    strknt = 0.0_wp

    Do irgd=1,rigid%n_types
       rgdtyp=rigid%list(0,irgd)

       If (rigid%frozen(0,rgdtyp) == 0) Then
          lrgd=rigid%list(-1,irgd)

          tmp=rigid%weight(0,rgdtyp)*Real(rigid%index_local(0,irgd),wp)/Real(lrgd,wp)

          strknt(1) = strknt(1) + tmp*rigid%vxx(irgd)*rigid%vxx(irgd)
          strknt(2) = strknt(2) + tmp*rigid%vxx(irgd)*rigid%vyy(irgd)
          strknt(3) = strknt(3) + tmp*rigid%vxx(irgd)*rigid%vzz(irgd)
          strknt(5) = strknt(5) + tmp*rigid%vyy(irgd)*rigid%vyy(irgd)
          strknt(6) = strknt(6) + tmp*rigid%vyy(irgd)*rigid%vzz(irgd)
          strknt(9) = strknt(9) + tmp*rigid%vzz(irgd)*rigid%vzz(irgd)
       End If
    End Do

    Call gsum(comm,strknt)

! Symmetrise

    strknt(4) = strknt(2)
    strknt(7) = strknt(3)
    strknt(8) = strknt(6)

  End Subroutine kinstrest

  

  Subroutine getvom(vom,vxx,vyy,vzz,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate system centre of mass momentum
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent(   Out ) :: vom(1:3)
    Real( Kind = wp ), Intent( In    ) :: vxx(:),vyy(:),vzz(:)
    Type( configuration_type ), Intent( InOut ) :: config
    Type(comms_type), Intent ( InOut ) :: comm

    Integer                 :: i

! total system mass not including frozen mass

    If (config%newjob_totmas) Then
      config%newjob_totmas = .false.

! For all unfrozen, free particles

      config%totmas = 0.0_wp
       Do i=1,config%natms
         If (config%lfrzn(i) == 0) config%totmas = config%totmas + config%weight(i)
       End Do

       Call gsum(comm,config%totmas)
     End If

     vom = 0.0_wp

     If (.not.config%lvom) Return

! For all unfrozen, free particles

     Do i=1,config%natms
       If (config%lfrzn(i) == 0) Then
         vom(1) = vom(1) + config%weight(i)*vxx(i)
         vom(2) = vom(2) + config%weight(i)*vyy(i)
         vom(3) = vom(3) + config%weight(i)*vzz(i)
       End If
     End Do

     Call gsum(comm,vom)
     If (config%totmas >= zero_plus) vom = vom/config%totmas

   End Subroutine getvom

   Subroutine getvom_rgd(vom,vxx,vyy,vzz,rigid,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to calculate system centre of mass momentum
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2013
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     Real( Kind = wp ), Intent(   Out ) :: vom(1:3)
     Real( Kind = wp ), Intent( In    ) :: vxx(:),vyy(:),vzz(:)
     Type( rigid_bodies_type ), Intent( In    ) :: rigid
     Type( configuration_type ), Intent( InOut ) :: config
     Type(comms_type), Intent ( InOut ) :: comm

     Integer                 :: i,j,irgd,lrgd,rgdtyp
     Real( Kind = wp )       :: tmp

! total system mass not including frozen mass

     If (config%newjob_totmas_r) Then
       config%newjob_totmas_r = .false.

       config%totmas_r = 0.0_wp

! For all unfrozen, free particles

       Do j=1,config%nfree
         i=config%lstfre(j)

         If (config%lfrzn(i) == 0) config%totmas_r = config%totmas_r + config%weight(i)
       End Do

! For all RBs without any frozen particles
! These with some frozen only turn and
! have zero net momentum enforced!

       Do irgd=1,rigid%n_types
         rgdtyp=rigid%list(0,irgd)

         lrgd=rigid%list(-1,irgd)
         If (rigid%frozen(0,rgdtyp) == 0) &
           config%totmas_r = config%totmas_r + rigid%weight(0,rgdtyp)*Real(rigid%index_local(0,irgd),wp)/Real(lrgd,wp)
       End Do

       Call gsum(comm,config%totmas_r)
     End If

     vom = 0.0_wp

     If (.not.config%lvom) Return

! For all unfrozen, free particles

     Do j=1,config%nfree
       i=config%lstfre(j)

       If (config%lfrzn(i) == 0) Then
         tmp=config%weight(i)
         vom(1) = vom(1) + tmp*vxx(i)
         vom(2) = vom(2) + tmp*vyy(i)
         vom(3) = vom(3) + tmp*vzz(i)
       End If
     End Do

     Do irgd=1,rigid%n_types
       rgdtyp=rigid%list(0,irgd)

! For all RBs without any frozen particles

       lrgd=rigid%list(-1,irgd)
       If (rigid%frozen(0,rgdtyp) == 0) Then
         tmp=rigid%weight(0,rgdtyp)*Real(rigid%index_local(0,irgd),wp)/Real(lrgd,wp)

         vom(1) = vom(1) + tmp*rigid%vxx(irgd)
         vom(2) = vom(2) + tmp*rigid%vyy(irgd)
         vom(3) = vom(3) + tmp*rigid%vzz(irgd)
       End If
     End Do

     Call gsum(comm,vom)
     If (config%totmas_r >= zero_plus) vom = vom/config%totmas_r

   End Subroutine getvom_rgd


   Subroutine cap_forces(temp,config,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to limit the absolute magnitude of forces.
! To be used in equilibration period only!!!
!
! copyright - daresbury laboratory
! author    - i.t.todorov july 2004
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     Real( Kind = wp ), Intent( In    ) :: temp
     Type(comms_type),  Intent( InOut ) :: comm
     Type( configuration_type ),  Intent( InOut ) :: config

     Integer           :: i
     Real( Kind = wp ) :: fmax2,fmod,scale,fcom(1:3)

    If (config%newjob_meg) Then
      config%newjob_meg = .false.

! Get the total number of non-frozen&non-massless particles

      config%meg=0.0_wp
      Do i=1,config%natms
        If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) &
          config%meg=config%meg+1.0_wp
      End Do
      Call gsum(comm,config%meg)
    End If

! maximum force permitted

    fmax2 = (boltz*temp*config%fmax)**2

! cap forces and conserve linear momentum
! for non-frozen&non-massless particles

    fcom = 0.0_wp
    Do i=1,config%natms
      If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
        fmod = config%parts(i)%fxx**2 + config%parts(i)%fyy**2 + config%parts(i)%fzz**2

        If (fmod > fmax2) Then
          scale  = Sqrt(fmax2/fmod)
          config%parts(i)%fxx = config%parts(i)%fxx*scale
          config%parts(i)%fyy = config%parts(i)%fyy*scale
          config%parts(i)%fzz = config%parts(i)%fzz*scale
        End If

! accumulate forces - to check on momentum conservation

        fcom(1) = fcom(1) + config%parts(i)%fxx
        fcom(2) = fcom(2) + config%parts(i)%fyy
        fcom(3) = fcom(3) + config%parts(i)%fzz
      End If
    End Do

! ensure net forces sum to zero

    Call gsum(comm,fcom)
    fcom = fcom/config%meg

! conserve momentum

    Do i=1,config%natms
      If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
        config%parts(i)%fxx = config%parts(i)%fxx - fcom(1)
        config%parts(i)%fyy = config%parts(i)%fyy - fcom(2)
        config%parts(i)%fzz = config%parts(i)%fzz - fcom(3)
      End If
    End Do

  End Subroutine cap_forces

End Module kinetics
