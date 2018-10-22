Module external_field

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global external field variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov july 2004
  ! refactoring:
  !           - a.m.elena march-october 2018
  !           - j.madge march-october 2018
  !           - a.b.g.chalk march-october 2018
  !           - i.scivetti march-october 2018
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use comms,   Only : comms_type,gcheck,gsum
  Use constants,   Only : twopi
  Use configuration,  Only : configuration_type
  Use particle, Only : corePart
  Use kinetics, Only : getcom_mol
  Use rigid_bodies, Only : rigid_bodies_type
  Use errors_warnings, Only : error
  Use numerics, Only : local_index,images
  Use shared_units, Only : update_shared_units
  Use statistics, Only : stats_type
  Use core_shell, Only : core_shell_type,SHELL_ADIABATIC
  Use rdfs, Only : rdf_type,usr_collect,usr_compute
  Use domains, Only : domains_type
  Implicit None

  Private

  ! External field type keys
  !> Null, no field
  Integer( Kind = wi ), Parameter, Public :: FIELD_NULL = 0
  !> Electric field
  Integer( Kind = wi ), Parameter, Public :: FIELD_ELECTRIC = 1
  !> Oscillating shear, orthorhombic box, $F_{x} = a * \cos (b*2*\pi*z/L)$
  Integer( Kind = wi ), Parameter, Public :: FIELD_SHEAR_OSCILLATING = 2
  !> Continuous shear, 2D perioidic box
  Integer( Kind = wi ), Parameter, Public :: FIELD_SHEAR_CONTINUOUS = 3
  !> Gravitational field
  Integer( Kind = wi ), Parameter, Public :: FIELD_GRAVITATIONAL = 4
  !> Magnetic field
  Integer( Kind = wi ), Parameter, Public :: FIELD_MAGNETIC = 5
  !> Containing sphere, $r^{-n}$ potential
  Integer( Kind = wi ), Parameter, Public :: FIELD_SPHERE = 6
  !> Repulsive wall (harmonic) starting at z0
  Integer( Kind = wi ), Parameter, Public :: FIELD_WALL = 7
  !> Piston wall pushing along the X=bxc direction
  Integer( Kind = wi ), Parameter, Public :: FIELD_WALL_PISTON = 8
  !> zres external field. Restrain molecule z-position (push in)
  Integer( Kind = wi ), Parameter, Public :: FIELD_ZRES = 9
  !> zres- external field (push out)
  Integer( Kind = wi ), Parameter, Public :: FIELD_ZRES_MINUS = 10
  !> zres+ external field (pull in)
  Integer( Kind = wi ), Parameter, Public :: FIELD_ZRES_PLUS = 11
  !> Oscillating electric field
  Integer( Kind = wi ), Parameter, Public :: FIELD_ELECTRIC_OSCILLATING = 12
  !> Umbrella potential harmonic restraint (pull in)
  Integer( Kind = wi ), Parameter, Public :: FIELD_UMBRELLA = 13

  !> Type to hold external field data
  Type, Public :: external_field_type
    Private

    !> Type of external field
    Integer( Kind = wi ), Public :: key = FIELD_NULL

    !> Mass?
    Real( Kind = wp ) :: mass = 0.0_wp

    !> Conversion factor to handle magnetic and electric units
    Real( Kind = wp ), Public :: conv_fact

    !> Field parameters
    Real( Kind = wp ), Allocatable, Public :: param(:)
    !> Number of parameters
    Integer( Kind = wi ), Public :: max_param
    Logical :: newjob


  Contains
    Private

    Procedure, Public :: init => allocate_external_field_arrays
    Final :: cleanup
  End Type external_field_type

  Public :: external_field_apply, external_field_correct

Contains

  !> Allocate and initialise the arrays of the external field type
  Subroutine allocate_external_field_arrays(T)
  Class( external_field_type) :: T

    Integer :: fail

    fail = 0

    Allocate (T%param(T%max_param), stat=fail)
    If (fail > 0) Call error(1019)

    T%param = 0.0_wp
  End Subroutine allocate_external_field_arrays

  Subroutine cleanup(T)
    Type( external_field_type ) :: T

    If (Allocated(T%param)) Then
      Deallocate(T%param)
    End If
  End subroutine cleanup

  Subroutine external_field_apply(time,leql,nsteql,nstep,cshell,stats,rdf, &
      ext_field,rigid,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for application of an external field
    !
    ! Note: Only one field at a time is allowed
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov november 2016
    ! amended   - i.t.todorov september 2017 :: zres, zrs+ and zrs- fields
    !                                           gsum stats%engfld and stats%virfld
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,           Intent( In    ) :: leql
    Integer,           Intent( In    ) :: nsteql,nstep
    Real( Kind = wp ), Intent( In    ) :: time ! for oscillating fields
    Type( stats_type ), Intent( Inout ) :: stats
    Type( core_shell_type ), Intent( Inout ) :: cshell
    Type( rdf_type ), Intent( InOut ) :: rdf
    Type( external_field_type ), Intent( InOut ) :: ext_field
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( Inout ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config

    Logical           :: safe,l1,l2
    Integer           :: i,j,ia,ib,ic,id,fail(1:2), &
      irgd,jrgd,lrgd,rgdtyp
    Real( Kind = wp ) :: gamma,rrr,rz,zdif,vxt,vyt,vzt,tmp,rtmp(1:2), &
      x(1:1),y(1:1),z(1:1),cmm(0:3),cm2(0:3)

    Integer,           Allocatable :: lstopt(:,:)
    Real( Kind = wp ), Allocatable :: oxt(:),oyt(:),ozt(:)
    Character( Len = 256 ) :: message

    ! energy and virial accumulators

    stats%engfld=0.0_wp
    stats%virfld=0.0_wp

    If (ext_field%key == FIELD_ELECTRIC) Then

      ! electric field: ext_field%param(1-3) are field components

      Do i=1,config%natms
        If (config%lfrzn(i) == 0) Then
          config%parts(i)%fxx=config%parts(i)%fxx + config%parts(i)%chge*ext_field%param(1)
          config%parts(i)%fyy=config%parts(i)%fyy + config%parts(i)%chge*ext_field%param(2)
          config%parts(i)%fzz=config%parts(i)%fzz + config%parts(i)%chge*ext_field%param(3)
        End If
      End Do

    Else If (ext_field%key == FIELD_SHEAR_OSCILLATING) Then

      ! oscillating shear: orthorhombic box:  Fx=a*Cos(b.2.pi.z/L)

      If (config%imcon /= 1 .and. config%imcon /= 2) Return

      rz=twopi/config%cell(9)

      Do i=1,config%natms
        If (config%lfrzn(i) == 0) config%parts(i)%fxx=config%parts(i)%fxx +&
          ext_field%param(1)*Cos(ext_field%param(2)*config%parts(i)%zzz*rz)
      End Do

    Else If (ext_field%key == FIELD_SHEAR_CONTINUOUS) Then

      ! continuous shear of walls : 2D periodic box (imcon=6)

      If (config%imcon /= 6) Return

      ! shear rate=ext_field%param(1) angstrom per ps for non-frozen
      ! and non-weightless atoms at Abs(z) > ext_field%param(2)

      If (rigid%total > 0) Then

        ! FPs
        Do j=1,config%nfree
          i=config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. Abs(config%parts(i)%zzz) > ext_field%param(2)) &
            config%vxx(i)=0.5_wp*Sign(ext_field%param(1),config%parts(i)%zzz)
        End Do

        ! RBs

        Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

          ! For all good RBs

          lrgd=rigid%list(-1,irgd)
          If (rigid%frozen(0,rgdtyp) == 0) Then
            x=rigid%xxx(irgd) ; y=rigid%yyy(irgd) ; z=rigid%zzz(irgd)
            Call images(config%imcon,config%cell,1,x,y,z)

            rz=z(1)
            If (Abs(rz) > ext_field%param(2)) Then
              tmp=0.5_wp*Sign(ext_field%param(1),rz)
              vxt=tmp-rigid%vxx(irgd)

              rigid%vxx(irgd)=tmp
              Do jrgd=1,lrgd
                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                If (i <= config%natms) config%vxx(i)=config%vxx(i)+vxt
              End Do
            End If
          End If
        End Do

      Else

        Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. Abs(config%parts(i)%zzz) > ext_field%param(2)) &
            config%vxx(i)=0.5_wp*Sign(ext_field%param(1),config%parts(i)%zzz)
        End Do

      End If

    Else If (ext_field%key == FIELD_GRAVITATIONAL) Then

      ! gravitational field: field components given by ext_field%param(1-3)

      If (cshell%keyshl == SHELL_ADIABATIC) Then

        Do i=1,config%natms
          If (config%lfrzn(i) == 0) Then
            config%parts(i)%fxx=config%parts(i)%fxx + ext_field%param(1)*config%weight(i)
            config%parts(i)%fyy=config%parts(i)%fyy + ext_field%param(2)*config%weight(i)
            config%parts(i)%fzz=config%parts(i)%fzz + ext_field%param(3)*config%weight(i)
          End If
        End Do

      Else

        fail=0
        Allocate (lstopt(1:2,1:cshell%mxshl),                       Stat=fail(1))
        Allocate (oxt(1:config%mxatms),oyt(1:config%mxatms),ozt(1:config%mxatms), Stat=fail(2))
        If (Any(fail > 0)) Then
          Write(message,'(a)') 'external_field_apply allocation failure'
          Call error(0,message)
        End If

        Do i=1,cshell%ntshl
          lstopt(1,i)=local_index(cshell%listshl(1,i),config%nlast,config%lsi,config%lsa)
          lstopt(2,i)=local_index(cshell%listshl(2,i),config%nlast,config%lsi,config%lsa)
        End Do

        Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            oxt(i)=ext_field%param(1)*config%weight(i)
            oyt(i)=ext_field%param(2)*config%weight(i)
            ozt(i)=ext_field%param(3)*config%weight(i)
          Else ! for the sake of massless sites of RBs
            oxt(i)=0.0_wp
            oyt(i)=0.0_wp
            ozt(i)=0.0_wp
          End If
        End Do

        If (cshell%lshmv_shl) Then
          Call update_shared_units(config,cshell%lishp_shl, &
            cshell%lashp_shl,oxt,oyt,ozt,domain,comm)
        End If

        ! Transfer cores' forces to shells

        Do i=1,cshell%ntshl
          ia=lstopt(1,i)
          ib=lstopt(2,i)
          If (ia > 0 .and. (ib > 0 .and. ib <= config%natms)) Then
            oxt(ib)=oxt(ia)
            oyt(ib)=oyt(ia)
            ozt(ib)=ozt(ia)
          End If
        End Do

        Do i=1,config%natms
          If (config%lfrzn(i) == 0) Then
            config%parts(i)%fxx=config%parts(i)%fxx + oxt(i)
            config%parts(i)%fyy=config%parts(i)%fyy + oyt(i)
            config%parts(i)%fzz=config%parts(i)%fzz + ozt(i)
          End If
        End Do

        Deallocate (lstopt,      Stat=fail(1))
        Deallocate (oxt,oyt,ozt, Stat=fail(2))
        If (Any(fail > 0)) Then
          Write(message,'(a)') 'external_field_apply deallocation failure'
          Call error(0,message)
        End If

      End If

    Else If (ext_field%key == FIELD_MAGNETIC) Then

      ! magnetic field: field components given by ext_field%param(1-3)

      If (cshell%keyshl == SHELL_ADIABATIC) Then

        Do i=1,config%natms
          If (config%lfrzn(i) == 0) Then
            config%parts(i)%fxx=config%parts(i)%fxx + &
              (config%vyy(i)*ext_field%param(3)-config%vzz(i)*ext_field%param(2))*config%parts(i)%chge
            config%parts(i)%fyy=config%parts(i)%fyy + &
              (config%vzz(i)*ext_field%param(1)-config%vxx(i)*ext_field%param(3))*config%parts(i)%chge
            config%parts(i)%fzz=config%parts(i)%fzz + &
              (config%vxx(i)*ext_field%param(2)-config%vyy(i)*ext_field%param(1))*config%parts(i)%chge
          End If
        End Do

      Else

        fail=0
        Allocate (lstopt(1:2,1:cshell%mxshl),                       Stat=fail(1))
        Allocate (oxt(1:config%mxatms),oyt(1:config%mxatms),ozt(1:config%mxatms), Stat=fail(2))
        If (Any(fail > 0)) Then
          Write(message,'(a)') 'external_field_apply allocation failure'
          Call error(0,message)
        End If

        Do i=1,cshell%ntshl
          lstopt(1,i)=local_index(cshell%listshl(1,i),config%nlast,config%lsi,config%lsa)
          lstopt(2,i)=local_index(cshell%listshl(2,i),config%nlast,config%lsi,config%lsa)
        End Do

        ! cores' velocities

        Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp) Then
            oxt(i)=config%vxx(i)
            oyt(i)=config%vyy(i)
            ozt(i)=config%vzz(i)
          Else
            oxt(i)=0.0_wp
            oyt(i)=0.0_wp
            ozt(i)=0.0_wp
          End If
        End Do

        If (cshell%lshmv_shl) Then
          Call update_shared_units(config,cshell%lishp_shl,cshell%lashp_shl,oxt,oyt,ozt,domain,comm)
        End If

        ! Transfer cores' velocities to shells

        Do i=1,cshell%ntshl
          ia=lstopt(1,i)
          ib=lstopt(2,i)
          If (ia > 0 .and. (ib > 0 .and. ib <= config%natms)) Then
            oxt(ib)=oxt(ia)
            oyt(ib)=oyt(ia)
            ozt(ib)=ozt(ia)
          End If
        End Do

        Do i=1,config%natms
          If (config%lfrzn(i) == 0) Then
            config%parts(i)%fxx=config%parts(i)%fxx + &
              (oyt(i)*ext_field%param(3)-ozt(i)*ext_field%param(2))*config%parts(i)%chge
            config%parts(i)%fyy=config%parts(i)%fyy + &
              (ozt(i)*ext_field%param(1)-oxt(i)*ext_field%param(3))*config%parts(i)%chge
            config%parts(i)%fzz=config%parts(i)%fzz + &
              (oxt(i)*ext_field%param(2)-oyt(i)*ext_field%param(1))*config%parts(i)%chge
          End If
        End Do

        Deallocate (lstopt,      Stat=fail(1))
        Deallocate (oxt,oyt,ozt, Stat=fail(2))
        If (Any(fail > 0)) Then
          Write(message,'(a,i0)') 'external_field_apply deallocation failure'
          Call error(0,message)
        End If

      End If

    Else If (ext_field%key == FIELD_SPHERE) Then

      ! containing sphere : r^(-n) potential

      Do i=1,config%natms
        If (config%lfrzn(i) == 0) Then
          rrr=Sqrt(config%parts(i)%xxx**2+config%parts(i)%yyy**2+config%parts(i)%zzz**2)
          If (rrr > ext_field%param(4)) Then
            rrr=ext_field%param(2)-rrr
            If (rrr < 0.0_wp) rrr=0.1_wp

            gamma=ext_field%param(1)*rrr**(-ext_field%param(3))
            stats%engfld=stats%engfld + gamma
            gamma=-ext_field%param(3)*gamma/(rrr*rrr)

            config%parts(i)%fxx=config%parts(i)%fxx + gamma*config%parts(i)%xxx
            config%parts(i)%fyy=config%parts(i)%fyy + gamma*config%parts(i)%yyy
            config%parts(i)%fzz=config%parts(i)%fzz + gamma*config%parts(i)%zzz
          End If
        End If
      End Do

      stats%virfld=-9.0_wp*stats%engfld

    Else If (ext_field%key == FIELD_WALL) Then

      ! repulsive wall (harmonic) starting at z0

      Do i=1,config%natms
        If (config%lfrzn(i) == 0 .and. ext_field%param(3)*config%parts(i)%zzz > ext_field%param(3)*ext_field%param(2)) Then
          zdif=config%parts(i)%zzz-ext_field%param(2)
          gamma=-ext_field%param(1)*zdif

          config%parts(i)%fzz=config%parts(i)%fzz + gamma
          stats%engfld=stats%engfld - gamma*zdif
        End If
      End Do

      stats%engfld=0.5_wp*stats%engfld

    Else If (ext_field%key == FIELD_WALL_PISTON) Then

      ! xpist - piston wall pushing down along the X=bxc direction
      ! ext_field%param(1) is the first atom of the layer of molecules (membrane) to be pushed
      ! ext_field%param(2) is the last atom of the layer of molecules (membrane) to be pushed
      ! ext_field%param(3) is the pressure applied to the layer of molecules (membrane) in the
      ! +X=bxc direction - i.e. left to right.  The layer plane is defined as _|_ bxc

      If (config%imcon /= 1 .and. config%imcon /= 2) Return

      ia = Nint(ext_field%param(1))
      ib = Nint(ext_field%param(2))

      If (ext_field%newjob) Then
        ext_field%newjob=.false.

        !        ext_field%mass=0.0_wp ! defined and initialise in external_field_module
        safe=.true.
        Do i=1,config%natms
          If ((config%ltg(i) >= ia .and. config%ltg(i) <= ib) .and. config%lfrzn(i) > 0) Then
            safe=.false.
          Else
            ext_field%mass=ext_field%mass+config%weight(i)
          End If
        End Do

        Call gsum(comm,ext_field%mass)
        Call gcheck(comm,safe)
        If (.not.safe) Call error(456)
      End If

      rtmp=0.0_wp  ! average velocity and force per atom in x direction of the piston
      Do i=1,config%natms ! preserve momentum and velocity in the direction of the push
        If (config%ltg(i) >= ia .and. config%ltg(i) <= ib) Then
          rtmp(1)=rtmp(1)+config%weight(i)*config%vxx(i) ; rtmp(2)=rtmp(2)+config%parts(i)%fxx
          config%vyy(i) = 0.0_wp                  ; config%parts(i)%fyy = 0.0_wp
          config%vzz(i) = 0.0_wp                  ; config%parts(i)%fzz = 0.0_wp
        End If
      End Do
      Call gsum(comm,rtmp) ! net velocity and force to ensure solid wall behaviour

      rtmp(1)=rtmp(1)/ext_field%mass             ! averaged velocity per particle
      rtmp(2)=(rtmp(2)+ext_field%param(3))/ext_field%mass ! averaged acceleration of the slab

      Do i=1,config%natms
        If (config%ltg(i) >= ia .and. config%ltg(i) <= ib) Then
          stats%engfld=config%weight(i)*(rtmp(1)-config%vxx(i))**2 ! must change E_kin to reflect solidity
          config%vxx(i)=rtmp(1)
          config%parts(i)%fxx=rtmp(2)*config%weight(i) ! force per particle
        End If
      End Do

      stats%engfld=0.5_wp*stats%engfld

    Else If (ext_field%key == FIELD_ZRES) Then

      ! zres external field: restrain molecule z-position (pull in)
      ! ext_field%param(1) is the index of first atom of restrained molecule
      ! ext_field%param(2) is the index of last atom of restrained molecule
      ! ext_field%param(3) is the restraining constant
      ! ext_field%param(4) is z-min (min limit in z-direction)
      ! ext_field%param(5) is z-max (max limit in z-direction)
      ! where ext_field%param(4) < ext_field%param(5)

      ia = Nint(ext_field%param(1))
      ib = Nint(ext_field%param(2))

      ! Get molecule's weight and COM

      Call getcom_mol(config,ia,ib,cmm,comm)

      ! Apply force corrections

      Do i=1,config%natms
        If (config%ltg(i) >= ia .and. config%ltg(i) <= ib .and. config%lfrzn(i) == 0) Then
          If (cmm(3) < ext_field%param(4) .or. cmm(3) > ext_field%param(5)) Then
            If (cmm(3) < ext_field%param(4)) zdif = cmm(3) - ext_field%param(4)
            If (cmm(3) > ext_field%param(5)) zdif = cmm(3) - ext_field%param(5)

            gamma=-ext_field%param(3)*zdif*config%weight(i)/cmm(0)
            config%parts(i)%fzz=config%parts(i)%fzz + gamma
            stats%engfld=stats%engfld - gamma*zdif
          End If
        End If
      End Do

      stats%engfld=0.5_wp*stats%engfld

    Else If (ext_field%key == FIELD_ZRES_MINUS) Then

      ! extension to zrs- external field (push out)
      ! ext_field%param(1) is the index of first atom of the water molecules to be restrained
      ! ext_field%param(2) is the index of last atom of the water molecules to be restrained
      ! ext_field%param(3) is the restraining constant
      ! ext_field%param(4) is z-min (min limit in z-direction)
      ! ext_field%param(5) is z-max (max limit in z-direction)
      ! where ext_field%param(4) < ext_field%param(5)
      ! This will keep water away from the membrane region, and allow
      ! the DMPC to slowly fill in the gaps before water molecules are let loose.

      ia = Nint(ext_field%param(1))
      ib = Nint(ext_field%param(2))
      Do i=1,config%natms
        If (config%ltg(i) >= ia .and. config%ltg(i) <= ib .and. config%lfrzn(i) == 0) Then
          If (config%parts(i)%zzz > ext_field%param(4) .and. config%parts(i)%zzz < ext_field%param(5)) Then
            tmp = ext_field%param(5) + ext_field%param(4)

            If (config%parts(i)%zzz <  tmp) Then
              zdif = config%parts(i)%zzz - ext_field%param(4)
            Else
              zdif = config%parts(i)%zzz - ext_field%param(5)
            End If

            gamma=-ext_field%param(3)*zdif
            config%parts(i)%fzz=config%parts(i)%fzz + gamma
            stats%engfld=stats%engfld - gamma*zdif
          End If
        End If
      End Do

      stats%engfld=0.5_wp*stats%engfld

    Else If (ext_field%key == FIELD_ZRES_PLUS) Then

      ! extension to zrs+ external field (pull in)
      ! ext_field%param(1) is the index of first atom of the water molecules to be restrained
      ! ext_field%param(2) is the index of last atom of the water molecules to be restrained
      ! ext_field%param(3) is the restraining constant
      ! ext_field%param(4) is z-min (min limit in z-direction)
      ! ext_field%param(5) is z-max (max limit in z-direction)
      ! where ext_field%param(4) < ext_field%param(5)
      ! This will keep water within this region.

      ia = Nint(ext_field%param(1))
      ib = Nint(ext_field%param(2))
      Do i=1,config%natms
        If (config%ltg(i) >= ia .and. config%ltg(i) <= ib .and. config%lfrzn(i) == 0) Then
          If (config%parts(i)%zzz < ext_field%param(4) .or. config%parts(i)%zzz > ext_field%param(5)) Then
            If (config%parts(i)%zzz < ext_field%param(4)) zdif = config%parts(i)%zzz - ext_field%param(4)
            If (config%parts(i)%zzz > ext_field%param(5)) zdif = config%parts(i)%zzz - ext_field%param(5)

            gamma=-ext_field%param(3)*zdif
            config%parts(i)%fzz=config%parts(i)%fzz + gamma
            stats%engfld=stats%engfld - gamma*zdif
          End If
        End If
      End Do

      stats%engfld=0.5_wp*stats%engfld

    Else If (ext_field%key == FIELD_ELECTRIC_OSCILLATING) Then

      ! extension to oscillating electric field: ext_field%param(1-3) are field components
      ! ext_field%param(4) is the oscillating frequency defined in ps^-1!

      tmp=Sin(time*ext_field%param(4)*twopi)
      Do i=1,config%natms
        If (config%lfrzn(i) == 0) Then
          config%parts(i)%fxx=config%parts(i)%fxx + config%parts(i)%chge*ext_field%param(1)*tmp
          config%parts(i)%fyy=config%parts(i)%fyy + config%parts(i)%chge*ext_field%param(2)*tmp
          config%parts(i)%fzz=config%parts(i)%fzz + config%parts(i)%chge*ext_field%param(3)*tmp
        End If
      End Do

    Else If (ext_field%key == FIELD_UMBRELLA) Then

      ! uphr external field: umbrella potential harmonic restraint (pull in)
      ! ext_field%param(1) is the index of first atom of the first atom group/molecule
      ! ext_field%param(2) is the index of last atom of the first atom group/molecule
      ! ext_field%param(3) is the index of first atom of the second atom group/molecule
      ! ext_field%param(4) is the index of last atom of the second atom group/molecule
      ! ext_field%param(5) is the umbrella force constant
      ! ext_field%param(6) is the com separation at the minimum of umbrella potential

      ia = Nint(ext_field%param(1))
      ib = Nint(ext_field%param(2))

      ! Get first molecule's weight and COM

      Call getcom_mol(config,ia,ib,cmm,comm)

      ic = Nint(ext_field%param(3))
      id = Nint(ext_field%param(4))

      ! Get second molecule's weight and COM

      Call getcom_mol(config,ic,id,cm2,comm)

      ! Apply PBC to COMS vector

      x(1)=cm2(1)-cmm(1)
      y(1)=cm2(2)-cmm(2)
      z(1)=cm2(3)-cmm(3)

      Call images(config%imcon,config%cell,1,x,y,z)

      rrr=Sqrt(x(1)**2+y(1)**2+z(1)**2)

      ! accumulate RDF for the 2 COMs every 50 timesteps and
      ! refresh USRDAT every 500 timesteps

      If ((.not.leql) .or. nstep >= nsteql) Then
        If (Mod(nstep,50)  == 0) Call usr_collect(rrr,rdf)
        If (Mod(nstep,500) == 0) Call usr_compute(rdf,config,comm)
      End If

      ! get force magnitude

      zdif  =rrr-ext_field%param(6)
      stats%engfld=-ext_field%param(5)*zdif/rrr

      ! Apply force corrections

      Do i=1,config%natms
        l1=(config%ltg(i) >= ia .and. config%ltg(i) <= ib)
        l2=(config%ltg(i) >= ic .and. config%ltg(i) <= id)

        If ((l1 .or. l2) .and. config%lfrzn(i) == 0) Then
          tmp=stats%engfld*config%weight(i)

          If (l2) Then
            tmp=tmp/cm2(0)
          Else ! If (l2) - no crossover
            tmp=-tmp/cmm(0)
          End If

          config%parts(i)%fxx=config%parts(i)%fxx + tmp*x(1)
          config%parts(i)%fyy=config%parts(i)%fyy + tmp*y(1)
          config%parts(i)%fzz=config%parts(i)%fzz + tmp*z(1)
        End If
      End Do

      stats%virfld=-ext_field%param(5)*zdif*rrr
      stats%engfld=0.5_wp*ext_field%param(5)*zdif**2

    Else

      ! unidentified field potential error exit

      Call error(454)

    End If

    ! sum up energy and virial contributions

    rtmp(1) = stats%engfld
    rtmp(2) = stats%virfld

    Call gsum(comm,rtmp)

    stats%engfld = rtmp(1)
    stats%virfld = rtmp(2)
  End Subroutine external_field_apply

  Subroutine external_field_correct(engfld,ext_field,rigid,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for correcting an external field application
    !
    ! Note: Only one field at a time is allowed
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2015
    ! amnded    - i.t.todorov september 2015 : gsum stats%engfld
    ! refactoring:
    !           - a.m.elena march-october 2018
    !           - j.madge march-october 2018
    !           - a.b.g.chalk march-october 2018
    !           - i.scivetti march-october 2018
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ), Intent(   Out ) :: engfld
    Type( external_field_type ), Intent( In    ) :: ext_field
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( comms_type ), Intent( InOut ) :: comm
    Type( configuration_type ),               Intent( InOut ) :: config

    Integer           :: i,j,ia,ib, irgd,jrgd,lrgd,rgdtyp
    Real( Kind = wp ) :: rz,vxt,tmp,rtmp(1:2), &
      x(1:1),y(1:1),z(1:1)

    If (ext_field%key == FIELD_SHEAR_CONTINUOUS) Then

      ! continuous shear of walls : 2D periodic box (imcon=6)

      If (config%imcon /= 6) Return

      ! shear rate=ext_field%param(1) angstrom per ps for non-frozen
      ! and non-weightless atoms at Abs(z) > ext_field%param(2)

      If (rigid%total > 0) Then

        ! FPs
        Do j=1,config%nfree
          i=config%lstfre(j)

          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. Abs(config%parts(i)%zzz) > ext_field%param(2)) &
            config%vxx(i)=0.5_wp*Sign(ext_field%param(1),config%parts(i)%zzz)
        End Do

        ! RBs

        Do irgd=1,rigid%n_types
          rgdtyp=rigid%list(0,irgd)

          ! For all good RBs

          lrgd=rigid%list(-1,irgd)
          If (rigid%frozen(0,rgdtyp) == 0) Then
            x=rigid%xxx(irgd) ; y=rigid%yyy(irgd) ; z=rigid%zzz(irgd)
            Call images(config%imcon,config%cell,1,x,y,z)

            rz=z(1)
            If (Abs(rz) > ext_field%param(2)) Then
              tmp=0.5_wp*Sign(ext_field%param(1),rz)
              vxt=tmp-rigid%vxx(irgd)

              rigid%vxx(irgd)=tmp
              Do jrgd=1,lrgd
                i=rigid%index_local(jrgd,irgd) ! local index of particle/site

                If (i <= config%natms) config%vxx(i)=config%vxx(i)+vxt
              End Do
            End If
          End If
        End Do

      Else

        Do i=1,config%natms
          If (config%lfrzn(i) == 0 .and. config%weight(i) > 1.0e-6_wp .and. Abs(config%parts(i)%zzz) > ext_field%param(2)) &
            config%vxx(i)=0.5_wp*Sign(ext_field%param(1),config%parts(i)%zzz)
        End Do

      End If

    Else If (ext_field%key == FIELD_WALL_PISTON) Then

      engfld = 0.0_wp

      ! xpist - piston wall pushing down along the X=bxc direction
      ! ext_field%param(1) is the first atom of the layer of molecules (membrane) to be pushed
      ! ext_field%param(2) is the last atom of the layer of molecules (membrane) to be pushed
      ! ext_field%param(3) is the pressure applied to the layer of molecules (membrane) in the
      ! +X=bxc direction - i.e. left to right.  The layer plane is defined as _|_ bxc

      If (config%imcon /= 1 .and. config%imcon /= 2) Return

      ia = Nint(ext_field%param(1))
      ib = Nint(ext_field%param(2))

      rtmp=0.0_wp  ! average velocity and force per atom in x direction of the piston
      Do i=1,config%natms ! preserve momentum and velocity in the direction of the push
        If (config%ltg(i) >= ia .and. config%ltg(i) <= ib) Then
          rtmp(1)=rtmp(1)+config%weight(i)*config%vxx(i) ; rtmp(2)=rtmp(2)+config%parts(i)%fxx
          config%vyy(i) = 0.0_wp                  ; config%parts(i)%fyy = 0.0_wp
          config%vzz(i) = 0.0_wp                  ; config%parts(i)%fzz = 0.0_wp
        End If
      End Do
      Call gsum(comm,rtmp) ! net velocity and force to ensure solid wall behaviour

      rtmp(1)=rtmp(1)/ext_field%mass             ! averaged velocity per particle
      rtmp(2)=(rtmp(2)+ext_field%param(3))/ext_field%mass ! averaged acceleration of the slab

      Do i=1,config%natms
        If (config%ltg(i) >= ia .and. config%ltg(i) <= ib) Then
          engfld=config%weight(i)*(rtmp(1)-config%vxx(i))**2 ! must change E_kin to reflect solidity
          config%vxx(i)=rtmp(1)
          config%parts(i)%fxx=rtmp(2)*config%weight(i) ! force per particle
        End If
      End Do

      engfld=0.5_wp*engfld
      Call gsum(comm,engfld)
    End If
  End Subroutine external_field_correct
End Module external_field
