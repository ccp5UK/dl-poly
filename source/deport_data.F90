Module deport_data
  Use kinds,            Only : wp
  Use comms,            Only : comms_type,gcheck,wp_mpi, Deport_tag, &
    Export_tag, MetLdExp_tag, ExpMplRM_tag, &
    PassUnit_tag,gsend,gwait,girecv,gmax,gsum
  Use constants,        Only : half_plus, half_minus 
  Use domains, Only : domains_type
  Use configuration, Only : configuration_type
  Use site, Only : site_type
  Use flow_control, Only : flow_type
  Use rigid_bodies, Only : rigid_bodies_type
  Use tethers,      Only : tethers_type
  Use bonds,        Only : bonds_type
  Use angles,       Only : angles_type
  Use dihedrals,    Only : dihedrals_type
  Use inversions, Only : inversions_type
  Use statistics, Only : stats_type
  Use minimise,     Only : minimise_type
  Use ewald,               Only : ewald_type
  Use mpole ,              Only : mpole_type,POLARISATION_CHARMM
  Use msd, Only : msd_type
  Use greenkubo,    Only : greenkubo_type
  Use core_shell,   Only : core_shell_type 
  Use constraints,  Only : constraints_type 
  Use errors_warnings, Only : error, warning
  Use mpoles_container, Only : rotate_mpoles, rotate_mpoles_d
  Use numerics, Only : local_index,dcell,invert,shellsort2,pbcshift
  Use pmf, Only : pmf_units_set, pmf_type
  Use build_book, Only : compress_book_intra
  Use shared_units, Only : pass_shared_units, tag_legend
  Use thermostat, Only : thermostat_type
  Use neighbours, Only : neighbours_type
  Use kim, Only : kim_type
  Implicit None

  Public :: deport_atomic_data, export_atomic_data

Contains

  Subroutine deport_atomic_data(mdir,lbook,lmsd,cshell,cons,pmf,stats,ewld,thermo,&
      green,bond,angle,dihedral,inversion,tether,neigh,minim,mpoles,rigid,domain, &
      config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to deport atomic and topological data of particles
    ! leaving this domain
    !
    ! NOTE: When executing on one node we need not get here at all!
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith & i.t.todorov december 2016
    ! contrib   - i.j.bush february 2014
    ! contrib   - m.a.seaton june 2014
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    Integer,            Intent( In    ) :: mdir
    Logical,            Intent( In    ) :: lbook
    Logical,            Intent( In    ) :: lmsd
    Type( pmf_type) , Intent( InOut ) :: pmf
    Type( constraints_type) , Intent( InOut ) :: cons
    Type( stats_type ), Intent( InOut ) :: stats
    Type( ewald_type ), Intent( InOut ) :: ewld
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( greenkubo_type ), Intent( InOut ) :: green
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( tethers_type ), Intent( InOut ) :: tether
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( minimise_type ), Intent( InOut ) :: minim
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( configuration_type ), Intent( InOut ) :: config
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut ) :: comm

    Logical           :: safe,lsx,lsy,lsz,lex,ley,lez,lwrap, &
      stay,safe1,check
    Integer           :: fail(1:3),iblock,jdnode,kdnode,        &
      imove,jmove,kmove,keep,                &
      i,j,k,l,ii,jj,kk,ll,                   &
      jxyz,ix,iy,iz,kx,ky,kz,                &
      newatm,iatm,jatm,katm,latm,matm,natm,  &
      jshels,kshels,jconst,kconst,jpmf,kpmf, &
      jrigid,krigid,jteths,kteths,           &
      jbonds,kbonds,jangle,kangle,           &
      jdihed,kdihed,jinver,kinver
    Real( Kind = wp ) :: uuu,vvv,www,xadd,yadd,zadd

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer
    Integer,           Dimension( : ), Allocatable :: lrgd
    Integer,           Dimension( : ), Allocatable :: ind_on,ind_off
    Integer,           Dimension( : ), Allocatable :: i1pmf,i2pmf

    Character( Len = 256 ) :: message
    fail=0
    Allocate (buffer(1:domain%mxbfdp),                   Stat=fail(1))
    Allocate (lrgd(-1:Max(rigid%max_list,rigid%max_rigid)),         Stat=fail(2))
    Allocate (ind_on(0:config%mxatms),ind_off(0:config%mxatms), Stat=fail(3))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'deport_atomic_data allocation failure 1'
      Call error(0,message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock=domain%mxbfdp/2

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! ls.  - wrap-around +1 in . direction (domain on the left MD cell border)
    ! le.  - wrap-around -1 in . direction (domain on the right MD cell border)
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0 ; ky = 0 ; kz = 0
    lsx = .false. ; lex = .false.
    lsy = .false. ; ley = .false.
    lsz = .false. ; lez = .false.
    If      (mdir == -1) Then ! Direction -x
      kx  = 1
      jxyz= 1
      lsx = (domain%idx == 0)

      jdnode = domain%map(1)
      kdnode = domain%map(2)
    Else If (mdir ==  1) Then ! Direction +x
      kx  = 1
      jxyz= 2
      lex = (domain%idx == domain%nx-1)

      jdnode = domain%map(2)
      kdnode = domain%map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky  = 1
      jxyz= 10
      lsy = (domain%idy == 0)

      jdnode = domain%map(3)
      kdnode = domain%map(4)
    Else If (mdir ==  2) Then ! Direction +y
      ky  = 1
      jxyz= 20
      ley = (domain%idy == domain%ny-1)

      jdnode = domain%map(4)
      kdnode = domain%map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz  = 1
      jxyz= 100
      lsz = (domain%idz == 0)

      jdnode = domain%map(5)
      kdnode = domain%map(6)
    Else If (mdir ==  3) Then ! Direction +z
      kz  = 1
      jxyz= 200
      lez = (domain%idz == domain%nz-1)

      jdnode = domain%map(6)
      kdnode = domain%map(5)
    Else
      Call error(42)
    End If

    ! Calculate PBC shift vector due to possible wrap around

    uuu=0.0_wp ; If (lsx) uuu=+1.0_wp ; If (lex) uuu=-1.0_wp
    vvv=0.0_wp ; If (lsy) vvv=+1.0_wp ; If (ley) vvv=-1.0_wp
    www=0.0_wp ; If (lsz) www=+1.0_wp ; If (lez) www=-1.0_wp

    lwrap = (Abs(uuu)+Abs(vvv)+Abs(www) > 0.5_wp)

    If (lwrap) Then
      xadd = config%cell(1)*uuu+config%cell(4)*vvv+config%cell(7)*www
      yadd = config%cell(2)*uuu+config%cell(5)*vvv+config%cell(8)*www
      zadd = config%cell(3)*uuu+config%cell(6)*vvv+config%cell(9)*www
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! buffer(1) and buffer(iblock+1) contain the actual number of
    ! particles to get transferred, imove and jmove are the lengths of
    ! the buffers

    imove=1
    jmove=1

    ! Initialise numbers of staying on and leaving off particles

    ind_on(0)=0 ; ind_off(0)=0

    ! Initialise array overflow flags

    safe=.true.
    safe1=.true.

    ! Initialise RB data set for preventing duplication

    lrgd(-1:rigid%n_types)=0

    ! LOOP OVER ALL PARTICLES ON THIS NODE

    Do i=1,config%natms
      stay=.false. ! the particle is assumed to be leaving

      ! If the particle is no longer scheduled to leave
      ! this domain in any direction

      If (config%ixyz(i) == 0) Then
        stay=.true.
      Else ! If (config%ixyz(i) > 0) Then ! Get the necessary halo indices
        ix=Mod(config%ixyz(i),10)           ! [0,1,2]
        iy=Mod(config%ixyz(i)-ix,100)       ! [0,10,20]
        iz=Mod(config%ixyz(i)-(ix+iy),1000) ! [0,100,200]

        ! Filter the move index for the selected direction

        j=ix*kx+iy*ky+iz*kz

        ! If the particle is scheduled to leave in the selected
        ! direction then reduce its move index, otherwise tag it as staying

        If (j == jxyz) Then
          config%ixyz(i)=config%ixyz(i)-jxyz
        Else
          stay=.true.
        End If
      End If

      If (stay) Then ! staying, keep it

        ii = ind_on(0)+1
        ind_on(0)  = ii
        ind_on(ii) = i

      Else ! not staying = leaving, pack all

        ii = ind_off(0)+1
        ind_off(0)  = ii
        ind_off(ii) = i

        ! If safe to proceed

        If (imove+18 <= iblock) Then

          ! pack positions and apply possible PBC shift for the receiver

          If (.not.lwrap) Then
            buffer(imove+1)=config%parts(i)%xxx
            buffer(imove+2)=config%parts(i)%yyy
            buffer(imove+3)=config%parts(i)%zzz
          Else
            buffer(imove+1)=config%parts(i)%xxx+xadd
            buffer(imove+2)=config%parts(i)%yyy+yadd
            buffer(imove+3)=config%parts(i)%zzz+zadd
          End If

          ! pack velocities

          buffer(imove+4)=config%vxx(i)
          buffer(imove+5)=config%vyy(i)
          buffer(imove+6)=config%vzz(i)

          ! pack forces

          buffer(imove+7)=config%parts(i)%fxx
          buffer(imove+8)=config%parts(i)%fyy
          buffer(imove+9)=config%parts(i)%fzz

          ! pack config indexing, site and move indexing arrays

          buffer(imove+10)=Real(config%ltg(i),wp)
          buffer(imove+11)=Real(config%lsite(i),wp)
          buffer(imove+12)=Real(config%ixyz(i),wp)

          ! pack initial positions

          buffer(imove+13)=stats%xin(i)
          buffer(imove+14)=stats%yin(i)
          buffer(imove+15)=stats%zin(i)

          ! pack final displacements

          buffer(imove+16)=stats%xto(i)
          buffer(imove+17)=stats%yto(i)
          buffer(imove+18)=stats%zto(i)

        Else

          safe=.false.

        End If
        imove=imove+18

        ! pack Langevin forces arrays

        If (thermo%l_langevin) Then
          If (imove+3 <= iblock) Then
            buffer(imove+1)=thermo%fxl(i)
            buffer(imove+2)=thermo%fyl(i)
            buffer(imove+3)=thermo%fzl(i)
          Else
            safe=.false.
          End If
          imove=imove+3
        End If

        ! pack minimisation arrays

        If (minim%transport) Then
          If (imove+3 <= iblock) Then
            buffer(imove+1)=minim%oxx(i)
            buffer(imove+2)=minim%oyy(i)
            buffer(imove+3)=minim%ozz(i)
          Else
            safe=.false.
          End If
          imove=imove+3
        End If

        ! pack frozen-frozen k-space SPME forces arrays

        If (ewld%lf_cp) Then
          If (imove+3 <= iblock) Then
            buffer(imove+1)=ewld%ffx(i)
            buffer(imove+2)=ewld%ffy(i)
            buffer(imove+3)=ewld%ffz(i)
          Else
            safe=.false.
          End If
          imove=imove+3
        End If

        ! pack k-space SPME forces arrays

        If (ewld%l_cp) Then
          If (imove+3 <= iblock) Then
            buffer(imove+1)=ewld%fcx(i)
            buffer(imove+2)=ewld%fcy(i)
            buffer(imove+3)=ewld%fcz(i)
          Else
            safe=.false.
          End If
          imove=imove+3
        End If

        ! pack initial velocities for VAF calculations

        If (green%samp > 0) Then
          If (imove+3 <= iblock) Then
            Do k=1,green%samp
              buffer(imove+1)=green%vxi(i,k)
              buffer(imove+2)=green%vyi(i,k)
              buffer(imove+3)=green%vzi(i,k)

              imove=imove+3
            End Do
          Else
            safe=.false.
          End If
        End If

        ! pack MSD arrays

        If (lmsd) Then
          If (imove+2*(6+stats%mxstak) <= iblock) Then
            jj=27+2*i
            buffer(imove+ 1)=stats%stpvl0(jj-1)
            buffer(imove+ 2)=stats%stpvl0(jj  )
            buffer(imove+ 3)=stats%stpval(jj-1)
            buffer(imove+ 4)=stats%stpval(jj  )
            buffer(imove+ 5)=stats%zumval(jj-1)
            buffer(imove+ 6)=stats%zumval(jj  )
            buffer(imove+ 7)=stats%ravval(jj-1)
            buffer(imove+ 8)=stats%ravval(jj  )
            buffer(imove+ 9)=stats%ssqval(jj-1)
            buffer(imove+10)=stats%ssqval(jj  )
            buffer(imove+11)=stats%sumval(jj-1)
            buffer(imove+12)=stats%sumval(jj  )
            Do kk=1,stats%mxstak
              l=2*kk   +12
              buffer(imove+l-1)=stats%stkval(kk,jj-1)
              buffer(imove+l  )=stats%stkval(kk,jj  )
            End Do
          Else
            safe=.false.
          End If
          imove=imove+2*(6+stats%mxstak)
        End If

        ! If intra-molecular entities exist in the system

        If (lbook) Then

          If (mpoles%max_mpoles > 0) Then ! pack topological array
            kk=mpoles%ltp(0,i)
            If (imove+1 <= iblock) Then
              imove=imove+1
              buffer(imove)=Real(kk,wp)
            Else
              imove=imove+1
              safe=.false.
            End If
            If (imove+kk <= iblock) Then
              Do k=1,kk
                imove=imove+1
                buffer(imove)=Real(mpoles%ltp(k,i),wp)
              End Do
            Else
              imove=imove+kk
              safe=.false.
            End If

            If (mpoles%key == POLARISATION_CHARMM) Then ! pack CHARMMing core-shell interactions array
              kk=mpoles%charmm(0,i)
              If (imove+1 <= iblock) Then
                imove=imove+1
                buffer(imove)=Real(kk,wp)
              Else
                imove=imove+1
                safe=.false.
              End If
              If (imove+kk <= iblock) Then
                Do k=1,kk
                  imove=imove+1
                  buffer(imove)=Real(mpoles%charmm(k,i),wp)
                End Do
              Else
                imove=imove+kk
                safe=.false.
              End If
            End If

          End If

          ! pack the exclusion list

          kk=neigh%list_excl(0,i)
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=Real(kk,wp)
          Else
            imove=imove+1
            safe=.false.
          End If
          If (imove+kk <= iblock) Then
            Do k=1,kk
              imove=imove+1
              buffer(imove)=Real(neigh%list_excl(k,i),wp)
            End Do
          Else
            imove=imove+kk
            safe=.false.
          End If

          ! the order of packing up intra bookkeeping arrays must be the same as
          ! their scanning in build_book_intra

          ! pack core-shell details

          jj=cshell%legshl(0,i) ; ii=Sign(1,jj) ; jj=Abs(jj)
          If (jj > 0) Then
            Do ll=1,jj
              If (imove+3 <= iblock) Then
                kk=cshell%legshl(ll,i)

                imove=imove+1
                buffer(imove)=Real(ii*cshell%listshl(0,kk),wp) ! negative for a shell particle

                Do k=1,2
                  imove=imove+1
                  buffer(imove)=Real(cshell%listshl(k,kk),wp)
                End Do
              Else
                imove=imove+3
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack constraint details

          jj=cons%legcon(0,i)
          If (jj > 0) Then
            Do ll=1,jj
              If (imove+3 <= iblock) Then
                kk=cons%legcon(ll,i)

                Do k=0,2
                  imove=imove+1
                  buffer(imove)=Real(cons%listcon(k,kk),wp)
                End Do
              Else
                imove=imove+3
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack PMF details

          jj=pmf%legpmf(0,i)
          If (jj > 0) Then
            Do ll=1,jj
              If (imove+pmf%mxtpmf(1)+pmf%mxtpmf(2)+2 <= iblock) Then
                kk=pmf%legpmf(ll,i)

                Do k=0,pmf%mxtpmf(1)
                  imove=imove+1
                  buffer(imove)=Real(pmf%listpmf(k,1,kk),wp)
                End Do

                Do k=0,pmf%mxtpmf(2)
                  imove=imove+1
                  buffer(imove)=Real(pmf%listpmf(k,2,kk),wp)
                End Do
              Else
                imove=imove+pmf%mxtpmf(1)+pmf%mxtpmf(2)+2
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack RB details

          !           jj=rigid%legend(0,i)
          !           If (jj > 0) Then
          !              l=12
          !              Do ll=1,jj
          !                 If (imove+rigid%max_list+l <= iblock) Then
          !                    kk=rigid%legend(ll,i)
          !
          !                    Do k=-1,rigid%list(-1,kk)
          !                       imove=imove+1
          !                       buffer(imove)=Real(rigid%list(k,kk),wp)
          !                    End Do
          !
          !                    buffer(imove+1)=rigid%q0(kk)
          !                    buffer(imove+2)=rigid%q1(kk)
          !                    buffer(imove+3)=rigid%q2(kk)
          !                    buffer(imove+4)=rigid%q3(kk)
          !                    imove=imove+4
          !
          !                    buffer(imove+1)=rigid%vxx(kk)
          !                    buffer(imove+2)=rigid%vyy(kk)
          !                    buffer(imove+3)=rigid%vzz(kk)
          !                    imove=imove+3
          !
          !                    buffer(imove+1)=rigid%oxx(kk)
          !                    buffer(imove+2)=rigid%oyy(kk)
          !                    buffer(imove+3)=rigid%ozz(kk)
          !                    imove=imove+3
          !                 Else
          !                    imove=imove+rigid%max_list+l
          !                    safe=.false.
          !                 End If
          !              End Do
          !           End If
          !           If (imove+1 <= iblock) Then
          !              imove=imove+1
          !              buffer(imove)=0.0_wp
          !           Else
          !              imove=imove+1
          !              safe=.false.
          !           End If

          jj=rigid%legend(0,i)
          If (jj > 0) Then
            l=12
            Do ll=1,jj
              If (imove+rigid%max_list+l <= iblock) Then
                kk=rigid%legend(ll,i)

                Do k=-1,1
                  imove=imove+1
                  buffer(imove)=Real(rigid%list(k,kk),wp)
                End Do

                ! Bypass if the data has already been sent

                If (lrgd(kk) == 0) Then
                  Do k=2,rigid%list(-1,kk)
                    imove=imove+1
                    buffer(imove)=Real(rigid%list(k,kk),wp)
                  End Do

                  buffer(imove+1)=rigid%q0(kk)
                  buffer(imove+2)=rigid%q1(kk)
                  buffer(imove+3)=rigid%q2(kk)
                  buffer(imove+4)=rigid%q3(kk)
                  imove=imove+4

                  buffer(imove+1)=rigid%vxx(kk)
                  buffer(imove+2)=rigid%vyy(kk)
                  buffer(imove+3)=rigid%vzz(kk)
                  imove=imove+3

                  buffer(imove+1)=rigid%oxx(kk)
                  buffer(imove+2)=rigid%oyy(kk)
                  buffer(imove+3)=rigid%ozz(kk)
                  imove=imove+3

                  lrgd(kk)=1
                End If
              Else
                imove=imove+rigid%max_list+l
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack tether details

          jj=tether%legtet(0,i)
          If (jj > 0) Then
            Do ll=1,jj
              If (imove+2 <= iblock) Then
                kk=tether%legtet(ll,i)

                Do k=0,1
                  imove=imove+1
                  buffer(imove)=Real(tether%listtet(k,kk),wp)
                End Do
              Else
                imove=imove+2
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack bond details

          jj=bond%legend(0,i)
          If (jj > 0) Then
            Do ll=1,jj
              If (imove+3 <= iblock) Then
                kk=bond%legend(ll,i)

                Do k=0,2
                  imove=imove+1
                  buffer(imove)=Real(bond%list(k,kk),wp)
                End Do
              Else
                imove=imove+3
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack valence angle details

          jj=angle%legend(0,i)
          If (jj > 0) Then
            Do ll=1,jj
              If (imove+4 <= iblock) Then
                kk=angle%legend(ll,i)

                Do k=0,3
                  imove=imove+1
                  buffer(imove)=Real(angle%list(k,kk),wp)
                End Do
              Else
                imove=imove+4
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack dihedral angle details

          jj=dihedral%legend(0,i)
          If (jj > 0) Then
            If (.not.dihedral%l_core_shell) Then ! dihedrals only have 4 members
              l=4
            Else                  ! dihedrals have 4+2 tracked members
              l=6
            End If
            Do ll=1,jj
              If (imove+l+1 <= iblock) Then
                kk=dihedral%legend(ll,i)

                Do k=0,l
                  imove=imove+1
                  buffer(imove)=Real(dihedral%list(k,kk),wp)
                End Do
              Else
                imove=imove+l+1
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

          ! pack inversion angle details

          jj=inversion%legend(0,i)
          If (jj > 0) Then
            Do ll=1,jj
              If (imove+5 <= iblock) Then
                kk=inversion%legend(ll,i)

                Do k=0,4
                  imove=imove+1
                  buffer(imove)=Real(inversion%list(k,kk),wp)
                End Do
              Else
                imove=imove+5
                safe=.false.
              End If
            End Do
          End If
          If (imove+1 <= iblock) Then
            imove=imove+1
            buffer(imove)=0.0_wp
          Else
            imove=imove+1
            safe=.false.
          End If

        End If

      End If
    End Do

    ! Check for array bound overflow (have arrays coped with outgoing data)

    Call gcheck(comm,safe)
    If (.not.safe) Call error(43)

    ! Restack arrays for leaving off particles with staying on ones
    ! Thanks to Victor Gamayunov

    k=ind_on(0)
    l=ind_off(0)
    Do ii=1,l
      keep=ind_off(ii) ; If (keep > ind_on(k)) Exit ! Thanks to Alin Elena
      i   =ind_on(k-ii+1)

      config%parts(keep)%xxx=config%parts(i)%xxx
      config%parts(keep)%yyy=config%parts(i)%yyy
      config%parts(keep)%zzz=config%parts(i)%zzz

      config%vxx(keep)=config%vxx(i)
      config%vyy(keep)=config%vyy(i)
      config%vzz(keep)=config%vzz(i)

      config%parts(keep)%fxx=config%parts(i)%fxx
      config%parts(keep)%fyy=config%parts(i)%fyy
      config%parts(keep)%fzz=config%parts(i)%fzz

      config%ltg(keep)=config%ltg(i)
      config%lsite(keep)=config%lsite(i)
      config%ixyz(keep)=config%ixyz(i)

      stats%xin(keep)=stats%xin(i)
      stats%yin(keep)=stats%yin(i)
      stats%zin(keep)=stats%zin(i)

      stats%xto(keep)=stats%xto(i)
      stats%yto(keep)=stats%yto(i)
      stats%zto(keep)=stats%zto(i)

      If (thermo%l_langevin) Then
        thermo%fxl(keep)=thermo%fxl(i)
        thermo%fyl(keep)=thermo%fyl(i)
        thermo%fzl(keep)=thermo%fzl(i)
      End If

      If (minim%transport) Then
        minim%oxx(keep)=minim%oxx(i)
        minim%oyy(keep)=minim%oyy(i)
        minim%ozz(keep)=minim%ozz(i)
      End If

      If (ewld%lf_cp) Then
        ewld%ffx(keep)=ewld%ffx(i)
        ewld%ffy(keep)=ewld%ffy(i)
        ewld%ffz(keep)=ewld%ffz(i)
      End If

      If (ewld%l_cp) Then
        ewld%fcx(keep)=ewld%fcx(i)
        ewld%fcy(keep)=ewld%fcy(i)
        ewld%fcz(keep)=ewld%fcz(i)
      End If

      If (green%samp > 0) Then
        green%vxi(keep,1:green%samp)=green%vxi(i,1:green%samp)
        green%vyi(keep,1:green%samp)=green%vyi(i,1:green%samp)
        green%vzi(keep,1:green%samp)=green%vzi(i,1:green%samp)
      End If

      If (lmsd) Then
        jj=27+2*i
        j =27+2*keep
        stats%stpvl0(j-1)=stats%stpvl0(jj-1)
        stats%stpvl0(j  )=stats%stpvl0(jj  )
        stats%stpval(j-1)=stats%stpval(jj-1)
        stats%stpval(j  )=stats%stpval(jj  )
        stats%zumval(j-1)=stats%zumval(jj-1)
        stats%zumval(j  )=stats%zumval(jj  )
        stats%ravval(j-1)=stats%ravval(jj-1)
        stats%ravval(j  )=stats%ravval(jj  )
        stats%ssqval(j-1)=stats%ssqval(jj-1)
        stats%ssqval(j  )=stats%ssqval(jj  )
        stats%sumval(j-1)=stats%sumval(jj-1)
        stats%sumval(j  )=stats%sumval(jj  )
        Do kk=1,stats%mxstak
          stats%stkval(kk,j-1)=stats%stkval(kk,jj-1)
          stats%stkval(kk,j  )=stats%stkval(kk,jj  )
        End Do
      End If

      If (lbook) Then
        If (mpoles%max_mpoles > 0) Then
          mpoles%ltp(:,keep)=mpoles%ltp(:,i)
          If (mpoles%key == POLARISATION_CHARMM) mpoles%charmm(:,keep)=mpoles%charmm(:,i)
        End If

        neigh%list_excl(:,keep)=neigh%list_excl(:,i)

        cshell%legshl(:,keep)=cshell%legshl(:,i)

        cons%legcon(:,keep)=cons%legcon(:,i)
        pmf%legpmf(:,keep)=pmf%legpmf(:,i)

        rigid%legend(:,keep)=rigid%legend(:,i)

        tether%legtet(:,keep)=tether%legtet(:,i)

        bond%legend(:,keep)=bond%legend(:,i)
        angle%legend(:,keep)=angle%legend(:,i)
        dihedral%legend(:,keep)=dihedral%legend(:,i)
        inversion%legend(:,keep)=inversion%legend(:,i)
      End If
    End Do
    keep=k ! How many particles are to be kept

    ! record of number of atoms for transfer

    buffer(1)=Real(config%natms-keep,wp)

    ! exchange information on buffer sizes

    Call girecv(comm,jmove,kdnode,Deport_tag)
    Call gsend(comm,imove,jdnode,Deport_tag)
    Call gwait(comm)

    ! exchange buffers between nodes (this is a MUST)

    If (jmove > 0) Then
      Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,Deport_tag)
    End If
    If (imove > 0) Then
      Call gsend(comm,buffer(1:imove),jdnode,Deport_tag)
    End If
    If (jmove > 0) Call gwait(comm)

    ! check arrays can cope with incoming atom numbers

    kmove=iblock+1
    jmove=Nint(buffer(kmove))

    config%natms=keep+jmove

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe=(config%natms <= config%mxatms)
    Call gcheck(comm,safe)
    If (.not.safe) Call error(44)

    Deallocate (ind_on,ind_off,                        Stat=fail(1))
    Allocate   (i1pmf(1:pmf%mxtpmf(1)),i2pmf(1:pmf%mxtpmf(2)), Stat=fail(2))
    If (Any(fail(1:2) > 0)) Then
      Write(message,'(a)') 'deport_atomic_data de/allocation failure'
      Call error(0,message)
    End If

    ! load transferred data

    Do i=1,jmove
      newatm=i+keep

      ! unpack positions

      config%parts(newatm)%xxx=buffer(kmove+1)
      config%parts(newatm)%yyy=buffer(kmove+2)
      config%parts(newatm)%zzz=buffer(kmove+3)

      ! unpack velocities

      config%vxx(newatm)=buffer(kmove+4)
      config%vyy(newatm)=buffer(kmove+5)
      config%vzz(newatm)=buffer(kmove+6)

      ! unpack forces

      config%parts(newatm)%fxx=buffer(kmove+7)
      config%parts(newatm)%fyy=buffer(kmove+8)
      config%parts(newatm)%fzz=buffer(kmove+9)

      ! unpack config indexing, site and move indexing arrays

      config%ltg(newatm)=Nint(buffer(kmove+10))
      config%lsite(newatm)=Nint(buffer(kmove+11))
      config%ixyz(newatm)=Nint(buffer(kmove+12))

      ! unpack initial positions arrays

      stats%xin(newatm)=buffer(kmove+13)
      stats%yin(newatm)=buffer(kmove+14)
      stats%zin(newatm)=buffer(kmove+15)

      ! unpack initial positions arrays

      stats%xto(newatm)=buffer(kmove+16)
      stats%yto(newatm)=buffer(kmove+17)
      stats%zto(newatm)=buffer(kmove+18)

      kmove=kmove+18

      ! unpack Langevin forces arrays

      If (thermo%l_langevin) Then
        thermo%fxl(newatm)=buffer(kmove+1)
        thermo%fyl(newatm)=buffer(kmove+2)
        thermo%fzl(newatm)=buffer(kmove+3)

        kmove=kmove+3
      End If

      ! unpack minimisation arrays

      If (minim%transport) Then
        minim%oxx(newatm)=buffer(kmove+1)
        minim%oyy(newatm)=buffer(kmove+2)
        minim%ozz(newatm)=buffer(kmove+3)

        kmove=kmove+3
      End If

      ! unpack frozen-frozen k-space SPME forces arrays

      If (ewld%lf_cp) Then
        ewld%ffx(newatm)=buffer(kmove+1)
        ewld%ffy(newatm)=buffer(kmove+2)
        ewld%ffz(newatm)=buffer(kmove+3)

        kmove=kmove+3
      End If

      ! unpack k-space SPME forces arrays

      If (ewld%l_cp) Then
        ewld%fcx(newatm)=buffer(kmove+1)
        ewld%fcy(newatm)=buffer(kmove+2)
        ewld%fcz(newatm)=buffer(kmove+3)

        kmove=kmove+3
      End If

      ! unpack initial velocities for VAF calculations

      If (green%samp > 0) Then
        Do k=1,green%samp
          green%vxi(newatm,k)=buffer(kmove+1)
          green%vyi(newatm,k)=buffer(kmove+2)
          green%vzi(newatm,k)=buffer(kmove+3)

          kmove=kmove+3
        End Do
      End If

      ! unpack MSD arrays

      If (lmsd) Then
        jj=27+2*newatm
        stats%stpvl0(jj-1)=buffer(kmove+1 )
        stats%stpvl0(jj  )=buffer(kmove+2 )
        stats%stpval(jj-1)=buffer(kmove+3 )
        stats%stpval(jj  )=buffer(kmove+4 )
        stats%zumval(jj-1)=buffer(kmove+5 )
        stats%zumval(jj  )=buffer(kmove+6 )
        stats%ravval(jj-1)=buffer(kmove+7 )
        stats%ravval(jj  )=buffer(kmove+8 )
        stats%ssqval(jj-1)=buffer(kmove+9 )
        stats%ssqval(jj  )=buffer(kmove+10)
        stats%sumval(jj-1)=buffer(kmove+11)
        stats%sumval(jj  )=buffer(kmove+12)
        Do kk=1,stats%mxstak
          l=2*kk                +12
          stats%stkval(kk,jj-1)=buffer(kmove+l-1)
          stats%stkval(kk,jj  )=buffer(kmove+l  )
        End Do

        kmove=kmove+2*(6+stats%mxstak)
      End If

      If (lbook) Then

        If (mpoles%max_mpoles > 0) Then ! unpack topological array
          kmove=kmove+1
          kk=Nint(buffer(kmove))
          mpoles%ltp(0,newatm)=kk
          Do k=1,kk
            kmove=kmove+1
            mpoles%ltp(k,newatm)=Nint(buffer(kmove))
          End Do
          mpoles%ltp(kk+1:neigh%max_exclude,newatm)=0

          If (mpoles%key == POLARISATION_CHARMM) Then ! unpack CHARMMing core-shell interactions array
            kmove=kmove+1
            kk=Nint(buffer(kmove))
            mpoles%charmm(0,newatm)=kk
            Do k=1,kk
              kmove=kmove+1
              mpoles%charmm(k,newatm)=Nint(buffer(kmove))
            End Do
            mpoles%charmm(kk+1:neigh%max_exclude,newatm)=0
          End If
        End If

        ! unpack the exclusion neigh%list

        kmove=kmove+1
        kk=Nint(buffer(kmove))
        neigh%list_excl(0,newatm)=kk
        Do k=1,kk
          kmove=kmove+1
          neigh%list_excl(k,newatm)=Nint(buffer(kmove))
        End Do
        neigh%list_excl(kk+1:neigh%max_exclude,newatm)=0

        ! the order of unpacking intra bookkeeping arrays must be the same as
        ! the order of their scanning in build_book_intra in order to rebuild
        ! correctly the new neigh%list arrays and create new legend arrays

        ! set initial intra counters

        jshels=cshell%ntshl

        jconst=cons%ntcons
        jpmf  =pmf%ntpmf

        jrigid=rigid%n_types

        jteths=tether%ntteth

        jbonds=bond%n_types
        jangle=angle%n_types
        jdihed=dihedral%n_types
        jinver=inversion%n_types

        ! unpack core-shell details

        cshell%legshl(:,newatm) = 0
        Do While (Abs(buffer(kmove+1)) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1)) ; ll=Sign(1,jj) ; jj=Abs(jj)
          iatm=Nint(buffer(kmove+2)) ! ll=1
          jatm=Nint(buffer(kmove+3)) ! ll=-1
          kmove=kmove+3

          ! check if core-shell unit already specified

          kshels=0
          check=.true.
          Do While (check .and. kshels < Min(jshels,cshell%mxshl))
            kshels=kshels+1
            check=.not.( jj   == cshell%listshl(0,kshels) .and. & ! core-shell units don't intersect
              iatm == cshell%listshl(1,kshels))        ! .and. jatm == cshell%listshl(2,kshels) )
          End Do

          ! add new core-shell unit

          If (check) Then
            jshels=jshels+1

            If (jshels <= cshell%mxshl) Then
              cshell%listshl(0,jshels)=jj
              cshell%listshl(1,jshels)=iatm
              cshell%listshl(2,jshels)=jatm

              Call tag_legend(safe1,newatm,ll*jshels,cshell%legshl,cshell%mxfshl)
            Else
              safe=.false.
              Write(message,'(a)') "too many core-shell units"
              Call warning(message)

            End If
          Else
            Call tag_legend(safe1,newatm,ll*kshels,cshell%legshl,cshell%mxfshl)
          End If
        End Do
        kmove=kmove+1

        ! unpack bond constraint details

        cons%legcon(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1))
          iatm=Nint(buffer(kmove+2))
          jatm=Nint(buffer(kmove+3))
          kmove=kmove+3

          ! check if constraint already specified

          kconst=0
          check=.true.
          Do While (check .and. kconst < Min(jconst,cons%mxcons))
            kconst=kconst+1
            check=.not.( jj   == cons%listcon(0,kconst) .and. &
              iatm == cons%listcon(1,kconst) .and. &
              jatm == cons%listcon(2,kconst) )
          End Do

          ! insert new constraint unit

          If (check) Then
            jconst=jconst+1

            If (jconst <= cons%mxcons) Then
              cons%listcon(0,jconst)=jj
              cons%listcon(1,jconst)=iatm
              cons%listcon(2,jconst)=jatm

              Call tag_legend(safe1,newatm,jconst,cons%legcon,cons%mxfcon)
            Else
              safe=.false.
              Call warning('too many constraint units')
            End If
          Else
            Call tag_legend(safe1,newatm,kconst,cons%legcon,cons%mxfcon)
          End If
        End Do
        kmove=kmove+1

        ! unpack PMF details

        pmf%legpmf(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1)) ! PMF global identifier
          kmove=kmove+1
          Do k=1,pmf%mxtpmf(1)
            kmove=kmove+1
            i1pmf(k)=Nint(buffer(kmove))
          End Do

          kmove=kmove+1            ! omit PMF units presence identifier
          Do k=1,pmf%mxtpmf(2)         ! and deal with it in pmf_units_set
            kmove=kmove+1
            i2pmf(k)=Nint(buffer(kmove))
          End Do

          ! check if PMF already specified

          kpmf=0
          check=.true.
          Do While (check .and. kpmf < Min(jpmf,pmf%mxpmf))
            kpmf=kpmf+1
            check=.not.(jj == pmf%listpmf(0,1,kpmf))
          End Do

          ! insert new PMF constraint

          If (check) Then
            jpmf=jpmf+1

            If (jpmf <= pmf%mxpmf) Then
              pmf%listpmf(0,1,jpmf)=jj
              Do k=1,pmf%mxtpmf(1)
                pmf%listpmf(k,1,jpmf)=i1pmf(k)
              End Do

              ! PMF units presence identifier holds zero temporarily
              ! it's dealt with in pmf_units_set

              pmf%listpmf(0,2,jpmf)=0
              Do k=1,pmf%mxtpmf(2)
                pmf%listpmf(k,2,jpmf)=i2pmf(k)
              End Do

              Call tag_legend(safe1,newatm,jpmf,pmf%legpmf,pmf%mxfpmf)
            Else
              safe=.false.
              Call warning('too many PMF units')
            End If
          Else
            Call tag_legend(safe1,newatm,kpmf,pmf%legpmf,pmf%mxfpmf)
          End If
        End Do
        kmove=kmove+1

        ! unpack RB details

        ! Initialise RB data set for preventing duplication

        lrgd(-1:rigid%max_list)=0

        !        rigid%legend(:,newatm) = 0
        !        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
        !           lrgd(-1)=Nint(buffer(kmove+1))
        !           kmove=kmove+1
        !           Do k=0,lrgd(-1)
        !              kmove=kmove+1
        !              lrgd(k)=Nint(buffer(kmove))
        !           End Do
        !
        !! check if RB already specified
        !
        !           krigid=0
        !           check=.true.
        !           Do While (check .and. krigid < Min(jrigid,rigid%max_rigid))
        !              krigid=krigid+1
        !              check=( lrgd( 0) == rigid%list( 0,krigid) .and. & ! Type
        !                      lrgd(-1) == rigid%list(-1,krigid) .and. & ! Size
        !                      lrgd( 1) == rigid%list( 1,krigid) )       ! Global ID of just the first
        !!                                                            ! member as RBs don't intersect
        !              check=.not.check
        !           End Do
        !
        !! insert new RB unit
        !
        !           If (check) Then
        !              jrigid=jrigid+1
        !
        !              If (jrigid <= rigid%max_rigid) Then
        !                 Do k=-1,lrgd(-1)
        !                    rigid%list(k,jrigid)=lrgd(k)
        !                 End Do
        !
        !                 Call tag_legend(safe1,newatm,jrigid,rigid%legend,rigid%max_frozen)
        !
        !                 rigid%q0(jrigid)=buffer(kmove+1)
        !                 rigid%q1(jrigid)=buffer(kmove+2)
        !                 rigid%q2(jrigid)=buffer(kmove+3)
        !                 rigid%q3(jrigid)=buffer(kmove+4)
        !                 kmove=kmove+4
        !
        !                 rigid%vxx(jrigid)=buffer(kmove+1)
        !                 rigid%vyy(jrigid)=buffer(kmove+2)
        !                 rigid%vzz(jrigid)=buffer(kmove+3)
        !                 kmove=kmove+3
        !
        !                 rigid%oxx(jrigid)=buffer(kmove+1)
        !                 rigid%oyy(jrigid)=buffer(kmove+2)
        !                 rigid%ozz(jrigid)=buffer(kmove+3)
        !                 kmove=kmove+3
        !              Else
        !                 safe=.false.
        !              End If
        !           Else
        !              l=10
        !              kmove=kmove+l ! Compensate for the 'l' extra unread buffers
        !              Call tag_legend(safe1,newatm,krigid,rigid%legend,rigid%max_frozen)
        !           End If
        !        End Do
        !        kmove=kmove+1

        rigid%legend(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          lrgd(-1)=Nint(buffer(kmove+1))
          kmove=kmove+1
          Do k=0,2
            kmove=kmove+1
            lrgd(k)=Nint(buffer(kmove))
          End Do

          ! check if RB already specified

          krigid=0
          check=.true.
          Do While (check .and. krigid < Min(jrigid,rigid%max_rigid))
            krigid=krigid+1
            check=( lrgd( 0) == rigid%list( 0,krigid) .and. & ! Type
              lrgd(-1) == rigid%list(-1,krigid) .and. & ! Size
              lrgd( 1) == rigid%list( 1,krigid) )       ! Global ID of just the first
            !                                                            ! member as RBs don't intersect
            check=.not.check
          End Do

          ! insert new RB unit

          If (check) Then
            Do k=3,lrgd(-1)
              kmove=kmove+1
              lrgd(k)=Nint(buffer(kmove))
            End Do

            jrigid=jrigid+1

            If (jrigid <= rigid%max_rigid) Then
              Do k=-1,lrgd(-1)
                rigid%list(k,jrigid)=lrgd(k)
              End Do

              Call tag_legend(safe1,newatm,jrigid,rigid%legend,rigid%max_frozen)

              rigid%q0(jrigid)=buffer(kmove+1)
              rigid%q1(jrigid)=buffer(kmove+2)
              rigid%q2(jrigid)=buffer(kmove+3)
              rigid%q3(jrigid)=buffer(kmove+4)
              kmove=kmove+4

              rigid%vxx(jrigid)=buffer(kmove+1)
              rigid%vyy(jrigid)=buffer(kmove+2)
              rigid%vzz(jrigid)=buffer(kmove+3)
              kmove=kmove+3

              rigid%oxx(jrigid)=buffer(kmove+1)
              rigid%oyy(jrigid)=buffer(kmove+2)
              rigid%ozz(jrigid)=buffer(kmove+3)
              kmove=kmove+3
            Else
              safe=.false.
              Call warning('too many rigid body units on node')
            End If
          Else
            If (lrgd(2) == 0) Then  ! Unduplication: Details have already been sent -
              kmove=kmove-1        ! back up once for the zero read at end of data-pack
            Else                    ! Data already present - jump over it as it's not needed
              l=10                 ! Compensate for the 'l' extra unread buffers
              kmove=kmove+l+lrgd(-1)-2
            End If
            Call tag_legend(safe1,newatm,krigid,rigid%legend,rigid%max_frozen)
          End If
        End Do
        kmove=kmove+1

        ! unpack tether details

        tether%legtet(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1))
          iatm=Nint(buffer(kmove+2))
          kmove=kmove+2

          ! check if tether unit already specified

          kteths=0
          check=.true.
          Do While (check .and. kteths < Min(jteths,tether%mxteth))
            kteths=kteths+1
            check=.not.( jj   == tether%listtet(0,kteths) .and. &
              iatm == tether%listtet(1,kteths) )
          End Do

          ! add new tether unit

          If (check) Then
            jteths=jteths+1

            If (jteths <= tether%mxteth) Then
              tether%listtet(0,jteths)=jj
              tether%listtet(1,jteths)=iatm

              Call tag_legend(safe1,newatm,jteths,tether%legtet,tether%mxftet)
            Else
              safe=.false.
              Call warning('too many tether units')
            End If
          Else
            Call tag_legend(safe1,newatm,kteths,tether%legtet,tether%mxftet)
          End If
        End Do
        kmove=kmove+1

        ! unpack bond details

        bond%legend(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1))
          iatm=Nint(buffer(kmove+2))
          jatm=Nint(buffer(kmove+3))
          kmove=kmove+3

          ! check if bond already specified

          kbonds=0
          check=.true.
          Do While (check .and. kbonds < Min(jbonds,bond%max_bonds))
            kbonds=kbonds+1
            check=.not.( jj   == bond%list(0,kbonds) .and. &
              iatm == bond%list(1,kbonds) .and. &
              jatm == bond%list(2,kbonds) )
          End Do

          ! insert new bond details

          If (check) Then
            jbonds=jbonds+1

            If (jbonds <= bond%max_bonds) Then
              bond%list(0,jbonds)=jj
              bond%list(1,jbonds)=iatm
              bond%list(2,jbonds)=jatm

              Call tag_legend(safe1,newatm,jbonds,bond%legend,bond%max_legend)
            Else
              safe=.false.
              Call warning('too many bond units')
            End If
          Else
            Call tag_legend(safe1,newatm,kbonds,bond%legend,bond%max_legend)
          End If
        End Do
        kmove=kmove+1

        ! unpack valence angle details

        angle%legend(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1))
          iatm=Nint(buffer(kmove+2))
          jatm=Nint(buffer(kmove+3))
          katm=Nint(buffer(kmove+4))
          kmove=kmove+4

          ! check if angle already specified

          kangle=0
          check=.true.
          Do While (check .and. kangle < Min(jangle,angle%max_angles))
            kangle=kangle+1
            check=.not.( jj   == angle%list(0,kangle) .and. &
              iatm == angle%list(1,kangle) .and. &
              jatm == angle%list(2,kangle) .and. &
              katm == angle%list(3,kangle) )
          End Do

          ! insert new angle details

          If (check) Then
            jangle=jangle+1

            If (jangle <= angle%max_angles) Then
              angle%list(0,jangle)=jj
              angle%list(1,jangle)=iatm
              angle%list(2,jangle)=jatm
              angle%list(3,jangle)=katm

              Call tag_legend(safe1,newatm,jangle,angle%legend,angle%max_legend)
            Else
              safe=.false.
              Call warning('too many angle units')
            End If
          Else
            Call tag_legend(safe1,newatm,kangle,angle%legend,angle%max_legend)
          End If
        End Do
        kmove=kmove+1

        ! unpack dihedral angle details

        dihedral%legend(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1))
          iatm=Nint(buffer(kmove+2))
          jatm=Nint(buffer(kmove+3))
          katm=Nint(buffer(kmove+4))
          latm=Nint(buffer(kmove+5))
          If (dihedral%l_core_shell) Then
            matm=Nint(buffer(kmove+6))
            natm=Nint(buffer(kmove+7))
            kmove=kmove+7
          Else
            kmove=kmove+5
          End If

          ! check if dihedral already specified

          kdihed=0
          check=.true.
          Do While (check .and. kdihed < Min(jdihed,dihedral%max_angles))
            kdihed=kdihed+1
            check=.not.( jj   == dihedral%list(0,kdihed) .and. &
              iatm == dihedral%list(1,kdihed) .and. &
              jatm == dihedral%list(2,kdihed) .and. &
              katm == dihedral%list(3,kdihed) .and. &
              latm == dihedral%list(4,kdihed) )
          End Do

          ! add new dihedral details

          If (check) Then
            jdihed=jdihed+1

            If (jdihed <= dihedral%max_angles) Then
              dihedral%list(0,jdihed)=jj
              dihedral%list(1,jdihed)=iatm
              dihedral%list(2,jdihed)=jatm
              dihedral%list(3,jdihed)=katm
              dihedral%list(4,jdihed)=latm
              If (dihedral%l_core_shell) Then
                dihedral%list(5,jdihed)=matm
                dihedral%list(6,jdihed)=natm
              End If

              Call tag_legend(safe1,newatm,jdihed,dihedral%legend,dihedral%max_legend)
            Else
              safe=.false.
              Call warning('too many dihedral units')
            End If
          Else
            Call tag_legend(safe1,newatm,kdihed,dihedral%legend,dihedral%max_legend)
          End If
        End Do
        kmove=kmove+1

        ! unpack inversion angle details

        inversion%legend(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
          jj=Nint(buffer(kmove+1))
          iatm=Nint(buffer(kmove+2))
          jatm=Nint(buffer(kmove+3))
          katm=Nint(buffer(kmove+4))
          latm=Nint(buffer(kmove+5))
          kmove=kmove+5

          ! check if inversion already specified

          kinver=0
          check=.true.
          Do While (check .and. kinver < Min(jinver,inversion%max_angles))
            kinver=kinver+1
            check=.not.( jj   == inversion%list(0,kinver) .and. &
              iatm == inversion%list(1,kinver) .and. &
              jatm == inversion%list(2,kinver) .and. &
              katm == inversion%list(3,kinver) .and. &
              latm == inversion%list(4,kinver) )
          End Do

          ! add new inversion details

          If (check) Then
            jinver=jinver+1

            If (jinver <= inversion%max_angles) Then
              inversion%list(0,jinver)=jj
              inversion%list(1,jinver)=iatm
              inversion%list(2,jinver)=jatm
              inversion%list(3,jinver)=katm
              inversion%list(4,jinver)=latm

              Call tag_legend(safe1,newatm,jinver,inversion%legend,inversion%max_legend)
            Else
              safe=.false.
              Call warning('too many inversion units')
            End If
          Else
            Call tag_legend(safe1,newatm,kinver,inversion%legend,inversion%max_legend)
          End If
        End Do
        kmove=kmove+1

        ! redefine intra counters

        cshell%ntshl =jshels

        cons%ntcons=jconst
        pmf%ntpmf =jpmf

        rigid%n_types =jrigid

        tether%ntteth=jteths

        bond%n_types=jbonds
        angle%n_types=jangle
        dihedral%n_types=jdihed
        inversion%n_types =jinver

      End If
    End Do

    ! check error flags

    Call gcheck(comm,safe)
    If (.not.safe) Call error(113)
    Call gcheck(comm,safe1)
    If (.not.safe1) Call error(114)

    Deallocate (buffer,      Stat=fail(1))
    Deallocate (lrgd,        Stat=fail(2))
    Deallocate (i1pmf,i2pmf, Stat=fail(3))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'deport_atomic_data deallocation failure'
      Call error(0,message)
    End If

  End Subroutine deport_atomic_data

  Subroutine export_atomic_data(mdir,domain,config,kim_data,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export atomic data in domain boundary regions
    ! for halo formation
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov december 2016
    ! contrib   - i.j.bush february 2016
    ! contrib   - h.a.boateng february 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer,            Intent( In    ) :: mdir
    Type( domains_type ), Intent( In    ) :: domain
    Type( kim_type ), Intent( InOut ) :: kim_data
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type ), Intent( InOut ) :: comm

    Logical           :: safe,lsx,lsy,lsz,lex,ley,lez,lwrap
    Integer           :: fail,iadd,limit,iblock,            &
      i,j,k,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
      jdnode,kdnode,imove,jmove,kmove,itmp
    Real( Kind = wp ) :: uuu,vvv,www,xadd,yadd,zadd

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    Character ( Len = 256 )   ::  message

    ! Number of transported quantities per particle

    iadd=6
    fail=0
    limit=iadd*domain%mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
    Allocate (buffer(1:limit), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'export_atomic_data allocation failure'
      Call error(0,message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock=limit/Merge(2,1,comm%mxnode > 1)

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! kxyz - corrected halo reduction factor particles haloing both +&- sides
    ! ls.  - wrap-around +1 in . direction (domain on the left MD cell border)
    ! le.  - wrap-around -1 in . direction (domain on the right MD cell border)
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0 ; ky = 0 ; kz = 0
    lsx = .false. ; lex = .false.
    lsy = .false. ; ley = .false.
    lsz = .false. ; lez = .false.
    If      (mdir == -1) Then ! Direction -x
      kx  = 1
      jxyz= 1
      kxyz= 3
      lsx = (domain%idx == 0)

      jdnode = domain%map(1)
      kdnode = domain%map(2)
    Else If (mdir ==  1) Then ! Direction +x
      kx  = 1
      jxyz= 2
      kxyz= 3
      lex = (domain%idx == domain%nx-1)

      jdnode = domain%map(2)
      kdnode = domain%map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky  = 1
      jxyz= 10
      kxyz= 30
      lsy = (domain%idy == 0)

      jdnode = domain%map(3)
      kdnode = domain%map(4)
    Else If (mdir ==  2) Then ! Direction +y
      ky  = 1
      jxyz= 20
      kxyz= 30
      ley = (domain%idy == domain%ny-1)

      jdnode = domain%map(4)
      kdnode = domain%map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz  = 1
      jxyz= 100
      kxyz= 300
      lsz = (domain%idz == 0)

      jdnode = domain%map(5)
      kdnode = domain%map(6)
    Else If (mdir ==  3) Then ! Direction +z
      kz  = 1
      jxyz= 200
      kxyz= 300
      lez = (domain%idz == domain%nz-1)

      jdnode = domain%map(6)
      kdnode = domain%map(5)
    Else
      Call error(46)
    End If

    ! Calculate PBC shift vector due to possible wrap around

    uuu=0.0_wp ; If (lsx) uuu=+1.0_wp ; If (lex) uuu=-1.0_wp
    vvv=0.0_wp ; If (lsy) vvv=+1.0_wp ; If (ley) vvv=-1.0_wp
    www=0.0_wp ; If (lsz) www=+1.0_wp ; If (lez) www=-1.0_wp

    lwrap = (Abs(uuu)+Abs(vvv)+Abs(www) > 0.5_wp)

    If (lwrap) Then
      xadd = config%cell(1)*uuu+config%cell(4)*vvv+config%cell(7)*www
      yadd = config%cell(2)*uuu+config%cell(5)*vvv+config%cell(8)*www
      zadd = config%cell(3)*uuu+config%cell(6)*vvv+config%cell(9)*www
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed

    imove=0
    jmove=0

    ! Initialise array overflow flags

    safe=.true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE

    Do i=1,config%nlast

      ! If the particle is within the remaining 'inverted halo' of this domain

      If (config%ixyz(i) > 0) Then

        ! Get the necessary halo indices

        ix=Mod(config%ixyz(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(config%ixyz(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(config%ixyz(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

        ! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

        ! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

          ! If safe to proceed

          If ((imove+iadd) <= iblock) Then

            ! pack positions and apply possible PBC shift for the receiver

            If (.not.lwrap) Then
              buffer(imove+1)=config%parts(i)%xxx
              buffer(imove+2)=config%parts(i)%yyy
              buffer(imove+3)=config%parts(i)%zzz
            Else
              buffer(imove+1)=config%parts(i)%xxx+xadd
              buffer(imove+2)=config%parts(i)%yyy+yadd
              buffer(imove+3)=config%parts(i)%zzz+zadd
            End If

            ! pack config indexing, site and remaining halo indexing arrays

            buffer(imove+4)=Real(config%ltg(i),wp)
            buffer(imove+5)=Real(config%lsite(i),wp)

            ! Use the corrected halo reduction factor when the particle is halo to both +&- sides

            buffer(imove+iadd)=Real(config%ixyz(i)-Merge(jxyz,kxyz,j == jxyz),wp)

          Else

            safe=.false.

          End If
          imove=imove+iadd

        End If

      End If

    End Do

    ! Check for array bound overflow (have arrays coped with outgoing data)

    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=Merge(2,1,comm%mxnode > 1)*imove
      Call gmax(comm,itmp)
      Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
      Call error(54)
    End If

    ! exchange information on buffer sizes

    If (comm%mxnode > 1) Then
      Call girecv(comm,jmove,kdnode,Export_tag)
      Call gsend(comm,imove,jdnode,Export_tag)
      Call gwait(comm)
    Else
      jmove=imove
    End If

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe=((config%nlast+jmove/iadd) <= config%mxatms)
    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=config%nlast+jmove/iadd
      Call gmax(comm,itmp)
      Call warning(160,Real(itmp,wp),Real(config%mxatms,wp),0.0_wp)
      Call error(56)
    End If

    ! exchange buffers between nodes (this is a MUST)

    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,Export_tag)
      End If
      If (imove > 0) Then
        Call gsend(comm,buffer(1:imove),jdnode,Export_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! openKIM halo indicators
    If (kim_data%active) Then
      i = Merge(2*mdir, -2*mdir-1, mdir>0)
      Call kim_data%kcomms%set(i, imove/iadd, config%nlast+1, config%nlast+jmove/iadd)
    End If


    ! load transferred data

    j=Merge(iblock,0,comm%mxnode > 1)
    Do i=1,jmove/iadd
      config%nlast=config%nlast+1

      ! unpack positions

      config%parts(config%nlast)%xxx=buffer(j+1)
      config%parts(config%nlast)%yyy=buffer(j+2)
      config%parts(config%nlast)%zzz=buffer(j+3)

      ! unpack config indexing, site and halo indexing arrays

      config%ltg(config%nlast)  =Nint(buffer(j+4))
      config%lsite(config%nlast)=Nint(buffer(j+5))

      ! unpack remaining halo indexing

      config%ixyz(config%nlast) =Nint(buffer(j+iadd))

      j=j+iadd
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'export_atomic_data deallocation failure'
      Call error(0,message)
    End If

  End Subroutine export_atomic_data

  Subroutine export_atomic_positions(mdir,mlast,ixyz0,domain,config,kim_data,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export atomic positions in domain boundary regions
    ! for halo refresh
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent( In    ) :: mdir,ixyz0(:)
    Integer, Intent( InOut ) :: mlast
    Type( domains_type ), Intent( In    ) :: domain
    Type( kim_type ), Intent( InOut ) :: kim_data
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type ), Intent (InOut) :: comm


    Logical           :: safe,lsx,lsy,lsz,lex,ley,lez,lwrap
    Integer           :: fail,iadd,limit,iblock,     &
      i,j,jxyz,ix,iy,iz,kx,ky,kz, &
      jdnode,kdnode,imove,jmove,itmp
    Real( Kind = wp ) :: uuu,vvv,www,xadd,yadd,zadd

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    Character ( Len = 256 )  :: message

    ! Number of transported quantities per particle

    iadd=3

    fail=0 ; limit=iadd*domain%mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
    Allocate (buffer(1:limit), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'export_atomic_positions allocation failure'
      Call error(0,message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock=limit/Merge(2,1,comm%mxnode > 1)

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! ls.  - wrap-around +1 in . direction (domain on the left MD cell border)
    ! le.  - wrap-around -1 in . direction (domain on the right MD cell border)
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0 ; ky = 0 ; kz = 0
    lsx = .false. ; lex = .false.
    lsy = .false. ; ley = .false.
    lsz = .false. ; lez = .false.
    If      (mdir == -1) Then ! Direction -x
      kx  = 1
      jxyz= 1
      lsx = (domain%idx == 0)

      jdnode = domain%map(1)
      kdnode = domain%map(2)
    Else If (mdir ==  1) Then ! Direction +x
      kx  = 1
      jxyz= 2
      lex = (domain%idx == domain%nx-1)

      jdnode = domain%map(2)
      kdnode = domain%map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky  = 1
      jxyz= 10
      lsy = (domain%idy == 0)

      jdnode = domain%map(3)
      kdnode = domain%map(4)
    Else If (mdir ==  2) Then ! Direction +y
      ky  = 1
      jxyz= 20
      ley = (domain%idy == domain%ny-1)

      jdnode = domain%map(4)
      kdnode = domain%map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz  = 1
      jxyz= 100
      lsz = (domain%idz == 0)

      jdnode = domain%map(5)
      kdnode = domain%map(6)
    Else If (mdir ==  3) Then ! Direction +z
      kz  = 1
      jxyz= 200
      lez = (domain%idz == domain%nz-1)

      jdnode = domain%map(6)
      kdnode = domain%map(5)
    Else
      Call error(46)
    End If

    ! Calculate PBC shift vector due to possible wrap around

    uuu=0.0_wp ; If (lsx) uuu=+1.0_wp ; If (lex) uuu=-1.0_wp
    vvv=0.0_wp ; If (lsy) vvv=+1.0_wp ; If (ley) vvv=-1.0_wp
    www=0.0_wp ; If (lsz) www=+1.0_wp ; If (lez) www=-1.0_wp

    lwrap = (Abs(uuu)+Abs(vvv)+Abs(www) > 0.5_wp)

    If (lwrap) Then
      xadd = config%cell(1)*uuu+config%cell(4)*vvv+config%cell(7)*www
      yadd = config%cell(2)*uuu+config%cell(5)*vvv+config%cell(8)*www
      zadd = config%cell(3)*uuu+config%cell(6)*vvv+config%cell(9)*www
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed

    imove=0
    jmove=0

    ! Initialise array overflow flags

    safe=.true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE

    Do i=1,mlast

      ! If the particle is within the remaining 'inverted halo' of this domain

      If (ixyz0(i) > 0) Then

        ! Get the necessary halo indices

        ix=Mod(ixyz0(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz0(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz0(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

        ! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

        ! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

          ! If safe to proceed

          If ((imove+iadd) <= iblock) Then

            ! pack positions and apply possible PBC shift for the receiver

            If (.not.lwrap) Then
              buffer(imove+1)=config%parts(i)%xxx
              buffer(imove+2)=config%parts(i)%yyy
              buffer(imove+3)=config%parts(i)%zzz
            Else
              buffer(imove+1)=config%parts(i)%xxx+xadd
              buffer(imove+2)=config%parts(i)%yyy+yadd
              buffer(imove+3)=config%parts(i)%zzz+zadd
            End If

          Else

            safe=.false.

          End If
          imove=imove+iadd

        End If

      End If

    End Do

    ! Check for array bound overflow (have arrays coped with outgoing data)

    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=Merge(2,1,comm%mxnode > 1)*imove
      Call gmax(comm,itmp)
      Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
      Call error(54)
    End If

    ! exchange information on buffer sizes

    If (comm%mxnode > 1) Then
      Call girecv(comm,jmove,kdnode,Export_tag)
      Call gsend(comm,imove,jdnode,Export_tag)
      Call gwait(comm)
    Else
      jmove=imove
    End If

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe=((mlast+jmove/iadd) <= config%mxatms)
    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=mlast+jmove/iadd
      Call gmax(comm,itmp)
      Call warning(160,Real(itmp,wp),Real(config%mxatms,wp),0.0_wp)
      Call error(56)
    End If

    ! exchange buffers between nodes (this is a MUST)

    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,Export_tag)
      End If
      If (imove > 0 ) Then
        Call gsend(comm,buffer(1:imove),jdnode,Export_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! openKIM halo indicators
    If (kim_data%active) Then
      i = Merge(2*mdir, -2*mdir-1, mdir>0)
      Call kim_data%kcomms%set(i, imove/iadd, mlast+1, mlast+jmove/iadd)
    End If

    ! load transferred data

    j=Merge(iblock,0,comm%mxnode > 1)
    Do i=1,jmove/iadd
      mlast=mlast+1

      ! unpack positions

      config%parts(mlast)%xxx=buffer(j+1)
      config%parts(mlast)%yyy=buffer(j+2)
      config%parts(mlast)%zzz=buffer(j+3)

      j=j+iadd
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'export_atomic_positions deallocation failure'
      Call error(0,message)
    End If

  End Subroutine export_atomic_positions

  Subroutine mpoles_rotmat_export(mdir,mlast,mxatms,ixyz0,mpoles,domain,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to export multipoles rotation and infinitesimally
    ! rotated matrices in the halo
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer, Intent( In    ) :: mdir,mxatms
    Integer, Intent( InOut ) :: mlast,ixyz0(:)
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut ) :: comm


    Logical           :: safe
    Integer           :: fail,iadd,limit,iblock,          &
      i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
      idl1,idl2,idl3,idl4,             &
      jdnode,kdnode,imove,jmove,itmp

    Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

    Character  ( Len =  256 )  ::  message

    ! Number of transported quantities per particle

    iadd=4*mpoles%max_mpoles+1
    idl1=mpoles%max_mpoles ; idl2=2*mpoles%max_mpoles ; idl3=3*mpoles%max_mpoles ; idl4=4*mpoles%max_mpoles

    fail=0 ; limit=iadd*domain%mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
    Allocate (buffer(1:limit), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'export_atomic_positions allocation failure'
      Call error(0,message)
    End If

    ! Set buffer limit (half for outgoing data - half for incoming)

    iblock=limit/Merge(2,1,comm%mxnode > 1)

    ! DIRECTION SETTINGS INITIALISATION

    ! define the neighbouring domains as sending and receiving with
    ! respect to the direction (mdir)
    ! k.   - direction selection factor
    ! jxyz - halo reduction factor
    ! ls.  - wrap-around +1 in . direction (domain on the left MD cell border)
    ! le.  - wrap-around -1 in . direction (domain on the right MD cell border)
    ! jdnode - destination (send to), kdnode - source (receive from)

    kx = 0 ; ky = 0 ; kz = 0
    If      (mdir == -1) Then ! Direction -x
      kx  = 1
      jxyz= 1
      kxyz= 3

      jdnode = domain%map(1)
      kdnode = domain%map(2)
    Else If (mdir ==  1) Then ! Direction +x
      kx  = 1
      jxyz= 2
      kxyz= 3

      jdnode = domain%map(2)
      kdnode = domain%map(1)
    Else If (mdir == -2) Then ! Direction -y
      ky  = 1
      jxyz= 10
      kxyz= 30

      jdnode = domain%map(3)
      kdnode = domain%map(4)
    Else If (mdir ==  2) Then ! Direction +y
      ky  = 1
      jxyz= 20
      kxyz= 30

      jdnode = domain%map(4)
      kdnode = domain%map(3)
    Else If (mdir == -3) Then ! Direction -z
      kz  = 1
      jxyz= 100
      kxyz= 300

      jdnode = domain%map(5)
      kdnode = domain%map(6)
    Else If (mdir ==  3) Then ! Direction +z
      kz  = 1
      jxyz= 200
      kxyz= 300

      jdnode = domain%map(6)
      kdnode = domain%map(5)
    Else
      Call error(176)
    End If

    ! Initialise counters for length of sending and receiving buffers
    ! imove and jmove are the actual number of particles to get haloed

    imove=0
    jmove=0

    ! Initialise array overflow flags

    safe=.true.

    ! LOOP OVER ALL PARTICLES ON THIS NODE

    Do i=1,mlast

      ! If the particle is within the remaining 'inverted halo' of this domain

      If (ixyz0(i) > 0) Then

        ! Get the necessary halo indices

        ix=Mod(ixyz0(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz0(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz0(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

        ! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

        ! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

          ! If safe to proceed

          If ((imove+iadd) <= iblock) Then

            ! pack rotation matrices and infinitely rotated matrices

            buffer(imove+1        : imove + idl1) = mpoles%global_frame(:,i)
            buffer(imove+1 + idl1 : imove + idl2) = mpoles%rotation_x(:,i)
            buffer(imove+1 + idl2 : imove + idl3) = mpoles%rotation_y(:,i)
            buffer(imove+1 + idl3 : imove + idl4) = mpoles%rotation_z(:,i)

            ! Use the corrected halo reduction factor when the particle is halo to both +&- sides

            buffer(imove+iadd)=Real(ixyz0(i)-Merge(jxyz,kxyz,j == jxyz),wp)

          Else

            safe=.false.

          End If

          imove=imove+iadd

        End If

      End If

    End Do

    ! Check for array bound overflow (have arrays coped with outgoing data)

    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=Merge(2,1,comm%mxnode > 1)*imove
      Call gmax(comm,itmp)
      Call warning(150,Real(itmp,wp),Real(limit,wp),0.0_wp)
      Call error(178)
    End If

    ! exchange information on buffer sizes

    If (comm%mxnode > 1) Then
      Call girecv(comm,jmove,kdnode,ExpMplRM_tag)
      Call gsend(comm,imove,jdnode,ExpMplRM_tag)
      Call gwait(comm)
    Else
      jmove=imove
    End If

    ! Check for array bound overflow (can arrays cope with incoming data)

    safe=((mlast+jmove/iadd) <= mxatms)
    Call gcheck(comm,safe)
    If (.not.safe) Then
      itmp=mlast+jmove/iadd
      Call gmax(comm,itmp)
      Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
      Call error(180)
    End If

    ! exchange buffers between nodes (this is a MUST)

    If (comm%mxnode > 1) Then
      If (jmove > 0) Then
        Call girecv(comm,buffer(iblock+1:iblock+jmove),kdnode,ExpMplRM_tag)
      End If
      If (imove > 0) Then
        Call gsend(comm,buffer(1:imove),jdnode,ExpMplRM_tag)
      End If
      If (jmove > 0) Call gwait(comm)
    End If

    ! load transferred data

    j=Merge(iblock,0,comm%mxnode > 1)
    Do i=1,jmove/iadd
      mlast=mlast+1

      ! unpack rotation matrices and infinitesimal rotation matrices

      mpoles%global_frame(:,mlast) = buffer(j+1        : j + idl1)
      mpoles%rotation_x(:,mlast) = buffer(j+1 + idl1 : j + idl2)
      mpoles%rotation_y(:,mlast) = buffer(j+1 + idl2 : j + idl3)
      mpoles%rotation_z(:,mlast) = buffer(j+1 + idl3 : j + idl4)

      ixyz0(mlast)=Nint(buffer(j+iadd))
      j=j+iadd
    End Do

    Deallocate (buffer, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'mpoles_rotmat_export deallocation failure'
      Call error(0,message)
    End If

  End Subroutine mpoles_rotmat_export

  Subroutine mpoles_rotmat_set_halo(mpoles,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to set and infinitesimally rotate native multipoles
    ! rotation matrices and then arrange exchange of halo data between
    ! neighbouring domains/nodes
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov february 2016
    ! contrib   - h.a.boateng february 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( domains_type ), Intent( In    ) :: domain
    Type( configuration_type ), Intent( InOut ) :: config
    Type( comms_type), Intent( InOut ) :: comm

    Logical :: safe
    Integer :: fail,i,mlast

    Integer, Allocatable :: ixyz0(:)

    Character ( Len = 256 )  ::  message

    Do i=1,config%natms
      mpoles%flg(i)=0
      If (mpoles%max_order < 3) Then
        Call rotate_mpoles_d(i,mpoles,config,comm)
      Else
        Call rotate_mpoles(i,mpoles,config,comm)
      End If
    End Do

    ! Communicate the matrices in the halo

    fail = 0
    Allocate (ixyz0(1:config%mxatms), Stat = fail)
    If (fail > 0) Then
      Write(message,'(a)') 'mpoles_rotmat_set_halo allocation failure'
      Call error(0,message)
    End If
    ixyz0(1:config%nlast) = config%ixyz(1:config%nlast)

    ! No halo, start with domain only particles

    mlast=config%natms

    ! exchange atom data in -/+ x directions

    Call mpoles_rotmat_export(-1,mlast,config%mxatms,ixyz0,mpoles,domain,comm)
    Call mpoles_rotmat_export( 1,mlast,config%mxatms,ixyz0,mpoles,domain,comm)

    ! exchange atom data in -/+ y directions

    Call mpoles_rotmat_export(-2,mlast,config%mxatms,ixyz0,mpoles,domain,comm)
    Call mpoles_rotmat_export( 2,mlast,config%mxatms,ixyz0,mpoles,domain,comm)

    ! exchange atom data in -/+ z directions

    Call mpoles_rotmat_export(-3,mlast,config%mxatms,ixyz0,mpoles,domain,comm)
    Call mpoles_rotmat_export( 3,mlast,config%mxatms,ixyz0,mpoles,domain,comm)

    ! check atom totals after data transfer

    safe=(mlast == config%nlast)
    Call gcheck(comm,safe)
    If (.not.safe) Call error(174)

    Deallocate (ixyz0, Stat = fail)
    If (fail > 0) Then
      Write(message,'(a)') 'mpoles_rotmat_set_halo deallocation failure'
      Call error(0,message)
    End If

  End Subroutine mpoles_rotmat_set_halo

  Subroutine relocate_particles(dvar,cutoff_extended,lbook,lmsd,megatm,flow,cshell,cons, &
      pmf,stats,ewld,thermo,green,bond,angle,dihedral,inversion,tether, &
      neigh,sites,minim,mpoles,rigid,domain,config,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 routine to arrange relocation of data between neighbouring
    ! domains/nodes after positions updates
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith august 1998
    ! amended   - i.t.todorov february 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Real( Kind = wp ), Intent( In    ) :: dvar,cutoff_extended
    Logical,           Intent( In    ) :: lbook
    Logical,           Intent( In    ) :: lmsd
    Integer,           Intent( In    ) :: megatm
    Type( flow_type), Intent( InOut ) :: flow
    Type( pmf_type), Intent( InOut ) :: pmf
    Type( core_shell_type), Intent( InOut ) :: cshell
    Type( constraints_type), Intent( InOut ) :: cons
    Type( stats_type ), Intent( InOut ) :: stats
    Type( ewald_type ), Intent( InOut ) :: ewld
    Type( thermostat_type ), Intent( InOut ) :: thermo
    Type( greenkubo_type ), Intent( InOut ) :: green
    Type( bonds_type ), Intent( InOut ) :: bond
    Type( angles_type ), Intent( InOut ) :: angle
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    Type( inversions_type ), Intent( InOut ) :: inversion
    Type( tethers_type ), Intent( InOut ) :: tether
    Type( neighbours_type ), Intent( InOut ) :: neigh
    Type( site_type ), Intent( In    ) :: sites
    Type( minimise_type ), Intent( InOut ) :: minim
    Type( mpole_type ), Intent( InOut ) :: mpoles
    Type( rigid_bodies_type ), Intent( InOut ) :: rigid
    Type( configuration_type ), Intent( InOut ) :: config
    Type( domains_type ), Intent( In    ) :: domain
    Type( comms_type ), Intent( InOut ) :: comm
    Real( Kind = wp ) :: cut

    Logical           :: safe(1:9)
    Integer           :: i,nlimit,ipx,ipy,ipz,itmp(1:9),jtmp(1:9)
    Real( Kind = wp ) :: big(1:3),det,celprp(1:10),rcell(1:9),x,y,z,tmp

    Character( Len = 256 ) :: message

    ! Define cut

    cut=cutoff_extended+1.0e-6_wp

    ! rescale mock cell vectors for non-periodic system

    If (config%imcon == 0 .or. config%imcon == 6) Then

      ! find maximum x,y,z positions

      big=0.0_wp

      Do i =1,config%natms
        big(1)=Max(big(1),Abs(config%parts(i)%xxx))
        big(2)=Max(big(2),Abs(config%parts(i)%yyy))
        big(3)=Max(big(3),Abs(config%parts(i)%zzz))
      End Do

      Call gmax(comm,big)

      If (config%imcon == 0) Then

        config%cell(1)=Max(2.0_wp*big(1)+cut,3.0_wp*cut,config%cell(1))
        config%cell(5)=Max(2.0_wp*big(2)+cut,3.0_wp*cut,config%cell(5))
        config%cell(9)=Max(2.0_wp*big(3)+cut,3.0_wp*cut,config%cell(9))

        config%cell(2)=0.0_wp
        config%cell(3)=0.0_wp
        config%cell(4)=0.0_wp
        config%cell(6)=0.0_wp
        config%cell(7)=0.0_wp
        config%cell(8)=0.0_wp

      Else If (config%imcon == 6) Then

        config%cell(9)=Max(2.0_wp*big(3)+cut,3.0_wp*cut,config%cell(9))

      End If

    End If

    ! Get the dimensional properties of the MD cell

    Call dcell(config%cell,celprp)

    ! Only when we have more than one domain we need real relocation

    If (comm%mxnode > 1) Then

      Call invert(config%cell,rcell,det)

      ! Convert atomic positions from MD centred
      ! Cartesian coordinates to reduced space ones.
      ! Populate the move (former halo) indicator array.
      ! Here we assume no particle has moved more than
      ! a link-cell width (> cutoff_extended) in any direction
      ! since last call to relocate as we don't test this!!!

      config%ixyz(1:config%natms)=0 ! Initialise move (former halo) indicator
      Do i=1,config%natms
        x=rcell(1)*config%parts(i)%xxx+rcell(4)*config%parts(i)%yyy+rcell(7)*config%parts(i)%zzz
        y=rcell(2)*config%parts(i)%xxx+rcell(5)*config%parts(i)%yyy+rcell(8)*config%parts(i)%zzz
        z=rcell(3)*config%parts(i)%xxx+rcell(6)*config%parts(i)%yyy+rcell(9)*config%parts(i)%zzz

        ! assign domain coordinates (call for errors)

        ipx=Int((x+0.5_wp)*domain%nx_real)
        ipy=Int((y+0.5_wp)*domain%ny_real)
        ipz=Int((z+0.5_wp)*domain%nz_real)

        If (domain%idx == 0) Then
          If (x < -half_plus) config%ixyz(i)=config%ixyz(i)+1
        Else
          If (ipx < domain%idx) config%ixyz(i)=config%ixyz(i)+1
        End If
        If (domain%idx == domain%nx-1) Then
          If (x >= half_minus) config%ixyz(i)=config%ixyz(i)+2
        Else
          If (ipx > domain%idx) config%ixyz(i)=config%ixyz(i)+2
        End If

        If (domain%idy == 0) Then
          If (y < -half_plus) config%ixyz(i)=config%ixyz(i)+10
        Else
          If (ipy < domain%idy) config%ixyz(i)=config%ixyz(i)+10
        End If
        If (domain%idy == domain%ny-1) Then
          If (y >= half_minus) config%ixyz(i)=config%ixyz(i)+20
        Else
          If (ipy > domain%idy) config%ixyz(i)=config%ixyz(i)+20
        End If

        If (domain%idz == 0) Then
          If (z < -half_plus) config%ixyz(i)=config%ixyz(i)+100
        Else
          If (ipz < domain%idz) config%ixyz(i)=config%ixyz(i)+100
        End If
        If (domain%idz == domain%nz-1) Then
          If (z >= half_minus) config%ixyz(i)=config%ixyz(i)+200
        Else
          If (ipz > domain%idz) config%ixyz(i)=config%ixyz(i)+200
        End If
      End Do

      ! exchange atom data in -/+ x directions

      Call deport_atomic_data(-1,lbook,lmsd,cshell,cons,pmf,stats,ewld,thermo, &
        green,bond,angle,dihedral,inversion,tether,neigh,minim,mpoles,rigid, &
        domain,config,comm)
      Call deport_atomic_data( 1,lbook,lmsd,cshell,cons,pmf,stats,ewld,thermo, &
        green,bond,angle,dihedral,inversion,tether,neigh,minim,mpoles,rigid, &
        domain,config,comm)

      ! exchange atom data in -/+ y directions

      Call deport_atomic_data(-2,lbook,lmsd,cshell,cons,pmf,stats,ewld,thermo, &
        green,bond,angle,dihedral,inversion,tether,neigh,minim,mpoles,rigid, &
        domain,config,comm)
      Call deport_atomic_data( 2,lbook,lmsd,cshell,cons,pmf,stats,ewld,thermo, &
        green,bond,angle,dihedral,inversion,tether,neigh,minim,mpoles,rigid, &
        domain,config,comm)

      ! exchange atom data in -/+ z directions

      Call deport_atomic_data(-3,lbook,lmsd,cshell,cons,pmf,stats,ewld,thermo, &
        green,bond,angle,dihedral,inversion,tether,neigh,minim,mpoles,rigid, &
        domain,config,comm)
      Call deport_atomic_data( 3,lbook,lmsd,cshell,cons,pmf,stats,ewld,thermo, &
        green,bond,angle,dihedral,inversion,tether,neigh,minim,mpoles,rigid, &
        domain,config,comm)

      ! check system for loss of atoms

      safe(1)=(All(config%ixyz(1:config%natms) == 0)) ; Call gcheck(comm,safe(1),"enforce")
      nlimit=config%natms ; Call gsum(comm,nlimit)
      If ((.not.safe(1)) .or. nlimit /= megatm) Call error(58)

      ! reassign atom properties

      Do i=1,config%natms
        config%atmnam(i)=sites%site_name(config%lsite(i))
        config%ltype(i)=sites%type_site(config%lsite(i))
        config%parts(i)%chge=sites%charge_site(config%lsite(i))
        config%weight(i)=sites%weight_site(config%lsite(i))
        config%lfrzn(i)=sites%freeze_site(config%lsite(i))
        config%lfree(i)=sites%free_site(config%lsite(i))
      End Do

      If (lbook) Then
        safe=.true. ! Initialise safety flag

        ! Change nlast and refresh record global atom indices for local sorting
        ! since particles may have been relocated across domains as the
        ! the old halo is invalid and a new one is not set yet.  This allows for
        ! local_index search over natms in pmf_units_set and compress_book_intra.
        ! Otherwise, everywhere else in the code, the search is over nlast as
        ! domain only indices are caught by the condition (1 >= index <= natms)!!!

        config%nlast=config%natms
        Do i=1,config%nlast
          config%lsi(i)=i
          config%lsa(i)=config%ltg(i)
        End Do
        Call shellsort2(config%nlast,config%lsi,config%lsa)

        ! Check safety of working arrays for all active bookkeeping arrays

        If (cshell%megshl > 0) safe(1)=(cshell%ntshl  <= cshell%mxshl )
        If (cons%m_con  > 0) safe(2)=(cons%ntcons <= cons%mxcons)
        If (pmf%megpmf > 0) safe(3)=(pmf%ntpmf  <= pmf%mxpmf )
        If (rigid%on) safe(4)=(rigid%n_types  <= rigid%max_rigid )
        If (tether%total > 0) safe(5)=(tether%ntteth <= tether%mxteth)
        If (bond%total > 0) safe(6)=(bond%n_types <= bond%max_bonds)
        If (angle%total > 0) safe(7)=(angle%n_types <= angle%max_angles)
        If (dihedral%total > 0) safe(8)=(dihedral%n_types <= dihedral%max_angles)
        If (inversion%total > 0) safe(9)=(inversion%n_types  <= inversion%max_angles )

        Call gcheck(comm,safe)

        If (Any(.not.safe)) Then
          itmp(1)=cshell%ntshl  ; jtmp(1)=cshell%mxshl
          itmp(2)=cons%ntcons ; jtmp(2)=cons%mxcons
          itmp(3)=pmf%ntpmf  ; jtmp(3)=pmf%mxpmf
          itmp(4)=rigid%n_types  ; jtmp(4)=rigid%max_rigid
          itmp(5)=tether%ntteth ; jtmp(5)=tether%mxteth
          itmp(6)=bond%n_types ; jtmp(6)=bond%max_bonds
          itmp(7)=angle%n_types ; jtmp(7)=angle%max_angles
          itmp(8)=dihedral%n_types ; jtmp(8)=dihedral%max_angles
          itmp(9)=inversion%n_types  ; jtmp(9)=inversion%max_angles

          Call gmax(comm,itmp(1:9))

          tmp=1.0_wp
          Do i=1,9
            tmp=Max(tmp,1.0_wp+Real(itmp(i),wp)/Real(Max(1,jtmp(i)),wp))
          End Do

          Write(message,'(a,i0)') 'estimated densvar value for passing this stage safely is : ', &
            Nint((dvar*tmp-1.0_wp)*100.0_wp+0.5_wp)
          Call warning(message,.true.)

        End If

        If (.not.safe(1)) Call error( 59)
        If (.not.safe(2)) Call error( 41)
        If (.not.safe(3)) Call error(488)
        If (.not.safe(4)) Call error(640)
        If (.not.safe(5)) Call error( 63)
        If (.not.safe(6)) Call error( 31)
        If (.not.safe(7)) Call error( 51)
        If (.not.safe(8)) Call error( 61)
        If (.not.safe(9)) Call error( 77)

        ! Update shared core-shell, constraint, PMF and RB units

        If (cshell%megshl > 0) Call pass_shared_units &
          (config,cshell%mxshl, Lbound(cshell%listshl,Dim=1),Ubound(cshell%listshl,Dim=1),&
          cshell%ntshl, cshell%listshl,cshell%mxfshl,&
          cshell%legshl,cshell%lshmv_shl,cshell%lishp_shl,cshell%lashp_shl, &
          flow%oldjob_shared_units,domain,comm,&
          rigid%q0,rigid%q1,rigid%q2,rigid%q3,rigid%vxx,rigid%vyy,rigid%vzz, &
          rigid%oxx,rigid%oyy,rigid%ozz)

        If (cons%m_con  > 0) Call pass_shared_units &
          (config,cons%mxcons,Lbound(cons%listcon,Dim=1),Ubound(cons%listcon,Dim=1),&
          cons%ntcons,cons%listcon,cons%mxfcon,cons%legcon,&
          cons%lshmv_con,cons%lishp_con,cons%lashp_con,flow%oldjob_shared_units,domain,comm,&
          rigid%q0,rigid%q1,rigid%q2,rigid%q3,rigid%vxx,rigid%vyy,rigid%vzz, &
          rigid%oxx,rigid%oyy,rigid%ozz)

        If (pmf%megpmf > 0) Call pmf_units_set(pmf,config,comm)

        If (rigid%on) Call pass_shared_units &
          (config,rigid%max_rigid, Lbound(rigid%list,Dim=1),Ubound(rigid%list,Dim=1),rigid%n_types, &
          rigid%list,rigid%max_frozen,rigid%legend,rigid%share,rigid%list_shared, &
          rigid%map_shared,flow%oldjob_shared_units,domain,comm,&
          rigid%q0,rigid%q1,rigid%q2,rigid%q3,rigid%vxx,rigid%vyy,rigid%vzz, &
          rigid%oxx,rigid%oyy,rigid%ozz)

        ! Compress the rest of the bookkeeping arrays if needed

        If (tether%total > 0) Call compress_book_intra &
          (tether%mxteth,tether%ntteth,Ubound(tether%listtet,Dim=1),&
          tether%listtet,tether%mxftet,tether%legtet,cons,config,comm)
        If (bond%total > 0) Then
          Call compress_book_intra(bond%max_bonds,bond%n_types, &
            Ubound(bond%list,Dim=1),bond%list,bond%max_legend,bond%legend,cons,config,comm)
        End If
        If (angle%total > 0) Then
          Call compress_book_intra(angle%max_angles,angle%n_types, &
            Ubound(angle%list,Dim=1),angle%list,angle%max_legend,angle%legend,cons,config,comm)
        End If
        If (dihedral%total > 0) Then
          Call compress_book_intra(dihedral%max_angles,dihedral%n_types, &
            Ubound(dihedral%list,Dim=1),dihedral%list,dihedral%max_legend,dihedral%legend,cons,config,comm)
        End If
        If (inversion%total > 0) Then
          Call compress_book_intra(inversion%max_angles,inversion%n_types, &
            Ubound(inversion%list,Dim=1),inversion%list,inversion%max_legend,inversion%legend,cons,config,comm)
        End If

      End If

    Else

      ! Restore periodic boundaries (re-bound > re-wrap)

      Call pbcshift(config%imcon,config%cell,config%natms,config%parts)

    End If

    ! Halt program if potential cutoff exceeds cell width

    If (cutoff_extended > Min(celprp(7),celprp(8),celprp(9))/2.0_wp) Call error(95)

  End Subroutine relocate_particles


End Module deport_data
