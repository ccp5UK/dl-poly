Module core_shell

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module for defining global core-shell interaction variables
  ! and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov december 2016
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds,           Only : wp
  Use comms,           Only : comms_type,gsync,gsum,gcheck,gmax
  Use setup,           Only : nrite,boltz,engunit,output,mxatms,mxatdm,zero_plus,&
    mxtmls,mxproc,mxlshp
  Use configuration,   Only : imcon,cell,natms,nlast,lsi,lsa,xxx,yyy,zzz,fxx,fyy,fzz, &
    weight,vxx,vyy,vzz,lfrzn,freeze_atoms
  Use parse,           Only : strip_blanks,lower_case
  Use shared_units,    Only : update_shared_units
  Use numerics,        Only : local_index,images
  Use errors_warnings, Only : error,warning,info
  Use statistics, Only : stats_type

  Implicit None
  Private
  Integer, Parameter, Public :: SHELL_ADIABATIC = 1
  InTeger, Parameter, Public :: SHELL_RELAXED = 2
  Type, Public :: core_shell_type
    Private
    Integer, Public :: mxshl,mxtshl,mxfshl,megshl,mxlshp
    Logical,                        Public :: lshmv_shl = .false.

    Integer,                        Public :: ntshl  = 0 , &
      ntshl1 = 0 , &
      ntshl2 = 0 , keyshl

    Real( Kind = wp ),              Public :: smax = 0.0_wp


    Integer,           Allocatable, Public :: numshl(:)
    Integer,           Allocatable, Public :: lstshl(:,:),listshl(:,:),legshl(:,:)
    Integer,           Allocatable, Public :: lishp_shl(:),lashp_shl(:)

    Real( Kind = wp ), Allocatable, Public :: prmshl(:,:)
  Contains
    Private
    Procedure, Public :: init => allocate_core_shell_arrays
    Procedure, Public :: deallocate_core_shell_tmp_arrays
    Final :: deallocate_core_shell_arrays
  End Type core_shell_type

  Public :: core_shell_forces
  Public :: core_shell_kinetic
  Public :: core_shell_on_top
  Public :: core_shell_quench
  Public :: core_shell_relax

Contains

  Subroutine allocate_core_shell_arrays(T,mxatdm,mxtmls,mxlshp,mxproc)
   Class( core_shell_type ) :: T
     Integer, Intent(In) :: mxatdm,mxtmls,mxlshp,mxproc

    Integer, Dimension( 1:6 ) :: fail

    fail = 0

    Allocate (T%numshl(1:mxtmls),                        Stat = fail(1))
    Allocate (T%lstshl(1:2,1:T%mxtshl),                    Stat = fail(2))
    Allocate (T%listshl(0:2,1:T%mxshl),                    Stat = fail(3))
    Allocate (T%legshl(0:T%mxfshl,1:mxatdm),               Stat = fail(4))
    Allocate (T%lishp_shl(1:mxlshp),T%lashp_shl(1:mxproc), Stat = fail(5))
    Allocate (T%prmshl(1:2,1:T%mxtshl),                    Stat = fail(6))

    If (Any(fail > 0)) Call error(1005)

    T%numshl  = 0
    T%lstshl  = 0
    T%listshl = 0
    T%legshl  = 0

    T%lishp_shl = 0 ; T%lashp_shl = 0

    T%prmshl  = 0.0_wp

  End Subroutine allocate_core_shell_arrays

  Subroutine deallocate_core_shell_tmp_arrays(T)
   Class( core_shell_type ) :: T

    Integer :: fail(2)

    fail = 0

    If (Allocated(T%numshl)) Deallocate (T%numshl, Stat = fail(1))
    If (Allocated(T%lstshl)) Deallocate (T%lstshl, Stat = fail(2))

    If (Any(fail > 0)) Call error(1030)

  End Subroutine deallocate_core_shell_tmp_arrays

  Subroutine deallocate_core_shell_arrays(T)
    Type( core_shell_type ) :: T

    Integer :: fail(7)

    fail = 0

    If (Allocated(T%numshl)) Deallocate (T%numshl, Stat = fail(1))
    If (Allocated(T%lstshl)) Deallocate (T%lstshl, Stat = fail(2))
    If (Allocated(T%listshl)) Deallocate (T%listshl, Stat = fail(3))
    If (Allocated(T%legshl)) Deallocate (T%legshl, Stat = fail(4))
    If (Allocated(T%lishp_shl)) Deallocate (T%lishp_shl, Stat = fail(5))
    If (Allocated(T%lashp_shl)) Deallocate (T%lashp_shl, Stat = fail(6))
    If (Allocated(T%prmshl)) Deallocate (T%prmshl, Stat = fail(7))
    If (Any(fail > 0)) Call error(1030)
  End Subroutine deallocate_core_shell_arrays


  Subroutine core_shell_forces(cshell,stat,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating core-shell model spring energy
    ! and force terms
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov june 2016
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( stats_type ), Intent( InOut ) :: stat
    Type( comms_type ),                  Intent( InOut ) :: comm

    Logical           :: safe
    Integer           :: fail(1:2),i,j,ia,ib,kk
    Real( Kind = wp ) :: rabsq,fx,fy,fz,gamma,omega,r_4_fac, &
      strs1,strs2,strs3,strs5,strs6,strs9,buffer(1:2)

    Logical,           Allocatable :: lunsafe(:)
    Integer,           Allocatable :: lstopt(:,:)
    Real( Kind = wp ), Allocatable :: xdab(:),ydab(:),zdab(:)
    Character( Len = 256 ) :: message

    fail=0
    Allocate (lunsafe(1:cshell%mxshl),lstopt(0:2,1:cshell%mxshl),      Stat=fail(1))
    Allocate (xdab(1:cshell%mxshl),ydab(1:cshell%mxshl),zdab(1:cshell%mxshl), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_forces allocation failure'
      Call error(0,message)
    End If

    r_4_fac = 1.0_wp/24.0_wp ! aharmonic shell coefficient = 1/(4!)

    ! calculate core-shell separation vectors

    Do i=1,cshell%ntshl
      lunsafe(i)=.false.

      ! indices of atoms in a core-shell

      ia=local_index(cshell%listshl(1,i),nlast,lsi,lsa) ; lstopt(1,i)=ia
      ib=local_index(cshell%listshl(2,i),nlast,lsi,lsa) ; lstopt(2,i)=ib

      lstopt(0,i)=0
      If (ia > 0 .and. ib > 0) Then ! Tag
        If (ia <= natms .or. ib <= natms) Then
          lstopt(0,i)=1
        End If
      Else                          ! Detect uncompressed unit
        If ( ((ia > 0 .and. ia <= natms) .or.   &
          (ib > 0 .and. ib <= natms)) .and. &
          (ia == 0 .or. ib == 0) ) lunsafe(i)=.true.
      End If

      ! components of bond vector

      If (lstopt(0,i) > 0) Then
        xdab(i)=xxx(ia)-xxx(ib)
        ydab(i)=yyy(ia)-yyy(ib)
        zdab(i)=zzz(ia)-zzz(ib)
      Else ! (DEBUG)
        xdab(i)=0.0_wp
        ydab(i)=0.0_wp
        zdab(i)=0.0_wp
      End If
    End Do

    ! Check for uncompressed units

    safe = .not. Any(lunsafe(1:cshell%ntshl))
    Call gcheck(comm,safe)
    If (.not.safe) Then
      Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
          Do i=1,cshell%ntshl
            If (lunsafe(i)) Then
              Write(message,'(a,2(i10,a))')     &
                'global unit number', cshell%listshl(0,i), &
                ' , with a head particle number', cshell%listshl(1,i),   &
                ' contributes towards next error'
              Call warning(message)
            End If
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(100)
    End If

    ! periodic boundary condition

    Call images(imcon,cell,cshell%ntshl,xdab,ydab,zdab)

    ! initialise stress tensor accumulators

    strs1=0.0_wp
    strs2=0.0_wp
    strs3=0.0_wp
    strs5=0.0_wp
    strs6=0.0_wp
    strs9=0.0_wp

    ! zero core-shell energy and virial accumulators

    stat%engshl=0.0_wp
    stat%virshl=0.0_wp

    ! loop over all specified core-shell units

    Do i=1,cshell%ntshl
      If (lstopt(0,i) > 0) Then

        ! indices of atoms in a core-shell

        ia=lstopt(1,i)
        ib=lstopt(2,i)

        ! define components of bond vector

        rabsq = xdab(i)**2+ydab(i)**2+zdab(i)**2

        ! index of potential function parameters

        kk=cshell%listshl(0,i)

        ! calculate scalar constant terms using spring potential function
        ! and the parameters in array prmshl

        omega=(0.5_wp*cshell%prmshl(1,kk)+r_4_fac*cshell%prmshl(2,kk)*rabsq)*rabsq
        gamma=cshell%prmshl(1,kk)+cshell%prmshl(2,kk)*rabsq

        ! calculate forces

        fx = -gamma*xdab(i)
        fy = -gamma*ydab(i)
        fz = -gamma*zdab(i)

        If (ia <= natms) Then

          fxx(ia)=fxx(ia)+fx
          fyy(ia)=fyy(ia)+fy
          fzz(ia)=fzz(ia)+fz

          ! calculate core-shell unit energy

          stat%engshl=stat%engshl+omega
          stat%virshl=stat%virshl+gamma*rabsq

          ! calculate stress tensor

          strs1 = strs1 + xdab(i)*fx
          strs2 = strs2 + xdab(i)*fy
          strs3 = strs3 + xdab(i)*fz
          strs5 = strs5 + ydab(i)*fy
          strs6 = strs6 + ydab(i)*fz
          strs9 = strs9 + zdab(i)*fz

        End If

        If (ib <= natms) Then

          fxx(ib)=fxx(ib)-fx
          fyy(ib)=fyy(ib)-fy
          fzz(ib)=fzz(ib)-fz

        End If

      End If
    End Do

    ! complete stress tensor

    stat%stress(1) = stat%stress(1) + strs1
    stat%stress(2) = stat%stress(2) + strs2
    stat%stress(3) = stat%stress(3) + strs3
    stat%stress(4) = stat%stress(4) + strs2
    stat%stress(5) = stat%stress(5) + strs5
    stat%stress(6) = stat%stress(6) + strs6
    stat%stress(7) = stat%stress(7) + strs3
    stat%stress(8) = stat%stress(8) + strs6
    stat%stress(9) = stat%stress(9) + strs9

    ! sum contributions to potential and virial

    If (comm%mxnode > 1) Then
      buffer(1)=stat%engshl
      buffer(2)=stat%virshl
      Call gsum(comm,buffer(1:2))
      stat%engshl=buffer(1)
      stat%virshl=buffer(2)
    End If

    Deallocate (lunsafe,lstopt, Stat=fail(1))
    Deallocate (xdab,ydab,zdab, Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_forces deallocation failure'
      Call error(0,message)
    End If

  End Subroutine core_shell_forces

  Subroutine core_shell_kinetic(shlke,cshell,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating the internal kinetic energy of
    ! core-shell units in the shell polarisation model
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov january 2008
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real( Kind = wp ),  Intent(   Out ) :: shlke
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Type( comms_type ), Intent( InOut ) :: comm

    Integer           :: i,j,k
    Real( Kind = wp ) :: rmu,rvx,rvy,rvz

    ! gather velocities of shared particles

    If (cshell%lshmv_shl) Then
      Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,vxx,vyy,vzz,comm)
    End If

    ! initialise energy

    shlke=0.0_wp

    ! loop over all specified core-shell pairs

    Do k=1,cshell%ntshl

      ! indices of atoms involved

      i=local_index(cshell%listshl(1,k),nlast,lsi,lsa)
      j=local_index(cshell%listshl(2,k),nlast,lsi,lsa)

      ! for all native and natively shared core-shell units

      If ((i > 0 .and. j > 0) .and. (i <= natms .or. j <= natms)) Then

        ! calculate reduced mass

        rmu=(weight(i)*weight(j))/(weight(i)+weight(j))

        ! frozen particles' velocities are zero
        ! calculate constraint vector normal

        rvx=vxx(j)-vxx(i)
        rvy=vyy(j)-vyy(i)
        rvz=vzz(j)-vzz(i)

        ! calculate core-shell internal kinetic energy

        If (i <= natms .and. j <= natms) Then

          ! for native core-shell units - full contribution

          shlke=shlke+0.5_wp*rmu*(rvx*rvx+rvy*rvy+rvz*rvz)

        Else

          ! If ( (i <= natms .and. j > natms) .or. (j <= natms .and. i > natms) ) Then
          ! for shared core-shell units - halved contribution

          shlke=shlke+0.25_wp*rmu*(rvx*rvx+rvy*rvy+rvz*rvz)

        End If

      End If

    End Do

    ! global sum of core-shell internal kinetic energy

    Call gsum(comm,shlke)

  End Subroutine core_shell_kinetic

  Subroutine core_shell_on_top(cshell,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for positioning shells on top of their cores
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov august 2010
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    Type( comms_type ),  Intent( InOut ) :: comm
    Type( core_shell_type ) , Intent( InOut ) :: cshell
    Logical :: safe
    Integer :: fail,i,j,ia,ib

    Logical, Allocatable :: lunsafe(:)
    Character( Len = 256 ) :: message

    fail=0
    Allocate (lunsafe(1:cshell%mxshl), Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'core_shell_on_top allocation failure'
      Call error(0,message)
    End If


    ! Coincide shells with their cores

    Do i=1,cshell%ntshl
      lunsafe(i)=.false.

      ! indices of atoms in a core-shell

      ia=local_index(cshell%listshl(1,i),nlast,lsi,lsa) ! This is a core
      ib=local_index(cshell%listshl(2,i),nlast,lsi,lsa) ! This is a shell

      ! For every shell in the domain get the coordinates of its
      ! corresponding core (which must be in the domain+hello
      ! area by construction, if not go to a controlled termination)

      If (ib > 0 .and. ib <= natms .and. ia > 0) Then
        xxx(ib)=xxx(ia)
        yyy(ib)=yyy(ia)
        zzz(ib)=zzz(ia)
      End If

      ! Detect uncompressed unit

      If ( ((ia > 0 .and. ia <= natms) .or.   &
        (ib > 0 .and. ib <= natms)) .and. &
        (ia == 0 .or. ib == 0) ) lunsafe(i)=.true.
    End Do

    ! Check for uncompressed units

    safe = .not. Any(lunsafe(1:cshell%ntshl))
    Call gcheck(comm,safe)
    If (.not.safe) Then
      Do j=0,comm%mxnode-1
        If (comm%idnode == j) Then
          Do i=1,cshell%ntshl
            If (lunsafe(i)) Then
              Write(message,'(a,2(i10,a))')     &
                'global unit number', cshell%listshl(0,i), &
                ' , with a head particle number', cshell%listshl(1,i),   &
                ' contributes towards next error'
              Call warning(message)
            End If
          End Do
        End If
        Call gsync(comm)
      End Do
      Call error(100)
    End If

    Deallocate (lunsafe, Stat=fail)
    If (fail > 0) Then
      Write(message,'(a)') 'core_shell_on_top deallocation failure'
      Call error(0,message)
    End If

  End Subroutine core_shell_on_top

  Subroutine core_shell_quench(safe,temp,cshell,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for quenching the internal bond energies of ions
    ! defined by shell model
    !
    ! copyright - daresbury laboratory
    ! author    - w.smith august 1999
    ! amended   - i.t.todorov november 2012
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,           Intent(   Out ) :: safe
    Type( core_shell_type ), Intent( InOut ) :: cshell
    Real( Kind = wp ), Intent( In    ) :: temp
    Type( comms_type), Intent( InOut ) :: comm

    Logical           :: safek
    Integer           :: ia,ib,k
    Real( Kind = wp ) :: dvx,dvy,dvz,pke,rmu,scl,tke,tmx,tmy,tmz

    ! Initialise safe flag

    safe=.true.

    ! gather velocities of shared particles

    If (cshell%lshmv_shl) Then
      Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,vxx,vyy,vzz,comm)
    End If

    ! permitted core-shell internal kinetic energy

    pke=boltz*temp*1.0e-4_wp

    ! amend core and shell velocities preserving total momentum

    Do k=1,cshell%ntshl
      ia=local_index(cshell%listshl(1,k),nlast,lsi,lsa)
      ib=local_index(cshell%listshl(2,k),nlast,lsi,lsa)

      If ((ia > 0 .and. ib > 0) .and. (ia <= natms .or. ib <= natms)) Then
        rmu=(weight(ia)*weight(ib))/(weight(ia)+weight(ib))

        ! frozen atoms have zero velocity (no worries)

        dvx=vxx(ib)-vxx(ia)
        dvy=vyy(ib)-vyy(ia)
        dvz=vzz(ib)-vzz(ia)

        tke=rmu*(dvx*dvx+dvy*dvy+dvz*dvz)

        safek = .not.( tke > pke .and. (tke-pke)/pke > 1.0e-3_wp )
        If (.not.safek) Then
          scl=Sqrt(pke/tke)

          tmx=weight(ia)*vxx(ia)+weight(ib)*vxx(ib)
          tmy=weight(ia)*vyy(ia)+weight(ib)*vyy(ib)
          tmz=weight(ia)*vzz(ia)+weight(ib)*vzz(ib)

          ! no corrections for frozen cores

          If (lfrzn(ia) == 0) Then
            vxx(ia)=tmx/(weight(ia)+weight(ib))-scl*rmu*dvx/weight(ia)
            vyy(ia)=tmy/(weight(ia)+weight(ib))-scl*rmu*dvy/weight(ia)
            vzz(ia)=tmz/(weight(ia)+weight(ib))-scl*rmu*dvz/weight(ia)

            vxx(ib)=tmx/(weight(ia)+weight(ib))+scl*rmu*dvx/weight(ib)
            vyy(ib)=tmy/(weight(ia)+weight(ib))+scl*rmu*dvy/weight(ib)
            vzz(ib)=tmz/(weight(ia)+weight(ib))+scl*rmu*dvz/weight(ib)
          Else
            vxx(ib)=tmx/(weight(ia)+weight(ib))+scl*rmu*dvx/weight(ib)
            vyy(ib)=tmy/(weight(ia)+weight(ib))+scl*rmu*dvy/weight(ib)
            vzz(ib)=tmz/(weight(ia)+weight(ib))+scl*rmu*dvz/weight(ib)
          End If

          safe=.false.
        End If
      End If
    End Do

    Call gcheck(comm,safe)

  End Subroutine core_shell_quench

  Subroutine core_shell_relax(l_str,relaxed,rdf_collect,rlx_tol,stpcfg,cshell,stat,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for relaxing shells to zero force using conjugate
    ! gradient method
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov & w.smith august 2014
    ! contrib   - a.m.elena february 2017
    ! contrib   - i.t.todorov february 2017
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Logical,            Intent( In    ) :: l_str
    Logical,            Intent( InOut ) :: relaxed,rdf_collect
    Real( Kind = wp ),  Intent( In    ) :: rlx_tol(1:2),stpcfg
    Type(core_shell_type), Intent( InOut ) :: cshell
    Type( stats_type ), Intent( InOut ) :: stat
    Type( comms_type ), Intent( InOut ) :: comm

    Logical,           Save :: newjob = .true. , l_rdf
    Integer,           Save :: keyopt
    Integer                 :: fail(1:2),i,ia,ib,jshl
    Real( Kind = wp ), Save :: grad_tol,eng_tol,dist_tol(1:2),   &
      step,eng,eng0,eng1,eng2,          &
      grad,grad0,grad1,grad2,onorm,sgn, &
      stride,gamma,x(1),y(1),z(1),fff(0:3)

    ! OUTPUT existence

    Logical               :: l_out
    Character( Len = 10 ) :: c_out

    ! Optimisation iteration and convergence limits

    Integer,      Parameter :: mxpass = 100

    Real( Kind = wp ), Save :: grad_pass

    Integer,           Allocatable       :: lstopt(:,:),lst_sh(:)
    Real( Kind = wp ), Allocatable       :: fxt(:),fyt(:),fzt(:)
    Real( Kind = wp ), Allocatable, Save :: oxt(:),oyt(:),ozt(:)
    Character( Len = 256 ) :: message

    fail=0
    Allocate (lstopt(1:2,1:cshell%mxshl),lst_sh(1:mxatms),      Stat=fail(1))
    Allocate (fxt(1:mxatms),fyt(1:mxatms),fzt(1:mxatms), Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_relax allocation failure'
      Call error(0,message)
    End If

    If (newjob) Then
      newjob = .false.

      ! At start no optimisation has been attempted yet

      keyopt = 0

      ! Passage accumulators are initialised in core_shell
      ! stat%passshl(1) - cycles counter
      ! stat%passshl(2) - access counter
      ! stat%passshl(3) - average cycles
      ! stat%passshl(4) - minimum cycles
      ! stat%passshl(5) - maximum cycles

      Allocate (oxt(1:cshell%mxshl),oyt(1:cshell%mxshl),ozt(1:cshell%mxshl), Stat=fail(1))
      If (fail(1) > 0) Then
        Write(message,'(a)') 'core_shell_relax allocation failure SAVE'
        Call error(0,message)
      End If
    End If

    ! Step length for relaxation

    If (rlx_tol(2) > zero_plus) Then

      ! Optionally specified

      step=rlx_tol(2)

    Else

      ! default if unspecified

      step=0.5_wp/cshell%smax

    End If

    If (keyopt == 0) Then

      ! No relaxation is yet attempted

      relaxed=.false.

      ! Minimum needed for a pass for this minimisation cycle

      grad_pass = Huge(1.0_wp)

      ! Avoid rdf calculation redundancy

      l_rdf=rdf_collect
      If (rdf_collect) rdf_collect=.false.

      ! Print header

      If (l_str) Then
        Write(message,'(a,3x,a,6x,a,11x,a,8x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
          'Relaxing shells to cores:','pass','eng_tot','grad_tol','dis_tol','dcs_max','tol=',rlx_tol(1),'step=',step
        Call info(message,.true.)
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)
      End If
    End If

    ! gather new forces on shared shells

    If (cshell%lshmv_shl) Then
      Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,fxx,fyy,fzz,comm)
    End If

    ! Load shell forces on cores (cores don't move during the shell relaxation)

    fxt(1:cshell%mxshl)=0.0_wp ; fyt(1:cshell%mxshl)=0.0_wp ; fzt(1:cshell%mxshl)=0.0_wp
    jshl=0
    Do i=1,cshell%ntshl
      ia=local_index(cshell%listshl(1,i),nlast,lsi,lsa)
      ib=local_index(cshell%listshl(2,i),nlast,lsi,lsa)
      If (ia > 0 .and. ia <= natms) Then ! THERE IS AN ib>0 FOR SURE
        jshl=jshl+1
        fxt(jshl)=fxx(ib)
        fyt(jshl)=fyy(ib)
        fzt(jshl)=fzz(ib)
      End If
      lstopt(1,i)=ia
      lstopt(2,i)=ib
    End Do

    ! Current configuration energy

    eng=stpcfg

    ! Initialise/get eng_tol & verify relaxed condition

    eng_tol=0.0_wp
    If (keyopt > 0) Then
      eng_tol=Abs(1.0_wp-eng2/eng)
    End If
    ! Current gradient (modulus of the total force on shells)

    grad=0.0_wp
    Do i=1,jshl
      grad=grad+fxt(i)**2+fyt(i)**2+fzt(i)**2
    End Do
    Call gsum(comm,grad)
    grad=Sqrt(grad)

    ! Get grad_tol & verify relaxed condition

    grad_tol=grad/Real(cshell%megshl,wp)
    relaxed=(grad_tol < rlx_tol(1))

    ! Initialise dist_tol

    dist_tol=0.0_wp

    ! CHECK FOR CONVERGENCE

    If (.not.relaxed) Then

      ! Increment main passage counter

      stat%passshl(1)=stat%passshl(1)+1.0_wp

      ! Minimum for passing

      grad_pass = Min(grad_pass,grad_tol)

      ! If in mxpass iterations we are not there, give up but
      ! allow for ten-fold boost in iteration cycle length
      ! for the very first MD step

      If (Nint(stat%passshl(2)) == 0) Then
        If (Nint(stat%passshl(1)) >= 10*mxpass) Then
          Call warning(330,rlx_tol(1),grad_pass,0.0_wp)
          Call error(474)
        End If
      Else
        If (Nint(stat%passshl(1)) >= mxpass) Then
          Call warning(330,rlx_tol(1),grad_pass,0.0_wp)
          Call error(474)
        End If
      End If

    Else

      Go To 100

    End If

    If      (keyopt == 0) Then

      ! Original configuration energy

      eng0=eng
      eng2=eng

      ! Original gradient (modulus of the total force on shells)

      onorm=grad
      grad0=grad
      grad2=grad

      ! Set original search direction

      oxt=0.0_wp ; oyt=0.0_wp ; ozt=0.0_wp
      Do i=1,jshl
        oxt(i)=fxt(i)
        oyt(i)=fyt(i)
        ozt(i)=fzt(i)
      End Do

      keyopt=1
      sgn=1.0_wp
      stride=sgn*step

    Else If (keyopt == 1) Then

      ! Line search along chosen direction

      eng1=eng0
      eng2=eng

      grad1=grad2
      grad2=0.0_wp
      Do i=1,jshl
        grad2=grad2+oxt(i)*fxt(i)+oyt(i)*fyt(i)+ozt(i)*fzt(i)
      End Do
      Call gsum(comm,grad2)
      grad2=sgn*grad2/onorm

      ! Linear extrapolation to minimum

      If (grad2 < 0.0_wp) Then ! BACK UP FROM THIS DIRECTION
        keyopt=2
        stride=sgn*step*grad2/(grad1-grad2)
      Else                     ! CARRY ON IN THIS DIRECTION
        stride=sgn*step
      End If

    Else If (keyopt == 2) Then

      ! Construct conjugate search vector

      eng1=eng2
      eng2=eng

      gamma=(grad/grad0)**2
      grad0=grad
      grad2=0.0_wp
      onorm=0.0_wp
      Do i=1,jshl
        oxt(i)=fxt(i)+gamma*oxt(i)
        oyt(i)=fyt(i)+gamma*oyt(i)
        ozt(i)=fzt(i)+gamma*ozt(i)

        onorm=onorm+oxt(i)**2+oyt(i)**2+ozt(i)**2
        grad2=grad2+oxt(i)*fxt(i)+oyt(i)*fyt(i)+ozt(i)*fzt(i)
      End Do
      Call gsum(comm,onorm)
      onorm=Sqrt(onorm)
      Call gsum(comm,grad2)
      grad2=grad2/onorm
      sgn=Sign(1.0_wp,grad2)
      grad2=sgn*grad2

      keyopt=1
      stride=sgn*step

    End If

    ! Load original shell forces on their cores in DD representation

    fxt(1:mxatdm)=0.0_wp ; fyt(1:mxatdm)=0.0_wp ; fzt(1:mxatdm)=0.0_wp
    jshl=0
    Do i=1,cshell%ntshl
      ia=lstopt(1,i)
      If (ia > 0 .and. ia <= natms) Then
        jshl=jshl+1
        fxt(ia)=oxt(jshl)
        fyt(ia)=oyt(jshl)
        fzt(ia)=ozt(jshl)
      End If
    End Do

    ! Exchange original shell forces on shared cores across domains

    If (cshell%lshmv_shl) Then
      Call update_shared_units(natms,nlast,lsi,lsa,cshell%lishp_shl,cshell%lashp_shl,fxt,fyt,fzt,comm)
    End If

    ! Move shells accordingly to their new positions

    Do i=1,cshell%ntshl
      ia=lstopt(1,i)
      ib=lstopt(2,i)
      If (ia > 0 .and. (ib > 0 .and. ib <= natms)) Then
        xxx(ib)=xxx(ib)+stride*fxt(ia)
        yyy(ib)=yyy(ib)+stride*fyt(ia)
        zzz(ib)=zzz(ib)+stride*fzt(ia)
        dist_tol(1)=Max(dist_tol(1),fxt(ia)**2+fyt(ia)**2+fzt(ia)**2) ! - shell move
        x(1)=xxx(ib)-xxx(ia) ; y(1)=yyy(ib)-yyy(ia) ; z(1)=zzz(ib)-zzz(ia)
        Call images(imcon,cell,1,x,y,z)
        dist_tol(2)=Max(dist_tol(2),x(1)**2+y(1)**2+z(1)**2) ! - core-shell separation
      End If
    End Do
    dist_tol=Sqrt(dist_tol)
    dist_tol(1)=dist_tol(1)*Abs(stride)
    Call gmax(comm,dist_tol)

    ! Fit headers in and Close and Open OUTPUT at every 25th print-out

    i=Nint(stat%passshl(1))
    If (l_str) Then
      Write(message,'(1x,i31,1x,1p,2e18.8,4x,f7.4,4x,f7.4,12x,e18.8)') i-1,stpcfg/engunit,grad_tol,dist_tol(1),dist_tol(2),eng_tol
      Call info(message,.true.)
      If (Mod(i,25) == 0) Then
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)
        Write(message,'(1x,a,3x,a,6x,a,11x,a,9x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
          'Relaxing shells to cores:','pass','eng_tot','grad_tol','ds_tol','dcs_max','tol=',rlx_tol(1),'step=',step
        Call info(message,.true.)
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)

        If (comm%idnode == 0) Then
          Inquire(File=Trim(output), Exist=l_out, Position=c_out)
          Call strip_blanks(c_out)
          Call lower_case(c_out)
          If (l_out .and. c_out(1:6) == 'append') Then
            Close(Unit=nrite)
            Open(Unit=nrite, File=Trim(output), Position='append')
          End If
        End If
      End If
    End If

    100 Continue

    If (relaxed) Then

      ! Final printout

      i=Nint(stat%passshl(1))
      If (.not.l_str) Then
        Write(message,'(a,4x,a,6x,a,11x,a,8x,a,4x,a,6x,a,1p,e11.4,3x,a,e11.4)') &
          'Relaxed shells to cores:','pass','eng_tot','grad_tol','dis_tol','dcs_max','tol=',rlx_tol(1),'step=',step
        Call info(message,.true.)
        Write(message,"(1x,130('-'))")
        Call info(message,.true.)
      End If
      Write(message,'(1x,i31,1x,1p,2e18.8,4x,f7.4,4x,f7.4,12x,e18.8)') &
        i-1,stpcfg/engunit,grad_tol,dist_tol(1),dist_tol(2),eng_tol
      Call info(message,.true.)
      Write(message,"(1x,130('-'))")
      Call info(message,.true.)

      ! Collect passage statistics

      stat%passshl(3)=stat%passshl(2)*stat%passshl(3)
      stat%passshl(2)=stat%passshl(2)+1.0_wp
      stat%passshl(3)=stat%passshl(3)/stat%passshl(2)+stat%passshl(1)/stat%passshl(2)
      stat%passshl(4)=Min(stat%passshl(1),stat%passshl(4))
      stat%passshl(5)=Max(stat%passshl(1),stat%passshl(5))

      ! Rewind keyopt and main passage counter

      keyopt =0
      stat%passshl(1)=0.0_wp

      ! Resume rdf calculations

      If (l_rdf) rdf_collect=l_rdf

      ! Zero shells' velocities and forces and redistribute
      ! the residual force to the rest of the system to prevent
      ! COM force generation

      lst_sh(1:natms)=0
      fff(0)=Real(natms,wp)
      fff(1:3)=0.0_wp
      Do i=1,cshell%ntshl
        ib=lstopt(2,i)
        If (ib > 0) Then
          lst_sh(ib)=1
          fff(0)=fff(0)-1.0_wp
          fff(1)=fff(1)+fxx(ib) ; fxx(ib)=0.0_wp ; vxx(ib)=0.0_wp
          fff(2)=fff(2)+fyy(ib) ; fyy(ib)=0.0_wp ; vyy(ib)=0.0_wp
          fff(3)=fff(3)+fzz(ib) ; fzz(ib)=0.0_wp ; vzz(ib)=0.0_wp
        End If
      End Do
      Call gsum(comm,fff)
      fff(1:3)=fff(1:3)/fff(0)
      Do i=1,natms
        If (lst_sh(i) == 0) Then
          fxx(i)=fxx(i)+fff(1)
          fyy(i)=fyy(i)+fff(2)
          fzz(i)=fzz(i)+fff(3)
        End If
      End Do

      ! Frozen atoms option

      Call freeze_atoms()
    End If

    Deallocate (lstopt,lst_sh, Stat=fail(1))
    Deallocate (fxt,fyt,fzt,   Stat=fail(2))
    If (Any(fail > 0)) Then
      Write(message,'(a)') 'core_shell_relax deallocation failure'
      Call error(0,message)
    End If

  End Subroutine core_shell_relax
End module core_shell
