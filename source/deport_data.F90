Module deport_data 

  Use kinds,            Only : wp
  Use comms,            Only : comms_type,gcheck,wp_mpi, Deport_tag, &
                               Export_tag, MetLdExp_tag, ExpMplRM_tag, &
                               PassUnit_tag,gsend,gwait,girecv

  Use setup
  Use domains

  Use configuration

  Use pmf,          Only : ntpmf,listpmf,legpmf

  Use rigid_bodies, Only : ntrgd,listrgd,legrgd, &
                                  q0,q1,q2,q3,          &
                                  rgdvxx,rgdvyy,rgdvzz, &
                                  rgdoxx,rgdoyy,rgdozz, &
                                  legrgd,lshmv_rgd,lishp_rgd,lashp_rgd

  Use tethers,      Only : ntteth,listtet,legtet

  Use bonds,        Only : ntbond,listbnd,legbnd
  Use angles,       Only : ntangl,listang,legang
  Use dihedrals,    Only : ntdihd,listdih,legdih,lx_dih
  Use inversions,   Only : ntinv,listinv,leginv

  Use statistics

  Use minimise,     Only : l_x,oxx,oyy,ozz
  Use langevin,     Only : fxl,fyl,fzl

  Use ewald,               Only : ewald_type
  Use mpole ,              Only : keyind,ltpatm,lchatm, &
                                 induce,indipx,indipy,indipz,rsdx,rsdy,rsdz, &
                                 mplgfr,mprotx,mproty,mprotz, mplflg

  Use msd
  Use greenkubo,    Only : vxi,vyi,vzi,vafsamp

  Use kim,    Only : kimim,idhalo


  Use core_shell,   Only : ntshl, listshl,legshl,lshmv_shl,lishp_shl,lashp_shl

  Use constraints,  Only : ntcons,listcon,legcon,lshmv_con,lishp_con,lashp_con

  Use errors_warnings, Only : error, warning
  Use mpoles_container, Only : rotate_mpoles, rotate_mpoles_d
  Use numerics, Only : local_index
  Use pmf, Only : pmf_units_set
  Use build_book, Only : compress_book_intra
  Use shared_units, Only : pass_shared_units, tag_legend
  Use thermostat, Only : thermostat_type
  Implicit None
  
  Public :: deport_atomic_data, export_atomic_data
  
  Contains 
  

Subroutine deport_atomic_data(mdir,lbook,ewld,thermo,comm)

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

  Logical,            Intent( In    ) :: lbook
  Integer,            Intent( In    ) :: mdir
  Type( ewald_type ), Intent( InOut ) :: ewld
  Type( thermostat_type ), Intent( In    ) :: thermo
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
  Allocate (buffer(1:mxbfdp),                   Stat=fail(1))
  Allocate (lrgd(-1:Max(mxlrgd,mxrgd)),         Stat=fail(2))
  Allocate (ind_on(0:mxatms),ind_off(0:mxatms), Stat=fail(3))
  If (Any(fail > 0)) Then
     Write(message,'(a)') 'deport_atomic_data allocation failure 1'
     Call error(0,message)
  End If

! Set buffer limit (half for outgoing data - half for incoming)

  iblock=mxbfdp/2

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
     lsx = (idx == 0)

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     lex = (idx == nprx-1)

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     lsy = (idy == 0)

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     ley = (idy == npry-1)

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     lsz = (idz == 0)

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     lez = (idz == nprz-1)

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(42)
  End If

! Calculate PBC shift vector due to possible wrap around

  uuu=0.0_wp ; If (lsx) uuu=+1.0_wp ; If (lex) uuu=-1.0_wp
  vvv=0.0_wp ; If (lsy) vvv=+1.0_wp ; If (ley) vvv=-1.0_wp
  www=0.0_wp ; If (lsz) www=+1.0_wp ; If (lez) www=-1.0_wp

  lwrap = (Abs(uuu)+Abs(vvv)+Abs(www) > 0.5_wp)

  If (lwrap) Then
     xadd = cell(1)*uuu+cell(4)*vvv+cell(7)*www
     yadd = cell(2)*uuu+cell(5)*vvv+cell(8)*www
     zadd = cell(3)*uuu+cell(6)*vvv+cell(9)*www
  End If

! Initialise counters for length of sending and receiving buffers
! buffer(1) and buffer(iblock+1) contain the actual number of
! particles to get transfered, imove and jmove are the lengths of
! the buffers

  imove=1
  jmove=1

! Initialise numbers of staying on and leaving off particles

  ind_on(0)=0 ; ind_off(0)=0

! Initialise array overflow flags

  safe=.true.
  safe1=.true.

! Initialise RB data set for preventing duplication

  lrgd(-1:ntrgd)=0

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,natms
     stay=.false. ! the particle is assumed to be leaving

! If the particle is no longer scheduled to leave
! this domain in any direction

     If (ixyz(i) == 0) Then
        stay=.true.
     Else ! If (ixyz(i) > 0) Then ! Get the necessary halo indices
        ix=Mod(ixyz(i),10)           ! [0,1,2]
        iy=Mod(ixyz(i)-ix,100)       ! [0,10,20]
        iz=Mod(ixyz(i)-(ix+iy),1000) ! [0,100,200]

! Filter the move index for the selected direction

        j=ix*kx+iy*ky+iz*kz

! If the particle is scheduled to leave in the selected
! direction then reduce its move index, otherwise tag it as staying

        If (j == jxyz) Then
           ixyz(i)=ixyz(i)-jxyz
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
              buffer(imove+1)=xxx(i)
              buffer(imove+2)=yyy(i)
              buffer(imove+3)=zzz(i)
           Else
              buffer(imove+1)=xxx(i)+xadd
              buffer(imove+2)=yyy(i)+yadd
              buffer(imove+3)=zzz(i)+zadd
           End If

! pack velocities

           buffer(imove+4)=vxx(i)
           buffer(imove+5)=vyy(i)
           buffer(imove+6)=vzz(i)

! pack forces

           buffer(imove+7)=fxx(i)
           buffer(imove+8)=fyy(i)
           buffer(imove+9)=fzz(i)

! pack config indexing, site and move indexing arrays

           buffer(imove+10)=Real(ltg(i),wp)
           buffer(imove+11)=Real(lsite(i),wp)
           buffer(imove+12)=Real(ixyz(i),wp)

! pack initial positions

           buffer(imove+13)=xin(i)
           buffer(imove+14)=yin(i)
           buffer(imove+15)=zin(i)

! pack final displacements

           buffer(imove+16)=xto(i)
           buffer(imove+17)=yto(i)
           buffer(imove+18)=zto(i)

        Else

           safe=.false.

        End If
        imove=imove+18

! pack Langevin forces arrays

        If (thermo%l_langevin) Then
           If (imove+3 <= iblock) Then
              buffer(imove+1)=fxl(i)
              buffer(imove+2)=fyl(i)
              buffer(imove+3)=fzl(i)
           Else
              safe=.false.
           End If
           imove=imove+3
        End If

! pack minimisation arrays

        If (l_x) Then
           If (imove+3 <= iblock) Then
              buffer(imove+1)=oxx(i)
              buffer(imove+2)=oyy(i)
              buffer(imove+3)=ozz(i)
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

        If (vafsamp > 0) Then
           If (imove+3 <= iblock) Then
              Do k=1,vafsamp
                buffer(imove+1)=vxi(i,k)
                buffer(imove+2)=vyi(i,k)
                buffer(imove+3)=vzi(i,k)

                imove=imove+3
              End Do
           Else
              safe=.false.
           End If
        End If

! pack MSD arrays

        If (l_msd) Then
           If (imove+2*(6+mxstak) <= iblock) Then
              jj=27+2*i
              buffer(imove+ 1)=stpvl0(jj-1)
              buffer(imove+ 2)=stpvl0(jj  )
              buffer(imove+ 3)=stpval(jj-1)
              buffer(imove+ 4)=stpval(jj  )
              buffer(imove+ 5)=zumval(jj-1)
              buffer(imove+ 6)=zumval(jj  )
              buffer(imove+ 7)=ravval(jj-1)
              buffer(imove+ 8)=ravval(jj  )
              buffer(imove+ 9)=ssqval(jj-1)
              buffer(imove+10)=ssqval(jj  )
              buffer(imove+11)=sumval(jj-1)
              buffer(imove+12)=sumval(jj  )
              Do kk=1,mxstak
                 l=2*kk   +12
                 buffer(imove+l-1)=stkval(kk,jj-1)
                 buffer(imove+l  )=stkval(kk,jj  )
              End Do
           Else
              safe=.false.
           End If
           imove=imove+2*(6+mxstak)
        End If

! If intra-molecular entities exist in the system

        If (lbook) Then

           If (mximpl > 0) Then ! pack topological array
              kk=ltpatm(0,i)
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
                    buffer(imove)=Real(ltpatm(k,i),wp)
                 End Do
              Else
                 imove=imove+kk
                 safe=.false.
              End If

              If (keyind == 1) Then ! pack CHARMMing core-shell interactions array
                 kk=lchatm(0,i)
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
                       buffer(imove)=Real(lchatm(k,i),wp)
                    End Do
                 Else
                    imove=imove+kk
                    safe=.false.
                 End If
              End If

           End If

! pack the exclusion list

           kk=lexatm(0,i)
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
                 buffer(imove)=Real(lexatm(k,i),wp)
              End Do
           Else
              imove=imove+kk
              safe=.false.
           End If

! the order of packing up intra bookkeeping arrays must be the same as
! their scanning in build_book_intra

! pack core-shell details

           jj=legshl(0,i) ; ii=Sign(1,jj) ; jj=Abs(jj)
           If (jj > 0) Then
              Do ll=1,jj
                 If (imove+3 <= iblock) Then
                    kk=legshl(ll,i)

                    imove=imove+1
                    buffer(imove)=Real(ii*listshl(0,kk),wp) ! negative for a shell particle

                    Do k=1,2
                       imove=imove+1
                       buffer(imove)=Real(listshl(k,kk),wp)
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

           jj=legcon(0,i)
           If (jj > 0) Then
              Do ll=1,jj
                 If (imove+3 <= iblock) Then
                    kk=legcon(ll,i)

                    Do k=0,2
                       imove=imove+1
                       buffer(imove)=Real(listcon(k,kk),wp)
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

           jj=legpmf(0,i)
           If (jj > 0) Then
              Do ll=1,jj
                 If (imove+mxtpmf(1)+mxtpmf(2)+2 <= iblock) Then
                    kk=legpmf(ll,i)

                    Do k=0,mxtpmf(1)
                       imove=imove+1
                       buffer(imove)=Real(listpmf(k,1,kk),wp)
                    End Do

                    Do k=0,mxtpmf(2)
                       imove=imove+1
                       buffer(imove)=Real(listpmf(k,2,kk),wp)
                    End Do
                 Else
                    imove=imove+mxtpmf(1)+mxtpmf(2)+2
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

!           jj=legrgd(0,i)
!           If (jj > 0) Then
!              l=12
!              Do ll=1,jj
!                 If (imove+mxlrgd+l <= iblock) Then
!                    kk=legrgd(ll,i)
!
!                    Do k=-1,listrgd(-1,kk)
!                       imove=imove+1
!                       buffer(imove)=Real(listrgd(k,kk),wp)
!                    End Do
!
!                    buffer(imove+1)=q0(kk)
!                    buffer(imove+2)=q1(kk)
!                    buffer(imove+3)=q2(kk)
!                    buffer(imove+4)=q3(kk)
!                    imove=imove+4
!
!                    buffer(imove+1)=rgdvxx(kk)
!                    buffer(imove+2)=rgdvyy(kk)
!                    buffer(imove+3)=rgdvzz(kk)
!                    imove=imove+3
!
!                    buffer(imove+1)=rgdoxx(kk)
!                    buffer(imove+2)=rgdoyy(kk)
!                    buffer(imove+3)=rgdozz(kk)
!                    imove=imove+3
!                 Else
!                    imove=imove+mxlrgd+l
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

           jj=legrgd(0,i)
           If (jj > 0) Then
              l=12
              Do ll=1,jj
                 If (imove+mxlrgd+l <= iblock) Then
                    kk=legrgd(ll,i)

                    Do k=-1,1
                       imove=imove+1
                       buffer(imove)=Real(listrgd(k,kk),wp)
                    End Do

! Bypass if the data has already been sent

                    If (lrgd(kk) == 0) Then
                       Do k=2,listrgd(-1,kk)
                          imove=imove+1
                          buffer(imove)=Real(listrgd(k,kk),wp)
                       End Do

                       buffer(imove+1)=q0(kk)
                       buffer(imove+2)=q1(kk)
                       buffer(imove+3)=q2(kk)
                       buffer(imove+4)=q3(kk)
                       imove=imove+4

                       buffer(imove+1)=rgdvxx(kk)
                       buffer(imove+2)=rgdvyy(kk)
                       buffer(imove+3)=rgdvzz(kk)
                       imove=imove+3

                       buffer(imove+1)=rgdoxx(kk)
                       buffer(imove+2)=rgdoyy(kk)
                       buffer(imove+3)=rgdozz(kk)
                       imove=imove+3

                       lrgd(kk)=1
                    End If
                 Else
                    imove=imove+mxlrgd+l
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

           jj=legtet(0,i)
           If (jj > 0) Then
              Do ll=1,jj
                 If (imove+2 <= iblock) Then
                    kk=legtet(ll,i)

                    Do k=0,1
                       imove=imove+1
                       buffer(imove)=Real(listtet(k,kk),wp)
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

           jj=legbnd(0,i)
           If (jj > 0) Then
              Do ll=1,jj
                 If (imove+3 <= iblock) Then
                    kk=legbnd(ll,i)

                    Do k=0,2
                       imove=imove+1
                       buffer(imove)=Real(listbnd(k,kk),wp)
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

           jj=legang(0,i)
           If (jj > 0) Then
              Do ll=1,jj
                 If (imove+4 <= iblock) Then
                    kk=legang(ll,i)

                    Do k=0,3
                       imove=imove+1
                       buffer(imove)=Real(listang(k,kk),wp)
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

           jj=legdih(0,i)
           If (jj > 0) Then
              If (.not.lx_dih) Then ! dihedrals only have 4 members
                 l=4
              Else                  ! dihedrals have 4+2 tracked members
                 l=6
              End If
              Do ll=1,jj
                 If (imove+l+1 <= iblock) Then
                    kk=legdih(ll,i)

                    Do k=0,l
                       imove=imove+1
                       buffer(imove)=Real(listdih(k,kk),wp)
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

           jj=leginv(0,i)
           If (jj > 0) Then
              Do ll=1,jj
                 If (imove+5 <= iblock) Then
                    kk=leginv(ll,i)

                    Do k=0,4
                       imove=imove+1
                       buffer(imove)=Real(listinv(k,kk),wp)
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

     xxx(keep)=xxx(i)
     yyy(keep)=yyy(i)
     zzz(keep)=zzz(i)

     vxx(keep)=vxx(i)
     vyy(keep)=vyy(i)
     vzz(keep)=vzz(i)

     fxx(keep)=fxx(i)
     fyy(keep)=fyy(i)
     fzz(keep)=fzz(i)

     ltg(keep)=ltg(i)
     lsite(keep)=lsite(i)
     ixyz(keep)=ixyz(i)

     xin(keep)=xin(i)
     yin(keep)=yin(i)
     zin(keep)=zin(i)

     xto(keep)=xto(i)
     yto(keep)=yto(i)
     zto(keep)=zto(i)

     If (thermo%l_langevin) Then
        fxl(keep)=fxl(i)
        fyl(keep)=fyl(i)
        fzl(keep)=fzl(i)
     End If

     If (l_x) Then
        oxx(keep)=oxx(i)
        oyy(keep)=oyy(i)
        ozz(keep)=ozz(i)
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

     If (vafsamp > 0) Then
        vxi(keep,1:vafsamp)=vxi(i,1:vafsamp)
        vyi(keep,1:vafsamp)=vyi(i,1:vafsamp)
        vzi(keep,1:vafsamp)=vzi(i,1:vafsamp)
     End If

     If (l_msd) Then
        jj=27+2*i
        j =27+2*keep
        stpvl0(j-1)=stpvl0(jj-1)
        stpvl0(j  )=stpvl0(jj  )
        stpval(j-1)=stpval(jj-1)
        stpval(j  )=stpval(jj  )
        zumval(j-1)=zumval(jj-1)
        zumval(j  )=zumval(jj  )
        ravval(j-1)=ravval(jj-1)
        ravval(j  )=ravval(jj  )
        ssqval(j-1)=ssqval(jj-1)
        ssqval(j  )=ssqval(jj  )
        sumval(j-1)=sumval(jj-1)
        sumval(j  )=sumval(jj  )
        Do kk=1,mxstak
           stkval(kk,j-1)=stkval(kk,jj-1)
           stkval(kk,j  )=stkval(kk,jj  )
        End Do
     End If

     If (lbook) Then
        If (mximpl > 0) Then
           ltpatm(:,keep)=ltpatm(:,i)
           If (keyind == 1) lchatm(:,keep)=lchatm(:,i)
        End If

        lexatm(:,keep)=lexatm(:,i)

        legshl(:,keep)=legshl(:,i)

        legcon(:,keep)=legcon(:,i)
        legpmf(:,keep)=legpmf(:,i)

        legrgd(:,keep)=legrgd(:,i)

        legtet(:,keep)=legtet(:,i)

        legbnd(:,keep)=legbnd(:,i)
        legang(:,keep)=legang(:,i)
        legdih(:,keep)=legdih(:,i)
        leginv(:,keep)=leginv(:,i)
     End If
  End Do
  keep=k ! How many particles are to be kept

! record of number of atoms for transfer

  buffer(1)=Real(natms-keep,wp)

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

  natms=keep+jmove

! Check for array bound overflow (can arrays cope with incoming data)

  safe=(natms <= mxatms)
  Call gcheck(comm,safe)
  If (.not.safe) Call error(44)

  Deallocate (ind_on,ind_off,                        Stat=fail(1))
  Allocate   (i1pmf(1:mxtpmf(1)),i2pmf(1:mxtpmf(2)), Stat=fail(2))
  If (Any(fail(1:2) > 0)) Then
     Write(message,'(a)') 'deport_atomic_data de/allocation failure'
     Call error(0,message)
  End If

! load transferred data

  Do i=1,jmove
     newatm=i+keep

! unpack positions

     xxx(newatm)=buffer(kmove+1)
     yyy(newatm)=buffer(kmove+2)
     zzz(newatm)=buffer(kmove+3)

! unpack velocities

     vxx(newatm)=buffer(kmove+4)
     vyy(newatm)=buffer(kmove+5)
     vzz(newatm)=buffer(kmove+6)

! unpack forces

     fxx(newatm)=buffer(kmove+7)
     fyy(newatm)=buffer(kmove+8)
     fzz(newatm)=buffer(kmove+9)

! unpack config indexing, site and move indexing arrays

     ltg(newatm)=Nint(buffer(kmove+10))
     lsite(newatm)=Nint(buffer(kmove+11))
     ixyz(newatm)=Nint(buffer(kmove+12))

! unpack initial positions arrays

     xin(newatm)=buffer(kmove+13)
     yin(newatm)=buffer(kmove+14)
     zin(newatm)=buffer(kmove+15)

! unpack initial positions arrays

     xto(newatm)=buffer(kmove+16)
     yto(newatm)=buffer(kmove+17)
     zto(newatm)=buffer(kmove+18)

     kmove=kmove+18

! unpack Langevin forces arrays

     If (thermo%l_langevin) Then
        fxl(newatm)=buffer(kmove+1)
        fyl(newatm)=buffer(kmove+2)
        fzl(newatm)=buffer(kmove+3)

        kmove=kmove+3
     End If

! unpack minimisation arrays

     If (l_x) Then
        oxx(newatm)=buffer(kmove+1)
        oyy(newatm)=buffer(kmove+2)
        ozz(newatm)=buffer(kmove+3)

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

     If (vafsamp > 0) Then
        Do k=1,vafsamp
          vxi(newatm,k)=buffer(kmove+1)
          vyi(newatm,k)=buffer(kmove+2)
          vzi(newatm,k)=buffer(kmove+3)

          kmove=kmove+3
        End Do
     End If

! unpack MSD arrays

     If (l_msd) Then
        jj=27+2*newatm
        stpvl0(jj-1)=buffer(kmove+1 )
        stpvl0(jj  )=buffer(kmove+2 )
        stpval(jj-1)=buffer(kmove+3 )
        stpval(jj  )=buffer(kmove+4 )
        zumval(jj-1)=buffer(kmove+5 )
        zumval(jj  )=buffer(kmove+6 )
        ravval(jj-1)=buffer(kmove+7 )
        ravval(jj  )=buffer(kmove+8 )
        ssqval(jj-1)=buffer(kmove+9 )
        ssqval(jj  )=buffer(kmove+10)
        sumval(jj-1)=buffer(kmove+11)
        sumval(jj  )=buffer(kmove+12)
        Do kk=1,mxstak
           l=2*kk                +12
           stkval(kk,jj-1)=buffer(kmove+l-1)
           stkval(kk,jj  )=buffer(kmove+l  )
        End Do

        kmove=kmove+2*(6+mxstak)
     End If

     If (lbook) Then

        If (mximpl > 0) Then ! unpack topological array
           kmove=kmove+1
           kk=Nint(buffer(kmove))
           ltpatm(0,newatm)=kk
           Do k=1,kk
              kmove=kmove+1
              ltpatm(k,newatm)=Nint(buffer(kmove))
           End Do
           ltpatm(kk+1:mxexcl,newatm)=0

           If (keyind == 1) Then ! unpack CHARMMing core-shell interactions array
              kmove=kmove+1
              kk=Nint(buffer(kmove))
              lchatm(0,newatm)=kk
              Do k=1,kk
                 kmove=kmove+1
                 lchatm(k,newatm)=Nint(buffer(kmove))
              End Do
              lchatm(kk+1:mxexcl,newatm)=0
           End If
        End If

! unpack the exclusion list

        kmove=kmove+1
        kk=Nint(buffer(kmove))
        lexatm(0,newatm)=kk
        Do k=1,kk
           kmove=kmove+1
           lexatm(k,newatm)=Nint(buffer(kmove))
        End Do
        lexatm(kk+1:mxexcl,newatm)=0

! the order of unpacking intra bookkeeping arrays must be the same as
! the order of their scanning in build_book_intra in order to rebuild
! correctly the new list arrays and create new legend arrays

! set initial intra counters

        jshels=ntshl

        jconst=ntcons
        jpmf  =ntpmf

        jrigid=ntrgd

        jteths=ntteth

        jbonds=ntbond
        jangle=ntangl
        jdihed=ntdihd
        jinver=ntinv

! unpack core-shell details

        legshl(:,newatm) = 0
        Do While (Abs(buffer(kmove+1)) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1)) ; ll=Sign(1,jj) ; jj=Abs(jj)
           iatm=Nint(buffer(kmove+2)) ! ll=1
           jatm=Nint(buffer(kmove+3)) ! ll=-1
           kmove=kmove+3

! check if core-shell unit already specified

           kshels=0
           check=.true.
           Do While (check .and. kshels < Min(jshels,mxshl))
              kshels=kshels+1
              check=.not.( jj   == listshl(0,kshels) .and. & ! core-shell units don't intersect
                           iatm == listshl(1,kshels))        ! .and. jatm == listshl(2,kshels) )
           End Do

! add new core-shell unit

           If (check) Then
              jshels=jshels+1

              If (jshels <= mxshl) Then
                 listshl(0,jshels)=jj
                 listshl(1,jshels)=iatm
                 listshl(2,jshels)=jatm

                 Call tag_legend(safe1,newatm,ll*jshels,legshl,mxfshl)
              Else
                 safe=.false.
                 Write(message,'(a)') "too many core-shell units"
                 Call warning(message)

              End If
           Else
              Call tag_legend(safe1,newatm,ll*kshels,legshl,mxfshl)
           End If
        End Do
        kmove=kmove+1

! unpack bond constraint details

        legcon(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           kmove=kmove+3

! check if constraint already specified

           kconst=0
           check=.true.
           Do While (check .and. kconst < Min(jconst,mxcons))
              kconst=kconst+1
              check=.not.( jj   == listcon(0,kconst) .and. &
                           iatm == listcon(1,kconst) .and. &
                           jatm == listcon(2,kconst) )
           End Do

! insert new constraint unit

           If (check) Then
              jconst=jconst+1

              If (jconst <= mxcons) Then
                 listcon(0,jconst)=jj
                 listcon(1,jconst)=iatm
                 listcon(2,jconst)=jatm

                 Call tag_legend(safe1,newatm,jconst,legcon,mxfcon)
              Else
                 safe=.false.
                 Call warning('too many constraint units')
              End If
           Else
              Call tag_legend(safe1,newatm,kconst,legcon,mxfcon)
           End If
        End Do
        kmove=kmove+1

! unpack PMF details

        legpmf(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1)) ! PMF global identifier
           kmove=kmove+1
           Do k=1,mxtpmf(1)
              kmove=kmove+1
              i1pmf(k)=Nint(buffer(kmove))
           End Do

           kmove=kmove+1            ! omit PMF units presence identifier
           Do k=1,mxtpmf(2)         ! and deal with it in pmf_units_set
              kmove=kmove+1
              i2pmf(k)=Nint(buffer(kmove))
           End Do

! check if PMF already specified

           kpmf=0
           check=.true.
           Do While (check .and. kpmf < Min(jpmf,mxpmf))
              kpmf=kpmf+1
              check=.not.(jj == listpmf(0,1,kpmf))
           End Do

! insert new PMF constraint

           If (check) Then
              jpmf=jpmf+1

              If (jpmf <= mxpmf) Then
                 listpmf(0,1,jpmf)=jj
                 Do k=1,mxtpmf(1)
                    listpmf(k,1,jpmf)=i1pmf(k)
                 End Do

! PMF units presence identifier holds zero temporarily
! it's dealt with in pmf_units_set

                 listpmf(0,2,jpmf)=0
                 Do k=1,mxtpmf(2)
                    listpmf(k,2,jpmf)=i2pmf(k)
                 End Do

                 Call tag_legend(safe1,newatm,jpmf,legpmf,mxfpmf)
              Else
                 safe=.false.
                 Call warning('too many PMF units')
              End If
           Else
              Call tag_legend(safe1,newatm,kpmf,legpmf,mxfpmf)
           End If
        End Do
        kmove=kmove+1

! unpack RB details

! Initialise RB data set for preventing duplication

        lrgd(-1:mxlrgd)=0

!        legrgd(:,newatm) = 0
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
!           Do While (check .and. krigid < Min(jrigid,mxrgd))
!              krigid=krigid+1
!              check=( lrgd( 0) == listrgd( 0,krigid) .and. & ! Type
!                      lrgd(-1) == listrgd(-1,krigid) .and. & ! Size
!                      lrgd( 1) == listrgd( 1,krigid) )       ! Global ID of just the first
!!                                                            ! member as RBs don't intersect
!              check=.not.check
!           End Do
!
!! insert new RB unit
!
!           If (check) Then
!              jrigid=jrigid+1
!
!              If (jrigid <= mxrgd) Then
!                 Do k=-1,lrgd(-1)
!                    listrgd(k,jrigid)=lrgd(k)
!                 End Do
!
!                 Call tag_legend(safe1,newatm,jrigid,legrgd,mxfrgd)
!
!                 q0(jrigid)=buffer(kmove+1)
!                 q1(jrigid)=buffer(kmove+2)
!                 q2(jrigid)=buffer(kmove+3)
!                 q3(jrigid)=buffer(kmove+4)
!                 kmove=kmove+4
!
!                 rgdvxx(jrigid)=buffer(kmove+1)
!                 rgdvyy(jrigid)=buffer(kmove+2)
!                 rgdvzz(jrigid)=buffer(kmove+3)
!                 kmove=kmove+3
!
!                 rgdoxx(jrigid)=buffer(kmove+1)
!                 rgdoyy(jrigid)=buffer(kmove+2)
!                 rgdozz(jrigid)=buffer(kmove+3)
!                 kmove=kmove+3
!              Else
!                 safe=.false.
!              End If
!           Else
!              l=10
!              kmove=kmove+l ! Compensate for the 'l' extra unread buffers
!              Call tag_legend(safe1,newatm,krigid,legrgd,mxfrgd)
!           End If
!        End Do
!        kmove=kmove+1

        legrgd(:,newatm) = 0
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
           Do While (check .and. krigid < Min(jrigid,mxrgd))
              krigid=krigid+1
              check=( lrgd( 0) == listrgd( 0,krigid) .and. & ! Type
                      lrgd(-1) == listrgd(-1,krigid) .and. & ! Size
                      lrgd( 1) == listrgd( 1,krigid) )       ! Global ID of just the first
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

              If (jrigid <= mxrgd) Then
                 Do k=-1,lrgd(-1)
                    listrgd(k,jrigid)=lrgd(k)
                 End Do

                 Call tag_legend(safe1,newatm,jrigid,legrgd,mxfrgd)

                 q0(jrigid)=buffer(kmove+1)
                 q1(jrigid)=buffer(kmove+2)
                 q2(jrigid)=buffer(kmove+3)
                 q3(jrigid)=buffer(kmove+4)
                 kmove=kmove+4

                 rgdvxx(jrigid)=buffer(kmove+1)
                 rgdvyy(jrigid)=buffer(kmove+2)
                 rgdvzz(jrigid)=buffer(kmove+3)
                 kmove=kmove+3

                 rgdoxx(jrigid)=buffer(kmove+1)
                 rgdoyy(jrigid)=buffer(kmove+2)
                 rgdozz(jrigid)=buffer(kmove+3)
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
              Call tag_legend(safe1,newatm,krigid,legrgd,mxfrgd)
           End If
        End Do
        kmove=kmove+1

! unpack tether details

        legtet(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           kmove=kmove+2

! check if tether unit already specified

           kteths=0
           check=.true.
           Do While (check .and. kteths < Min(jteths,mxteth))
              kteths=kteths+1
              check=.not.( jj   == listtet(0,kteths) .and. &
                           iatm == listtet(1,kteths) )
           End Do

! add new tether unit

           If (check) Then
              jteths=jteths+1

              If (jteths <= mxteth) Then
                 listtet(0,jteths)=jj
                 listtet(1,jteths)=iatm

                 Call tag_legend(safe1,newatm,jteths,legtet,mxftet)
              Else
                 safe=.false.
                 Call warning('too many tether units')
              End If
           Else
              Call tag_legend(safe1,newatm,kteths,legtet,mxftet)
           End If
        End Do
        kmove=kmove+1

! unpack bond details

        legbnd(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           kmove=kmove+3

! check if bond already specified

           kbonds=0
           check=.true.
           Do While (check .and. kbonds < Min(jbonds,mxbond))
              kbonds=kbonds+1
              check=.not.( jj   == listbnd(0,kbonds) .and. &
                           iatm == listbnd(1,kbonds) .and. &
                           jatm == listbnd(2,kbonds) )
           End Do

! insert new bond details

           If (check) Then
              jbonds=jbonds+1

              If (jbonds <= mxbond) Then
                 listbnd(0,jbonds)=jj
                 listbnd(1,jbonds)=iatm
                 listbnd(2,jbonds)=jatm

                 Call tag_legend(safe1,newatm,jbonds,legbnd,mxfbnd)
              Else
                 safe=.false.
                 Call warning('too many bond units')
              End If
           Else
              Call tag_legend(safe1,newatm,kbonds,legbnd,mxfbnd)
           End If
        End Do
        kmove=kmove+1

! unpack valence angle details

        legang(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           katm=Nint(buffer(kmove+4))
           kmove=kmove+4

! check if angle already specified

           kangle=0
           check=.true.
           Do While (check .and. kangle < Min(jangle,mxangl))
              kangle=kangle+1
              check=.not.( jj   == listang(0,kangle) .and. &
                           iatm == listang(1,kangle) .and. &
                           jatm == listang(2,kangle) .and. &
                           katm == listang(3,kangle) )
           End Do

! insert new angle details

           If (check) Then
              jangle=jangle+1

              If (jangle <= mxangl) Then
                 listang(0,jangle)=jj
                 listang(1,jangle)=iatm
                 listang(2,jangle)=jatm
                 listang(3,jangle)=katm

                 Call tag_legend(safe1,newatm,jangle,legang,mxfang)
              Else
                 safe=.false.
                 Call warning('too many angle units')
              End If
           Else
              Call tag_legend(safe1,newatm,kangle,legang,mxfang)
           End If
        End Do
        kmove=kmove+1

! unpack dihedral angle details

        legdih(:,newatm) = 0
        Do While (buffer(kmove+1) > 0.0_wp .and. safe)
           jj=Nint(buffer(kmove+1))
           iatm=Nint(buffer(kmove+2))
           jatm=Nint(buffer(kmove+3))
           katm=Nint(buffer(kmove+4))
           latm=Nint(buffer(kmove+5))
           If (lx_dih) Then
              matm=Nint(buffer(kmove+6))
              natm=Nint(buffer(kmove+7))
              kmove=kmove+7
           Else
              kmove=kmove+5
           End If

! check if dihedral already specified

           kdihed=0
           check=.true.
           Do While (check .and. kdihed < Min(jdihed,mxdihd))
              kdihed=kdihed+1
              check=.not.( jj   == listdih(0,kdihed) .and. &
                           iatm == listdih(1,kdihed) .and. &
                           jatm == listdih(2,kdihed) .and. &
                           katm == listdih(3,kdihed) .and. &
                           latm == listdih(4,kdihed) )
           End Do

! add new dihedral details

           If (check) Then
              jdihed=jdihed+1

              If (jdihed <= mxdihd) Then
                 listdih(0,jdihed)=jj
                 listdih(1,jdihed)=iatm
                 listdih(2,jdihed)=jatm
                 listdih(3,jdihed)=katm
                 listdih(4,jdihed)=latm
                 If (lx_dih) Then
                    listdih(5,jdihed)=matm
                    listdih(6,jdihed)=natm
                 End If

                 Call tag_legend(safe1,newatm,jdihed,legdih,mxfdih)
              Else
                 safe=.false.
                 Call warning('too many dihedral units')
              End If
           Else
              Call tag_legend(safe1,newatm,kdihed,legdih,mxfdih)
           End If
        End Do
        kmove=kmove+1

! unpack inversion angle details

        leginv(:,newatm) = 0
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
           Do While (check .and. kinver < Min(jinver,mxinv))
              kinver=kinver+1
              check=.not.( jj   == listinv(0,kinver) .and. &
                           iatm == listinv(1,kinver) .and. &
                           jatm == listinv(2,kinver) .and. &
                           katm == listinv(3,kinver) .and. &
                           latm == listinv(4,kinver) )
           End Do

! add new inversion details

           If (check) Then
              jinver=jinver+1

              If (jinver <= mxinv) Then
                 listinv(0,jinver)=jj
                 listinv(1,jinver)=iatm
                 listinv(2,jinver)=jatm
                 listinv(3,jinver)=katm
                 listinv(4,jinver)=latm

                 Call tag_legend(safe1,newatm,jinver,leginv,mxfinv)
              Else
                 safe=.false.
                 Call warning('too many inversion units')
              End If
           Else
              Call tag_legend(safe1,newatm,kinver,leginv,mxfinv)
           End If
        End Do
        kmove=kmove+1

! redefine intra counters

        ntshl =jshels

        ntcons=jconst
        ntpmf =jpmf

        ntrgd =jrigid

        ntteth=jteths

        ntbond=jbonds
        ntangl=jangle
        ntdihd=jdihed
        ntinv =jinver

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

Subroutine export_atomic_data(mdir,comm)

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
  Type( comms_type ), Intent( InOut ) :: comm 

  Logical           :: safe,lsx,lsy,lsz,lex,ley,lez,lwrap
  Integer           :: fail,iadd,limit,iblock,            &
                       i,j,k,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
                       jdnode,kdnode,imove,jmove,kmove,itmp
  Real( Kind = wp ) :: uuu,vvv,www,xadd,yadd,zadd

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  Character ( Len = 256 )   ::  message

! Number of transported quantities per particle

  iadd=6+Merge(6,0,induce)

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
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
     lsx = (idx == 0)

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     kxyz= 3
     lex = (idx == nprx-1)

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     kxyz= 30
     lsy = (idy == 0)

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     kxyz= 30
     ley = (idy == npry-1)

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     kxyz= 300
     lsz = (idz == 0)

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     kxyz= 300
     lez = (idz == nprz-1)

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(46)
  End If

! Calculate PBC shift vector due to possible wrap around

  uuu=0.0_wp ; If (lsx) uuu=+1.0_wp ; If (lex) uuu=-1.0_wp
  vvv=0.0_wp ; If (lsy) vvv=+1.0_wp ; If (ley) vvv=-1.0_wp
  www=0.0_wp ; If (lsz) www=+1.0_wp ; If (lez) www=-1.0_wp

  lwrap = (Abs(uuu)+Abs(vvv)+Abs(www) > 0.5_wp)

  If (lwrap) Then
     xadd = cell(1)*uuu+cell(4)*vvv+cell(7)*www
     yadd = cell(2)*uuu+cell(5)*vvv+cell(8)*www
     zadd = cell(3)*uuu+cell(6)*vvv+cell(9)*www
  End If

! Initialise counters for length of sending and receiving buffers
! imove and jmove are the actual number of particles to get haloed

  imove=0
  jmove=0

! Initialise array overflow flags

  safe=.true.

! LOOP OVER ALL PARTICLES ON THIS NODE

  Do i=1,nlast

! If the particle is within the remaining 'inverted halo' of this domain

     If (ixyz(i) > 0) Then

! Get the necessary halo indices

        ix=Mod(ixyz(i),10)           ! [0,1,2,3=1+2]
        iy=Mod(ixyz(i)-ix,100)       ! [0,10,20,30=10+20]
        iz=Mod(ixyz(i)-(ix+iy),1000) ! [0,100,200,300=100+200]

! Filter the halo index for the selected direction

        j=ix*kx+iy*ky+iz*kz

! If the particle is within the correct halo for the selected direction

        If (j == jxyz .or. (j > jxyz .and. Mod(j,3) == 0)) Then

! If safe to proceed

           If ((imove+iadd) <= iblock) Then

! pack positions and apply possible PBC shift for the receiver

              If (.not.lwrap) Then
                 buffer(imove+1)=xxx(i)
                 buffer(imove+2)=yyy(i)
                 buffer(imove+3)=zzz(i)
              Else
                 buffer(imove+1)=xxx(i)+xadd
                 buffer(imove+2)=yyy(i)+yadd
                 buffer(imove+3)=zzz(i)+zadd
              End If

! pack config indexing, site and remaining halo indexing arrays

              buffer(imove+4)=Real(ltg(i),wp)
              buffer(imove+5)=Real(lsite(i),wp)

! pack induced dipoles data

              If (induce) Then
                 kmove=imove+5

                 buffer(kmove+1)=indipx(i)
                 buffer(kmove+2)=indipy(i)
                 buffer(kmove+3)=indipz(i)
                 buffer(kmove+4)=rsdx(i)
                 buffer(kmove+5)=rsdy(i)
                 buffer(kmove+6)=rsdz(i)
              End If

! Use the corrected halo reduction factor when the particle is halo to both +&- sides

              buffer(imove+iadd)=Real(ixyz(i)-Merge(jxyz,kxyz,j == jxyz),wp)

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

  safe=((nlast+jmove/iadd) <= mxatms)
  Call gcheck(comm,safe)
  If (.not.safe) Then
     itmp=nlast+jmove/iadd
     Call gmax(comm,itmp)
     Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
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

  If (kimim /= ' ') Then
     i = Abs(2*mdir)+Sign(mdir,1) ! Merge( 2*mdir , -2*mdir-1 , mdir > 0 )
     idhalo(0,i)=imove/iadd       ! atoms to send
     idhalo(1,i)=nlast+1          ! first atom to receive
     idhalo(2,i)=nlast+jmove/iadd ! last atom to receive
  End If

! load transferred data

  j=Merge(iblock,0,comm%mxnode > 1)
  Do i=1,jmove/iadd
     nlast=nlast+1

! unpack positions

     xxx(nlast)=buffer(j+1)
     yyy(nlast)=buffer(j+2)
     zzz(nlast)=buffer(j+3)

! unpack config indexing, site and halo indexing arrays

     ltg(nlast)  =Nint(buffer(j+4))
     lsite(nlast)=Nint(buffer(j+5))

! unpack induced dipoles data

     If (induce) Then
        k=j+5

        indipx(nlast)=buffer(k+1)
        indipy(nlast)=buffer(k+2)
        indipz(nlast)=buffer(k+3)
        rsdx(nlast)  =buffer(k+4)
        rsdy(nlast)  =buffer(k+5)
        rsdz(nlast)  =buffer(k+6)
     End If

! unpack remaining halo indexing

     ixyz(nlast) =Nint(buffer(j+iadd))

     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'export_atomic_data deallocation failure'
     Call error(0,message)
  End If

End Subroutine export_atomic_data

Subroutine export_atomic_positions(mdir,mlast,ixyz0,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export atomic positions in domain boundary regions
! for halo refresh
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer, Intent( In    ) :: mdir,ixyz0(1:mxatms)
  Integer, Intent( InOut ) :: mlast
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

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
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
     lsx = (idx == 0)

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     lex = (idx == nprx-1)

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     lsy = (idy == 0)

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     ley = (idy == npry-1)

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     lsz = (idz == 0)

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     lez = (idz == nprz-1)

     jdnode = map(6)
     kdnode = map(5)
  Else
     Call error(46)
  End If

! Calculate PBC shift vector due to possible wrap around

  uuu=0.0_wp ; If (lsx) uuu=+1.0_wp ; If (lex) uuu=-1.0_wp
  vvv=0.0_wp ; If (lsy) vvv=+1.0_wp ; If (ley) vvv=-1.0_wp
  www=0.0_wp ; If (lsz) www=+1.0_wp ; If (lez) www=-1.0_wp

  lwrap = (Abs(uuu)+Abs(vvv)+Abs(www) > 0.5_wp)

  If (lwrap) Then
     xadd = cell(1)*uuu+cell(4)*vvv+cell(7)*www
     yadd = cell(2)*uuu+cell(5)*vvv+cell(8)*www
     zadd = cell(3)*uuu+cell(6)*vvv+cell(9)*www
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
                 buffer(imove+1)=xxx(i)
                 buffer(imove+2)=yyy(i)
                 buffer(imove+3)=zzz(i)
              Else
                 buffer(imove+1)=xxx(i)+xadd
                 buffer(imove+2)=yyy(i)+yadd
                 buffer(imove+3)=zzz(i)+zadd
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

  safe=((mlast+jmove/iadd) <= mxatms)
  Call gcheck(comm,safe)
  If (.not.safe) Then
     itmp=mlast+jmove/iadd
     Call gmax(comm,itmp)
     Call warning(160,Real(itmp,wp),Real(mxatms,wp),0.0_wp)
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

  If (kimim /= ' ') Then
     i = Abs(2*mdir)+Sign(mdir,1) ! Merge( 2*mdir , -2*mdir-1 , mdir > 0 )
     idhalo(0,i)=imove/iadd       ! atoms to send
     idhalo(1,i)=mlast+1          ! first atom to receive
     idhalo(2,i)=mlast+jmove/iadd ! last atom to receive
  End If

! load transferred data

  j=Merge(iblock,0,comm%mxnode > 1)
  Do i=1,jmove/iadd
     mlast=mlast+1

! unpack positions

     xxx(mlast)=buffer(j+1)
     yyy(mlast)=buffer(j+2)
     zzz(mlast)=buffer(j+3)

     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'export_atomic_positions deallocation failure'
     Call error(0,message)
  End If

End Subroutine export_atomic_positions



Subroutine mpoles_rotmat_export(mdir,mlast,ixyz0,comm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 routine to export multipoles rotation and infinitesimally
! rotated matrices in the halo
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Integer, Intent( In    ) :: mdir
  Integer, Intent( InOut ) :: mlast,ixyz0(1:mxatms)
  Type( comms_type ), Intent( InOut ) :: comm


  Logical           :: safe
  Integer           :: fail,iadd,limit,iblock,          &
                       i,j,jxyz,kxyz,ix,iy,iz,kx,ky,kz, &
                       idl1,idl2,idl3,idl4,             &
                       jdnode,kdnode,imove,jmove,itmp

  Real( Kind = wp ), Dimension( : ), Allocatable :: buffer

  Character  ( Len =  256 )  ::  message

! Number of transported quantities per particle

  iadd=4*mximpl+1
  idl1=mximpl ; idl2=2*mximpl ; idl3=3*mximpl ; idl4=4*mximpl

  fail=0 ; limit=iadd*mxbfxp ! limit=Merge(1,2,mxnode > 1)*iblock*iadd
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

     jdnode = map(1)
     kdnode = map(2)
  Else If (mdir ==  1) Then ! Direction +x
     kx  = 1
     jxyz= 2
     kxyz= 3

     jdnode = map(2)
     kdnode = map(1)
  Else If (mdir == -2) Then ! Direction -y
     ky  = 1
     jxyz= 10
     kxyz= 30

     jdnode = map(3)
     kdnode = map(4)
  Else If (mdir ==  2) Then ! Direction +y
     ky  = 1
     jxyz= 20
     kxyz= 30

     jdnode = map(4)
     kdnode = map(3)
  Else If (mdir == -3) Then ! Direction -z
     kz  = 1
     jxyz= 100
     kxyz= 300

     jdnode = map(5)
     kdnode = map(6)
  Else If (mdir ==  3) Then ! Direction +z
     kz  = 1
     jxyz= 200
     kxyz= 300

     jdnode = map(6)
     kdnode = map(5)
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

              buffer(imove+1        : imove + idl1) = mplgfr(:,i)
              buffer(imove+1 + idl1 : imove + idl2) = mprotx(:,i)
              buffer(imove+1 + idl2 : imove + idl3) = mproty(:,i)
              buffer(imove+1 + idl3 : imove + idl4) = mprotz(:,i)

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

     mplgfr(:,mlast) = buffer(j+1        : j + idl1)
     mprotx(:,mlast) = buffer(j+1 + idl1 : j + idl2)
     mproty(:,mlast) = buffer(j+1 + idl2 : j + idl3)
     mprotz(:,mlast) = buffer(j+1 + idl3 : j + idl4)

     ixyz0(mlast)=Nint(buffer(j+iadd))
     j=j+iadd
  End Do

  Deallocate (buffer, Stat=fail)
  If (fail > 0) Then
     Write(message,'(a)') 'mpoles_rotmat_export deallocation failure'
     Call error(0,message)
  End If

End Subroutine mpoles_rotmat_export


Subroutine mpoles_rotmat_set_halo(comm)

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

  Type( comms_type), Intent( InOut ) :: comm

  Logical :: safe
  Integer :: fail,i,mlast

  Integer, Allocatable :: ixyz0(:)

  Character ( Len = 256 )  ::  message

  Do i=1,natms
     mplflg(i)=0
     If (mxompl < 3) Then
        Call rotate_mpoles_d(i,comm)
     Else
        Call rotate_mpoles(i,comm)
     End If
  End Do

! Communicate the matrices in the halo

  fail = 0
  Allocate (ixyz0(1:mxatms), Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'mpoles_rotmat_set_halo allocation failure'
     Call error(0,message)
  End If
  ixyz0(1:nlast) = ixyz(1:nlast)

! No halo, start with domain only particles

  mlast=natms

! exchange atom data in -/+ x directions

  Call mpoles_rotmat_export(-1,mlast,ixyz0,comm)
  Call mpoles_rotmat_export( 1,mlast,ixyz0,comm)

! exchange atom data in -/+ y directions

  Call mpoles_rotmat_export(-2,mlast,ixyz0,comm)
  Call mpoles_rotmat_export( 2,mlast,ixyz0,comm)

! exchange atom data in -/+ z directions

  Call mpoles_rotmat_export(-3,mlast,ixyz0,comm)
  Call mpoles_rotmat_export( 3,mlast,ixyz0,comm)

! check atom totals after data transfer

  safe=(mlast == nlast)
  Call gcheck(comm,safe)
  If (.not.safe) Call error(174)

  Deallocate (ixyz0, Stat = fail)
  If (fail > 0) Then
     Write(message,'(a)') 'mpoles_rotmat_set_halo deallocation failure'
     Call error(0,message)
  End If

End Subroutine mpoles_rotmat_set_halo

Subroutine relocate_particles       &
           (dvar,rlnk,lbook,megatm, &
           megshl,m_con,megpmf,     &
           m_rgd,megtet,            &
           megbnd,megang,megdih,    &
           meginv,ewld,thermo,comm)

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


  Logical,           Intent( In    ) :: lbook
  Integer,           Intent( In    ) :: megatm,              &
                                        megshl,m_con,megpmf, &
                                        m_rgd,megtet,        &
                                        megbnd,megang,megdih,meginv
  Real( Kind = wp ), Intent( In    ) :: dvar,rlnk
  Type( ewald_type ), Intent( InOut ) :: ewld
  Type( thermostat_type ), Intent( In    ) :: thermo
  Type( comms_type ), Intent( InOut ) :: comm
  Real( Kind = wp ), Save :: cut

  Logical           :: safe(1:9)
  Integer           :: i,nlimit,ipx,ipy,ipz,itmp(1:9),jtmp(1:9)
  Real( Kind = wp ) :: big(1:3),det,celprp(1:10),rcell(1:9),x,y,z,tmp

  Character( Len = 256 ) :: message

! Define cut

  cut=rlnk+1.0e-6_wp

! rescale mock cell vectors for non-periodic system

  If (imcon == 0 .or. imcon == 6) Then

! find maximum x,y,z positions

     big=0.0_wp

     Do i =1,natms
        big(1)=Max(big(1),Abs(xxx(i)))
        big(2)=Max(big(2),Abs(yyy(i)))
        big(3)=Max(big(3),Abs(zzz(i)))
     End Do

     Call gmax(comm,big)

     If (imcon == 0) Then

        cell(1)=Max(2.0_wp*big(1)+cut,3.0_wp*cut,cell(1))
        cell(5)=Max(2.0_wp*big(2)+cut,3.0_wp*cut,cell(5))
        cell(9)=Max(2.0_wp*big(3)+cut,3.0_wp*cut,cell(9))

        cell(2)=0.0_wp
        cell(3)=0.0_wp
        cell(4)=0.0_wp
        cell(6)=0.0_wp
        cell(7)=0.0_wp
        cell(8)=0.0_wp

     Else If (imcon == 6) Then

        cell(9)=Max(2.0_wp*big(3)+cut,3.0_wp*cut,cell(9))

     End If

  End If

! Get the dimensional properties of the MD cell

  Call dcell(cell,celprp)

! Only when we have more than one domain we need real relocation

  If (comm%mxnode > 1) Then

     Call invert(cell,rcell,det)

! Convert atomic positions from MD centred
! Cartesian coordinates to reduced space ones.
! Populate the move (former halo) indicator array.
! Here we assume no particle has moved more than
! a link-cell width (> rlnk) in any direction
! since last call to relocate as we don't test this!!!

     ixyz(1:natms)=0 ! Initialise move (former halo) indicator
     Do i=1,natms
        x=rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i)
        y=rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i)
        z=rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i)

! assign domain coordinates (call for errors)

        ipx=Int((x+0.5_wp)*nprx_r)
        ipy=Int((y+0.5_wp)*npry_r)
        ipz=Int((z+0.5_wp)*nprz_r)

        If (idx == 0) Then
           If (x < -half_plus) ixyz(i)=ixyz(i)+1
        Else
           If (ipx < idx) ixyz(i)=ixyz(i)+1
        End If
        If (idx == nprx-1) Then
           If (x >= half_minus) ixyz(i)=ixyz(i)+2
        Else
           If (ipx > idx) ixyz(i)=ixyz(i)+2
        End If

        If (idy == 0) Then
           If (y < -half_plus) ixyz(i)=ixyz(i)+10
        Else
           If (ipy < idy) ixyz(i)=ixyz(i)+10
        End If
        If (idy == npry-1) Then
           If (y >= half_minus) ixyz(i)=ixyz(i)+20
        Else
           If (ipy > idy) ixyz(i)=ixyz(i)+20
        End If

        If (idz == 0) Then
           If (z < -half_plus) ixyz(i)=ixyz(i)+100
        Else
           If (ipz < idz) ixyz(i)=ixyz(i)+100
        End If
        If (idz == nprz-1) Then
           If (z >= half_minus) ixyz(i)=ixyz(i)+200
        Else
           If (ipz > idz) ixyz(i)=ixyz(i)+200
        End If
     End Do

! exchange atom data in -/+ x directions

     Call deport_atomic_data(-1,lbook,ewld,thermo,comm)
     Call deport_atomic_data( 1,lbook,ewld,thermo,comm)

! exchange atom data in -/+ y directions

     Call deport_atomic_data(-2,lbook,ewld,thermo,comm)
     Call deport_atomic_data( 2,lbook,ewld,thermo,comm)

! exchange atom data in -/+ z directions

     Call deport_atomic_data(-3,lbook,ewld,thermo,comm)
     Call deport_atomic_data( 3,lbook,ewld,thermo,comm)

! check system for loss of atoms

     safe(1)=(All(ixyz(1:natms) == 0)) ; Call gcheck(comm,safe(1),"enforce")
     nlimit=natms ; Call gsum(comm,nlimit)
     If ((.not.safe(1)) .or. nlimit /= megatm) Call error(58)

! reassign atom properties

     Do i=1,natms
        atmnam(i)=sitnam(lsite(i))
        ltype(i)=typsit(lsite(i))
        chge(i)=chgsit(lsite(i))
        weight(i)=wgtsit(lsite(i))
        lfrzn(i)=frzsit(lsite(i))
        lfree(i)=fresit(lsite(i))
     End Do

     If (lbook) Then
        safe=.true. ! Initialise safety flag

! Change nlast and refresh record global atom indices for local sorting
! since particles may have been relocated across domains as the
! the old halo is invalid and a new one is not set yet.  This allows for
! local_index search over natms in pmf_units_set and compress_book_intra.
! Otherwise, everywhere else in the code, the search is over nlast as
! domain only indices are caught by the condition (1 >= index <= natms)!!!

        nlast=natms
        Do i=1,nlast
           lsi(i)=i
           lsa(i)=ltg(i)
        End Do
        Call shellsort2(nlast,lsi,lsa)

! Check safety of working arrays for all active bookkeeping arrays

        If (megshl > 0) safe(1)=(ntshl  <= mxshl )
        If (m_con  > 0) safe(2)=(ntcons <= mxcons)
        If (megpmf > 0) safe(3)=(ntpmf  <= mxpmf )
        If (m_rgd  > 0) safe(4)=(ntrgd  <= mxrgd )
        If (megtet > 0) safe(5)=(ntteth <= mxteth)
        If (megbnd > 0) safe(6)=(ntbond <= mxbond)
        If (megang > 0) safe(7)=(ntangl <= mxangl)
        If (megdih > 0) safe(8)=(ntdihd <= mxdihd)
        If (meginv > 0) safe(9)=(ntinv  <= mxinv )

        Call gcheck(comm,safe)

        If (Any(.not.safe)) Then
           itmp(1)=ntshl  ; jtmp(1)=mxshl
           itmp(2)=ntcons ; jtmp(2)=mxcons
           itmp(3)=ntpmf  ; jtmp(3)=mxpmf
           itmp(4)=ntrgd  ; jtmp(4)=mxrgd
           itmp(5)=ntteth ; jtmp(5)=mxteth
           itmp(6)=ntbond ; jtmp(6)=mxbond
           itmp(7)=ntangl ; jtmp(7)=mxangl
           itmp(8)=ntdihd ; jtmp(8)=mxdihd
           itmp(9)=ntinv  ; jtmp(9)=mxinv

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

        If (megshl > 0) Call pass_shared_units &
     (mxshl, Lbound(listshl,Dim=1),Ubound(listshl,Dim=1),ntshl, listshl,mxfshl,legshl,lshmv_shl,lishp_shl,lashp_shl,comm,&
     q0,q1,q2,q3,rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz)

        If (m_con  > 0) Call pass_shared_units &
     (mxcons,Lbound(listcon,Dim=1),Ubound(listcon,Dim=1),ntcons,listcon,mxfcon,legcon,lshmv_con,lishp_con,lashp_con,comm,&
     q0,q1,q2,q3,rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz)

        If (megpmf > 0) Call pmf_units_set(comm)

        If (m_rgd  > 0) Call pass_shared_units &
     (mxrgd, Lbound(listrgd,Dim=1),Ubound(listrgd,Dim=1),ntrgd, listrgd,mxfrgd,legrgd,lshmv_rgd,lishp_rgd,lashp_rgd,comm,&
   q0,q1,q2,q3,rgdvxx,rgdvyy,rgdvzz,rgdoxx,rgdoyy,rgdozz   )

! Compress the rest of the bookkeeping arrays if needed

        If (megtet > 0) Call compress_book_intra &
           (mxteth,ntteth,Ubound(listtet,Dim=1),listtet,mxftet,legtet,comm)

        If (megbnd > 0) Call compress_book_intra &
           (mxbond,ntbond,Ubound(listbnd,Dim=1),listbnd,mxfbnd,legbnd,comm)
        If (megang > 0) Call compress_book_intra &
           (mxangl,ntangl,Ubound(listang,Dim=1),listang,mxfang,legang,comm)
        If (megdih > 0) Call compress_book_intra &
           (mxdihd,ntdihd,Ubound(listdih,Dim=1),listdih,mxfdih,legdih,comm)
        If (meginv > 0) Call compress_book_intra &
           (mxinv,ntinv,  Ubound(listinv,Dim=1),listinv,mxfinv,leginv,comm)

     End If

  Else

! Restore periodic boundaries (re-bound > re-wrap)

     Call pbcshift(imcon,cell,natms,xxx,yyy,zzz)

  End If

! Halt program if potential cutoff exceeds cell width

  If (rlnk > Min(celprp(7),celprp(8),celprp(9))/2.0_wp) Call error(95)

End Subroutine relocate_particles


End Module deport_data
