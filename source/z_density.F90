Module z_density

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_4 module declaring global z-density variables and arrays
  !
  ! copyright - daresbury laboratory
  ! author    - i.t.todorov march 2014
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, Only : wp,wi
  Use comms,  Only : comms_type,gsum
  Use constants,  Only : nzdndt
  Use site, Only : site_type
  Use configuration, Only : configuration_type
  Use particle,        Only : corePart
  Use errors_warnings, Only : error,info
  Implicit None

  Private

  !> Type containing z density data
  Type, Public :: z_density_type
    Private

    !> Collection switch
    Logical, Public :: l_collect
    !> Printing switch
    Logical, Public :: l_print
    !> Number of configurations sampled
    Integer( Kind = wi ), Public :: n_samples = 0
    !> Collection frequency in steps
    Integer( Kind = wi ), Public :: frequency
    !> Bin width
    Real( Kind = wp ), Public :: bin_width
    !> z density
    Real( Kind = wp ), Allocatable, Public :: density(:,:)
  Contains
    Private

    Final :: cleanup
  End Type z_density_type

  Public :: allocate_z_density_arrays, z_density_collect, z_density_compute

Contains

  Subroutine allocate_z_density_arrays(zdensity,max_grid_rdf,mxatyp)
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Integer( Kind = wi ), Intent( In    ) :: max_grid_rdf
    Integer( Kind = wi ), Intent( In    ) :: mxatyp

    Integer :: fail

    fail = 0

    Allocate (zdensity%density(1:max_grid_rdf,1:mxatyp), Stat = fail)

    If (fail > 0) Call error(1016)

    zdensity%density = 0.0_wp

  End Subroutine allocate_z_density_arrays

  Subroutine z_density_collect(max_grid_rdf,zdensity,config)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for accumulating statistic for z-density profile
    !
    ! copyright - daresbury laboratory
    ! author    - i.t.todorov march 2014
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer( Kind = wi ), Intent( In    ) :: max_grid_rdf
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Type( configuration_type ), Intent( InOut ) :: config 

    Integer           :: i,k,l
    Real( Kind = wp ) :: zlen,zleno2,rdelr

    ! accumulator

    zdensity%n_samples=zdensity%n_samples+1

    ! length of cell in z direction

    zlen=Abs(config%cell(3))+Abs(config%cell(6))+Abs(config%cell(9))

    ! half of z length

    zleno2=zlen*0.5_wp

    ! grid interval for density profiles

    rdelr=Real(max_grid_rdf,wp)/zlen

    ! set up atom iatm type and accumulate statistic

    Do i=1,config%natms
      k=config%ltype(i)

      l=Min(1+Int((config%parts(i)%zzz+zleno2)*rdelr),max_grid_rdf)

      zdensity%density(l,k)=zdensity%density(l,k) + 1.0_wp
    End Do

  End Subroutine z_density_collect


  Subroutine z_density_compute(config,max_grid_rdf,zdensity,sites,comm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! dl_poly_4 subroutine for calculating z-density profile from
    ! accumulated data
    !
    ! copyright - daresbury laboratory
    ! author    - t.forester march 1994
    ! amended   - i.t.todorov march 2014
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Integer( Kind = wi ), Intent( In    ) :: max_grid_rdf
    Type( configuration_type ), Intent( InOut ) :: config
    Type( z_density_type ), Intent( InOut ) :: zdensity
    Type( site_type ), Intent( In    ) :: sites
    Type( comms_type ), Intent( InOut ) :: comm

    Integer           :: j,k
    Real( Kind = wp ) :: zlen,delr,dvolz,factor,rho,rho1,rrr,sum,sum1
    Character( Len = 256 ) :: messages(2)

    Write(messages(1),'(a)') 'z density profiles:'
    Write(messages(2),'(2x,a,i8,a)') 'calculated using ',zdensity%n_samples,' configurations'
    Call info(messages,2,.true.)

    ! open Z density file and write headers

    If (comm%idnode == 0) Then
      Open(Unit=nzdndt, File='ZDNDAT', Status='replace')
      Write(nzdndt,'(a)') config%cfgname
      Write(nzdndt,'(2i10)') sites%ntype_atom,max_grid_rdf
    End If

    ! length of cell in z direction

    zlen=Abs(config%cell(3))+Abs(config%cell(6))+Abs(config%cell(9))

    ! grid interval for density profiles

    delr=zlen/Real(max_grid_rdf,wp)

    ! volume of z-strip

    dvolz=(config%volm/zlen)*delr

    ! normalisation factor

    zdensity%n_samples=Max(zdensity%n_samples,1)
    factor=1.0_wp/(Real(zdensity%n_samples,wp)*dvolz)

    ! for every species

    Do k=1,sites%ntype_atom
      Write(messages(1),'(2x,a,a8)') 'rho(r): ',sites%unique_atom(k)
      Write(messages(2),'(9x,a1,6x,a3,9x,a4)') 'r','rho','n(r)'
      Call info(messages,2,.true.)
      If (comm%idnode == 0) Then
        Write(nzdndt,'(a8)') sites%unique_atom(k)
      End If

      ! global sum of data on all nodes

      Call gsum(comm,zdensity%density(1:max_grid_rdf,k))

      ! running integration of z-density

      sum=0.0_wp

      ! loop over distances

      Do j=1,max_grid_rdf
        rrr=(Real(j,wp)-0.5_wp)*delr - zlen*0.5_wp
        rho=zdensity%density(j,k)*factor
        sum=sum + rho*dvolz

        ! null it if < 1.0e-6_wp

        If (rho < 1.0e-6_wp) Then
          rho1 = 0.0_wp
        Else
          rho1 = rho
        End If

        If (sum < 1.0e-6_wp) Then
          sum1 = 0.0_wp
        Else
          sum1 = sum
        End If

        ! print out information

        Write(messages(1),'(2x,f10.4,1p,2e14.6)') rrr,rho1,sum1
        Call info(messages(1),.true.)
        If (comm%idnode == 0) Then
          Write(nzdndt,"(1p,2e14.6)") rrr,rho
        End If
      End Do
    End Do

    If (comm%idnode == 0) Close(Unit=nzdndt)

  End Subroutine z_density_compute

  Subroutine cleanup(zdensity)
    Type( z_density_type ) :: zdensity

    If (Allocated(zdensity%density)) Then
      Deallocate(zdensity%density)
    End If
  End Subroutine cleanup
End Module z_density
