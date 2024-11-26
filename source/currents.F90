Module currents
  Use comms,           Only: comms_type,&
                             gsum
  Use configuration,   Only: configuration_type
  Use constants,       Only: czero
  Use errors_warnings, Only: error
  Use filename,        Only: file_type
  Use kinds,           Only: wp
  Use site,            Only: site_type

  Implicit None

  Type, Public :: current_type
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_5 type for calculating momentum currents, energy currents,
  ! and density from supplied KPOINTS.
  !
  ! Refer to the "Collective properties" chapter of
  ! Balucani, U. and Zoppi, M., 1995. Dynamics of the liquid state
  ! (Vol. 10). Clarendon Press.
  !
  ! author    - a.m.elena 2019
  ! re-write  - h.l.devereux 2024
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !> e^{ik\cdot r_n(t)}
    Complex(Kind=wp), Allocatable :: density_jlk(:, :, :)
    !> (v_n(t)\cdot \hat{k})\hat{k}e^{ik\cdot r_n(t)}
    Complex(Kind=wp), Allocatable :: longitudinal_jlk(:, :, :)
    !> (v_n(t)-(v_n(t)\cdot \hat{k})\hat{k})e^{ik\cdot r_n(t)}
    Complex(Kind=wp), Allocatable :: transverse_jlk(:, :, :)
    !> 0.5 * [mv_n(t)^2 + \sum_{j}U_{nj}(t)]
    Complex(Kind=wp), Allocatable :: energy_density_jlk(:, :, :)
    !> 0.5[e_n v_n - 0.5(v_n(t)^b+v_j(t)^b)(r_nj(t)^a r_nj(t)^b/|r_nj(t)|^2)|r_{nj}(t)|U'(r_{nj}(t))(1-exp(ik\cdot r_n(t))/ik\cdot r_n(t))]
    Complex(Kind=wp), Allocatable :: energy_jlk(:, :, :)
    Complex(Kind=wp), Allocatable :: stress_jlk(:, :, :)
    Integer                       :: nkpoints, lag, yaml_indent=6
    Integer                       :: file_handle = -2
    Logical                       :: on = .false., &
                                     io_yaml = .false., &
                                     k_energy_stress_current_on = .false.

  Contains
    Private
    Procedure, Public  :: init
    Procedure, Public  :: compute
    Procedure, Private :: write_yaml_current_block
    Final              :: cleanup
  End Type

Contains

  Subroutine init(T, nk, lag, fcurrent, comm, config, types, io_yaml)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_5 subroutine to setup currents and CURRENTS file.
  !
  ! author    - a.m.elena 2019
  ! re-write  - h.l.devereux 2024
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(current_type)                               :: T
    Integer,                            Intent(In   ) :: nk, lag, types
    Type(file_type),                    Intent(InOut) :: fcurrent
    Type(comms_type),                   Intent(In   ) :: comm
    Type(configuration_type),           Intent(In   ) :: config
    Logical,                  Optional, Intent(In   ) :: io_yaml


    Allocate (T%density_jlk(nk, 1, types))
    Allocate (T%longitudinal_jlk(nk, 3, types))
    Allocate (T%transverse_jlk(nk, 3, types))
    Allocate (T%energy_density_jlk(nk, 3, types))
    If (T%k_energy_stress_current_on) Then
      Allocate(T%energy_jlk(nk, 3, types))
      Allocate(T%stress_jlk(nk, 6, types))
    End If
    T%nkpoints = nk
    T%lag = lag
    If (Present(io_yaml)) Then 
      T%io_yaml = io_yaml
    Else 
      T%io_yaml = .false.
    End If
    If (comm%idnode == 0) Then
      Open (Newunit=fcurrent%unit_no, File=fcurrent%filename, Status='unknown', Action="Write")
      If (T%io_yaml) Then
        Write (fcurrent%unit_no,'(a)') "%YAML 1.2"
        Write (fcurrent%unit_no,'(a)') "---"
        Write (fcurrent%unit_no, '(a,a)') "title: ", Trim(config%cfgname)
        Write (fcurrent%unit_no, '(a)') "timesteps: "
      End If
      T%file_handle = fcurrent%unit_no
    End If
  End Subroutine init

  Subroutine write_yaml_current_block(T, sites, current, name, last)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_5 subroutine to write a current block in yaml format.
  !
  ! author  - h.l.devereux 2024
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(current_type),           Intent(InOut) :: T
    Type(site_type),               Intent(In   ) :: sites
    Complex(Kind=wp), Allocatable, Intent(In   ) :: current(:, :, :)
    Character(Len=*),              Intent(In   ) :: name
    Logical,                       Intent(In   ) :: last

    Integer          :: atype
    Character(Len=2) :: ending

    If (.not. last) Then
      ending = "},"
    Else
      ending = "} "
    End If

    Write (T%file_handle, '(a)') Repeat(" ", T%yaml_indent)//name//": {"
    Do atype = 1, sites%mxatyp
      Write (T%file_handle, '(a, *(g16.8, :, ","))', advance="no") &
        Repeat(" ", T%yaml_indent*3)//Trim(sites%unique_atom(atype))//": [", current(:, :, atype)
      If (atype < sites%mxatyp) Then
        Write (T%file_handle, '(a)') "],"
      Else
        Write (T%file_handle, '(a)') "]"
        Write (T%file_handle, '(a)') Repeat(" ", T%yaml_indent*2)//Trim(ending)
      End If
    End Do
  End Subroutine write_yaml_current_block

  Subroutine compute(T, config, time, comm, sites, pp_energy, pp_virial, pp_cur_stress)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! dl_poly_5 subroutine to calculate the currents. Note pp_energy and
  ! pp_virial are computed with forces.
  !
  ! author    - a.m.elena 2019
  ! re-write  - h.l.devereux 2024
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(current_type),      Intent(InOut)           :: T
    Type(configuration_type), Intent(In   )           :: config
    Real(Kind=wp),            Intent(In   )           :: time
    Type(comms_type),         Intent(InOut)           :: comm
    Type(site_type),          Intent(In   )           :: sites
    Real(Kind=wp),            Intent(In   )           :: pp_energy(:)
    Complex(Kind=wp),         Intent(In   ), Optional :: pp_virial(:, :, :)
    Complex(Kind=wp),         Intent(In   ), Optional :: pp_cur_stress(:, :, :)

    Integer          :: i, k, atype
    Real(Kind=wp)    :: velocity(1:3), khat(1:3), tmp(1:3), eng(1:3)
    Complex(Kind=wp) :: dens

    If (T%k_energy_stress_current_on) Then
      If ((.not. Present(pp_virial)) .or. (.not. Present(pp_cur_stress))) Then
        Call error(0, "No per-particle arrays passed to compute currents")
      End If
    End If

    T%density_jlk = czero
    T%longitudinal_jlk = czero
    T%transverse_jlk = czero
    T%energy_density_jlk = czero
    If (T%k_energy_stress_current_on) Then
      T%energy_jlk = czero
      T%stress_jlk = czero
    End If

    Do k = 1, config%k%n
      khat = config%k%r(:, k) / norm2(config%k%r(:, k))
      Do i = 1, config%natms
        atype = config%ltype(i)
        dens = Exp(Cmplx(0.0_wp, &
          Dot_product(config%k%r(:, k), [config%parts(i)%xxx, config%parts(i)%yyy, config%parts(i)%zzz]), wp))
        velocity = [config%vxx(i), config%vyy(i), config%vzz(i)]
        tmp = Dot_product(velocity, khat)*khat

        T%density_jlk(k, 1, atype) = T%density_jlk(k, 1, atype) + dens
        T%longitudinal_jlk(k, :, atype) = T%longitudinal_jlk(k, :, atype) + tmp*dens
        T%transverse_jlk(k, :, atype) = T%transverse_jlk(k, :, atype) + (velocity-tmp)*dens
        
        eng = 0.5_wp*Dot_product(velocity, velocity)+pp_energy(i)
        T%energy_density_jlk(k, :, atype) = T%energy_density_jlk(k, :, atype) + eng*dens
        If (T%k_energy_stress_current_on) Then
          T%energy_jlk(k, :, atype) = T%energy_jlk(k, :, atype) + &
            (eng*velocity+0.5_wp*pp_virial(i, k, :))*dens
          T%stress_jlk(k, :, atype) = T%stress_jlk(k, :, atype) + pp_cur_stress(i, k, :)*dens
        End If
      End Do
    End Do
    Call gsum(comm, T%density_jlk)
    Call gsum(comm, T%longitudinal_jlk)
    Call gsum(comm, T%transverse_jlk)
    Call gsum(comm, T%energy_density_jlk)
    If (T%k_energy_stress_current_on) Then
      Call gsum(comm, T%energy_jlk)
      Call gsum(comm, T%stress_jlk)
    End If

    If (comm%idnode == 0) Then
      If (T%io_yaml) Then 
        Write (T%file_handle, '(a, g16.8, a)') "  - { time: ", time, ","
        Call T%write_yaml_current_block(sites, T%density_jlk, "density", .false.)
        Call T%write_yaml_current_block(sites, T%longitudinal_jlk, "longitudinal", .false.)
        Call T%write_yaml_current_block(sites, T%transverse_jlk, "transverse", .false.)
        If (T%k_energy_stress_current_on) Then
          Call T%write_yaml_current_block(sites, T%energy_density_jlk, "energy_density", .false.)
          Call T%write_yaml_current_block(sites, T%stress_jlk, "stress", .false.)
          Call T%write_yaml_current_block(sites, T%energy_jlk, "energy", .true.)
        Else
          Call T%write_yaml_current_block(sites, T%energy_density_jlk, "energy_density", .true.)
        End If
        Write(T%file_handle, '(a)') "    }"
      Else
        Do atype = 1, sites%mxatyp
          Write (T%file_handle, '(g16.8,3a,*(g16.8, :, ","))') time, ", ", &
            Trim(sites%unique_atom(atype)), ", ", T%density_jlk(:, 1, atype)
          Write (T%file_handle, '(g16.8,3a,*(g16.8, :, ","))') time, ", ", &
            Trim(sites%unique_atom(atype)), ", ", T%longitudinal_jlk(:, :, atype)
          Write (T%file_handle, '(g16.8,3a,*(g16.8, :, ","))') time, ", ", &
            Trim(sites%unique_atom(atype)), ", ", T%transverse_jlk(:, :, atype)
          Write (T%file_handle, '(g16.8,3a,*(g16.8, :, ","))') time, ", ", &
            Trim(sites%unique_atom(atype)), ", ", T%energy_density_jlk(:, :, atype)
          If (T%k_energy_stress_current_on) Then
            Write (T%file_handle, '(g16.8,3a,*(g16.8, :, ","))') time, ", ", &
              Trim(sites%unique_atom(atype)), ", ", T%stress_jlk(:, :, atype)
            Write (T%file_handle, '(g16.8,3a,*(g16.8, :, ","))') time, ", ", &
              Trim(sites%unique_atom(atype)), ", ", T%energy_jlk(:, :, atype)
          End If
        End Do
      End If
    End If

  End Subroutine compute

  Subroutine cleanup(T)
    Type(current_type) :: T

    If (Allocated(T%density_jlk)) Deallocate (T%density_jlk)
    If (Allocated(T%longitudinal_jlk)) Deallocate (T%longitudinal_jlk)
    If (Allocated(T%transverse_jlk)) Deallocate (T%transverse_jlk)
    If (T%k_energy_stress_current_on) Then
        If (Allocated(T%energy_jlk)) Deallocate (T%energy_jlk)
        If (Allocated(T%stress_jlk)) Deallocate (T%stress_jlk)
    End If
  End Subroutine cleanup

End Module currents

