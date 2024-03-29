!> Module for run time analysis
Module analysis

  Use angles,        Only: angles_compute,&
                           angles_type
  Use bonds,         Only: bonds_compute,&
                           bonds_type
  Use comms,         Only: comms_type
  Use configuration, Only: configuration_type
  Use constants,     Only: zero_plus
  Use dihedrals,     Only: dihedrals_compute,&
                           dihedrals_type
  Use filename,      Only: file_type
  Use greenkubo,     Only: greenkubo_type,&
                           vaf_compute
  Use inversions,    Only: inversions_compute,&
                           inversions_type
  Use kinds,         Only: wi,&
                           wp
  Use rdfs,          Only: calculate_errors,&
                           calculate_errors_jackknife,&
                           rdf_compute,&
                           rdf_type,&
                           usr_compute
  Use site,          Only: site_type
  Use statistics,    Only: stats_type
  Use thermostat,    Only: thermostat_type
  Use z_density,     Only: z_density_compute,&
                           z_density_type

  Implicit None

  Private

  Public :: analysis_result

Contains

  !> Calculate and print final analysis
  Subroutine analysis_result( &
    rcut, thermo, &
    bond, angle, dihedral, inversion, stats, &
    green, zdensity, sites, rdf, config, files, comm)

    !> Cut off
    Real(Kind=wp),            Intent(In   ) :: rcut
    Type(thermostat_type),    Intent(In   ) :: thermo
    Type(bonds_type),         Intent(InOut) :: bond
    Type(angles_type),        Intent(InOut) :: angle
    Type(dihedrals_type),     Intent(InOut) :: dihedral
    Type(inversions_type),    Intent(InOut) :: inversion
    Type(stats_type),         Intent(In   ) :: stats
    Type(greenkubo_type),     Intent(In   ) :: green
    Type(z_density_type),     Intent(InOut) :: zdensity
    Type(site_type),          Intent(InOut) :: sites
    Type(rdf_type),           Intent(InOut) :: rdf
    Type(configuration_type), Intent(InOut) :: config
    Type(file_type),          Intent(InOut) :: files(:)
    Type(comms_type),         Intent(InOut) :: comm

    Integer(Kind=wi) :: i
    Real(Kind=wp)    :: avvol, temp

!> Ensemble key
!> Bonds data
!> Angles data
!> Dihedrals data
!> Inversion angles data
!> Statistics data
!> Greenkubo data
!> Z density data
!> Site data
!> RDF data
!> Config data
!> Comms
!> Average volume

    temp = stats%sumval(2)

    ! If still running in the pure equilibration regime - NO AVERAGES
    If (stats%numacc == 0) Then
      Return
    End If

    ! Average volume
    avvol = stats%sumval(19)
    ! Scale densities for average volume and average volume and config%cell
    Do i = 1, sites%ntype_atom
      sites%dens(i) = sites%dens(i) * (config%volm / avvol)
    End Do

    ! Redefine volume for analysis routines
    config%volm = avvol
    ! Redefine config%cell dimensions for analysis routines for npt/nst
    If (thermo%ensemble >= 20) Then
      Do i = 1, 9
        config%cell(i) = stats%sumval(36 + sites%ntype_atom + i)
      End Do
    End If

    ! Calculate and print radial distribution functions
    ! If block average errors, output that, else if jackknife errors output those, else just RDF.
    If (rdf%l_collect .and. rdf%l_print .and. rdf%n_configs > 0 .and. rdf%l_errors_block) Then
      Call calculate_errors(temp, rcut, sites, rdf, config, files, comm)
    End If
    If (rdf%l_collect .and. rdf%l_print .and. rdf%n_configs > 0 .and. rdf%l_errors_jack .and. .not. rdf%l_errors_block) Then
      Call calculate_errors_jackknife(temp, rcut, sites, rdf, config, files, comm)
    End If
    If (rdf%l_collect .and. rdf%l_print .and. rdf%n_configs > 0 .and. .not. (rdf%l_errors_block .or. rdf%l_errors_jack)) Then
      Call rdf_compute(stats%lpana, rcut, temp, sites, rdf, config, files, comm)
    End If
    If (rdf%n_configs_usr > 0) Call usr_compute(rdf, config, comm)

    ! calculate and print z-density profile
    If (zdensity%l_collect .and. zdensity%l_print .and. zdensity%n_samples > 0) Then
      Call z_density_compute(config, zdensity, sites, comm)
    End If

    ! calculate and print velocity autocorrelation function
    If (green%samp > 0 .and. green%l_print .and. green%vafcount > zero_plus) Then
      Call vaf_compute(thermo%tstep, sites%num_type_nf, sites%mxatyp, green)
    End If

    ! Calculate and print PDFs
    If (stats%lpana) Then
      If (bond%bin_pdf > 0 .and. bond%n_frames > 0) Then
        Call bonds_compute(temp, sites%unique_atom, bond, config, comm)
      End If
      If (angle%bin_adf > 0 .and. angle%n_frames > 0) Then
        Call angles_compute(temp, sites%unique_atom, angle, config, comm)
      End If
      If (dihedral%bin_adf > 0 .and. dihedral%n_frames > 0) Then
        Call dihedrals_compute(temp, sites%unique_atom, dihedral, config, comm)
      End If
      If (inversion%bin_adf > 0 .and. inversion%n_frames > 0) Then
        Call inversions_compute(temp, sites%unique_atom, inversion, config, comm)
      End If
    End If

  End Subroutine analysis_result
End Module analysis
