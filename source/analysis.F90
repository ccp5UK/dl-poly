!> Module for run time analysis
Module analysis
  Use kinds,         Only : wp,wi
  Use setup,         Only : zero_plus
  Use configuration, Only : volm,cell
  Use site,          Only : dens,ntpatm
  Use statistics,    Only : stats_type
  Use bonds,         Only : bonds_type,bonds_compute
  Use angles,        Only : angles_type,angles_compute
  Use dihedrals,     Only : dihedrals_type,dihedrals_compute
  Use inversions,    Only : ncfinv,mxginv1,inversions_compute
  Use greenkubo,     Only : greenkubo_type,vaf_compute
  Use rdfs,          Only : ncfrdf,ncfusr,l_errors_block,l_errors_jack,rdf_compute, &
                            usr_compute,calculate_errors,calculate_errors_jackknife
  Use z_density,     Only : ncfzdn,z_density_compute
  Use comms,         Only : comms_type
  Implicit None

  Private

  Public :: analysis_result

Contains

  !> Calculate and print final analysis
  Subroutine analysis_result(lrdf,lzdn,lpana,lprdf,lpzdn, &
                             nstep,tstep,rcut,temp,ensemble, &
                             bond,angle,dihedral,stats,green,comm)

    Logical, Intent( In    ) :: lrdf,lzdn,lpana,lprdf,lpzdn
    !> Number of simulation steps
    Integer( Kind = wi ), Intent( In    ) :: nstep
    !> Simulation tiime step (ps)
    Real( Kind = wp ), Intent( In    ) :: tstep
    !> Cut off
    Real( Kind = wp ), Intent( In    ) :: rcut
    !> Temperature
    Real( Kind = wp ), Intent( In    ) :: temp
    !> Ensemble key
    Integer( Kind = wi ), Intent( In    ) :: ensemble
    !> Bonds data
    Type( bonds_type ), Intent( InOut ) :: bond
    !> Angles data
    Type( angles_type ), Intent( InOut ) :: angle
    !> Dihedrals data
    Type( dihedrals_type ), Intent( InOut ) :: dihedral
    !> Statistics data
    Type( stats_type ), Intent( In    ) :: stats
    !> Greenkubo data
    Type( greenkubo_type ), Intent( In    ) :: green
    !> Comms
    Type( comms_type ), Intent( InOut ) :: comm

    !> Average volume
    Real( Kind = wp ) :: avvol
    Integer( Kind = wi ) :: i

    ! If still running in the pure equilibration regime - NO AVERAGES
    If (stats%numacc == 0) Then
      Return
    End If

    ! Average volume
    avvol = stats%sumval(19)
    ! Scale densities for average volume and average volume and cell
    Do i = 1,ntpatm
      dens(i) = dens(i)*(volm/avvol)
    End Do

    ! Redefine volume for analysis routines
    volm = avvol
    ! Redefine cell dimensions for analysis routines for npt/nst
    If (ensemble >= 20) Then
      Do i = 1, 9
        cell(i) = stats%sumval(36+ntpatm+i)
      End Do
    End If

    ! Calculate and print radial distribution functions
    ! If block average errors, output that, else if jackknife errors output those, else just RDF.
    If (lrdf .and. lprdf .and. ncfrdf > 0 .and. l_errors_block) Then
      Call calculate_errors(temp, rcut, nstep, comm)
    End If
    If (lrdf .and. lprdf .and. ncfrdf > 0 .and. l_errors_jack .and. .not. l_errors_block) Then
      Call calculate_errors_jackknife(temp, rcut, nstep, comm)
    End If
    If (lrdf .and. lprdf .and. ncfrdf > 0 .and. .not.(l_errors_block .or. l_errors_jack)) Then
      Call rdf_compute(lpana,rcut,temp,comm)
    End IF
    If (ncfusr > 0) Call usr_compute(comm)

    ! calculate and print z-density profile
    If (lzdn .and. lpzdn .and. ncfzdn > 0) Then
      Call z_density_compute(comm)
    End If

    ! calculate and print velocity autocorrelation function
    If (green%samp > 0 .and. green%l_print .and. green%vafcount > zero_plus) Then
      Call vaf_compute(tstep,green,comm)
    End If

    ! Calculate and print PDFs
    If (lpana) Then
      If (bond%bin_pdf > 0 .and. bond%n_frames > 0) Call bonds_compute(temp,bond,comm)
      If (angle%bin_adf > 0 .and. angle%n_frames > 0) Call angles_compute(temp,angle,comm)
      If (dihedral%bin_adf > 0 .and. dihedral%n_frames > 0) Call dihedrals_compute(temp,dihedral,comm)
      If (mxginv1 > 0 .and. ncfinv > 0) Call inversions_compute(temp,comm)
    End If
  End Subroutine analysis_result
End Module analysis
