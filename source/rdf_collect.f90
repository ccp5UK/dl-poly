Subroutine rdf_collect(iatm,rcut,rsqdf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for radial
! distribution functions
!
! Note: to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - t.forester march 1994
! amended   - i.t.todorov february 2012
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,      Only : mxlist,mxatms,mxgrdf,zero_plus
  Use site_module,       Only : numtyp
  Use config_module,     Only : natms,ltg,ltype,list
  Use statistics_module, Only : ntprdf,lstrdf,rdf

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: rsqdf

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: rcsq,rdelr

  Integer                 :: idi,jatm,ai,aj,keyrdf,kk,ll,m
  Real( Kind = wp )       :: rrr,rsq


  If (newjob) Then
     newjob = .false.

! set cutoff condition for pair forces and grid interval for rdf tables

     rcsq = rcut*rcut
     rdelr= Real(mxgrdf,wp)/rcut
  End If

! global identity of iatm

  idi=ltg(iatm)

! set up atom iatm type and exclude it if absent crystallographically

  ai=ltype(iatm)

  If (numtyp(ai) > zero_plus) Then

! start of primary loop for rdf accumulation

     Do m=1,list(0,iatm)

! atomic and type indices

        jatm=list(m,iatm)
        aj=ltype(jatm)

        If (numtyp(aj) > zero_plus .and. jatm <= natms .or. idi < ltg(jatm)) Then

! rdf function indices

           keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
           kk=lstrdf(keyrdf)

! only for valid interactions specified for a look up

           If (kk > 0 .and. kk <= ntprdf) Then

! apply truncation of potential

              rsq=rsqdf(m)

              If (rsq < rcsq) Then
                 rrr=Sqrt(rsq)
                 ll=Int(rrr*rdelr+0.999999_wp)

! accumulate correlation

                 rdf(ll,kk) = rdf(ll,kk) + 1.0_wp
              End If

           End If

        End If

     End Do

  End If

End Subroutine rdf_collect
