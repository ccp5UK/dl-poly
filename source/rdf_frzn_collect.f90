Subroutine rdf_frzn_collect(iatm,rcut,rsqdf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for radial
! distribution functions of frozen pairs
!
! Note: to be used as part of two_body_forces
!
! copyright - daresbury laboratory
! author    - i.t.todorov march 2014
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,  Only : mxlist,mxgrdf
  Use config_module, Only : natms,ltg,ltype,list
  Use rdf_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatm
  Real( Kind = wp ),                        Intent( In    ) :: rcut
  Real( Kind = wp ), Dimension( 1:mxlist ), Intent( In    ) :: rsqdf

  Logical,           Save :: newjob = .true.
  Real( Kind = wp ), Save :: rcsq,rdelr

  Integer                 :: limit,idi,jatm,ai,aj,keyrdf,kk,ll,m
  Real( Kind = wp )       :: rrr,rsq


  If (newjob) Then
     newjob = .false.

! set cutoff condition for pair forces and grid interval for rdf tables

     rcsq = rcut*rcut
     rdelr= Real(mxgrdf,wp)/rcut
  End If

! global identity and type of iatm

  idi=ltg(iatm)
  ai=ltype(iatm)

! Get list limit

  limit=list(-2,iatm)-list(-1,iatm)

! start of primary loop for rdf accumulation

  Do m=1,limit

! atomic and type indices

     jatm=list(list(-1,iatm)+m,iatm)
     aj=ltype(jatm)

     If (jatm <= natms .or. idi < ltg(jatm)) Then

! rdf function indices

        keyrdf=(Max(ai,aj)*(Max(ai,aj)-1))/2+Min(ai,aj)
        kk=lstrdf(keyrdf)

! only for valid interactions specified for a look up

        If (kk > 0 .and. kk <= ntprdf) Then

! apply truncation of potential

           rsq=rsqdf(m)

           If (rsq < rcsq) Then
              rrr=Sqrt(rsq)
              ll=Min(1+Int(rrr*rdelr),mxgrdf)

! accumulate correlation

              rdf(ll,kk) = rdf(ll,kk) + 1.0_wp
           End If

        End If

     End If

  End Do

End Subroutine rdf_frzn_collect
