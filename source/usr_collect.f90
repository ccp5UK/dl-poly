Subroutine usr_collect(rrt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for accumulating statistic for USR RDFs
!
! Note: to be used in external_field_apply
!
! copyright - daresbury laboratory
! author    - i.t.todorov & a.brukhno november 2016
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds_f90
  Use setup_module,  Only : mxgusr
  Use rdf_module,    Only : ncfusr,rusr,usr

  Implicit None

  Real( Kind = wp ), Intent( In    ) :: rrt

  Integer           :: ll
  Real( Kind = wp ) :: rdelr

! set cutoff condition for pair forces and grid interval for rdf tables

  rdelr= Real(mxgusr,wp)/rusr

  If (rrt < rusr) Then ! apply truncation of potential
     ll=Min(1+Int(rrt*rdelr),mxgusr)
     usr(ll) = usr(ll) + 1.0_wp ! accumulate correlation
     ncfusr = ncfusr + 1        ! Increment sample
  End If

End Subroutine usr_collect
