Subroutine pmf_tags(lstitr,indpmf,pxx,pyy,pzz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! dl_poly_4 subroutine for identifying, constructing and indexing PMF
! constraints' vectors for iterative (constraints) algorithms
!
! Note: must be used in conjunction with integration algorithms
!
! copyright - daresbury laboratory
! author    - i.t.todorov february 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Use kinds, only : wp
  Use setup_module

  Use configuration, Only : natms,nlast,lsi,lsa,lfrzn
  Use pmf_module,    Only : ntpmf,listpmf

  Implicit None

  Logical,           Intent( InOut ) :: lstitr(1:mxatms)
  Integer,           Intent(   Out ) :: indpmf(1:Max(mxtpmf(1),mxtpmf(2)),1:2,1:mxpmf)
  Real( Kind = wp ), Intent(   Out ) :: pxx(1:mxpmf),pyy(1:mxpmf),pzz(1:mxpmf)

  Integer :: ipmf,jpmf,j,k,local_index

! Loop over all local to this node PMF constraints, their units and
! save the indices of the members that are present on my domain (no halo)
! update lstitr

  Do ipmf=1,ntpmf

! Initialise indices

     indpmf(:,:,ipmf) = 0

     Do jpmf=1,2
        If (listpmf(0,2,ipmf) == jpmf .or. listpmf(0,2,ipmf) == 3) Then
           Do k=1,mxtpmf(jpmf)
              j=local_index(listpmf(k,jpmf,ipmf),nlast,lsi,lsa)
              indpmf(k,jpmf,ipmf)=j
              If (j > 0 .and. j <= natms) lstitr(j)=(lfrzn(j) == 0)
           End Do
        End If
     End Do

  End Do

! Get PMF units' COM vectors

  Call pmf_coms(indpmf,pxx,pyy,pzz)

End Subroutine pmf_tags
