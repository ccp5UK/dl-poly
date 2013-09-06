! DL_POLY_4 NVIDIA GPU & OpenMP Port
! Irish Centre for High-End Computing (ICHEC)
! http://www.ichec.ie
!
! Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
! collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
! W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
!
! Distributed under the same license that the original, unmodified,
! DL_POLY_4 is. You should have received these sources from the
! STFC Daresbury Laboratory.

! dl_poly_cuda_wrappers.f90: A set of wrapper functions
! to allow original FORTRAN90 sources to be portably called from C
! without modifying the original .f90 sources

Subroutine wrapper_f_invert(a,b,d) bind(c)
  Use kinds_f90

  Implicit None

  Real( Kind = wp ), Dimension( 1:9 ), Intent( In    ) :: a
  Real( Kind = wp ), Dimension( 1:9 ), Intent(   Out ) :: b
  Real( Kind = wp ),                   Intent(   Out ) :: d

  Call invert(a,b,d)
End Subroutine wrapper_f_invert

Subroutine wrapper_f_error(kode) bind(c)
  Implicit None

  Integer, Intent( In    ) :: kode

  Call error(kode)
End Subroutine wrapper_f_error

Subroutine wrapper_f_erfcgen(rcut,alpha,mxgrid,erc,fer) bind(c)
  Use kinds_f90

  Integer,                                  Intent( In    ) :: mxgrid
  Real( Kind = wp ),                        Intent( In    ) :: rcut,alpha
  Real( Kind = wp ), Dimension( 1:mxgrid ), Intent(   Out ) :: erc,fer

  Call erfcgen(rcut,alpha,mxgrid,erc,fer)
End Subroutine wrapper_f_erfcgen

Subroutine wrapper_f_two_body_forces_cuda_helper(&
     iatmb,iters,keyfce,rho,safe,&
     imcon,rcut,rvdw,alpha,epsq,stress,&
     engmet,virmet,engvdw,virvdw,engcpe_rl,vircpe_rl) bind(c)
  Use kinds_f90
  Use setup_module,      Only : mxatms

  Implicit None

  Integer,           Intent(In   ) :: iatmb,iters,imcon,keyfce

  Real( Kind = wp ), Intent(In   ) :: rcut,rvdw,alpha,epsq
  Real( Kind = wp),  Dimension(1:mxatms), Intent(In) :: rho
  Real( Kind = wp ), Intent(InOut) :: engmet,virmet,engvdw,virvdw,engcpe_rl,vircpe_rl
  Real( Kind = wp ), Dimension( 1:9 ),&
                     Intent(InOut) :: stress

  Logical,           Intent(InOut) :: safe

  Call two_body_forces_cuda_helper(iatmb,iters,keyfce,rho,safe,                   &
                                   imcon,rcut,rvdw,alpha,epsq,stress,             &
                                   engmet,virmet,engvdw,virvdw,engcpe_rl,vircpe_rl)
End Subroutine wrapper_f_two_body_forces_cuda_helper

Subroutine wrapper_f_link_cell_pairs_helper                   &
           (ibegin,iend,nix,niy,niz,                          &
           nlx,nly, nlx0e, nly0e, nlz0e, nlx1s, nly1s, nlz1s, &
           lct_start,at_list,ncells,nlp,nlp3,nsbcll,megfrz,rcsq) bind(c)
  Use kinds_f90
  Use setup_module

  Implicit None

  Integer, Dimension(1:mxatms), Intent(In) :: at_list
  Integer,                      Intent(In) :: ncells,nlp,nsbcll,nlp3,megfrz, &
                                              nlx,nly, nlx0e, nly0e, nlz0e,  &
                                              nlx1s, nly1s, nlz1s,           &
                                              ibegin, iend

  Integer, Dimension(1:ncells+1), Intent(In) :: lct_start
  Integer, Dimension(1:nlp3),     Intent(In) :: nix,niy,niz;

  Real (Kind=wp),                 Intent(In) :: rcsq

  Call link_cell_pairs_helper(ibegin,iend,nix,niy,niz,      &
                              nlx,nly, nlx0e, nly0e, nlz0e, &
                              nlx1s, nly1s, nlz1s,          &
                              lct_start,at_list,            &
                              ncells,nlp,nlp3,nsbcll,megfrz,rcsq)
End Subroutine wrapper_f_link_cell_pairs_helper

Subroutine wrapper_f_link_cell_pairs_remove_exclusions_helper(ibegin,iend) bind(c)
  Implicit None
  Integer, Intent(In) :: ibegin,iend

  Call link_cell_pairs_remove_exclusions_helper(ibegin,iend)
End Subroutine wrapper_f_link_cell_pairs_remove_exclusions_helper

Subroutine wrapper_f_spme_forces_helper                                 &
           (ibegin,iiters, rcell,tmp0,fff, ixdb,iydb,izdb, ixx,iyy,izz, &
           bspx,bspy,bspz, bsdx,bsdy,bsdz, qqc_domain, ixt,iyt,izt) bind(c)
  Use kinds_f90
  Use setup_module

  Implicit None

  Integer,                        Intent( In    ) :: ixt,iyt,izt, ixdb,iydb,izdb, &
                                                     ibegin,iiters
  Integer, Dimension( 1:mxatms ), Intent( In    ) :: ixx,iyy,izz

  Real( Kind = wp ),                               Intent( In    ) :: tmp0
  Real( Kind = wp ), Dimension( 1:9 ),             Intent( In    ) :: rcell
  Real( Kind = wp ), Dimension( 0:3 ),             Intent( InOut ) :: fff
  Real( Kind = wp ), Dimension( 1:mxspl,1:mxatms), Intent( In    ) :: bsdx,bsdy,bsdz, &
                                                                      bspx,bspy,bspz
  Real( Kind = wp ), Dimension( ixdb:ixt, iydb:iyt, izdb:izt ), Intent( In ) :: qqc_domain

  Call spme_forces_helper                                           &
       (ibegin,iiters, rcell,tmp0,fff, ixdb,iydb,izdb, ixx,iyy,izz, &
       bspx,bspy,bspz, bsdx,bsdy,bsdz, qqc_domain, ixt,iyt,izt)
End Subroutine wrapper_f_spme_forces_helper

Subroutine wrapper_f_ewald_spme_forces_ccarray_helper &
           (ibegin,niters,                            &
           block_x,block_y,block_z,ixb,iyb,izb,ixt,   &
           iyt,izt,ixx,iyy,izz,it,                    &
           bspx,bspy,bspz,qqc_local) bind(c)
  Use kinds_f90
  Use setup_module
  Use config_module,  Only : natms,chge

  Implicit None

  Integer, Intent(In) :: ibegin,niters,block_x,block_y,block_z,ixb,iyb,izb,ixt,iyt,izt
  Integer,              Dimension(1:mxatms),         Intent(In)    :: ixx,iyy,izz,it

  Real( Kind = wp ),    Dimension(1:mxspl,1:mxatms), Intent(In)    :: bspx,bspy,bspz
  Real( Kind = wp ),    Dimension(1:block_x,1:block_y,1:block_z), &
                                                     Intent(InOut) :: qqc_local

  Call ewald_spme_forces_ccarray_helper         &
       (ibegin,niters,                          &
       block_x,block_y,block_z,ixb,iyb,izb,ixt, &
       iyt,izt,ixx,iyy,izz,it,                  &
       bspx,bspy,bspz,qqc_local)

End Subroutine wrapper_f_ewald_spme_forces_ccarray_helper


Subroutine wrapper_f_ewald_spme_forces_ccarray_final_reduction &
           (block_x,block_y,block_z,qqcd,qqcs) bind(c)
  Use kinds_f90
  Implicit None

  Integer, Intent(In)    :: block_x,block_y,block_z
  Real( Kind = wp ),    Dimension(1:block_x,1:block_y,1:block_z), Intent(InOut) :: qqcd
  Real( Kind = wp ),    Dimension(1:block_x,1:block_y,1:block_z), Intent(In   ) :: qqcs

  Call ewald_spme_forces_ccarray_final_reduction(block_x,block_y,block_z,qqcd,qqcs)
End Subroutine wrapper_f_ewald_spme_forces_ccarray_final_reduction

Function wrapper_f_local_index(global_index,search_limit,rank,list) bind(c)
  Implicit None

  Integer,                   Intent( In    ) :: global_index,search_limit
  Integer, Dimension( 1:* ), Intent( In    ) :: rank,list

  Integer :: wrapper_f_local_index, local_index

  wrapper_f_local_index = local_index(global_index,search_limit,rank,list)

  Return
End Function wrapper_f_local_index


Subroutine wrapper_f_bspgen_helper &
           (iatmb,iatme,nospl,xxx,yyy,zzz,bspx,bspy,bspz,bsdx,bsdy,bsdz) bind(c)
  Use kinds_f90
  Use setup_module

  Implicit None

  Integer,                                  Intent( In    ) :: iatmb,iatme,nospl
  Real( Kind = wp ), Dimension( 1:mxatms ), Intent( In    ) :: xxx,yyy,zzz

  Real( Kind = wp ), Dimension( 1:mxspl , 1:mxatms ), Intent(   Out ) :: &
                                             bsdx,bsdy,bsdz,bspx,bspy,bspz

  Call bspgen_helper &
       (iatmb,iatme,nospl,xxx,yyy,zzz,bspx,bspy,bspz,bsdx,bsdy,bsdz)
End Subroutine wrapper_f_bspgen_helper
