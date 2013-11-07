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

Module dl_poly_cuda_module
Use kinds_f90

! Define interfaces for CUDA-port-related C functions

Interface

   ! metal_ld_compute
   Subroutine metal_ld_compute_cuda_initialise                  &
              (IsListOnline, mxatms,                            &
              natms, mxgrid, ntpmet, mxmet, mxatdm, mxlist,     &
              xxx, yyy, zzz, list, ltype, ltpmet, lstmet, vmet, &
              dmet, cell, rho) bind(c)
     Use kinds_f90

     Integer                                               :: IsListOnline, mxatms,  &
                                                              natms, mxgrid, ntpmet, &
                                                              mxmet, mxatdm, mxlist

     Integer,            Dimension( -2:mxlist,1:mxatdm )    :: list
     Integer,            Dimension( 1:mxatms )             :: ltype
     Integer,            Dimension( 1:mxmet )              :: ltpmet, lstmet

     Real ( Kind = wp ), Dimension( 1:9 )                  :: cell
     Real ( Kind = wp ), Dimension( 1:mxatms )             :: xxx, yyy, zzz, rho
     Real ( Kind = wp ), Dimension( 1:mxgrid,1:mxmet,1:2 ) :: vmet, dmet
   End Subroutine metal_ld_compute_cuda_initialise

   Subroutine metal_ld_compute_cuda_invoke() bind(c)
   End Subroutine metal_ld_compute_cuda_invoke

   Subroutine metal_ld_compute_cuda_finalise() bind(c)
   End Subroutine metal_ld_compute_cuda_finalise

   Subroutine dl_poly_cuda_allocate_pinned_default &
              (n, unitSize, cpqqq_local) bind(c)
     Use iso_c_binding
     Integer     :: n, unitSize
     Type(c_ptr) :: cpqqq_local
   End Subroutine dl_poly_cuda_allocate_pinned_default

   Subroutine dl_poly_cuda_deallocate_pinned_default &
              (cpqqq_local) bind(c)
     Use iso_c_binding
     Type(c_ptr) :: cpqqq_local
   End Subroutine dl_poly_cuda_deallocate_pinned_default


   ! spme_container
   Subroutine spme_container_cuda_bspgen_initialise &
              (unroll,natms,mxspl,nospl,xxx,yyy,zzz) bind(c)
     Use kinds_f90
     Integer                                 :: unroll, mxspl, natms, nospl
     Real( Kind = wp ), Dimension( 1:natms ) :: xxx, yyy, zzz
   End Subroutine

   Subroutine spme_container_cuda_bspgen_invoke(bspx,bspy,bspz,bsdx,bsdy,bsdz) bind(c)
     Use kinds_f90
     Use setup_module, Only : mxspl, mxatms
     Real( Kind = wp ), Dimension( 1:mxspl,1:mxatms) :: bspx, bspy, bspz, &
                                                        bsdx, bsdy, bsdz
   End Subroutine Spme_container_cuda_bspgen_invoke

   Subroutine spme_container_cuda_bspgen_finalise() bind(c)
   End Subroutine spme_container_cuda_bspgen_finalise

   Subroutine spme_container_cuda_bspgen_set_leave_bspdxyz_data_online &
              (truth) bind(c)
     Logical :: truth
   End Subroutine spme_container_cuda_bspgen_set_leave_bspdxyz_data_online


   Subroutine spme_container_cuda_bspgen_set_is_under_ewald_spme_forces &
        (truth, natms) bind(c)
     Logical :: truth
     Integer :: natms
   End Subroutine spme_container_cuda_bspgen_set_is_under_ewald_spme_forces


   Subroutine spme_container_cuda_bspgen_kill_bspdxyz_data() bind(c)
   End Subroutine spme_container_cuda_bspgen_kill_bspdxyz_data


   Subroutine ewald_spme_forces_cuda_ccarray_initialise                &
              (nlast,mxspl,mxatms,ixb,iyb,izb,ixt,iyt,izt,ixx,iyy,izz, &
              bspx,bspy,bspz,bsdx,bsdy,bsdz,it,chge,                   &
              block_x,block_y,block_z,qqc_local) bind(c)
     Use kinds_f90
     Integer :: nlast, mxspl, mxatms,     &
                ixb,iyb,izb, ixt,iyt,izt, &
                block_x, block_y, block_z
     Integer,           Dimension( 1:mxatms )                      :: ixx, iyy, izz, it
     Real( Kind = wp ), Dimension( 1:mxatms )                      :: chge
     Real( Kind = wp ), Dimension( 1:mxspl,1:mxatms )              :: bspx, bspy, bspz, &
                                                                      bsdx, bsdy, bsdz
     Real( Kind = wp ), Dimension( 1:block_x,1:block_y,1:block_z ) :: qqc_local
   End Subroutine ewald_spme_forces_cuda_ccarray_initialise

   Subroutine ewald_spme_forces_cuda_ccarray_invoke () bind(c)
   End Subroutine ewald_spme_forces_cuda_ccarray_invoke

   Subroutine ewald_spme_forces_cuda_ccarray_finalise () bind(c)
   End Subroutine ewald_spme_forces_cuda_ccarray_finalise


   ! spme_forces
   Subroutine spme_forces_cuda_initialise                 &
              (natms, mxspl, mxatms, kmaxa, kmaxb, kmaxc, &
              l_cp, tmp, rcell, qqc_domain, chge,         &
              bsdx, bsdy, bsdz, bspx, bspy, bspz,         &
              ixx, iyy, izz, ixdb, ixt, iydb, iyt, izdb, izt) bind(c)
     Use kinds_f90
     Logical :: l_cp
     Integer :: natms, mxspl, mxatms, kmaxa, kmaxb, kmaxc
     Integer, Dimension( 1:mxatms ) :: ixx, iyy, izz

     Real( Kind = wp ) :: tmp
     Real( Kind = wp ), Dimension( 1:9 )                        :: rcell
     Real( Kind = wp ), Dimension( 1:mxatms )                   :: chge
     Real( Kind = wp ), Dimension( 1:mxspl,1:mxatms )           :: bspx, bspy, bspz, &
                                                                   bsdx, bsdy, bsdz
     Real( Kind = wp ), Dimension( ixdb:ixt,iydb:iyt,izdb:izt ) :: qqc_domain
   End Subroutine spme_forces_cuda_initialise

   Subroutine spme_forces_cuda_invoke(fff, fxx, fyy, fzz, fcx, fcy, fcz) bind(c)
     Use kinds_f90
     Use setup_module, Only : mxatms
     Real( Kind = wp ), Dimension( 1:mxatms ) :: fxx, fyy, fzz, &
                                                 fcx, fcy, fcz
     Real( Kind = wp ), Dimension( 0:3 )      :: fff
   End Subroutine spme_forces_cuda_invoke

   Subroutine spme_forces_cuda_finalise() bind(c)
   End Subroutine spme_forces_cuda_finalise


   ! ewald_spme_forces cccharge
   Subroutine ewald_spme_forces_cuda_cccharge_initialise              &
              (block_x, block_y, block_z, kmaxa, kmaxb, kmaxc,        &
              qqq_local, bscx, bscy, bscz, index_x, index_y, index_z, &
              rcell, rcpct2, ralph, strs) bind(c)
     Use kinds_f90

     Integer :: block_x, block_y, block_z, kmaxa, kmaxb, kmaxc
     Integer,              Dimension(1:block_x) :: index_x
     Integer,              Dimension(1:block_y) :: index_y
     Integer,              Dimension(1:block_z) :: index_z

     Real   ( Kind = wp )                                             :: rcpct2, ralph
     Real   ( Kind = wp ), Dimension( 1:9 )                           :: rcell, strs
     Complex( Kind = wp ), Dimension( 1:kmaxa)                        :: bscx
     Complex( Kind = wp ), Dimension( 1:kmaxb)                        :: bscy
     Complex( Kind = wp ), Dimension( 1:kmaxc)                        :: bscz
     Complex( Kind = wp ), Dimension( 1:block_x,1:block_y,1:block_z ) :: qqq_local
   End Subroutine ewald_spme_forces_cuda_cccharge_initialise

   Subroutine ewald_spme_forces_cuda_cccharge_invoke() bind(c)
   End Subroutine ewald_spme_forces_cuda_cccharge_invoke

   Subroutine ewald_spme_forces_cuda_cccharge_finalise() bind(c)
   End Subroutine ewald_spme_forces_cuda_cccharge_finalise


   ! two_body_forces
   Subroutine two_body_forces_cuda_initialise                          &
              (cuda_ilo,natms,mxlist,mxatms,mxatdm,mxmet,              &
              ntpmet,ntpvdw,keyfce,imcon,                              &
              cell,xxx,yyy,zzz,list,ltype,ltpmet,lstmet,vmet,dmet,rho, &
              mxgrid, mxvdw, ltpvdw, lstvdw, ls_vdw, vvdw, gvdw, rvdw, ltg,    &
              chge, rcut, alpha, epsq) bind(c)
              !!Modified by BBG: Third line from above is changed. Rho is
              !removed. 
     Use kinds_f90
     Use iso_c_binding

     Logical( c_bool ), Value :: ls_vdw

     Integer :: cuda_ilo, natms, mxlist, mxatms, mxatdm, mxmet, mxgrid, &
                mxvdw, ntpmet, ntpvdw, keyfce, imcon
     Integer, Dimension( 1:mxvdw )                        :: lstvdw, ltpvdw
     Integer, Dimension( -2:mxlist,1:mxatdm )              :: list
     Integer, Dimension( 1:mxatms )                       :: ltype, ltg
     Integer, Dimension( 1:mxmet )                        :: ltpmet, lstmet

     Real( Kind = wp )                                    :: rvdw, rcut, alpha, epsq
     Real( Kind = wp ), Dimension( 1:mxgrid,1:mxvdw )     :: vvdw, gvdw
     Real( Kind = wp ), Dimension( 1:9 )                  :: cell
     Real( Kind = wp ), Dimension( 1:mxatms )             :: xxx, yyy, zzz, rho, chge
     Real( Kind = wp ), Dimension( 1:mxgrid,1:mxmet,1:2 ) :: vmet, dmet
   End Subroutine two_body_forces_cuda_initialise

   Subroutine two_body_forces_cuda_invoke             &
              (fxx, fyy, fzz, stress, engmet, virmet, &
              engvdw, virvdw, engcpe_rl, vircpe_rl, safe) bind(c)
     Use kinds_f90
     Use setup_module, Only : mxatms
     Logical :: safe
     Real( Kind = wp )                   :: engmet, virmet, engvdw, virvdw, &
                                            engcpe_rl, vircpe_rl
     Real( Kind = wp ), Dimension( 1:9 )      :: stress
     Real( Kind = wp ), Dimension( 1:mxatms ) :: fxx, fyy, fzz
   End Subroutine two_body_forces_cuda_invoke

   Subroutine two_body_forces_cuda_finalise() bind(c)
   End Subroutine two_body_forces_cuda_finalise


   ! link_cell_pairs
   Subroutine link_cell_pairs_cuda_initialise                        &
              (natms,mxatms,ncells,mxlist,mxatdm,nsbcll,nlp,nlx,nly, &
              nlx0e, nly0e, nlz0e, nlx1s, nly1s, nlz1s,              &
              xxx,yyy,zzz,lfrzn,at_list,lct_start,list,      &
              nlp3,nix,niy,niz,rcsq,lbook,mxexcl,lexatm,ltg,megfrz) bind(c)
     Use kinds_f90
     Logical :: lbook
     Integer :: natms,mxatms,ncells,mxlist,mxatdm,nsbcll,nlp,nlx,nly,   &
                nlx0e, nly0e, nlz0e, nlx1s, nly1s, nlz1s, nlp3,mxexcl,megfrz
     Integer,           Dimension( 1:nlp3 )           :: nix,niy,niz
     Integer,           Dimension( 1:mxatms )          :: at_list, lfrzn, ltg
     Integer,           Dimension( 0:ncells+1 )        :: lct_start
     Integer,           Dimension( -2:mxlist,1:mxatdm ) :: list
     Integer,           Dimension( 0:mxexcl,1:mxatdm ) :: lexatm

     Real( Kind = wp ), Dimension( 1:mxatms )          :: xxx, yyy, zzz
     Real( Kind = wp )                                 :: rcsq
   End Subroutine link_cell_pairs_cuda_initialise

   Subroutine link_cell_pairs_cuda_invoke() bind(c)
   End Subroutine link_cell_pairs_cuda_invoke

   Subroutine link_cell_pairs_cuda_finalise() bind(c)
   End Subroutine link_cell_pairs_cuda_finalise

   Subroutine link_cell_pairs_cuda_pull_lists(natms) bind(c)
     Integer, value :: natms
   End Subroutine link_cell_pairs_cuda_pull_lists

   Subroutine link_cell_pairs_cuda_invoke_sort_atoms() bind(c)
   End Subroutine link_cell_pairs_cuda_invoke_sort_atoms

   Subroutine link_cell_pairs_cuda_invoke_remove_exclusions() bind(c)
   End Subroutine link_cell_pairs_cuda_invoke_remove_exclusions


   ! constraints_shake interfaces (called in both _lfv and _vv versions)
   Subroutine constraints_shake_cuda_initialise                        &
              (ntcons, mxcons, mxatms, natms,imcon,                       &
              lsi,lsa,lishp_con,lashp_con,mop,mxbuff,nlast,               &
              lstopt, lfrzn, listcon, listot, prmcon, weight,             &
              dxx, dyy, dzz, dxt, dyt, dzt, dt2, tstep2, tolnce, cell,    &
              xxx, yyy, zzz, xxt, yyt, zzt, strcon, is_vv) bind(c)
     Use kinds_f90
     Use iso_c_binding
     Use setup_module, Only : mxlshp, mxproc, mxtcon

     Logical( c_bool ), Value :: is_vv

     Integer :: ntcons, mxcons, mxatms, mxbuff, nlast, natms, imcon
     Integer, Dimension( 1:mxlshp )     :: lishp_con
     Integer, Dimension( 1:mxproc )     :: lashp_con
     Integer, Dimension( 1:mxatms )     :: lsi, lsa, lfrzn, listot
     Integer, Dimension( 1:26 )         :: mop
     Integer, Dimension( 0:2,1:mxcons ) :: lstopt
     Integer, Dimension( 0:2,1:mxcons ) :: listcon

     Real( Kind = wp )                        :: tstep2, tolnce
     Real( Kind = wp ), Dimension( 1:mxtcon ) :: prmcon
     Real( Kind = wp ), Dimension( 1:9 )      :: cell, strcon
     Real( Kind = wp ), Dimension( 1:mxatms ) :: xxx,yyy,zzz,xxt,yyt,zzt,weight
     Real( Kind = wp ), Dimension( 1:mxcons ) :: dxx,dyy,dzz,dxt,dyt,dzt,dt2
   End Subroutine constraints_shake_cuda_initialise

   Subroutine constraints_shake_cuda_gather_hs_scatter_dv() bind(c)
   End Subroutine constraints_shake_cuda_gather_hs_scatter_dv

   Subroutine constraints_shake_cuda_gather_dv_scatter_hs() bind(c)
   End Subroutine constraints_shake_cuda_gather_dv_scatter_hs

   Subroutine constraints_shake_cuda_invoke_bh() bind(c)
   End Subroutine constraints_shake_cuda_invoke_bh

   Subroutine constraints_shake_cuda_invoke_correct_positions() bind(c)
   End Subroutine constraints_shake_cuda_invoke_correct_positions

   Subroutine constraints_shake_cuda_copy_local_positions_to_host() bind(c)
   End Subroutine constraints_shake_cuda_copy_local_positions_to_host

   Subroutine constraints_shake_cuda_invoke_strcon_extraction() bind(c)
   End Subroutine constraints_shake_cuda_invoke_strcon_extraction

   Subroutine constraints_shake_cuda_finalise() bind(c)
   End Subroutine constraints_shake_cuda_finalise

   Logical Function constraints_shake_cuda_has_recv_list_been_constructed() bind(c)
   End Function constraints_shake_cuda_has_recv_list_been_constructed

   Logical Function dl_poly_cuda_constraints_shake_ntcons_enough_to_offload(ntcons) bind(c)
     Integer :: ntcons
   End Function dl_poly_cuda_constraints_shake_ntcons_enough_to_offload

   Subroutine constraints_shake_cuda_commit_receive_list(rilist, rilistcnt) bind(c)
     Use setup_module, Only : mxatms
     Integer :: rilistcnt
     Integer, Dimension( 0:mxatms-1 ) :: rilist
   End Subroutine constraints_shake_cuda_commit_receive_list

   Subroutine constraints_shake_cuda_invoke_th(esig) bind(c)
     Use kinds_f90
     Real( Kind = wp ) :: esig
   End Subroutine constraints_shake_cuda_invoke_th


! Initialisation Routines
   Subroutine dl_poly_cuda_initialise1(true, false) bind(c)
     Logical :: true, false
   End Subroutine dl_poly_cuda_initialise1

   Subroutine dl_poly_cuda_initialise2(mxatms,mxcons) bind(c)
     Integer :: mxatms, mxcons
   End Subroutine dl_poly_cuda_initialise2

   Subroutine dl_poly_cuda_finalise() bind(c)
   End Subroutine dl_poly_cuda_finalise


! Device offload and query routines
   Subroutine dl_poly_cuda_offload_set(offload_link_cell_pairs,           &
                                       offload_link_cell_pairs_re,        &
                                       offload_tbforces,                  &
                                       offload_constraints_shake,         &
                                       offload_metal_ld_compute,          &
                                       offload_ewald_spme_forces,         &
                                       offload_spme_forces,               &
                                       offload_bspgen,                    &
                                       offload_ewald_spme_forces_ccarray, &
                                       offload_ewald_spme_forces_cccharge) bind(c)
     Use iso_c_binding
     Logical( c_bool ), Value :: offload_link_cell_pairs,           &
                                 offload_link_cell_pairs_re,        &
                                 offload_tbforces,                  &
                                 offload_constraints_shake,         &
                                 offload_metal_ld_compute,          &
                                 offload_ewald_spme_forces,         &
                                 offload_spme_forces,               &
                                 offload_bspgen,                    &
                                 offload_ewald_spme_forces_ccarray, &
                                 offload_ewald_spme_forces_cccharge

   End Subroutine dl_poly_cuda_offload_set

   Logical Function dl_poly_cuda_is_cuda_capable () bind(c)
   End Function dl_poly_cuda_is_cuda_capable

   Logical Function dl_poly_cuda_offload_tbforces() bind(c)
   End Function dl_poly_cuda_offload_tbforces

   Logical Function dl_poly_cuda_offload_link_cell_pairs() bind(c)
   End Function dl_poly_cuda_offload_link_cell_pairs

   Logical Function dl_poly_cuda_offload_link_cell_pairs_re() bind(c)
   End Function dl_poly_cuda_offload_link_cell_pairs_re

   Logical Function dl_poly_cuda_offload_constraints_shake() bind(c)
   End Function dl_poly_cuda_offload_constraints_shake

   Logical Function dl_poly_cuda_offload_metal_ld_compute() bind(c)
   End Function dl_poly_cuda_offload_metal_ld_compute

   Logical Function dl_poly_cuda_offload_ewald_spme_forces() bind(c)
   End Function dl_poly_cuda_offload_ewald_spme_forces

   Logical Function dl_poly_cuda_offload_ewald_spme_forces_ccarray() bind(c)
   End Function dl_poly_cuda_offload_ewald_spme_forces_ccarray

   Logical Function dl_poly_cuda_offload_ewald_spme_forces_cccharge() bind(c)
   End Function dl_poly_cuda_offload_ewald_spme_forces_cccharge

   Logical Function dl_poly_cuda_offload_spme_forces() bind(c)
   End Function dl_poly_cuda_offload_spme_forces

   Logical Function dl_poly_cuda_offload_bspgen() bind(c)
   End Function dl_poly_cuda_offload_bspgen


! Timing Interfaces:
   Subroutine dump_timings() bind(c)
   End Subroutine dump_timings

   Subroutine start_timing_total() bind(c)
   End Subroutine start_timing_total

   Subroutine stop_timing_total() bind(c)
   End Subroutine stop_timing_total

   Subroutine start_timing_any() bind(c)
   End Subroutine start_timing_any

   Subroutine stop_timing_any() bind(c)
   End Subroutine stop_timing_any

   Subroutine start_timing_angles_forces() bind(c)
   End Subroutine start_timing_angles_forces

   Subroutine stop_timing_angles_forces() bind(c)
   End Subroutine stop_timing_angles_forces

   Subroutine start_timing_bonds_forces() bind(c)
   End Subroutine start_timing_bonds_forces

   Subroutine stop_timing_bonds_forces() bind(c)
   End Subroutine stop_timing_bonds_forces

   Subroutine start_timing_bspgen() bind(c)
   End Subroutine start_timing_bspgen

   Subroutine stop_timing_bspgen() bind(c)
   End Subroutine stop_timing_bspgen

   Subroutine start_timing_bspgen_cuda_k1() bind(c)
   End Subroutine start_timing_bspgen_cuda_k1

   Subroutine start_timing_bspgen_cuda_write() bind(c)
   End Subroutine start_timing_bspgen_cuda_write

   Subroutine stop_timing_bspgen_cuda_write() bind(c)
   End Subroutine stop_timing_bspgen_cuda_write

   Subroutine start_timing_constraints_shake() bind(c)
   End Subroutine start_timing_constraints_shake

   Subroutine stop_timing_constraints_shake() bind(c)
   End Subroutine stop_timing_constraints_shake

   Subroutine start_timing_constraints_shake_cuda_initialise() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_initialise

   Subroutine stop_timing_constraints_shake_cuda_initialise() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_initialise

   Subroutine start_timing_constraints_shake_cuda_k1_th() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_k1_th

   Subroutine start_timing_constraints_shake_cuda_k1_bh() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_k1_bh

   Subroutine stop_timing_constraints_shake_cuda_k1_bh() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_k1_bh

   Subroutine start_timing_constraints_shake_cuda_ccforces_k3() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_ccforces_k3

   Subroutine stop_timing_constraints_shake_cuda_ccforces_k3() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_ccforces_k3

   Subroutine start_timing_constraints_shake_cuda_read() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_read

   Subroutine stop_timing_constraints_shake_cuda_read() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_read

   Subroutine start_timing_constraints_shake_cuda_write() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_write

   Subroutine stop_timing_constraints_shake_cuda_write() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_write

   Subroutine start_timing_constraints_shake_cuda_ccforces_finalise() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_ccforces_finalise

   Subroutine stop_timing_constraints_shake_cuda_ccforces_finalise() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_ccforces_finalise

   Subroutine start_timing_constraints_shake_cuda_install_red_struct() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_install_red_struct

   Subroutine stop_timing_constraints_shake_cuda_install_red_struct() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_install_red_struct

   Subroutine start_timing_constraints_shake_cuda_cmeib_k1() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_cmeib_k1

   Subroutine stop_timing_constraints_shake_cuda_cmeib_k1() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_cmeib_k1

   Subroutine start_timing_constraints_shake_cuda_invoke_correct_positions() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_invoke_correct_positions

   Subroutine stop_timing_constraints_shake_cuda_invoke_correct_positions() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_invoke_correct_positions

   Subroutine start_timing_constraints_shake_cuda_gather_dv_scatter_hs() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_gather_dv_scatter_hs

   Subroutine stop_timing_constraints_shake_cuda_gather_dv_scatter_hs() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_gather_dv_scatter_hs

   Subroutine start_timing_constraints_shake_cuda_gather_hs_scatter_dv() bind(c)
   End Subroutine start_timing_constraints_shake_cuda_gather_hs_scatter_dv

   Subroutine stop_timing_constraints_shake_cuda_gather_hs_scatter_dv() bind(c)
   End Subroutine stop_timing_constraints_shake_cuda_gather_hs_scatter_dv

   Subroutine start_timing_update_shared_units() bind(c)
   End Subroutine start_timing_update_shared_units

   Subroutine stop_timing_update_shared_units() bind(c)
   End Subroutine stop_timing_update_shared_units

   Subroutine start_timing_core_shell_forces() bind(c)
   End Subroutine start_timing_core_shell_forces

   Subroutine stop_timing_core_shell_forces() bind(c)
   End Subroutine stop_timing_core_shell_forces

   Subroutine start_timing_coul_cp_forces() bind(c)
   End Subroutine start_timing_coul_cp_forces

   Subroutine stop_timing_coul_cp_forces() bind(c)
   End Subroutine stop_timing_coul_cp_forces

   Subroutine start_timing_coul_dddp_forces() bind(c)
   End Subroutine start_timing_coul_dddp_forces

   Subroutine stop_timing_coul_dddp_forces() bind(c)
   End Subroutine stop_timing_coul_dddp_forces

   Subroutine start_timing_coul_fscp_forces() bind(c)
   End Subroutine start_timing_coul_fscp_forces

   Subroutine stop_timing_coul_fscp_forces() bind(c)
   End Subroutine stop_timing_coul_fscp_forces

   Subroutine start_timing_coul_rfp_forces() bind(c)
   End Subroutine start_timing_coul_rfp_forces

   Subroutine stop_timing_coul_rfp_forces() bind(c)
   End Subroutine stop_timing_coul_rfp_forces

   Subroutine start_timing_dlpfft3() bind(c)
   End Subroutine start_timing_dlpfft3

   Subroutine stop_timing_dlpfft3() bind(c)
   End Subroutine stop_timing_dlpfft3

   Subroutine start_timing_ewald_excl_forces() bind(c)
   End Subroutine start_timing_ewald_excl_forces

   Subroutine stop_timing_ewald_excl_forces() bind(c)
   End Subroutine stop_timing_ewald_excl_forces

   Subroutine start_timing_ewald_real_forces() bind(c)
   End Subroutine start_timing_ewald_real_forces

   Subroutine stop_timing_ewald_real_forces() bind(c)
   End Subroutine stop_timing_ewald_real_forces

   Subroutine start_timing_ewald_spme_forces() bind(c)
   End Subroutine start_timing_ewald_spme_forces

   Subroutine stop_timing_ewald_spme_forces() bind(c)
   End Subroutine stop_timing_ewald_spme_forces

   Subroutine start_timing_ewald_spme_forces_any() bind(c)
   End Subroutine start_timing_ewald_spme_forces_any

   Subroutine stop_timing_ewald_spme_forces_any() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_any

   Subroutine start_timing_ewald_spme_forces_cuda_cccharge_k1() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_cccharge_k1

   Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_k1() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_k1

   Subroutine start_timing_ewald_spme_forces_cuda_cccharge_read() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_cccharge_read

   Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_read() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_read

   Subroutine start_timing_ewald_spme_forces_cuda_cccharge_write() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_cccharge_write

   Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_write() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_write

   Subroutine start_timing_ewald_spme_forces_cuda_cccharge_finalise() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_cccharge_finalise

   Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_finalise() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_cccharge_finalise

   Subroutine start_timing_ewald_spme_forces_ccarray() bind(c)
   End Subroutine start_timing_ewald_spme_forces_ccarray

   Subroutine stop_timing_ewald_spme_forces_ccarray() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_ccarray

   Subroutine start_timing_ewald_spme_forces_cuda_ccarray_preprocess() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_ccarray_preprocess

   Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_preprocess() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_preprocess

   Subroutine start_timing_ewald_spme_forces_cuda_ccarray_k1() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_ccarray_k1

   Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_k1() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_k1

   Subroutine start_timing_ewald_spme_forces_cuda_ccarray_k2() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_ccarray_k2

   Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_k2() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_k2

   Subroutine start_timing_ewald_spme_forces_cuda_ccarray_read() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_ccarray_read

   Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_read() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_read

   Subroutine start_timing_ewald_spme_forces_cuda_ccarray_write() bind(c)
   End Subroutine start_timing_ewald_spme_forces_cuda_ccarray_write

   Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_write() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_cuda_ccarray_write

   Subroutine start_timing_ewald_spme_forces_ccharge() bind(c)
   End Subroutine start_timing_ewald_spme_forces_ccharge

   Subroutine stop_timing_ewald_spme_forces_ccharge() bind(c)
   End Subroutine stop_timing_ewald_spme_forces_ccharge

   Subroutine start_timing_fft() bind(c)
   End Subroutine start_timing_fft

   Subroutine stop_timing_fft() bind(c)
   End Subroutine stop_timing_fft

   Subroutine start_timing_images() bind(c)
   End Subroutine start_timing_images

   Subroutine stop_timing_images() bind(c)
   End Subroutine stop_timing_images

   Subroutine start_timing_four_body_forces() bind(c)
   End Subroutine start_timing_four_body_forces

   Subroutine stop_timing_four_body_forces() bind(c)
   End Subroutine stop_timing_four_body_forces

   Subroutine start_timing_inversions_forces() bind(c)
   End Subroutine start_timing_inversions_forces

   Subroutine stop_timing_inversions_forces() bind(c)
   End Subroutine stop_timing_inversions_forces

   Subroutine start_timing_metal_forces() bind(c)
   End Subroutine start_timing_metal_forces

   Subroutine stop_timing_metal_forces() bind(c)
   End Subroutine stop_timing_metal_forces

   Subroutine start_timing_metal_ld_compute() bind(c)
   End Subroutine start_timing_metal_ld_compute

   Subroutine stop_timing_metal_ld_compute() bind(c)
   End Subroutine stop_timing_metal_ld_compute

   Subroutine start_timing_metal_ld_compute_cuda_write() bind(c)
   End Subroutine start_timing_metal_ld_compute_cuda_write

   Subroutine stop_timing_metal_ld_compute_cuda_write() bind(c)
   End Subroutine stop_timing_metal_ld_compute_cuda_write

   Subroutine start_timing_metal_ld_compute_cuda_read() bind(c)
   End Subroutine start_timing_metal_ld_compute_cuda_read

   Subroutine stop_timing_metal_ld_compute_cuda_read() bind(c)
   End Subroutine stop_timing_metal_ld_compute_cuda_read

   Subroutine start_timing_metal_ld_compute_cuda_k0() bind(c)
   End Subroutine start_timing_metal_ld_compute_cuda_k0

   Subroutine stop_timing_metal_ld_compute_cuda_k0() bind(c)
   End Subroutine stop_timing_metal_ld_compute_cuda_k0

   Subroutine start_timing_metal_ld_compute_cuda_k1() bind(c)
   End Subroutine start_timing_metal_ld_compute_cuda_k1

   Subroutine stop_timing_metal_ld_compute_cuda_k1() bind(c)
   End Subroutine stop_timing_metal_ld_compute_cuda_k1

   Subroutine start_timing_metal_ld_compute_cuda_k2() bind(c)
   End Subroutine start_timing_metal_ld_compute_cuda_k2

   Subroutine stop_timing_metal_ld_compute_cuda_k2() bind(c)
   End Subroutine stop_timing_metal_ld_compute_cuda_k2

   Subroutine start_timing_metal_ld_compute_cuda_k3() bind(c)
   End Subroutine start_timing_metal_ld_compute_cuda_k3

   Subroutine stop_timing_metal_ld_compute_cuda_k3() bind(c)
   End Subroutine stop_timing_metal_ld_compute_cuda_k3

   Subroutine start_timing_spme_forces() bind(c)
   End Subroutine start_timing_spme_forces

   Subroutine stop_timing_spme_forces() bind(c)
   End Subroutine stop_timing_spme_forces

   Subroutine start_timing_spme_forces_k1() bind(c)
   End Subroutine start_timing_spme_forces_k1

   Subroutine stop_timing_spme_forces_k1() bind(c)
   End Subroutine stop_timing_spme_forces_k1

   Subroutine start_timing_spme_forces_read() bind(c)
   End Subroutine start_timing_spme_forces_read

   Subroutine stop_timing_spme_forces_read() bind(c)
   End Subroutine stop_timing_spme_forces_read

   Subroutine start_timing_spme_forces_write() bind(c)
   End Subroutine start_timing_spme_forces_write

   Subroutine stop_timing_spme_forces_write() bind(c)
   End Subroutine stop_timing_spme_forces_write

   Subroutine start_timing_spme_forces_finalise() bind(c)
   End Subroutine start_timing_spme_forces_finalise

   Subroutine stop_timing_spme_forces_finalise() bind(c)
   End Subroutine stop_timing_spme_forces_finalise

   Subroutine start_timing_tersoff_forces() bind(c)
   End Subroutine start_timing_tersoff_forces

   Subroutine stop_timing_tersoff_forces() bind(c)
   End Subroutine stop_timing_tersoff_forces

   Subroutine start_timing_three_body_forces() bind(c)
   End Subroutine start_timing_three_body_forces

   Subroutine stop_timing_three_body_forces() bind(c)
   End Subroutine stop_timing_three_body_forces

   Subroutine start_timing_three_body_forces_cuda_scan() bind(c)
   End Subroutine start_timing_three_body_forces_cuda_scan

   Subroutine stop_timing_three_body_forces_cuda_scan() bind(c)
   End Subroutine stop_timing_three_body_forces_cuda_scan

   Subroutine start_timing_three_body_forces_cuda_reformat() bind(c)
   End Subroutine start_timing_three_body_forces_cuda_reformat

   Subroutine stop_timing_three_body_forces_cuda_reformat() bind(c)
   End Subroutine stop_timing_three_body_forces_cuda_reformat

   Subroutine start_timing_three_body_forces_cuda_write() bind(c)
   End Subroutine start_timing_three_body_forces_cuda_write

   Subroutine stop_timing_three_body_forces_cuda_write() bind(c)
   End Subroutine stop_timing_three_body_forces_cuda_write

   Subroutine start_timing_three_body_forces_cuda_kernel() bind(c)
   End Subroutine start_timing_three_body_forces_cuda_kernel

   Subroutine stop_timing_three_body_forces_cuda_kernel() bind(c)
   End Subroutine stop_timing_three_body_forces_cuda_kernel

   Subroutine start_timing_three_body_forces_cuda_read() bind(c)
   End Subroutine start_timing_three_body_forces_cuda_read

   Subroutine stop_timing_three_body_forces_cuda_read() bind(c)
   End Subroutine stop_timing_three_body_forces_cuda_read

   Subroutine start_timing_three_body_forces_cuda_finalise() bind(c)
   End Subroutine start_timing_three_body_forces_cuda_finalise

   Subroutine stop_timing_three_body_forces_cuda_finalise() bind(c)
   End Subroutine stop_timing_three_body_forces_cuda_finalise

   Subroutine start_timing_two_body_forces() bind(c)
   End Subroutine start_timing_two_body_forces

   Subroutine stop_timing_two_body_forces() bind(c)
   End Subroutine stop_timing_two_body_forces

   Subroutine start_timing_two_body_forces_any() bind(c)
   End Subroutine start_timing_two_body_forces_any

   Subroutine stop_timing_two_body_forces_any() bind(c)
   End Subroutine stop_timing_two_body_forces_any

   Subroutine start_timing_two_body_forces_cuda_km1() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_km1

   Subroutine stop_timing_two_body_forces_cuda_km1() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_km1

   Subroutine start_timing_two_body_forces_cuda_k0() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_k0

   Subroutine stop_timing_two_body_forces_cuda_k0() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_k0

   Subroutine start_timing_two_body_forces_cuda_k1() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_k1

   Subroutine stop_timing_two_body_forces_cuda_k1() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_k1

   Subroutine start_timing_two_body_forces_cuda_k2() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_k2

   Subroutine stop_timing_two_body_forces_cuda_k2() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_k2

   Subroutine start_timing_two_body_forces_cuda_k2_reduce() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_k2_reduce

   Subroutine stop_timing_two_body_forces_cuda_k2_reduce() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_k2_reduce

   Subroutine start_timing_two_body_forces_cuda_k3() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_k3

   Subroutine stop_timing_two_body_forces_cuda_k3() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_k3

   Subroutine start_timing_two_body_forces_cuda_finalise() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_finalise

   Subroutine stop_timing_two_body_forces_cuda_finalise() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_finalise

   Subroutine start_timing_two_body_forces_cuda_read() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_read

   Subroutine stop_timing_two_body_forces_cuda_read() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_read

   Subroutine start_timing_two_body_forces_cuda_write() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_write

   Subroutine stop_timing_two_body_forces_cuda_write() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_write

   Subroutine start_timing_two_body_forces_cuda_reformat() bind(c)
   End Subroutine start_timing_two_body_forces_cuda_reformat

   Subroutine stop_timing_two_body_forces_cuda_reformat() bind(c)
   End Subroutine stop_timing_two_body_forces_cuda_reformat

   Subroutine start_timing_link_cell_pairs() bind(c)
   End Subroutine start_timing_link_cell_pairs

   Subroutine stop_timing_link_cell_pairs() bind(c)
   End Subroutine stop_timing_link_cell_pairs

   Subroutine start_timing_link_cell_pairs_any() bind(c)
   End Subroutine start_timing_link_cell_pairs_any

   Subroutine stop_timing_link_cell_pairs_any() bind(c)
   End Subroutine stop_timing_link_cell_pairs_any

   Subroutine start_timing_link_cell_pairs_cuda_k0() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_k0

   Subroutine stop_timing_link_cell_pairs_cuda_k0() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_k0

   Subroutine start_timing_link_cell_pairs_cuda_k1() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_k1

   Subroutine stop_timing_link_cell_pairs_cuda_k1() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_k1

   Subroutine start_timing_link_cell_pairs_cuda_k2() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_k2

   Subroutine stop_timing_link_cell_pairs_cuda_k2() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_k2

   Subroutine start_timing_link_cell_pairs_cuda_finalise() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_finalise

   Subroutine stop_timing_link_cell_pairs_cuda_finalise() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_finalise

   Subroutine start_timing_link_cell_pairs_cuda_read() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_read

   Subroutine stop_timing_link_cell_pairs_cuda_read() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_read

   Subroutine start_timing_link_cell_pairs_cuda_write() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_write

   Subroutine stop_timing_link_cell_pairs_cuda_write() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_write

   Subroutine start_timing_link_cell_pairs_cuda_remove_excluded() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_remove_excluded

   Subroutine stop_timing_link_cell_pairs_cuda_remove_excluded() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_remove_excluded

   Subroutine start_timing_link_cell_pairs_cuda_sort_atoms() bind(c)
   End Subroutine start_timing_link_cell_pairs_cuda_sort_atoms

   Subroutine stop_timing_link_cell_pairs_cuda_sort_atoms() bind(c)
   End Subroutine stop_timing_link_cell_pairs_cuda_sort_atoms

   Subroutine start_timing_link_cell_pairs_sparse_list_transfer() bind(c)
   End Subroutine start_timing_link_cell_pairs_sparse_list_transfer

   Subroutine stop_timing_link_cell_pairs_sparse_list_transfer() bind(c)
   End Subroutine stop_timing_link_cell_pairs_sparse_list_transfer

   Subroutine start_timing_link_cell_pairs_sparse_list_transfer_gxscan() bind(c)
   End Subroutine start_timing_link_cell_pairs_sparse_list_transfer_gxscan

   Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_gxscan() bind(c)
   End Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_gxscan

   Subroutine start_timing_link_cell_pairs_sparse_list_transfer_pack() bind(c)
   End Subroutine start_timing_link_cell_pairs_sparse_list_transfer_pack

   Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_pack() bind(c)
   End Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_pack

   Subroutine start_timing_link_cell_pairs_sparse_list_transfer_reorder() bind(c)
   End Subroutine start_timing_link_cell_pairs_sparse_list_transfer_reorder

   Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_reorder() bind(c)
   End Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_reorder

   Subroutine start_timing_link_cell_pairs_sparse_list_transfer_read() bind(c)
   End Subroutine start_timing_link_cell_pairs_sparse_list_transfer_read

   Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_read() bind(c)
   End Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_read

   Subroutine start_timing_link_cell_pairs_sparse_list_transfer_unpack() bind(c)
   End Subroutine start_timing_link_cell_pairs_sparse_list_transfer_unpack

   Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_unpack() bind(c)
   End Subroutine stop_timing_link_cell_pairs_sparse_list_transfer_unpack

   Subroutine start_timing_vdw_forces() bind(c)
   End Subroutine start_timing_vdw_forces

   Subroutine stop_timing_vdw_forces() bind(c)
   End Subroutine stop_timing_vdw_forces
End Interface

! Miscellaneous functions/subroutines needed for the CUDA port

Contains
  Function remap_lb_real_2d(array, lb1, lb2) result(ptr)
    Integer, intent(in) :: lb1, lb2
    Real( Kind = wp), dimension(lb1:,lb2:), intent(in), target :: array
    Real( Kind = wp), dimension(:,:), pointer :: ptr
    ptr => array
  End Function remap_lb_real_2d

  Function remap_lb_real_3d(array, lb1, lb2, lb3) result(ptr)
    Integer, intent(in) :: lb1, lb2, lb3
    Real( Kind = wp), dimension(lb1:,lb2:,lb3:), intent(in), target :: array
    Real( Kind = wp), dimension(:,:,:), pointer :: ptr
    ptr => array
  End Function remap_lb_real_3d

  Function remap_lb_complex_3d(array, lb1, lb2, lb3) result(ptr)
    Integer, intent(in) :: lb1, lb2, lb3
    Complex( Kind = wp), dimension(lb1:,lb2:,lb3:), intent(in), target :: array
    Complex( Kind = wp), dimension(:,:,:), pointer :: ptr
    ptr => array
  End Function remap_lb_complex_3d

  Function remap_lb_integer_2d(array, lb1, lb2) result(ptr)
    Integer, intent(in) :: lb1, lb2
    Integer( Kind = ip), dimension(lb1:,lb2:), intent(in), target :: array
    Integer( Kind = ip), dimension(:,:), pointer :: ptr
    ptr => array
  End Function remap_lb_integer_2d

! At initialisation, check certain parameters to see if offloading to
! GPU is possible (i.e. if the relevant functionality has been implemented).

  Subroutine dl_poly_cuda_check_offload_conditions(keyfce, imcon)
    Use iso_c_binding
    Use metal_module, Only : ld_met,ntpmet
    Use vdw_module  , Only : ld_vdw
    Use setup_module, Only : mxspl,nrite

    Implicit None

    Integer               :: keypot
    Integer, Intent( In ) :: keyfce, imcon

    Logical( c_bool ) :: offload_link_cell_pairs    = .true., & !Modified by BBG: Original true. 
                         offload_link_cell_pairs_re = .true., & !Modified by BBG: Original false.
                         offload_tbforces           = .false., & !Modified by BBG: Original true.
                         offload_constraints_shake  = .false., & !Modified by BBG: Original true.
                         offload_metal_ld_compute   = .false., & !Modified by BBG: Original true.
                         offload_ewald_spme_forces  = .true., & !Modified by BBG: Original true.
                         offload_spme_forces        = .true., & !Modified by BBG: Original true.
                         offload_bspgen             = .true., & !Modified by BBG: Original true.
                         offload_ewald_spme_forces_ccarray  = .true., & !Modified by BBG: Original true.
                         offload_ewald_spme_forces_cccharge = .true. !Modified by BBG: Original true.

    Write(nrite,*)

    ! Compute keypot variable
    Call metal_ld_compute_get_keypot(keypot)

    ! Link_cell_pairs:
    ! malysaght250112: offload_link_cell_pairs_re should be set to .false. until fix is implemented
    If ( .not. offload_link_cell_pairs_re ) Then
       Write(nrite,'(1x,a,/)') 'Disabling CUDA acceleration for link_cell_pairs_remove_exclusions (if applicable)'
    Endif

    ! Two-body forces:
    If (ld_vdw) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (two_body_forces): ld_vdw is true'
       offload_tbforces = .false.
    Endif

    If (ld_met) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (two_body_forces): ld_met is true'
       offload_tbforces = .false.
    End If

    If ( keyfce /= 0 .and. keyfce /= 2 ) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (two_body_forces): found unsupported keyfce value'
       Write(nrite,'(1x,a)') 'Only keyfce = 0,2 are supported'
       offload_tbforces = .false.
    Endif

    If ( imcon /= 2 .and. imcon /= 3 ) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (two_body_forces): found unsupported imcon value'
       Write(nrite,'(1x,a)') 'Only imcon = 2,3 are supported'
       offload_tbforces = .false.
    Endif

    If ( keypot == 0 .and. ntpmet > 0 ) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (two_body_forces): found unsupported keypot (0) and ntpmet (>0) values'
       offload_tbforces = .false.
    Endif

    If ( .not. offload_tbforces ) Then
       Write(nrite,'(1x,a,/)') 'Disabling CUDA acceleration for two_body_forces (if applicable)'
    Endif


    ! Constraints Shake
    If ( imcon /= 2 .and. imcon /= 3 ) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (constraints_shake): found unsupported imcon value'
       Write(nrite,'(1x,a)') 'Only imcon = 2,3 are supported'
       offload_constraints_shake = .false.
    Endif

    If ( .not. offload_constraints_shake ) Then
       Write(nrite,'(1x,a,/)') 'Disabling CUDA acceleration for constraints_shake (if applicable)'
    Endif


    ! Metal LD Compute
    If ( ld_met ) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (metal_ld_compute): ld_met is true'
       offload_metal_ld_compute = .false.
    End If

    If ( keypot == 0 ) Then
       Write(nrite,'(1x,a)') 'CUDA Port Warning (metal_ld_compute): found unsupported keypot value (0)'
       offload_metal_ld_compute = .false.
    Endif

    If ( .not. offload_metal_ld_compute ) Then
       Write(nrite,'(1x,a,/)') 'Disabling CUDA acceleration for metal_ld_compute (if applicable)'
    Endif

    ! Ewald_spme_forces (top-level)
    ! no unimplemented functionality at top level

    ! Ewald_spme_forces -> bspgen (spme_container.f90)
    ! Ewald_spme_forces -> spme_forces (ewald_spme_forces.f90)
    If ( mxspl /= 8 ) Then
       Write(nrite,'(1x,a  )') 'CUDA Port Warning (spme_forces): found non-default mxspl value (mxspl /= 8)'
       Write(nrite,'(1x,a,/)') 'Disabling CUDA acceleration for spme_forces (if applicable)'
       offload_spme_forces = .false.
    Endif

    ! Ewald_spme_forces -> ewald_spme_forces_ccharge

    ! Ewald_spme_forces -> ewald_spme_forces_ccarray

    ! bspgen (spme_container.f90) - Always offloaded

    ! Set some static truth variables so that C routines can
    ! check what will be offloaded (dl_poly_init_cu.cu)
    Call dl_poly_cuda_offload_set                &
             (offload_link_cell_pairs,           &
              offload_link_cell_pairs_re,        &
              offload_tbforces,                  &
              offload_constraints_shake,         &
              offload_metal_ld_compute,          &
              offload_ewald_spme_forces,         &
              offload_spme_forces,               &
              offload_bspgen,                    &
              offload_ewald_spme_forces_ccarray, &
              offload_ewald_spme_forces_cccharge)

  End Subroutine dl_poly_cuda_check_offload_conditions

End Module dl_poly_cuda_module


