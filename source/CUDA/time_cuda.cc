/*
* DL_POLY_4 NVIDIA GPU & OpenMP Port
* Irish Centre for High-End Computing (ICHEC)
* http://www.ichec.ie
*
* Developed by Christos Kartsaklis (christos.kartsaklis@ichec.ie) in
* collaboration with I.T. Todorov (i.t.todorov@dl.ac.uk) and
* W. Smith (w.smith@dl.ac.uk) from STFC Daresbury Laboratory.
*
* Distributed under the same license that the original, unmodified,
* DL_POLY_4 is. You should have received these sources from the
* STFC Daresbury Laboratory.
*/

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>

#define INCLUDED_BY_TIME_C 1
#include "dl_poly_cu.h"

std::list<int>* constraints_shake_cuda_get_recvlisttmp() {
  static std::list<int> sTheList;
  return (&sTheList);
}

typedef enum {
  angles_forces=0,
  any,
  bonds_forces,
  bspgen,
  bspgen_cuda_k1,
  bspgen_cuda_write,
  constraints_shake,
  constraints_shake_cuda_initialise,
  constraints_shake_cuda_k1_th,
  constraints_shake_cuda_k1_bh,
  constraints_shake_cuda_ccforces_k3,
  constraints_shake_cuda_read,
  constraints_shake_cuda_write,
  constraints_shake_cuda_ccforces_finalise,
  constraints_shake_cuda_install_red_struct,
  constraints_shake_cuda_cmeib_k1,
  constraints_shake_cuda_invoke_correct_positions,
  constraints_shake_cuda_gather_dv_scatter_hs,
  constraints_shake_cuda_gather_hs_scatter_dv,
  update_shared_units,
  core_shell_forces,
  coul_cp_forces,
  coul_dddp_forces,
  coul_fscp_forces,
  coul_rfp_forces,
  dlpfft3,
  dlpfft3_cuda_k1,
  dlpfft3_cuda_k2,
  dlpfft3_cuda_read,
  dlpfft3_cuda_write,
  ewald_excl_forces,
  ewald_real_forces,
  ewald_spme_forces,
  ewald_spme_forces_any,
  ewald_spme_forces_dcft3,
  ewald_spme_forces_cuda_cccharge_k1,
  ewald_spme_forces_cuda_cccharge_read,
  ewald_spme_forces_cuda_cccharge_write,
  ewald_spme_forces_cuda_cccharge_finalise,
  ewald_spme_forces_ccarray,
  ewald_spme_forces_cuda_ccarray_preprocess,
  ewald_spme_forces_cuda_ccarray_k1,
  ewald_spme_forces_cuda_ccarray_k2,
  ewald_spme_forces_cuda_ccarray_read,
  ewald_spme_forces_cuda_ccarray_write,
  ewald_spme_forces_ccharge,
  fft,
  four_body_forces,
  images,
  inversions_forces,
  link_cell_pairs,
  link_cell_pairs_any,
  link_cell_pairs_cuda_k0,
  link_cell_pairs_cuda_k1,
  link_cell_pairs_cuda_k2,
  link_cell_pairs_cuda_finalise,
  link_cell_pairs_cuda_read,
  link_cell_pairs_cuda_write,
  link_cell_pairs_cuda_remove_excluded,
  link_cell_pairs_cuda_sort_atoms,
  link_cell_pairs_sparse_list_transfer,
  link_cell_pairs_sparse_list_transfer_gxscan,
  link_cell_pairs_sparse_list_transfer_pack,
  link_cell_pairs_sparse_list_transfer_reorder,
  link_cell_pairs_sparse_list_transfer_read,
  link_cell_pairs_sparse_list_transfer_unpack,
  metal_forces,
  metal_ld_compute,
  metal_ld_compute_cuda_write,
  metal_ld_compute_cuda_read,
  metal_ld_compute_cuda_k0,
  metal_ld_compute_cuda_k1,
  metal_ld_compute_cuda_k2,
  metal_ld_compute_cuda_k3,
  spme_forces,
  spme_forces_k1,
  spme_forces_read,
  spme_forces_finalise,
  spme_forces_write,
  tersoff_forces,
  three_body_forces,
  three_body_forces_cuda_kernel,
  three_body_forces_cuda_read,
  three_body_forces_cuda_scan,
  three_body_forces_cuda_reformat,
  three_body_forces_cuda_write,
  three_body_forces_cuda_finalise,
  two_body_forces,
  two_body_forces_any,
  two_body_forces_cuda_kernel,
  two_body_forces_cuda_km1,
  two_body_forces_cuda_k0,
  two_body_forces_cuda_k1,
  two_body_forces_cuda_k2,
  two_body_forces_cuda_k2_reduce,
  two_body_forces_cuda_k3,
  two_body_forces_cuda_finalise,
  two_body_forces_cuda_read,
  two_body_forces_cuda_write,
  two_body_forces_cuda_reformat,
  vdw_forces,
  total,
  num_
} component_t;


struct timeval mStarted[num_];
double         mTotal[num_];
unsigned long long mCalls[num_];

struct {
  struct timeval mStarted;
  double mTimeSpent;
  long long mBytes;
} sIO[2]; // read/write


void start_timing(int aIndex) {
#if (!DISABLE_TIMING)
  static int lIsFirstCall = 1;
  if (lIsFirstCall) {
    int lI;
    for(lI=0 ; lI<num_ ; lI++) {
      mTotal[lI] = 0.0;
      mCalls[lI] = 0ULL;
    }
    lIsFirstCall = 0;
  }
  gettimeofday(&mStarted[aIndex], NULL);
  mCalls[aIndex]++;
#endif
}

void stop_timing(int aIndex) {
#if (!DISABLE_TIMING)
  struct timeval lTV;
  gettimeofday(&lTV, NULL);
  mTotal[aIndex] += secsfromtimeval(lTV) - secsfromtimeval(mStarted[aIndex]);
#endif
}


#define define_timer(NAME)\
  extern "C" void start_timing_##NAME##() { start_timing(NAME); }\
  extern "C" void stop_timing_##NAME##() { stop_timing(NAME); }

define_timer(total);
define_timer(any);
define_timer(angles_forces);
define_timer(bonds_forces);
define_timer(bspgen);
define_timer(bspgen_cuda_k1);
define_timer(bspgen_cuda_write);
define_timer(constraints_shake);
define_timer(constraints_shake_cuda_initialise);
define_timer(constraints_shake_cuda_k1_th);
define_timer(constraints_shake_cuda_k1_bh);
define_timer(constraints_shake_cuda_read);
define_timer(constraints_shake_cuda_write);
define_timer(constraints_shake_cuda_ccforces_finalise);
define_timer(constraints_shake_cuda_install_red_struct);
define_timer(constraints_shake_cuda_cmeib_k1);
define_timer(constraints_shake_cuda_ccforces_k3);
define_timer(constraints_shake_cuda_invoke_correct_positions);
define_timer(constraints_shake_cuda_gather_dv_scatter_hs);
define_timer(constraints_shake_cuda_gather_hs_scatter_dv);
define_timer(update_shared_units);
define_timer(core_shell_forces);
define_timer(coul_cp_forces);
define_timer(coul_dddp_forces);
define_timer(coul_fscp_forces);
define_timer(coul_rfp_forces);
define_timer(dlpfft3);
define_timer(dlpfft3_cuda_k1);
define_timer(dlpfft3_cuda_k2);
define_timer(ewald_spme_forces_cuda_ccarray_preprocess);
define_timer(dlpfft3_cuda_read);
define_timer(dlpfft3_cuda_write);
define_timer(ewald_excl_forces);
define_timer(ewald_real_forces);
define_timer(ewald_spme_forces);
define_timer(ewald_spme_forces_any);
define_timer(ewald_spme_forces_dcft3);
define_timer(ewald_spme_forces_cuda_cccharge_k1);
define_timer(ewald_spme_forces_cuda_cccharge_read);
define_timer(ewald_spme_forces_cuda_cccharge_write);
define_timer(ewald_spme_forces_cuda_cccharge_finalise);
define_timer(ewald_spme_forces_ccharge);
define_timer(ewald_spme_forces_ccarray);
define_timer(ewald_spme_forces_cuda_ccarray_k1);
define_timer(ewald_spme_forces_cuda_ccarray_k2);
define_timer(ewald_spme_forces_cuda_ccarray_read);
define_timer(ewald_spme_forces_cuda_ccarray_write);
define_timer(fft);
define_timer(four_body_forces);
define_timer(images);
define_timer(inversions_forces);
define_timer(link_cell_pairs);
define_timer(link_cell_pairs_any);
define_timer(link_cell_pairs_cuda_k0);
define_timer(link_cell_pairs_cuda_k1);
define_timer(link_cell_pairs_cuda_k2);
define_timer(link_cell_pairs_cuda_finalise);
define_timer(link_cell_pairs_cuda_read);
define_timer(link_cell_pairs_cuda_write);
define_timer(link_cell_pairs_cuda_remove_excluded);
define_timer(link_cell_pairs_cuda_sort_atoms);
define_timer(link_cell_pairs_sparse_list_transfer);
define_timer(link_cell_pairs_sparse_list_transfer_gxscan);
define_timer(link_cell_pairs_sparse_list_transfer_pack);
define_timer(link_cell_pairs_sparse_list_transfer_reorder);
define_timer(link_cell_pairs_sparse_list_transfer_read);
define_timer(link_cell_pairs_sparse_list_transfer_unpack);
define_timer(metal_forces);
define_timer(metal_ld_compute);
define_timer(metal_ld_compute_cuda_write);
define_timer(metal_ld_compute_cuda_read);
define_timer(metal_ld_compute_cuda_k0);
define_timer(metal_ld_compute_cuda_k1);
define_timer(metal_ld_compute_cuda_k2);
define_timer(metal_ld_compute_cuda_k3);
define_timer(spme_forces);
define_timer(spme_forces_k1);
define_timer(spme_forces_read);
define_timer(spme_forces_finalise);
define_timer(spme_forces_write);
define_timer(tersoff_forces);
define_timer(three_body_forces_cuda_kernel);
define_timer(three_body_forces);
define_timer(three_body_forces_cuda_reformat);
define_timer(three_body_forces_cuda_scan);
define_timer(three_body_forces_cuda_read);
define_timer(three_body_forces_cuda_write);
define_timer(three_body_forces_cuda_finalise);
define_timer(two_body_forces);
define_timer(two_body_forces_any);
define_timer(two_body_forces_cuda_kernel);
define_timer(two_body_forces_cuda_km1);
define_timer(two_body_forces_cuda_k0);
define_timer(two_body_forces_cuda_k1);
define_timer(two_body_forces_cuda_k2);
define_timer(two_body_forces_cuda_k2_reduce);
define_timer(two_body_forces_cuda_k3);
define_timer(two_body_forces_cuda_finalise);
define_timer(two_body_forces_cuda_read);
define_timer(two_body_forces_cuda_write);
define_timer(two_body_forces_cuda_reformat);
define_timer(vdw_forces);

extern "C" int dl_poly_cuda_process();

extern "C" void dump_timings() {
#if (!DISABLE_TIMING)
  stop_timing_total();
  printf("%3d: %35s : %4s %lf\n", dl_poly_cuda_process(), "", "", mTotal[total]);
#define dump_timer(NAME)\
  printf("%3d: %70s : %10llu %3d%% %06lf %06lf\n", dl_poly_cuda_process(), #NAME, mCalls[NAME],\
         (int)(100.0 * (mTotal[NAME] / mTotal[total])), mTotal[NAME], mCalls[NAME]>0 ? (mTotal[NAME]/mCalls[NAME]) : 0.0)
  dump_timer(angles_forces);
  dump_timer(any);
  dump_timer(bonds_forces);
  dump_timer(bspgen);
  dump_timer(bspgen_cuda_k1);
  dump_timer(bspgen_cuda_write);
  dump_timer(constraints_shake);
  dump_timer(constraints_shake_cuda_initialise);
  dump_timer(constraints_shake_cuda_k1_th);
  dump_timer(constraints_shake_cuda_k1_bh);
  dump_timer(constraints_shake_cuda_ccforces_k3);
  dump_timer(constraints_shake_cuda_read);
  dump_timer(constraints_shake_cuda_write);
  dump_timer(constraints_shake_cuda_ccforces_finalise);
  dump_timer(constraints_shake_cuda_install_red_struct);
  dump_timer(constraints_shake_cuda_cmeib_k1);
  dump_timer(constraints_shake_cuda_invoke_correct_positions);
  dump_timer(constraints_shake_cuda_gather_dv_scatter_hs);
  dump_timer(constraints_shake_cuda_gather_hs_scatter_dv);
  dump_timer(update_shared_units);
  dump_timer(core_shell_forces);
  dump_timer(coul_cp_forces);
  dump_timer(coul_dddp_forces);
  dump_timer(coul_fscp_forces);
  dump_timer(coul_rfp_forces);
  dump_timer(dlpfft3);
  dump_timer(ewald_excl_forces);
  dump_timer(ewald_real_forces);
  dump_timer(ewald_spme_forces);
  dump_timer(ewald_spme_forces_any);
  dump_timer(ewald_spme_forces_cuda_cccharge_k1);
  dump_timer(ewald_spme_forces_cuda_cccharge_read);
  dump_timer(ewald_spme_forces_cuda_cccharge_write);
  dump_timer(ewald_spme_forces_cuda_cccharge_finalise);
  dump_timer(ewald_spme_forces_ccarray);
  dump_timer(ewald_spme_forces_cuda_ccarray_preprocess);
  dump_timer(ewald_spme_forces_cuda_ccarray_k1);
  dump_timer(ewald_spme_forces_cuda_ccarray_k2);
  dump_timer(ewald_spme_forces_cuda_ccarray_read);
  dump_timer(ewald_spme_forces_cuda_ccarray_write);
  dump_timer(ewald_spme_forces_ccharge);
  dump_timer(fft);
  dump_timer(images);
  dump_timer(four_body_forces);
  dump_timer(inversions_forces);
  dump_timer(metal_forces);
  dump_timer(metal_ld_compute);
  dump_timer(metal_ld_compute_cuda_write);
  dump_timer(metal_ld_compute_cuda_read);
  dump_timer(metal_ld_compute_cuda_k0);
  dump_timer(metal_ld_compute_cuda_k1);
  dump_timer(metal_ld_compute_cuda_k2);
  dump_timer(metal_ld_compute_cuda_k3);
  dump_timer(spme_forces);
  dump_timer(spme_forces_k1);
  dump_timer(spme_forces_read);
  dump_timer(spme_forces_write);
  dump_timer(spme_forces_finalise);
  dump_timer(tersoff_forces);
  dump_timer(three_body_forces);
  dump_timer(three_body_forces_cuda_scan);
  dump_timer(three_body_forces_cuda_reformat);
  dump_timer(three_body_forces_cuda_write);
  dump_timer(three_body_forces_cuda_kernel);
  dump_timer(three_body_forces_cuda_read);
  dump_timer(three_body_forces_cuda_finalise);
  dump_timer(two_body_forces);
  dump_timer(two_body_forces_any);
  dump_timer(two_body_forces_cuda_km1);
  dump_timer(two_body_forces_cuda_k0);
  dump_timer(two_body_forces_cuda_k1);
  dump_timer(two_body_forces_cuda_k2);
  dump_timer(two_body_forces_cuda_k2_reduce);
  dump_timer(two_body_forces_cuda_k3);
  dump_timer(two_body_forces_cuda_finalise);
  dump_timer(two_body_forces_cuda_read);
  dump_timer(two_body_forces_cuda_write);
  dump_timer(two_body_forces_cuda_reformat);
  dump_timer(link_cell_pairs);
  dump_timer(link_cell_pairs_any);
  dump_timer(link_cell_pairs_cuda_k0);
  dump_timer(link_cell_pairs_cuda_k1);
  dump_timer(link_cell_pairs_cuda_k2);
  dump_timer(link_cell_pairs_cuda_finalise);
  dump_timer(link_cell_pairs_cuda_read);
  dump_timer(link_cell_pairs_cuda_write);
  dump_timer(link_cell_pairs_cuda_remove_excluded);
  dump_timer(link_cell_pairs_cuda_sort_atoms);
  dump_timer(link_cell_pairs_sparse_list_transfer);
  dump_timer(link_cell_pairs_sparse_list_transfer_gxscan);
  dump_timer(link_cell_pairs_sparse_list_transfer_pack);
  dump_timer(link_cell_pairs_sparse_list_transfer_reorder);
  dump_timer(link_cell_pairs_sparse_list_transfer_read);
  dump_timer(link_cell_pairs_sparse_list_transfer_unpack);
  dump_timer(vdw_forces);
  printf("%3d: bytes read: %12lld @ %6lldmbytes/sec\n", dl_poly_cuda_process(),
         sIO[0].mBytes, (long long)(sIO[0].mTimeSpent>0.0 ? (((double) sIO[0].mBytes)/1000000.0) / sIO[0].mTimeSpent : 0.0));
#endif
}

