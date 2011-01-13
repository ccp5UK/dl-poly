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

#ifndef _DL_POLY_CUDA_H
#define _DL_POLY_CUDA_H

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>


std::list<int>* constraints_shake_cuda_get_recvlisttmp();

#define secsfromtimeval(T)					\
  (((double) (T).tv_sec) + ((double) (T).tv_usec)/1000000.0)

/* This should match the precision of the real in the F90 dl_poly sources; it
 * is important as host<->device memory copies _will_not_ convert data.
 */
#if (!defined(CFG_DOUBLE_PRECISION))
#  define CFG_DOUBLE_PRECISION         1
#endif

#if (CFG_DOUBLE_PRECISION)
#  define CFG_DL_POLY_REAL             double
#  define real                         double
#  define real_tex                     int2
#else
#  define CFG_DL_POLY_REAL             float
#  define real                         float
#  define real_tex                     float
#endif

template<typename T_> struct complex {
  T_ mReal;
  T_ mImag;
};

#if (!defined(CFG_COMPUTE_MAJOR) || !defined(CFG_COMPUTE_MINOR))
#  error "either of CFG_COMPUTE_{MAJOR, MINOR} haven't been defined"
#endif

#if ((CFG_COMPUTE_MAJOR == 1) && (CFG_COMPUTE_MINOR < 3))
#  error "Compilation error: Compute capability (CFG_COMPUTE_MAJOR.CFG_COMPUTE_MINOR) \
must be set to 1.3 or greater in Makefile"
#endif

#define CFG_OVERLAP_WITH_HOST        1

/* For testing purposes only -- forces those who include dl_poly_cuda_common.cu
 * to use macros/functions that do not guarantee multiply-add separation.
 *
 * Update: for compute capability >=2.0, nvcc compiles to IEEE-compilant code
 * by default.
 */
#define CFG_DISABLE_ACCURATE_VERSIONS (CFG_COMPUTE_MAJOR>=2)

/* GPU Architecture Specification
 */
#define CFG_GRID_DIMENSION_MAX_SIZE                ((64*1024)-1)
#define CFG_SHARED_MEMORY_COMPILER_HEADROOM        (128)
#define CFG_SHARED_MEMORY_SIZE                     ((16*1024) - CFG_SHARED_MEMORY_COMPILER_HEADROOM)
#define CFG_WARP_SIZE                              32

// for compute >= 2.0
#define CFG_UNIFIED_ADDRESS_SPACE    (CFG_COMPUTE_MAJOR >= 2)

#if !CFG_UNIFIED_ADDRESS_SPACE
#  define DECLARE_DYNAMIC_SHARED(TYPE) extern __shared__ TYPE shared[]
#else
#  define DECLARE_DYNAMIC_SHARED(TYPE)\
  extern __shared__ TYPE _shared[];    \
  volatile TYPE * shared = _shared
#endif

/* If set to 1 and limits permit so, some kernels will employ 24bit
 * integer arithmetics as for certain architectures they have better
 * throughput.
 */
#define CFG_USE_24BIT_IAOPS          1

#define PARALLEL_BUILD               (1)

# if (PARALLEL_BUILD)
#  ifndef MPICH_IGNORE_CXX_SEEK
#    define MPICH_IGNORE_CXX_SEEK
#  endif
#  include <mpi.h>
# endif

template<typename T_> class host_pinned_buffer {
 private:
  T_ *mPointer;
  int mLength;
  unsigned int mFlags;
 public:
  inline host_pinned_buffer() { mPointer=NULL; mLength=0; mFlags=0; }

  inline void realloc(int aRetain, int aNewLength) {

    if (aNewLength>mLength) {
      if (aRetain) {
        CUDA_SAFE_CALL(cudaFreeHost(mPointer));
        mLength = aNewLength;
        CUDA_SAFE_CALL(cudaHostAlloc(&mPointer, mLength*sizeof(T_), mFlags));
      } else {
        T_ *lTmpPointer;
        CUDA_SAFE_CALL(cudaHostAlloc(&lTmpPointer, aNewLength*sizeof(T_), mFlags));
        memcpy(lTmpPointer, mPointer, mLength*sizeof(T_));
        CUDA_SAFE_CALL(cudaFreeHost(mPointer));
        mPointer = lTmpPointer;
        mLength = aNewLength;
      }
    }
  }
  inline ~host_pinned_buffer() {
    if (mLength>0) {
      CUDA_SAFE_CALL(cudaFreeHost(mPointer));
    }
  }
  inline T_* pointer() { return (mPointer); }
  inline int length() { return(mLength); }
};

extern "C" void dl_poly_cuda_reduce_constant_memory_by(const char *aSymbol);
extern "C" int dl_poly_cuda_process();

extern "C" int dl_poly_cuda_fortran_true();
extern "C" int dl_poly_cuda_fortran_false();
extern "C" int dl_poly_cuda_is_cuda_capable();
extern "C" int dl_poly_cuda_is_cuda_capable();
extern "C" int dl_poly_cuda_offload_tbforces();
extern "C" int dl_poly_cuda_offload_link_cell_pairs();
extern "C" int dl_poly_cuda_offload_ewald_spme_forces();
extern "C" int dl_poly_cuda_offload_spme_forces();
extern "C" int dl_poly_cuda_offload_metal_ld_compute();

extern "C" void* dl_poly_cuda_get_buffer_constraints_shake_ijso();
extern "C" void* dl_poly_cuda_get_buffer_constraints_shake_ijkpos();
extern "C" void* dl_poly_cuda_get_buffer_constraints_shake_prmcon_k();

template<typename T_> T_ dl_poly_cuda_getenv(const char *aName, T_ aDefaultValue);
template<> int dl_poly_cuda_getenv(const char *aName, int aDefault);


/* Splits the interval [aB,aE] into aN (at most) disjoint intervals based on the
 * percentages specified in aR. The result is placed in aOI; if an interval i
 * (0<=i<aN) has been created, then aOI[3*i]==1 (otherwise zero) and the interval
 * is [aOI[3*i+1], aOI[3*i+2]. If an interval has been created, and i>0 then it
 * will hold that aOI[3*i+1]==aOI[3*(i-1)+2]+1.
 */
void split_range(int aN, int aB, int aE, float *aR, int *aOI);

extern "C" int link_cell_pairs_cuda_is_in_valid_context();
extern "C" void* link_cell_pairs_cuda_get_list();
extern "C" void link_cell_pairs_cuda_pull_lists(int aLastIteration);
extern "C" void link_cell_pairs_cuda_push_lists(int aFirstIteration);
extern "C" int two_body_forces_cuda_overlaps_with_host();

extern "C" int* spme_container_cuda_bspgen_leave_bspdxyz_data_online();
extern "C" void spme_container_cuda_bspgen_grab_bspdxyz_data
                (void** aBSPX, void ** aBSPY, void **aBSPZ,
                 void** aBSDX, void ** aBSDY, void **aBSDZ);

/**
 * @param aTruth If ==.true. (fortran), then the bsp{x,y,z} data will remain
 * on the device until, spme_container_cuda_bspgen_kill_bspdxyz_data() is
 * invoked.
 */
extern "C" void spme_container_cuda_bspgen_set_leave_bspdxyz_data_online(int *aTruth);

/**
 * Deallocates the bsp{x,y,z} data that have been left online.
 */
extern "C" void spme_container_cuda_bspgen_kill_bspdxyz_data();

extern "C" int dl_poly_cuda_offload_tbforces();

/**
 * @param aTruth If ==.true. (fortran), then the transfer of the bsp{x,y,z} data
 * back to the host can be reduced subject to the needs of the host-part of the
 * ewald spme forces.
 * @param aNATMS Needed as the iteration space of spme_forces is 1:natms (rather
 * than 1:nlast for bspgen and ccarray).
 */
extern "C" void spme_container_cuda_bspgen_set_is_under_ewald_spme_forces
                (int *aTruth, int *aNATMS);


/* @return The percentage of iterations that the CUDA code will handle.
 */
extern "C" double spme_forces_cuda_percentage_of_iterations_offloaded_to_the_device();

/* @return The percentage of iterations that the CUDA code will handle.
 */
extern "C" double ewald_spme_forces_cuda_ccarray_percentage_of_iterations_offloaded_to_the_device();

/* Existing DL_POLY_4 FORTRAN Subroutines: These are wrapped so they can be portably called
   from CUDA (using iso_c_binding) without modifying the original f90 sources
 */
extern "C" void wrapper_f_invert(real*, real*, real*);
extern "C" void wrapper_f_error(int*);
extern "C" void wrapper_f_erfcgen(real*,real*,int*,real*,real*);
extern "C" void wrapper_f_two_body_forces_cuda_helper
                (int *aIATMB, int *aITERS, int *aKEYFCE, real *aRHO, int *aSAFE,
                 int *aIMCON, real *aRCUT, real *aRVDW, real *aALPHA, real *aEPSQ,
                 real *aSTRESS, real *ENGMET, real *aVIRMET, real *aENGVDW,
                 real *aVIRVDW, real *aENGCPE_RL, real *aVIRCPE_RL);
extern "C" void wrapper_f_link_cell_pairs_helper
                (int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,
                 int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,real*);
extern "C" void wrapper_f_link_cell_pairs_remove_exclusions_helper(int*,int*);
extern "C" void wrapper_f_spme_forces_helper
                (int*,int*,real*,real*,real*,int*,int*,int*,int*,int*,int*,
                 real*,real*,real*,real*,real*,real*,real*,
                 int*,int*,int*);
extern "C" void wrapper_f_ewald_spme_forces_ccarray_helper
                (int*,int*,int*,int*,int*,int*,int*,int*,int*,
                 int*,int*,int*,int*,int*,int*,
                 real*,real*,real*,real*);
extern "C" void wrapper_f_ewald_spme_forces_ccarray_final_reduction
                (int*,int*,int*,real*,real*);
extern "C" int wrapper_f_local_index(int*,int*,int*,int*);

/* Copies list(*,1)..list(*,aLastIteration) from device to host memory. Depending on
 * what the state of data is after link_cell_pairs_cuda_invoke_ completes, there are
 * three possible cases
 * (1) if the requested portion, or a larger one, has been computed locally, then the
 *     call returns immediately, otherwise
 * (2) copy what is missing.
 */
extern "C" void link_cell_pairs_sparseListTransfer_any(int aIATM_Begin, int aN);

/* Timing-specific; disabling will keep the symbols available for the
 * FORTRAN sources.
 */
#define DISABLE_TIMING 0
#if (DISABLE_TIMING)
#if (!INCLUDED_BY_TIME_C)
#define start_timing_two_body_forces_cuda_k1() do {} while(0)
#define stop_timing_two_body_forces_cuda_k1() do {} while(0)
#define start_timing_two_body_forces_cuda_k2() do {} while(0)
#define stop_timing_two_body_forces_cuda_k2() do {} while(0)
#define start_timing_two_body_forces_cuda_k2_reduce() do {} while (0)
#define stop_timing_two_body_forces_cuda_k2_reduce() do {} while (0)
#define start_timing_two_body_forces_cuda_k3() do {} while(0)
#define stop_timing_two_body_forces_cuda_k3() do {} while(0)
#define start_timing_two_body_forces_cuda_finalise() do {} while(0)
#define stop_timing_two_body_forces_cuda_finalise() do {} while(0)
#define start_timing_two_body_forces_cuda_read() do {} while(0)
#define stop_timing_two_body_forces_cuda_read() do {} while(0)
#define start_timing_two_body_forces_cuda_write() do {} while(0)
#define stop_timing_two_body_forces_cuda_write() do {} while(0)
#define start_timing_link_cell_pairs_cuda_k0() do {} while(0)
#define start_timing_link_cell_pairs_cuda_k1() do {} while(0)
#define start_timing_link_cell_pairs_cuda_k2() do {} while(0)
#define start_timing_link_cell_pairs_cuda_finalise() do {} while(0)
#define start_timing_link_cell_pairs_cuda_read() do {} while(0)
#define start_timing_link_cell_pairs_cuda_write() do {} while(0)
#define start_timing_link_cell_pairs_sparse_list_transfer_reorder() do {} while (0)
#define stop_timing_link_cell_pairs_sparse_list_transfer_reorder() do {} while (0)
#define start_timing_link_cell_pairs_sparse_list_transfer_gxscan() do {} while (0)
#define stop_timing_link_cell_pairs_sparse_list_transfer_gxscan() do {} while (0)
#define start_timing_link_cell_pairs_sparse_list_transfer_pack() do {} while (0)
#define stop_timing_link_cell_pairs_sparse_list_transfer_pack() do {} while (0)
#define start_timing_link_cell_pairs_sparse_list_transfer_read() do {} while (0)
#define stop_timing_link_cell_pairs_sparse_list_transfer_read() do {} while (0)
#define start_timing_link_cell_pairs_sparse_list_transfer_unpack() do {} while (0)
#define stop_timing_link_cell_pairs_sparse_list_transfer_unpack() do {} while (0)
#define start_timing_link_cell_pairs_sparse_list_transfer() do {} while (0)
#define stop_timing_link_cell_pairs_sparse_list_transfer() do {} while (0)
#define stop_timing_link_cell_pairs_cuda_k0() do {} while(0)
#define stop_timing_link_cell_pairs_cuda_k1() do {} while(0)
#define stop_timing_link_cell_pairs_cuda_k2() do {} while(0)
#define stop_timing_link_cell_pairs_cuda_finalise() do {} while(0)
#define stop_timing_link_cell_pairs_cuda_read() do {} while(0)
#define stop_timing_link_cell_pairs_cuda_write() do {} while(0)
#define start_timing_spme_forces_k1() do {} while (0)
#define stop_timing_spme_forces_k1() do {} while (0)
#define start_timing_spme_forces_read() do {} while (0)
#define stop_timing_spme_forces_read() do {} while (0)
#define start_timing_spme_forces_write() do {} while (0)
#define stop_timing_spme_forces_write() do {} while (0)
#define start_timing_spme_forces_finalise() do {} while (0)
#define stop_timing_spme_forces_finalise() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_cccharge_k1() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_cccharge_read() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_cccharge_write() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_cccharge_k1() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_cccharge_read() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_cccharge_write() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_cccharge_finalise() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_cccharge_finalise() do {} while (0)
#define start_timing_constraints_shake_cuda_initialise() do {} while (0)
#define stop_timing_constraints_shake_cuda_initialise() do {} while (0)
#define start_timing_constraints_shake_cuda_write() do {} while (0)
#define stop_timing_constraints_shake_cuda_write() do {} while (0)
#define start_timing_constraints_shake_cuda_k1_th() do {} while (0)
#define stop_timing_constraints_shake_cuda_k1_th() do {} while (0)
#define start_timing_constraints_shake_cuda_read() do {} while (0)
#define stop_timing_constraints_shake_cuda_read() do {} while (0)
#define start_timing_constraints_shake_cuda_k1_bh() do {} while (0)
#define stop_timing_constraints_shake_cuda_k1_bh() do {} while (0)
#define start_timing_constraints_shake_cuda_install_red_struct() do {} while (0)
#define stop_timing_constraints_shake_cuda_install_red_struct() do {} while (0)
#define start_timing_constraints_shake_cuda_gather_dv_scatter_hs() do {} while (0)
#define stop_timing_constraints_shake_cuda_gather_dv_scatter_hs() do {} while (0)
#define start_timing_constraints_shake_cuda_gather_hs_scatter_dv() do {} while (0)
#define stop_timing_constraints_shake_cuda_gather_hs_scatter_dv() do {} while (0)
#define start_timing_constraints_shake_cuda_invoke_correct_positions() do {} while (0)
#define stop_timing_constraints_shake_cuda_invoke_correct_positions() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_ccarray_k1() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_ccarray_k1() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_ccarray_k2() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_ccarray_k2() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_ccarray_read() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_ccarray_read() do {} while (0)
#define start_timing_ewald_spme_forces_cuda_ccarray_write() do {} while (0)
#define stop_timing_ewald_spme_forces_cuda_ccarray_write() do {} while (0)
#define start_timing_link_cell_pairs_cuda_sort_atoms() do {} while (0)
#define stop_timing_link_cell_pairs_cuda_sort_atoms() do {} while (0)
#define start_timing_bspgen_cuda_k1() do {} while (0)
#define stop_timing_bspgen_cuda_k1() do {} while (0)
#define start_timing_bspgen_cuda_write() do {} while (0)
#define stop_timing_bspgen_cuda_write() do {} while (0)
#define start_timing_metal_ld_compute_cuda_write() do {} while (0)
#define stop_timing_metal_ld_compute_cuda_write() do {} while (0)
#define start_timing_metal_ld_compute_cuda_read() do {} while (0)
#define stop_timing_metal_ld_compute_cuda_read() do {} while (0)
#define start_timing_metal_ld_compute_cuda_k0() do {} while (0)
#define stop_timing_metal_ld_compute_cuda_k0() do {} while (0)
#define start_timing_metal_ld_compute_cuda_k1() do {} while (0)
#define stop_timing_metal_ld_compute_cuda_k1() do {} while (0)
#define start_timing_metal_ld_compute_cuda_k2() do {} while (0)
#define stop_timing_metal_ld_compute_cuda_k2() do {} while (0)
#define start_timing_metal_ld_compute_cuda_k3() do {} while (0)
#define stop_timing_metal_ld_compute_cuda_k3() do {} while (0)
#endif
#else
extern "C" void start_timing_spme_forces_k1();
extern "C" void stop_timing_spme_forces_k1();
extern "C" void start_timing_spme_forces_read();
extern "C" void stop_timing_spme_forces_read();
extern "C" void start_timing_spme_forces_write();
extern "C" void stop_timing_spme_forces_write();
extern "C" void start_timing_spme_forces_finalise();
extern "C" void stop_timing_spme_forces_finalise();
extern "C" void start_timing_two_body_forces_cuda_k0();
extern "C" void stop_timing_two_body_forces_cuda_k0();
extern "C" void start_timing_two_body_forces_cuda_k1();
extern "C" void stop_timing_two_body_forces_cuda_k1();
extern "C" void start_timing_two_body_forces_cuda_k2();
extern "C" void stop_timing_two_body_forces_cuda_k2();
extern "C" void start_timing_two_body_forces_cuda_k2_reduce();
extern "C" void stop_timing_two_body_forces_cuda_k2_reduce();
extern "C" void start_timing_two_body_forces_cuda_k3();
extern "C" void stop_timing_two_body_forces_cuda_k3();
extern "C" void start_timing_two_body_forces_cuda_finalise();
extern "C" void stop_timing_two_body_forces_cuda_finalise();
extern "C" void start_timing_two_body_forces_cuda_read();
extern "C" void stop_timing_two_body_forces_cuda_read();
extern "C" void start_timing_two_body_forces_cuda_write();
extern "C" void stop_timing_two_body_forces_cuda_write();
extern "C" void start_timing_link_cell_pairs_cuda_k0();
extern "C" void start_timing_link_cell_pairs_cuda_k1();
extern "C" void start_timing_link_cell_pairs_cuda_k2();
extern "C" void start_timing_link_cell_pairs_cuda_finalise();
extern "C" void start_timing_link_cell_pairs_cuda_read();
extern "C" void start_timing_link_cell_pairs_cuda_write();
extern "C" void start_timing_link_cell_pairs_sparse_list_transfer_reorder();
extern "C" void stop_timing_link_cell_pairs_sparse_list_transfer_reorder();
extern "C" void start_timing_link_cell_pairs_sparse_list_transfer_gxscan();
extern "C" void stop_timing_link_cell_pairs_sparse_list_transfer_gxscan();
extern "C" void start_timing_link_cell_pairs_sparse_list_transfer_pack();
extern "C" void stop_timing_link_cell_pairs_sparse_list_transfer_pack();
extern "C" void start_timing_link_cell_pairs_sparse_list_transfer_read();
extern "C" void stop_timing_link_cell_pairs_sparse_list_transfer_read();
extern "C" void start_timing_link_cell_pairs_sparse_list_transfer_unpack();
extern "C" void stop_timing_link_cell_pairs_sparse_list_transfer_unpack();
extern "C" void start_timing_link_cell_pairs_sparse_list_transfer();
extern "C" void stop_timing_link_cell_pairs_sparse_list_transfer();
extern "C" void stop_timing_link_cell_pairs_cuda_k0();
extern "C" void stop_timing_link_cell_pairs_cuda_k1();
extern "C" void stop_timing_link_cell_pairs_cuda_k2();
extern "C" void stop_timing_link_cell_pairs_cuda_finalise();
extern "C" void stop_timing_link_cell_pairs_cuda_read();
extern "C" void stop_timing_link_cell_pairs_cuda_write();
extern "C" void start_timing_ewald_spme_forces_cuda_cccharge_k1();
extern "C" void start_timing_ewald_spme_forces_cuda_cccharge_read();
extern "C" void start_timing_ewald_spme_forces_cuda_cccharge_write();
extern "C" void stop_timing_ewald_spme_forces_cuda_cccharge_k1();
extern "C" void stop_timing_ewald_spme_forces_cuda_cccharge_read();
extern "C" void stop_timing_ewald_spme_forces_cuda_cccharge_write();
extern "C" void start_timing_ewald_spme_forces_cuda_cccharge_finalise();
extern "C" void stop_timing_ewald_spme_forces_cuda_cccharge_finalise();
extern "C" void start_timing_constraints_shake_cuda_initialise();
extern "C" void stop_timing_constraints_shake_cuda_initialise();
extern "C" void start_timing_constraints_shake_cuda_write();
extern "C" void stop_timing_constraints_shake_cuda_write();
extern "C" void start_timing_constraints_shake_cuda_k1_th();
extern "C" void stop_timing_constraints_shake_cuda_k1_th();
extern "C" void start_timing_constraints_shake_cuda_read();
extern "C" void stop_timing_constraints_shake_cuda_read();
extern "C" void start_timing_constraints_shake_cuda_k1_bh();
extern "C" void stop_timing_constraints_shake_cuda_k1_bh();
extern "C" void start_timing_constraints_shake_cuda_install_red_struct();
extern "C" void stop_timing_constraints_shake_cuda_install_red_struct();
extern "C" void start_timing_constraints_shake_cuda_gather_dv_scatter_hs();
extern "C" void stop_timing_constraints_shake_cuda_gather_dv_scatter_hs();
extern "C" void start_timing_constraints_shake_cuda_gather_hs_scatter_dv();
extern "C" void stop_timing_constraints_shake_cuda_gather_hs_scatter_dv();
extern "C" void start_timing_constraints_shake_cuda_invoke_correct_positions();
extern "C" void stop_timing_constraints_shake_cuda_invoke_correct_positions();
extern "C" void start_timing_ewald_spme_forces_cuda_ccarray_k1();
extern "C" void stop_timing_ewald_spme_forces_cuda_ccarray_k1();
extern "C" void start_timing_ewald_spme_forces_cuda_ccarray_k2();
extern "C" void stop_timing_ewald_spme_forces_cuda_ccarray_k2();
extern "C" void start_timing_ewald_spme_forces_cuda_ccarray_read();
extern "C" void stop_timing_ewald_spme_forces_cuda_ccarray_read();
extern "C" void start_timing_ewald_spme_forces_cuda_ccarray_write();
extern "C" void stop_timing_ewald_spme_forces_cuda_ccarray_write();
extern "C" void start_timing_link_cell_pairs_cuda_sort_atoms();
extern "C" void stop_timing_link_cell_pairs_cuda_sort_atoms();
extern "C" void start_timing_bspgen_cuda_k1();
extern "C" void stop_timing_bspgen_cuda_k1();
extern "C" void start_timing_bspgen_cuda_write();
extern "C" void stop_timing_bspgen_cuda_write();
extern "C" void start_timing_metal_ld_compute_cuda_write();
extern "C" void stop_timing_metal_ld_compute_cuda_write();
extern "C" void start_timing_metal_ld_compute_cuda_read();
extern "C" void stop_timing_metal_ld_compute_cuda_read();
extern "C" void start_timing_metal_ld_compute_cuda_k0();
extern "C" void stop_timing_metal_ld_compute_cuda_k0();
extern "C" void start_timing_metal_ld_compute_cuda_k1();
extern "C" void stop_timing_metal_ld_compute_cuda_k1();
extern "C" void start_timing_metal_ld_compute_cuda_k2();
extern "C" void stop_timing_metal_ld_compute_cuda_k2();
extern "C" void start_timing_metal_ld_compute_cuda_k3();
extern "C" void stop_timing_metal_ld_compute_cuda_k3();

#endif


#endif //!_DL_POLY_CUDA_H
