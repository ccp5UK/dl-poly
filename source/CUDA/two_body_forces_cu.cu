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

#include <assert.h>
#include <stdlib.h>
#include <cuda.h>
#include <cutil.h>


/* TODOs:
 * 01. Check those "fabs" call what PTX instructions they translate to --
 *     don't want any conversions in the way.
 * 02. metal_forces: "safe" is not propagated to host.
 * 03. there are redundant "images" calculations
 * 04. we're copying more to the device than needed -- must take the
 *     primary atom offset into account.
 */

#include "dl_poly_cu.h"
#define CFG_ACCURATE                                          1
#include "dl_poly_common_cu.cu"

/* If set to 1, then the code will constantly be monitoring itself
 * in order to decide on the most appropriate host/device and
 * device/device ratios.
 *  Env var: dlpolycuda_twobody_dynamic_ratios
 */
#define CFG_DYNAMIC_RATIOS                                    (CFG_OVERLAP_WITH_HOST)

/* This is the maximum number of iterations that the device can
 * be assigned. When overlapping is enabled, the host executes
 * this value minus what the device has been assigned.
 *   Env var: dlpolycuda_twobody_max_device_iterations
 */
#define CFG_MAX_DEVICE_ITERATIONS                             30000

/* The number of threadblocks assigned to K1.
 *   Env var: dlpolycuda_twobody_k1_grid_dims_y
 */

#define CFG_K1_GRID_DIMS_Y                                    900
#define CFG_K1_DEBUG                                          0

#if (!CFG_K1_DEBUG)
#  if (CFG_COMPUTE_MAJOR != 2)
#    define CFG_K1_THREADS_PER_BLOCK                            64
#  else
//   this value works well with fermis
#    define CFG_K1_THREADS_PER_BLOCK                            32
#  endif
#else
#  define CFG_K1_THREADS_PER_BLOCK                              1
#endif

#define CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE  (CFG_OVERLAP_WITH_HOST ? 0.7 :1.0)

/* The number of buffers dedicated for the parallelisation of K2. Due
 * to the compiler performing heavy unrolling, this option should be
 * tuned with care. It should be possible to set it quite high (>10)
 * if the k2 and k2_reduce threadblock sizes are shrinked.
 *   Env var: dlpolycuda_twobody_k2_unroll_slices
 * default value 4
 */
#define CFG_K2_UNROLL_SLICES                                  6

#define CFG_K2_RED_GRID_DIMS_Y                                300

#define CFG_K2_RED_GRID_DIMS                                  dim3(3,CFG_K2_RED_GRID_DIMS_Y,1)

/* The number of threads per block for the k2 kernel.
 *   Env var: dlpolycuda_twobody_k2_block_dims_x
 */
#define CFG_K2_BLOCK_DIMS_X                                   256
#define CFG_K2_BLOCK_DIMS                                     dim3(CFG_K2_BLOCK_DIMS_X,1,1)
#define CFG_K2_RED_BLOCK_DIMS_X                               512
#define CFG_K2_RED_BLOCK_DIMS                                 dim3(CFG_K2_RED_BLOCK_DIMS_X,1,1)
#define CFG_PERCENTAGE_OF_ITERATIONS_TO_BE_OVERLAPPED_WITH_K1 0.6
#define CFG_LSTLTPVDW_FETCH_FROM_CONSTANT_MEMORY              1
#define CFG_MXVDW_MAX_VALUE                                   195
#define CFG_K2_DEBUG                                          0
#define CFG_MAX_LIMIT                                         512


extern "C" double two_body_forces_cuda_piotodevice() {
  return (CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE);
}


extern "C" int two_body_forces_cuda_overlaps_with_host() {
  return (CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE!=1.0);
}

template<typename T_> struct constant_data {

  struct metal_forces {
    int mKEYPOT;
    int *mLTPMET, *mLSTMET;
    T_  *mVMET,   *mDMET;
    T_  *mRHO;
  }    mMET;

  // vdw_forces-specific:
  struct vdw_forces {
#if (CFG_LSTLTPVDW_FETCH_FROM_CONSTANT_MEMORY)
    int mLTPVDW[CFG_MXVDW_MAX_VALUE];
    int mLSTVDW[CFG_MXVDW_MAX_VALUE];
#else
    int *mLSTVDW;
    int *mLTPVDW;
#endif
    T_  *mVVDW;
    T_  *mGVDW;
    T_   mRDR;
    T_   mRCSQ;
    T_   mRVDW;
    bool mLS_VDW;
  }    mVDW;

  struct ewald_real_forces {
    T_  *mCHGE;
    T_   mRCSQ;
    T_   mDREWD;
    T_   mRDREWD;
    T_   mR4PIE0;
    T_   mEPSQ;
    T_   mRCUT;
    T_   mALPHA;
    T_  *mERC;
    T_  *mFER;
  }    mEWR;

  int  mMXMET;
  int  mNATMS;
  int  mMXLIST;
  int  mMXATMS;
  int  mMXGRID;
  T_   mCELL[9];
  T_   mRCELL[9]; // reciprocals of mCELL
  T_   mICELL[9]; // the following results from invert_ against the mCELL field.
  T_  *mXXX;
  T_  *mYYY;
  T_  *mZZZ;

  int *mLIST; // reminder: list has mxlist+1 rows
  int *mLTYPE;
  int *mLTG;
};

__device__ __constant__ constant_data<real> CONSTANT_DATA;
static constant_data<real>                  sCD;

template<typename T_> struct host_data {
  int mK2_Unroll_Slices;
  int mK2_Block_Dims_X;
  int mK1_Grid_Dims_Y;
  int mUnroll;
  int mDynamicRatios;

  int mIsListAlreadyOnline, mFreeList;
  int  mNTPMET, mNTPVDW, mKEYFCE, mIMCON;
  int * __restrict__ mLIST;
  T_  *mDEV_Scratchpad;
  T_  *mHST_PinnedScratchpad;
  int  mScratchpadSize;

  double mPercentageOfIterationsOffloadedToTheDevice;
  double mPercentageOfHostIterationsToBeOverlappedWithK1;
  int * mLTG;
  int * mLTYPE;
  T_ *mXXX, *mYYY, *mZZZ, *mCELL;

  struct metal_forces {
    T_ *mRHO;
  } mMET;

  struct vdw_forces {
    T_ mRVDW;
  } mVDW;

  T_ *mDEV_Spad0, *mHST_Spad0;
};
static host_data<real> sHD;

template<typename T_> __global__ void memsetDevicePointer(T_ *aPointer, int aN, T_ aValue) {
  int lOUT_Slot = threadIdx.x;
  while (lOUT_Slot < aN) {
    aPointer[lOUT_Slot] = aValue;
    lOUT_Slot += blockDim.x;
  }
}

texture<real_tex, 1, cudaReadModeElementType> TEXTURE_DATA_VDW_GVDW;
texture<real_tex, 1, cudaReadModeElementType> TEXTURE_DATA_VDW_VVDW;
texture<real_tex, 1, cudaReadModeElementType> TEXTURE_DATA_EWR_ERC;
texture<real_tex, 1, cudaReadModeElementType> TEXTURE_DATA_EWR_FER;
texture<real_tex, 1, cudaReadModeElementType> TEXTURE_DATA_XXX, TEXTURE_DATA_YYY, TEXTURE_DATA_ZZZ;

/* @param aIsListOnline If 1, then link_cell_pairs_cuda has left a copy of the list
 *   in device memory which we can use.
 */
extern "C" void two_body_forces_cuda_initialise(
      int *aIsListOnline,
      int *aNATMS, int *aMXLIST, int *aMXATMS, int *aMXATDM, int *aMXMET,
      int *aNTPMET, int *aNTPVDW, int *aKEYFCE, int *aIMCON,
      real *aCELL, real *aXXX, real *aYYY, real *aZZZ,
      int *aLIST, int *aLTYPE,
      int *aLTPMET, int *aLSTMET, real *aVMET, real *aDMET, real *aRHO,
      int *aMXGRID, int *aMXVDW, int *aLTPVDW, int *aLSTVDW, bool aLS_VDW,
      real *aVVDW, real *aGVDW, real *aRVDW, int *aLTG,
      real *aCHGE, real *aRCUT, real *aALPHA, real *aEPSQ) {


  sHD.mK1_Grid_Dims_Y   = dl_poly_cuda_getenv("dlpolycuda_twobody_k1_grid_dims_y",
                                              CFG_K1_GRID_DIMS_Y);
  sHD.mK2_Unroll_Slices = dl_poly_cuda_getenv("dlpolycuda_twobody_k2_unroll_slices",
                                              CFG_K2_UNROLL_SLICES);
  sHD.mK2_Block_Dims_X  = dl_poly_cuda_getenv("dlpolycuda_twobody_k2_block_dims_x",
                                              CFG_K2_BLOCK_DIMS_X);
  sHD.mUnroll           = dl_poly_cuda_getenv("dlpolycuda_twobody_max_device_iterations",
                                              CFG_MAX_DEVICE_ITERATIONS);
  sHD.mDynamicRatios    = dl_poly_cuda_getenv("dlpolycuda_twobody_dynamic_ratios",
                                              CFG_DYNAMIC_RATIOS);
  static int sSetRatios = 1;

  if (sSetRatios) {
    sHD.mPercentageOfIterationsOffloadedToTheDevice     =
      CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE;
    sHD.mPercentageOfHostIterationsToBeOverlappedWithK1 =
      CFG_PERCENTAGE_OF_ITERATIONS_TO_BE_OVERLAPPED_WITH_K1;

    if (sHD.mDynamicRatios) {
      sSetRatios = 0;
    }
  }

  if (sHD.mDynamicRatios) {
    if (sHD.mPercentageOfIterationsOffloadedToTheDevice==0.0 ||
        sHD.mPercentageOfIterationsOffloadedToTheDevice==1.0) {
      sHD.mDynamicRatios = 0;
      printf("%s::%s: warning: disabled host/device load-balancing\n", __FILE__, __FUNCTION__);
    }
  }

  sHD.mIsListAlreadyOnline = *aIsListOnline;
  sHD.mLIST                =  aLIST;
  sCD.mNATMS               = *aNATMS;
  sCD.mMXLIST              = *aMXLIST;
  sHD.mScratchpadSize      = (sHD.mUnroll*(512*3 + 3 + 12) + 12)*sizeof(real);

  CUDA_SAFE_CALL(cudaMalloc(&sHD.mDEV_Scratchpad, sHD.mScratchpadSize));
  CUDA_SAFE_CALL(cudaHostAlloc(&sHD.mHST_PinnedScratchpad,12*sHD.mUnroll*sizeof(real),0));

  sHD.mNTPMET = *aNTPMET;
  sHD.mNTPVDW = *aNTPVDW;
  sHD.mKEYFCE = *aKEYFCE;
  sHD.mIMCON  = *aIMCON;
  sHD.mLTG    =  aLTG;
  sHD.mLTYPE  =  aLTYPE;
  sHD.mXXX = aXXX;
  sHD.mYYY = aYYY;
  sHD.mZZZ = aZZZ;
  sHD.mCELL = aCELL;
  sHD.mVDW.mRVDW = *aRVDW;

  // this is what has been implemented:
  if (!((sHD.mKEYFCE==0 || sHD.mKEYFCE==2) &&
        (sHD.mIMCON==2  || sHD.mIMCON==3))) {
    printf("%s::%s: stub: can only handle keyfce==0,2,imcon==2,3; found %d,%d\n",
           __FILE__, __FUNCTION__, sHD.mKEYFCE, sHD.mIMCON);
    exit(-1);
  }

  CUDA_SAFE_CALL(cudaMalloc(&sHD.mDEV_Spad0, sHD.mK2_Unroll_Slices*(*aMXATMS)*sizeof(real)*3));
  CUDA_SAFE_CALL(cudaHostAlloc(&sHD.mHST_Spad0, (*aMXATMS)*sizeof(real)*3, 0));

  memsetDevicePointer<real><<<1, 512>>>(sHD.mDEV_Spad0, sHD.mK2_Unroll_Slices*3*(*aMXATMS), (real)0);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  cudaError_t lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);

  sCD.mVDW.mLS_VDW =  aLS_VDW;
  sHD.mVDW.mRVDW   = *aRVDW;
  sCD.mMXATMS      = *aMXATMS;
  sCD.mMXGRID      = *aMXGRID;
  sCD.mMXMET       = *aMXMET;

  for (int lI=0 ; lI<9 ; lI++) {
    real lCELL        = aCELL[lI];
    sCD.mCELL[lI]     = lCELL;
    sCD.mRCELL[lI] = ((real) 1) / lCELL;
  }
  real lDummy;
  wrapper_f_invert(sCD.mCELL, sCD.mICELL, &lDummy);

  const real lDLRPOT = (*aRVDW) / (real)(sCD.mMXGRID - 4);
  sCD.mVDW.mRDR      = ((real)1) / lDLRPOT;
  sCD.mVDW.mRCSQ     = (*aRVDW)*(*aRVDW);

  sCD.mEWR.mRCUT     = *aRCUT;
  sCD.mEWR.mDREWD    = sCD.mEWR.mRCUT / (real)(sCD.mMXGRID - 4);
  sCD.mEWR.mALPHA    = *aALPHA;
  sCD.mEWR.mR4PIE0   = (real)138935.4835;
  sCD.mEWR.mEPSQ     = *aEPSQ;
  sCD.mEWR.mRCSQ     = sCD.mEWR.mRCUT * sCD.mEWR.mRCUT;
  sCD.mEWR.mRDREWD   = ((real) 1) / sCD.mEWR.mDREWD;

  real *lERC = NULL, *lFER = NULL;

  int lMXATDM = *aMXATDM;
  int lMXVDW  = *aMXVDW;
  int lMXGRID = *aMXGRID;

  CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mXXX,         sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mYYY,         sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mZZZ,         sCD.mMXATMS*sizeof(real)));

  BIND_TEXTURE_REAL_1D(TEXTURE_DATA_XXX, sCD.mXXX, sCD.mMXATMS);
  BIND_TEXTURE_REAL_1D(TEXTURE_DATA_YYY, sCD.mYYY, sCD.mMXATMS);
  BIND_TEXTURE_REAL_1D(TEXTURE_DATA_ZZZ, sCD.mZZZ, sCD.mMXATMS);

  if (sHD.mNTPMET > 0) {
    /* This static part is concerned with the newjob clause in metal_forces,
     * which configures the 'jeypot' variable.
     */
    static int sHasKeypotBeenComputed = 0;
    if (!sHasKeypotBeenComputed) {
      int lKEYPOT = 0;
      for (int lL=1 ; lL<=sHD.mNTPMET ; lL++) {
        lKEYPOT = aLTPMET[lL-1];
        if (lL>1) {
          if (lKEYPOT != aLTPMET[(lL-1)]) {
            int lError = 92;
            wrapper_f_error(&lError);
          }
        }
      }
      sCD.mMET.mKEYPOT = lKEYPOT;
      sHasKeypotBeenComputed = 1;

      if (sCD.mMET.mKEYPOT==0) {
        // TODO -- easy to implement as code is already there but *has not been tested*
        printf("%s::%s: stub: keypot==0.\n", __FILE__, __FUNCTION__);
        exit(-1);
      }
    }
    sHD.mMET.mRHO = aRHO;

    CUDA_SAFE_CALL(cudaMalloc(&sCD.mMET.mLSTMET, sCD.mMXMET*sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mMET.mLTPMET, sCD.mMXMET*sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mMET.mVMET,   sCD.mMXGRID*sCD.mMXMET*2*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mMET.mDMET,   sCD.mMXGRID*sCD.mMXMET*2*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mMET.mRHO,    sCD.mMXATMS*sizeof(real)));
  }

  /* Use the list if it is already online:
   */
  if (sHD.mIsListAlreadyOnline) {
    sCD.mLIST = (int*) link_cell_pairs_cuda_get_list();
    sHD.mFreeList = 0;
  } else {
    // note: +512 is headroom space
    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mLIST, (512+lMXATDM*(1+sCD.mMXLIST))*sizeof(int)));
    sHD.mFreeList = 1;
  }

  CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mLTYPE,  sCD.mMXATMS*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mLTG,    sCD.mMXATMS*sizeof(int)));

  if (sHD.mNTPVDW > 0) {
    /* vdw_forces.f90-specific initialisations.
     */
#if (CFG_LSTLTPVDW_FETCH_FROM_CONSTANT_MEMORY)
    /* If l{st,pt}vdw are fetched from the constant memory, make sure that
     * there is enough room available.
     */
    if (CFG_MXVDW_MAX_VALUE<lMXVDW) {
      printf("%s::%s: cannot handle mxvdw=%d; Try setting"
             " CFG_LSTLTPVDW_FETCH_FROM_CONSTANT_MEMORY to zero"
	     " or increase CFG_MXVDW_MAX_VALUE\n",
	     __FILE__, __FUNCTION__, lMXVDW);
      exit(-1);
    }
#else
    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mVDW.mLSTVDW, lMXVDW*sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mVDW.mLTPVDW, lMXVDW*sizeof(int)));
#endif

    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mVDW.mVVDW,   sCD.mMXGRID*lMXVDW*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mVDW.mGVDW,   sCD.mMXGRID*lMXVDW*sizeof(real)));

    BIND_TEXTURE_REAL_1D(TEXTURE_DATA_VDW_VVDW, sCD.mVDW.mVVDW, sCD.mMXGRID*lMXVDW);

    BIND_TEXTURE_REAL_1D(TEXTURE_DATA_VDW_GVDW, sCD.mVDW.mGVDW, sCD.mMXGRID*lMXVDW);
  }

  if (sHD.mKEYFCE==2) {
    /* ewald_real_forces.f90-specific initialisations:
     */
    lERC = new real[sCD.mMXGRID];
    lFER = new real[sCD.mMXGRID];
    wrapper_f_erfcgen(aRCUT, aALPHA, aMXGRID, lERC, lFER);

    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mEWR.mCHGE,   sCD.mMXATMS*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mEWR.mERC,    sCD.mMXGRID*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mEWR.mFER,    sCD.mMXGRID*sizeof(real)));

    BIND_TEXTURE_REAL_1D(TEXTURE_DATA_EWR_ERC, sCD.mEWR.mERC, sCD.mMXGRID);
    BIND_TEXTURE_REAL_1D(TEXTURE_DATA_EWR_FER, sCD.mEWR.mFER, sCD.mMXGRID);
  }
  /* Copy the data over to the device.
   */
  start_timing_two_body_forces_cuda_write();
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mXXX, aXXX, sCD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mYYY, aYYY, sCD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mZZZ, aZZZ, sCD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));

  if (!sHD.mIsListAlreadyOnline) {
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mLIST, aLIST, lMXATDM*(1+sCD.mMXLIST)*sizeof(int), cudaMemcpyHostToDevice));
    sHD.mIsListAlreadyOnline = 1;
  }

  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLTYPE, aLTYPE, sCD.mMXATMS*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLTG,   aLTG,   sCD.mMXATMS*sizeof(int), cudaMemcpyHostToDevice));

  if (sHD.mNTPMET>0) {
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mMET.mLSTMET, aLSTMET, sCD.mMXMET*sizeof(int),
                              cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mMET.mLTPMET, aLTPMET, sCD.mMXMET*sizeof(int),
                              cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mMET.mVMET, aVMET, sCD.mMXGRID*sCD.mMXMET*2*sizeof(real),
                              cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mMET.mDMET, aDMET, sCD.mMXGRID*sCD.mMXMET*2*sizeof(real),
                              cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mMET.mRHO, aRHO, sCD.mMXATMS*sizeof(real),
                              cudaMemcpyHostToDevice));
  }

  if (sHD.mNTPVDW > 0) {
#if (CFG_LSTLTPVDW_FETCH_FROM_CONSTANT_MEMORY)
    for (int lI=0 ; lI<lMXVDW ; lI++) {
      sCD.mVDW.mLSTVDW[lI] = aLSTVDW[lI];
      sCD.mVDW.mLTPVDW[lI] = aLTPVDW[lI];
    }
#else
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mVDW.mLSTVDW, aLSTVDW, lMXVDW*sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mVDW.mLTPVDW, aLTPVDW, lMXVDW*sizeof(int), cudaMemcpyHostToDevice));
#endif
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mVDW.mVVDW,   aVVDW,    lMXGRID*lMXVDW*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mVDW.mGVDW,   aGVDW,    lMXGRID*lMXVDW*sizeof(real), cudaMemcpyHostToDevice));
  }

  if (sHD.mKEYFCE==2) {
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mEWR.mCHGE,   aCHGE,    sCD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mEWR.mERC,    lERC,     sCD.mMXGRID*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mEWR.mFER,    lFER,     sCD.mMXGRID*sizeof(real), cudaMemcpyHostToDevice));
    delete[] lERC;
    delete[] lFER;
  }

  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA, (void*)&sCD, sizeof(constant_data<real>)));
  stop_timing_two_body_forces_cuda_write();
}

extern "C" void two_body_forces_cuda_finalise() {
  CUDA_SAFE_CALL(cudaFreeHost(sHD.mHST_PinnedScratchpad));
  CUDA_SAFE_CALL(cudaFree(sHD.mDEV_Scratchpad));
  CUDA_SAFE_CALL(cudaFreeHost(sHD.mHST_Spad0));
  CUDA_SAFE_CALL(cudaFree(sHD.mDEV_Spad0));

  CUDA_SAFE_CALL(cudaFree(sCD.mXXX));
  CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_XXX));
  CUDA_SAFE_CALL(cudaFree(sCD.mYYY));
  CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_YYY));
  CUDA_SAFE_CALL(cudaFree(sCD.mZZZ));
  CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_ZZZ));

  if (sHD.mFreeList) {
    CUDA_SAFE_CALL(cudaFree(sCD.mLIST));
  }

  CUDA_SAFE_CALL(cudaFree(sCD.mLTYPE));
  CUDA_SAFE_CALL(cudaFree(sCD.mLTG));

  if (sHD.mNTPMET>0) {
    CUDA_SAFE_CALL(cudaFree(sCD.mMET.mLSTMET));
    CUDA_SAFE_CALL(cudaFree(sCD.mMET.mLTPMET));
    CUDA_SAFE_CALL(cudaFree(sCD.mMET.mVMET));
    CUDA_SAFE_CALL(cudaFree(sCD.mMET.mDMET));
    CUDA_SAFE_CALL(cudaFree(sCD.mMET.mRHO));
  }

  if (sHD.mNTPVDW > 0) {
#if (!CFG_LSTLTPVDW_FETCH_FROM_CONSTANT_MEMORY)
    CUDA_SAFE_CALL(cudaFree(sCD.mVDW.mLSTVDW));
    CUDA_SAFE_CALL(cudaFree(sCD.mVDW.mLTPVDW));
#endif
    CUDA_SAFE_CALL(cudaFree(sCD.mVDW.mVVDW));
    CUDA_SAFE_CALL(cudaFree(sCD.mVDW.mGVDW));
    CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_VDW_VVDW));
    CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_VDW_GVDW));
  }

  if (sHD.mKEYFCE==2) {
    CUDA_SAFE_CALL(cudaFree(sCD.mEWR.mCHGE));
    CUDA_SAFE_CALL(cudaFree(sCD.mEWR.mERC));
    CUDA_SAFE_CALL(cudaFree(sCD.mEWR.mFER));
    CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_EWR_ERC));
    CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_EWR_FER));
  }
}


template<typename T_> __device__ void obtain_iatm_specific(
  int aIATM, T_& aXXX_I, T_& aYYY_I, T_& aZZZ_I, int& aLIMIT, int& aIDI) {

  aLIMIT = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1), 0, aIATM);

  if (aLIMIT > 0) {
    aXXX_I = CONSTANT_DATA.mXXX[aIATM-1];
    aYYY_I = CONSTANT_DATA.mYYY[aIATM-1];
    aZZZ_I = CONSTANT_DATA.mZZZ[aIATM-1];
    aIDI   = CONSTANT_DATA.mLTG[aIATM-1];
  }
}


template<typename T_, int IMCON_> __device__ void obtain_jatm_specific(
  int aJATM, T_ aXXX_I, T_ aYYY_I, T_ aZZZ_I, T_& aXDF, T_& aYDF, T_& aZDF, T_& aRSQDF, int& aIDJ) {
  T_ lXXX_J = fetch_real_1d(TEXTURE_DATA_XXX, aJATM-1);//CONSTANT_DATA.mXXX[aJATM-1];
  T_ lYYY_J = fetch_real_1d(TEXTURE_DATA_YYY, aJATM-1);//CONSTANT_DATA.mYYY[aJATM-1];
  T_ lZZZ_J = fetch_real_1d(TEXTURE_DATA_ZZZ, aJATM-1);//CONSTANT_DATA.mZZZ[aJATM-1];

  aIDJ      = CONSTANT_DATA.mLTG[aJATM-1];
  aXDF      = aXXX_I - lXXX_J;
  aYDF      = aYYY_I - lYYY_J;
  aZDF      = aZZZ_I - lZZZ_J;

  IMAGES(CONSTANT_DATA, T_, IMCON_, aXDF, aYDF, aZDF);

  aRSQDF = addp3(aXDF,aXDF, aYDF,aYDF, aZDF,aZDF);
}

template<typename T_> __device__ T_* vmet(int aX, int aY, int aZ) {
  return(dev_f3d_address<T_,0>(CONSTANT_DATA.mMET.mVMET, 1,1,1,
			       CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX, aY, aZ));
}

template<typename T_> __device__ T_* dmet(int aX, int aY, int aZ) {
  return(dev_f3d_address<T_,0>(CONSTANT_DATA.mMET.mDMET, 1,1,1,
			       CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX, aY, aZ));
}

template<typename T_> __device__ T_ metal_forces_calc_gamma(T_ aPPx, T_ aT1, T_ aT2) {
  T_ lSelect_T1T2 = (aPPx < (T_)0) ? aT1 : aT2;
  T_ lPlusMinus1  = (aPPx < (T_)0) ? ((T_) 1) : ((T_) -1);

  T_ lGamma = madd(lSelect_T1T2, ((T_)0.5), mul(add(aT2,-aT1),add(aPPx, lPlusMinus1)));
  return (lGamma);
}

template<typename T_>
__device__ T_ metal_forces_calc_gamma(T_ *aBase, int aX, int aY, int aZ, T_ aPPx) {
  T_ lGVK0 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX+3, aY, aZ);
  T_ lGVK1 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX+4, aY, aZ);
  T_ lGVK2 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX+5, aY, aZ);

  T_ lT1 = madd(lGVK1, aPPx, (lGVK1 - lGVK0));
  T_ lT2 = madd(lGVK1, aPPx, (lGVK2 - lGVK1));
  return(metal_forces_calc_gamma(aPPx, lT1, lT2));
}

template<typename T_>
__device__ T_ metal_forces_calc_gamma_vmet(int aX, int aY, int aZ, T_ aPPx) {
  return(metal_forces_calc_gamma(CONSTANT_DATA.mMET.mVMET, aX, aY, aZ, aPPx));
}

template<typename T_>
__device__ T_ metal_forces_calc_gamma_dmet(int aX, int aY, int aZ, T_ aPPx) {
  return(metal_forces_calc_gamma(CONSTANT_DATA.mMET.mDMET, aX, aY, aZ, aPPx));
}

template<typename T_> __device__ void metal_forces_KEYPOTNeqZero
#if (CFG_UNIFIED_ADDRESS_SPACE)
    (int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
     volatile T_& aFX_I, volatile T_& aFY_I, volatile T_& aFZ_I,
     T_& aFX_J, T_& aFY_J, T_& aFZ_J,
     T_& aSTRS1, T_& aSTRS2, T_& aSTRS3, volatile T_& aSTRS5, volatile T_& aSTRS6, volatile T_& aSTRS9, T_& aENGMET, T_& aVIRMET){
#else
    (int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
     T_& aFX_I, T_& aFY_I, T_& aFZ_I,
     T_& aFX_J, T_& aFY_J, T_& aFZ_J,
     T_& aSTRS1, T_& aSTRS2, T_& aSTRS3, T_& aSTRS5, T_& aSTRS6, T_& aSTRS9, T_& aENGMET, T_& aVIRMET){
#endif

  int lAI          = CONSTANT_DATA.mLTYPE[aIATM-1];
  int lAJ          = CONSTANT_DATA.mLTYPE[aJATM-1];
  int lKEY         = max(lAI,lAJ) * (max(lAI,lAJ)-1) / 2 + min(lAI,lAJ);
  int lK0          = CONSTANT_DATA.mMET.mLSTMET[lKEY - 1];
  T_  lVMET_1_K0_1 = *vmet<T_>(1,lK0,1);
  T_  lVMET_3_K0_1 = *vmet<T_>(3,lK0,1);
  int lLTPMET_K0   = CONSTANT_DATA.mMET.mLTPMET[lK0-1];

  if (lLTPMET_K0>=0 && fabs(lVMET_1_K0_1)>(T_)0 && aRSQ<=(lVMET_3_K0_1*lVMET_3_K0_1)) {
    T_ lVMET_2_K0_1 = *vmet<T_>(2,lK0,1);
    T_ lVMET_4_K0_1 = *vmet<T_>(4,lK0,1);
    T_ lRDR         = ((T_) 1) / lVMET_4_K0_1;
    T_ lRRR         = sqrt(aRSQ) - lVMET_2_K0_1;
    int lL          = nint(lRRR*lRDR);
    T_ lPPP         = mul(lRRR,lRDR) - (T_)lL;
    T_ lGAMMA1      = metal_forces_calc_gamma_vmet(lL, lK0, 2, lPPP);
    T_ lGAMMA2      = metal_forces_calc_gamma_dmet(lL, lK0, 2, lPPP);
    int4 lXZ        = lAI>lAJ ? make_int4(1,2,2,2) : make_int4(2,2,1,2);


    T_ lGAMMA       = (lGAMMA1 - lGAMMA2*(CONSTANT_DATA.mMET.mRHO[aIATM-1]*(*dmet<T_>(lXZ.x,lK0,lXZ.y)) +
                                          CONSTANT_DATA.mMET.mRHO[aJATM-1]*(*dmet<T_>(lXZ.z,lK0,lXZ.w)))) / aRSQ;

    T_  lFX        = mul(lGAMMA, aXDF);
    T_  lFY        = mul(lGAMMA, aYDF);
    T_  lFZ        = mul(lGAMMA, aZDF);
    aFX_I         += lFX;
    aFY_I         += lFY;
    aFZ_I         += lFZ;

    if (aJATM <= CONSTANT_DATA.mNATMS) {
      aFX_J      -= lFX;
      aFY_J      -= lFY;
      aFZ_J      -= lFZ;
    }

    if (aJATM <= CONSTANT_DATA.mNATMS || aIDI < aIDJ) {
      T_ lVK0 = *vmet<T_>(lL+3,lK0,1);
      T_ lVK1 = *vmet<T_>(lL+4,lK0,1);
      T_ lVK2 = *vmet<T_>(lL+5,lK0,1);

      T_ lT1  = madd(lVK1, lPPP, (lVK1 - lVK0));
      T_ lT2  = madd(lVK1, lPPP, (lVK2 - lVK1));

      aENGMET += metal_forces_calc_gamma(lPPP, lT1, lT2);
      aVIRMET -= mul(lGAMMA, aRSQ);

      aSTRS1       += mul(aXDF, lFX);
      aSTRS2       += mul(aXDF, lFY);
      aSTRS3       += mul(aXDF, lFZ);
      aSTRS5       += mul(aYDF, lFY);
      aSTRS6       += mul(aYDF, lFZ);
      aSTRS9       += mul(aZDF, lFZ);
    }
  }
}

template<typename T_, int ISKEYPOTZERO_> __device__ inline void metal_forces
#if (CFG_COMPUTE_MAJOR >= 2)
    (int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
     volatile T_& aFX_I, volatile T_& aFY_I, volatile T_& aFZ_I,
     T_& aFX_J, T_& aFY_J, T_& aFZ_J,
     T_& aSTRS1, T_& aSTRS2, T_& aSTRS3, volatile T_& aSTRS5, volatile T_& aSTRS6, volatile T_& aSTRS9, T_& aENGMET, T_& aVIRMET){
#else
    (int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
     T_& aFX_I, T_& aFY_I, T_& aFZ_I,
     T_& aFX_J, T_& aFY_J, T_& aFZ_J,
     T_& aSTRS1, T_& aSTRS2, T_& aSTRS3, T_& aSTRS5, T_& aSTRS6, T_& aSTRS9, T_& aENGMET, T_& aVIRMET){
#endif

  if (!ISKEYPOTZERO_) {
    metal_forces_KEYPOTNeqZero(aIATM, aJATM, aIDI, aIDJ, aXDF, aYDF, aZDF, aRSQ,
                               aFX_I, aFY_I, aFZ_I, aFX_J, aFY_J, aFZ_J,
                               aSTRS1, aSTRS2, aSTRS3, aSTRS5, aSTRS6, aSTRS9, aENGMET, aVIRMET);
    return;
  }
}

template<typename T_> __device__ void vdw_forces(
#if (CFG_COMPUTE_MAJOR >= 2)
    int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
    volatile T_& aFX_I, volatile T_& aFY_I, volatile T_& aFZ_I,
    T_& aFX_J, T_& aFY_J, T_& aFZ_J,
    T_& aSTRS1, T_& aSTRS2, T_& aSTRS3,
    volatile T_& aSTRS5, volatile T_& aSTRS6, volatile T_& aSTRS9,
    volatile T_& aENGVDW, volatile T_& aVIRVDW){
#else
  int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
    T_& aFX_I, T_& aFY_I, T_& aFZ_I,
    T_& aFX_J, T_& aFY_J, T_& aFZ_J,
    T_& aSTRS1, T_& aSTRS2, T_& aSTRS3, T_& aSTRS5, T_& aSTRS6, T_& aSTRS9, T_& aENGVDW, T_& aVIRVDW){
#endif

  int lAI       = CONSTANT_DATA.mLTYPE[aIATM-1];
  int lAJ       = CONSTANT_DATA.mLTYPE[aJATM-1];
  int lKEY      = max(lAI,lAJ) * (max(lAI,lAJ)-1) / 2 + min(lAI,lAJ);
  int lK        = CONSTANT_DATA.mVDW.mLSTVDW[lKEY - 1];
  T_  lVVDW_1_K = *F2D_ADDRESS(CONSTANT_DATA.mVDW.mVVDW, 1, 1, CONSTANT_DATA.mMXGRID, 1, lK);
  int lLTPVDW_K = CONSTANT_DATA.mVDW.mLTPVDW[lK - 1];

  if (lLTPVDW_K>=0 && fabs(lVVDW_1_K) > (T_)0 && aRSQ < CONSTANT_DATA.mVDW.mRCSQ) {
    T_  lRRR           = sqrt(aRSQ);
    int lL             = (int) trunc(mul(lRRR, CONSTANT_DATA.mVDW.mRDR));
    T_  lGVDW_L0_K     = fetch_real_1d(TEXTURE_DATA_VDW_GVDW, (lK-1)*CONSTANT_DATA.mMXGRID + lL+0-1);
    T_  lGVDW_L1_K     = fetch_real_1d(TEXTURE_DATA_VDW_GVDW, (lK-1)*CONSTANT_DATA.mMXGRID + lL+1-1);
    T_  lGVDW_L2_K     = fetch_real_1d(TEXTURE_DATA_VDW_GVDW, (lK-1)*CONSTANT_DATA.mMXGRID + lL+2-1);
    T_  lGVDW_MXGRID_K = fetch_real_1d(TEXTURE_DATA_VDW_GVDW, (lK-1)*CONSTANT_DATA.mMXGRID + CONSTANT_DATA.mMXGRID-1);

    T_  lPPP       = (lRRR*CONSTANT_DATA.mVDW.mRDR) - (T_)lL;
    T_  lGT1       = madd(lGVDW_L0_K, (lGVDW_L1_K - lGVDW_L0_K), lPPP);
    T_  lGT2       = madd(lGVDW_L1_K, (lGVDW_L2_K - lGVDW_L1_K), (lPPP - (T_)1));
    T_  lGAMMA     = (madd(lGT1, (lGT2 - lGT1)*lPPP, ((T_) 0.5))) / aRSQ;

    T_   lRVDW     = CONSTANT_DATA.mVDW.mRVDW;
    bool lLS_VDW   = CONSTANT_DATA.mVDW.mLS_VDW;

    if (lLS_VDW){
      lGAMMA -= lGVDW_MXGRID_K/(mul(lRRR, lRVDW));
    }

    T_  lFX        = mul(lGAMMA, aXDF);
    T_  lFY        = mul(lGAMMA, aYDF);
    T_  lFZ        = mul(lGAMMA, aZDF);
    aFX_I         += lFX;
    aFY_I         += lFY;
    aFZ_I         += lFZ;

    if (aJATM <= CONSTANT_DATA.mNATMS) {
      aFX_J      -= lFX;
      aFY_J      -= lFY;
      aFZ_J      -= lFZ;
    }

    if (aJATM <= CONSTANT_DATA.mNATMS || aIDI < aIDJ) {
      T_ lVVDW_L0_K = fetch_real_1d(TEXTURE_DATA_VDW_VVDW, (lK-1)*CONSTANT_DATA.mMXGRID + lL+0-1);
      T_ lVVDW_L1_K = fetch_real_1d(TEXTURE_DATA_VDW_VVDW, (lK-1)*CONSTANT_DATA.mMXGRID + lL+1-1);
      T_ lVVDW_L2_K = fetch_real_1d(TEXTURE_DATA_VDW_VVDW, (lK-1)*CONSTANT_DATA.mMXGRID + lL+2-1);
      T_ lVT1       = madd(lVVDW_L0_K, (lVVDW_L1_K - lVVDW_L0_K), lPPP);
      T_ lVT2       = madd(lVVDW_L1_K, (lVVDW_L2_K - lVVDW_L1_K), (lPPP - (T_)1));

      aENGVDW      += madd(lVT1, (lVT2 - lVT1), mul(lPPP, (T_)0.5));

      //RN :: If (ls_vdw) engvdw = engvdw - gvdw(mxgrid,k)*(rrr/rvdw-1.0_wp) - vvdw(mxgrid,k)
      if (lLS_VDW){
        T_ lVVDW_MXGRID_K = fetch_real_1d(TEXTURE_DATA_VDW_VVDW,
                                          (lK-1)*CONSTANT_DATA.mMXGRID + CONSTANT_DATA.mMXGRID-1);
        aENGVDW    -= add(mul(lGVDW_MXGRID_K, (lRRR/lRVDW - (T_)1.0)), lVVDW_MXGRID_K);
      }

      aVIRVDW      += -mul(lGAMMA, aRSQ);
      aSTRS1       += mul(aXDF, lFX);
      aSTRS2       += mul(aXDF, lFY);
      aSTRS3       += mul(aXDF, lFZ);
      aSTRS5       += mul(aYDF, lFY);
      aSTRS6       += mul(aYDF, lFZ);
      aSTRS9       += mul(aZDF, lFZ);
    }
  }
}

template<typename T_> __device__ inline void erfcgen(int aI, T_& aERC, T_& aFER) {
  aERC = fetch_real_1d(TEXTURE_DATA_EWR_ERC, aI-1);
  aFER = fetch_real_1d(TEXTURE_DATA_EWR_FER, aI-1);
}


template<typename T_> __device__ void ewald_real_forces(
#if (CFG_COMPUTE_MAJOR >= 2)
    int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
    volatile T_& aFX_I, volatile T_& aFY_I, volatile T_& aFZ_I,
    T_& aFX_J, T_& aFY_J, T_& aFZ_J,
    T_& aSTRS1, T_& aSTRS2, T_& aSTRS3,
    volatile T_& aSTRS5, volatile T_& aSTRS6, volatile T_& aSTRS9,
    volatile T_& aENGCPE_RL, volatile T_& aVIRCPE_RL) {
#else
    int aIATM, int aJATM, int aIDI, int aIDJ, T_ aXDF, T_ aYDF, T_ aZDF, T_ aRSQ,
    T_& aFX_I, T_& aFY_I, T_& aFZ_I,
    T_& aFX_J, T_& aFY_J, T_& aFZ_J,
    T_& aSTRS1, T_& aSTRS2, T_& aSTRS3, T_& aSTRS5, T_& aSTRS6, T_& aSTRS9, T_& aENGCPE_RL, T_& aVIRCPE_RL) {
#endif

  T_ lCHGEA = CONSTANT_DATA.mEWR.mCHGE[aIATM - 1];

  if (fabs(lCHGEA) > (T_)0) {

     lCHGEA = mul(lCHGEA, CONSTANT_DATA.mEWR.mR4PIE0) / CONSTANT_DATA.mEWR.mEPSQ;
     T_ lCHGPRD = CONSTANT_DATA.mEWR.mCHGE[aJATM - 1];

     if (fabs(lCHGPRD) > (T_)0) {
        lCHGPRD = mul(lCHGPRD, lCHGEA);

        if (aRSQ < CONSTANT_DATA.mEWR.mRCSQ)   {
          T_  lRRR       = sqrt(aRSQ);
          int lK         = (int) trunc(mul(lRRR, CONSTANT_DATA.mEWR.mRDREWD));

          T_ lERC0, lERC1, lERC2, lFER0, lFER1, lFER2;
          erfcgen(lK+0, lERC0, lFER0);
          erfcgen(lK+1, lERC1, lFER1);
          erfcgen(lK+2, lERC2, lFER2);

          const T_  lPPP    = madd(-(T_)lK, lRRR, CONSTANT_DATA.mEWR.mRDREWD);
          const T_  lFT1    = madd(lFER0, (lFER1 - lFER0), lPPP);
          const T_  lFT2    = madd(lFER1, (lFER2 - lFER1), (lPPP - ((T_) 1.0)));
          const T_  lEGAMMA = mul(madd(lFT1, (lFT2 - lFT1), lPPP*(T_)0.5), lCHGPRD);
          T_ lFX               = mul(lEGAMMA, aXDF);
          T_ lFY               = mul(lEGAMMA, aYDF);
          T_ lFZ               = mul(lEGAMMA, aZDF);

          aFX_I        += lFX;
          aFY_I        += lFY;
          aFZ_I        += lFZ;
          if (aJATM <= CONSTANT_DATA.mNATMS) {
             aFX_J     -= lFX;
             aFY_J     -= lFY;
             aFZ_J     -= lFZ;
          }
          if (aJATM <= CONSTANT_DATA.mNATMS || aIDI < aIDJ) {
             const T_ lET1 = madd(lERC0, (lERC1 - lERC0),lPPP);
             const T_ lET2 = madd(lERC1, (lERC2 - lERC1), (lPPP - (T_)1.0));
             aENGCPE_RL   += mul(madd(lET1, (lET2 - lET1), mul(lPPP, (T_)0.5)), lCHGPRD);
             aVIRCPE_RL   += -mul(lEGAMMA,aRSQ);
             aSTRS1       += mul(aXDF, lFX);
             aSTRS2       += mul(aXDF, lFY);
             aSTRS3       += mul(aXDF, lFZ);
             aSTRS5       += mul(aYDF, lFY);
             aSTRS6       += mul(aYDF, lFZ);
             aSTRS9       += mul(aZDF, lFZ);
         }
      }
    }
  }
}


/* It is easy to run out of register space when many threads are invovled. Instead of
 * performing both computations at once, this is done in stages. If RELOAD_ is set to
 * 1, the force, stress, virial etc. updates will be added to what had been previously
 * computed. For instance, vdw_forces may run with RELOAD_==0, and ewalds_real_forces
 * after with RELOAD_==1.
 *
 * @param aStripeSizeInItems The number of T_ items to be stored on a per-list basis.
 * @param aOUT The destination buffer for the f{xx,yy,zz} updates against the primary
 *        atom in each list and the secondary atoms for each list. The organisation
 *        for each list, in T_ datums, is as follows; consequtive lists differ by
 *        aStripeSizeInItems*sizeof(T_) bytes; secondary atoms have a SoA layout:
 *         01: primary atom's fxx update
 *         02: primary atom's fyy update
 *         03: primary atom's fzz update
 *         next limit(0,i) datums: secondary atoms' fxx updates
 *         next limit(0,i) datums: secondary atoms' fyy updates
 *         next limit(0,i) datums: secondary atoms' fzz updates
 */
// for compute >=2.0 & cudaFuncSetCacheConfig
#define CFG_K1_CPP_MAGLED_NAME_FORMATTER\
  "_Z23two_body_forces_cuda_k1I%cLj%dELi%dELi%dELi%dELi%dELi%dEEviiPT_S1_i"

template<typename T_, unsigned int BX_, int IMCON_, int METFORCES_, int ISKEYPOTZERO_, int VDWFORCES_, int KEYFCE_>
__global__ void two_body_forces_cuda_k1(int aUnroll, int aIATM_Begin, T_ *aOUT, T_ *aOUT_Red, int aStripeSizeInItems) {
  DECLARE_DYNAMIC_SHARED(T_);

#if (CFG_COMPUTE_MAJOR >= 2)
#define lSTRS5 shared[3*BX_ + threadIdx.x]
#define lSTRS6 shared[4*BX_ + threadIdx.x]
#define lSTRS9 shared[5*BX_ + threadIdx.x]

  T_ lSTRS1=0, lSTRS2=0, lSTRS3=0;
  lSTRS5=0;
  lSTRS6=0;
  lSTRS9=0;

#define lENGVDW shared[6*BX_ + threadIdx.x]
#define lVIRVDW shared[7*BX_ + threadIdx.x]

  lENGVDW=0; lVIRVDW=0; T_ lENGMET=0;

#define lENGEWR shared[8*BX_ + threadIdx.x]
#define lVIREWR shared[9*BX_ + threadIdx.x]

  lENGEWR=0; lVIREWR=0; T_ lVIRMET=0;

#else

  T_ lSTRS1=0, lSTRS2=0, lSTRS3=0, lSTRS5=0, lSTRS6=0, lSTRS9=0;
  T_ lENGVDW=0, lVIRVDW=0, lENGMET=0;
  T_ lENGEWR=0, lVIREWR=0, lVIRMET=0;

#endif


  for (int lQ=blockIdx.y ; lQ<aUnroll ; lQ+=gridDim.y) {
    int lIATM = aIATM_Begin + lQ;

    int lLIMIT, lIDI;
    T_ lXXX_I, lYYY_I, lZZZ_I;
    obtain_iatm_specific<T_>(lIATM, lXXX_I, lYYY_I, lZZZ_I, lLIMIT, lIDI);

    if (lLIMIT>0) {

      /* Initialise everything to 0 so that those threads used for padding limit to
       * a power of 2 will not contribute data that may alter the reduction results.
       */
#if (CFG_COMPUTE_MAJOR >= 2)
#define lFX_I shared[0*BX_ + threadIdx.x]
#define lFY_I shared[1*BX_ + threadIdx.x]
#define lFZ_I shared[2*BX_ + threadIdx.x]

      lFX_I=0; lFY_I=0; lFZ_I=0;

#else

      T_ lFX_I=0, lFY_I=0, lFZ_I=0;

#endif

      T_ *lOUT    = aOUT + (3 + 3*aStripeSizeInItems)*lQ;
      T_ *lOUT_FJ = lOUT + 3;

      for (int lTIdx = threadIdx.x ; lTIdx<lLIMIT ; lTIdx+=BX_) {
        T_ lFX_J=0, lFY_J=0, lFZ_J=0;
        int lJATM = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1), 1+lTIdx, lIATM);
        int lIDJ;
        T_ lXDF, lYDF, lZDF, lRSQDF;

        obtain_jatm_specific<T_,IMCON_>(lJATM, lXXX_I, lYYY_I, lZZZ_I, lXDF, lYDF, lZDF, lRSQDF, lIDJ);

        if (METFORCES_) {
          metal_forces<T_,ISKEYPOTZERO_>(lIATM, lJATM, lIDI, lIDJ, lXDF, lYDF, lZDF, lRSQDF,
                                         lFX_I, lFY_I, lFZ_I, lFX_J, lFY_J, lFZ_J, lSTRS1,
                                         lSTRS2, lSTRS3, lSTRS5, lSTRS6, lSTRS9, lENGMET, lVIRMET);
        }

        if (VDWFORCES_) {
          vdw_forces<T_>(lIATM, lJATM, lIDI, lIDJ, lXDF, lYDF, lZDF, lRSQDF,
                         lFX_I, lFY_I, lFZ_I, lFX_J, lFY_J, lFZ_J, lSTRS1,
                         lSTRS2, lSTRS3, lSTRS5, lSTRS6, lSTRS9, lENGVDW, lVIRVDW);
        }

        if (KEYFCE_==2) {
          ewald_real_forces<T_>(lIATM, lJATM, lIDI, lIDJ, lXDF, lYDF, lZDF, lRSQDF,
                                lFX_I, lFY_I, lFZ_I, lFX_J, lFY_J, lFZ_J, lSTRS1,
                                lSTRS2, lSTRS3, lSTRS5, lSTRS6, lSTRS9, lENGEWR, lVIREWR);
        }
        // store the f{xx, yy, zz}(jatm) data:
        lOUT_FJ[0*CFG_MAX_LIMIT + lTIdx] = lFX_J;
        lOUT_FJ[1*CFG_MAX_LIMIT + lTIdx] = lFY_J;
        lOUT_FJ[2*CFG_MAX_LIMIT + lTIdx] = lFZ_J;
      }

#if (CFG_COMPUTE_MAJOR < 2)
      // in >=2.0, data are already there
      shared[ 0*BX_ + threadIdx.x] = lFX_I;
      shared[ 1*BX_ + threadIdx.x] = lFY_I;
      shared[ 2*BX_ + threadIdx.x] = lFZ_I;
#endif

      if (BX_>32) __syncthreads();

#if !CFG_UNIFIED_ADDRESS_SPACE
      psum<T_,BX_,3>();
#else
      psum_uas<T_,BX_,3>(shared);
#endif

      if (BX_>=3) {
        if (threadIdx.x<3) lOUT[threadIdx.x] = shared[threadIdx.x*BX_];
      } else {
        if (threadIdx.x==0) {
          lOUT[0] = shared[ 0*BX_];
          lOUT[1] = shared[ 1*BX_];
          lOUT[2] = shared[ 2*BX_];
        }
      }
      if (BX_>32) __syncthreads();
    }
  }

  shared[0*BX_ + threadIdx.x] = lSTRS1;
  shared[1*BX_ + threadIdx.x] = lSTRS2;
  shared[2*BX_ + threadIdx.x] = lSTRS3;
  shared[3*BX_ + threadIdx.x] = lSTRS5;
  shared[4*BX_ + threadIdx.x] = lSTRS6;
  shared[5*BX_ + threadIdx.x] = lSTRS9;

  if (VDWFORCES_) {
    shared[6*BX_ + threadIdx.x] = lENGVDW;
    shared[7*BX_ + threadIdx.x] = lVIRVDW;
  }

  if (KEYFCE_>0) {
    shared[8*BX_ + threadIdx.x] = lENGEWR;
    shared[9*BX_ + threadIdx.x] = lVIREWR;
  }

  __syncthreads();

#if !CFG_UNIFIED_ADDRESS_SPACE
  psum<T_,BX_,10>();
#else
  psum_uas<T_,BX_,10>(shared);
#endif

  if (threadIdx.x==0) {
    aOUT_Red[gridDim.y*0 + blockIdx.y]  = shared[ 0*BX_];
    aOUT_Red[gridDim.y*1 + blockIdx.y]  = shared[ 1*BX_];
    aOUT_Red[gridDim.y*2 + blockIdx.y]  = shared[ 2*BX_];
    aOUT_Red[gridDim.y*3 + blockIdx.y]  = shared[ 3*BX_];
    aOUT_Red[gridDim.y*4 + blockIdx.y]  = shared[ 4*BX_];
    aOUT_Red[gridDim.y*5 + blockIdx.y]  = shared[ 5*BX_];

    if (VDWFORCES_) {
      aOUT_Red[gridDim.y*8 + blockIdx.y] = shared[ 6*BX_];
      aOUT_Red[gridDim.y*9 + blockIdx.y] = shared[ 7*BX_];
    }

    if (KEYFCE_>0) {
      aOUT_Red[gridDim.y*10 + blockIdx.y] = shared[ 8*BX_];
      aOUT_Red[gridDim.y*11 + blockIdx.y] = shared[ 9*BX_];
    }
  }

  if (METFORCES_) {
    __syncthreads();
    shared[0*BX_ + threadIdx.x] = lENGMET;
    shared[1*BX_ + threadIdx.x] = lVIRMET;
    __syncthreads();

#if !CFG_UNIFIED_ADDRESS_SPACE
    psum<T_,BX_,2>();
#else
    psum_uas<T_,BX_,2>(shared);
#endif

    if (threadIdx.x==0) {
      aOUT_Red[gridDim.y*6 + blockIdx.y] = shared[0*BX_];
      aOUT_Red[gridDim.y*7 + blockIdx.y] = shared[1*BX_];
    }
  }
}


/* k3 reduces the strs*, energy and virial data.
 */
template<typename T_, unsigned int B_> __global__ void two_body_forces_cuda_k3(
    unsigned int aIterations, T_ *aIN, T_ *aOUT) {
  DECLARE_DYNAMIC_SHARED(T_);

  T_ lVAL = (T_)0;
  for (unsigned int lI=0 ; lI<aIterations ; lI+=blockDim.x) {
    unsigned int lMyIndex = lI + threadIdx.x;
    if (lMyIndex < aIterations) {
      lVAL += aIN[blockIdx.y*aIterations + lMyIndex];
    }
  }
  shared[threadIdx.x] = lVAL;
  __syncthreads();
#if !CFG_UNIFIED_ADDRESS_SPACE
  psum<T_, B_, 1>();
#else
  psum_uas<T_, B_, 1>(shared);
#endif

  if (threadIdx.x==0) {
    aOUT[blockIdx.y] = shared[0];
  }
}

template<typename T_, unsigned int BX_>
__global__ void two_body_forces_cuda_k2(int lIATM_Begin, int lUnroll, T_ *aIN, int lIN_StripeSizeInItems, T_ *aOUT) {

  T_ *lOUT = aOUT + blockIdx.x*CONSTANT_DATA.mMXATMS + blockIdx.y*3*CONSTANT_DATA.mMXATMS;

  int *lListBaseAddress = F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1), 0, lIATM_Begin+blockIdx.y);

  for (int lI=blockIdx.y ; lI<lUnroll ; lI+=gridDim.y, lListBaseAddress+=gridDim.y*(CONSTANT_DATA.mMXLIST+1)) {
    int   lLIMIT   = lListBaseAddress[0];
    real *lIN_Base = aIN + lI*lIN_StripeSizeInItems + 3 + blockIdx.x*CFG_MAX_LIMIT;

    if (BX_>32) __syncthreads(); // wait for gmem stores from prev round

    for (int lJ=threadIdx.x ; lJ<lLIMIT ; lJ+=BX_) {
      int lJATM = lListBaseAddress[lJ+1];
      T_ lV = lIN_Base[lJ];

      if (lV != (T_)0) {
        lOUT[lJATM-1] += lIN_Base[lJ];
      }
    }
  }

  if (BX_>32) __syncthreads();

  for (int lI=blockIdx.y+threadIdx.x*gridDim.y ; lI<lUnroll ; lI+=BX_*gridDim.y) {
    if ((*F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1), 0, lIATM_Begin+lI))>0) {
      T_ lV = aIN[lI*lIN_StripeSizeInItems + blockIdx.x];
      if (lV != (T_)0)
        lOUT[lIATM_Begin+lI-1] += lV;
    }
  }

}


/* This is here for debugging purposes; set CFG_K2_DEBUG to 1 to enable it.
 */
template<typename T_>
__global__ void two_body_forces_cuda_k2_debug
    (int lIATM_Begin, int lUnroll, T_ *aIN, int lIN_StripeSizeInItems, T_ *aOUT) {

  for (int lI=lIATM_Begin ; lI<=lIATM_Begin+lUnroll-1 ; lI++) {
    int lLIMIT = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1), 0, lI);

    if (threadIdx.x<lLIMIT) {
      int lJATM = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1), 1+threadIdx.x, lI);
      aOUT[blockIdx.x*CONSTANT_DATA.mMXATMS + lJATM -1] +=
          aIN[(lI-lIATM_Begin)*lIN_StripeSizeInItems + 3 + blockIdx.x*CFG_MAX_LIMIT + threadIdx.x];

      if (threadIdx.x==0) {
        aOUT[blockIdx.x*CONSTANT_DATA.mMXATMS + lI-1] += aIN[(lI-lIATM_Begin)*lIN_StripeSizeInItems + blockIdx.x];
      }
    }
    __syncthreads();
  }
}


template<typename T_, unsigned int GY_, unsigned int BX_, int SL_>
__global__ void two_body_forces_cuda_k2_reduce(int lUnroll, T_ *aOUT) {
  /* If only a single slice has been modified then that is the first one and
   * thus there's nothing for us to do.
   */
  if (SL_==1){
    return;
  }else{
    int lInterSliceGapInItems = 3*CONSTANT_DATA.mMXATMS;
    for (int lI=blockIdx.y*BX_+threadIdx.x ; lI<CONSTANT_DATA.mNATMS ; lI+=GY_*BX_) {
        T_ *lBase = aOUT + blockIdx.x*CONSTANT_DATA.mMXATMS;
        T_ lR = (T_)0;
        for (int lJ=0 ; lJ<SL_ ; lJ++) {
            lR += lBase[lI + lJ*lInterSliceGapInItems];
            if (lJ>0)
                lBase[lI + lJ*lInterSliceGapInItems]=(T_)0;
        }
        lBase[lI] = lR;
    }
  }
}

template<typename T_, unsigned int B_, int METFORCES_, int ISKEYPOTZERO_, int VDWFORCES_, int KEYFCE_>
inline void two_body_forces_cuda_k1_imcon_switch(dim3 aK1_GridDims, int aSharedMemorySize, int aUnroll,
                                                 int aIATM_Begin, T_ *aOUT, T_ *aOUT_Red, int aStripeSizeInItems) {

#if (CFG_COMPUTE_MAJOR >= 2)
  char lKernelMangledName[256];
  sprintf(lKernelMangledName, CFG_K1_CPP_MAGLED_NAME_FORMATTER,
          sizeof(T_)==8 ? 'd' : 'f', B_, sHD.mIMCON, METFORCES_, ISKEYPOTZERO_, VDWFORCES_, KEYFCE_);
  CUDA_SAFE_CALL(cudaFuncSetCacheConfig(lKernelMangledName, cudaFuncCachePreferL1));
#endif

#define INVOKE_IMCON(IMCON)\
  two_body_forces_cuda_k1<real, B_, IMCON, METFORCES_, ISKEYPOTZERO_, VDWFORCES_, KEYFCE_>\
      <<<aK1_GridDims, B_, aSharedMemorySize>>>(aUnroll, aIATM_Begin, aOUT, aOUT_Red, aStripeSizeInItems)

  switch (sHD.mIMCON) {
  case 2: { INVOKE_IMCON(2); break; }
  case 3: { INVOKE_IMCON(3); break; }
  default: { printf("%s::%s: (E) stub.\n",__FILE__,__FUNCTION__); exit(-1); }
  }
}

template<typename T_, unsigned int B_, int METFORCES_, int ISKEYPOTZERO_, int VDWFORCES_>
inline void two_body_forces_cuda_k1_imcon_keyfce_switch(dim3 aK1_GridDims, int aSharedMemorySize,
                                                        int aUnroll, int aIATM_Begin, T_ *aOUT,
                                                        T_ *aOUT_Red, int aStripeSizeInItems) {

#define INVOKE_KEYFCE(KEYFCE)\
  two_body_forces_cuda_k1_imcon_switch<real, B_, METFORCES_, ISKEYPOTZERO_, VDWFORCES_, KEYFCE>\
      (aK1_GridDims, aSharedMemorySize, aUnroll, aIATM_Begin, aOUT, aOUT_Red, aStripeSizeInItems)

  switch (sHD.mKEYFCE) {
    case 0: { INVOKE_KEYFCE(0); break; }
    case 2: { INVOKE_KEYFCE(2); break; }
    default: { printf("%s::%s: (E) stub.\n",__FILE__,__FUNCTION__); exit(-1); }
  }
}

template<typename T_, unsigned int B_, int METFORCES_, int ISKEYPOTZERO_>
inline void two_body_forces_cuda_k1_imcon_keyfce_vdwforces_switch(
    dim3 aK1_GridDims, int aSharedMemorySize, int aUnroll, int aIATM_Begin, T_ *aOUT, T_ *aOUT_Red, int aStripeSizeInItems) {

#define INVOKE_VDWFORCES(VDWFORCES)\
  two_body_forces_cuda_k1_imcon_keyfce_switch<real, B_, METFORCES_, ISKEYPOTZERO_, VDWFORCES>\
      (aK1_GridDims, aSharedMemorySize, aUnroll, aIATM_Begin, aOUT, aOUT_Red, aStripeSizeInItems)

  if (sHD.mNTPVDW>0) {
    INVOKE_VDWFORCES(1);
  } else {
    INVOKE_VDWFORCES(0);
  }
}

template<typename T_, unsigned int B_>
inline void two_body_forces_cuda_k1_imcon_keyfce_iskeypotzero_metforces_switch
    (dim3 aK1_GridDims, int aSharedMemorySize, int aUnroll, int aIATM_Begin,
     T_ *aOUT, T_ *aOUT_Red, int aStripeSizeInItems) {

#define INVOKE_ISKEYPOTZERO_METFORCES(METFORCES,ISKEYPOTZERO)					\
  two_body_forces_cuda_k1_imcon_keyfce_vdwforces_switch<real, B_, METFORCES, ISKEYPOTZERO>\
      (aK1_GridDims, aSharedMemorySize, aUnroll, aIATM_Begin, aOUT, aOUT_Red, aStripeSizeInItems)

  if (sHD.mNTPMET>0) {
    if (sCD.mMET.mKEYPOT==0) {
      INVOKE_ISKEYPOTZERO_METFORCES(1,1);
    } else {
      INVOKE_ISKEYPOTZERO_METFORCES(1,0);
    }
  } else {
    INVOKE_ISKEYPOTZERO_METFORCES(0,0);
  }
}
// shorter macro, for clarity
#define two_body_forces_cuda_k1_select two_body_forces_cuda_k1_imcon_keyfce_iskeypotzero_metforces_switch


template<typename T_>
inline void two_body_forces_cuda_k2_bx_switch
    (int aIATM_DevBeginC, int aUnroll_D, T_* aOUT_DEV, int aStripeSizeInItems) {

#define INVOKE_BX(BX)							\
  two_body_forces_cuda_k2<T_, BX><<<dim3(3,sHD.mK2_Unroll_Slices,1), dim3(BX,1,1)>>>	\
      (aIATM_DevBeginC, aUnroll_D, aOUT_DEV, aStripeSizeInItems, sHD.mDEV_Spad0)

  switch (sHD.mK2_Block_Dims_X) {
    case 512: { INVOKE_BX(512); break; }
    case 256: { INVOKE_BX(256); break; }
    case 128: { INVOKE_BX(128); break; }
    case  64: { INVOKE_BX( 64); break; }
    case  32: { INVOKE_BX( 32); break; }
    case  16: { INVOKE_BX( 16); break; }
    default: {
      printf("%s::%s (e) stub reached; sHD.mK2_Block_Dims_X==%d\n",
             __FILE__, __FUNCTION__, sHD.mK2_Block_Dims_X);
      exit(-1);
    }
  }
}

extern "C" void two_body_forces_cuda_invoke(
                  real *aFXX, real *aFYY, real *aFZZ,
                  real *aSTRS,
                  real *aENGMET, real *aVIRMET,
                  real *aENGVDW, real *aVIRVDW,
                  real *aENGCPE_RL, real *aVIRCPE_RL,
                  int *aSAFE) {

  /* The host will work on the first iterations, with the remaining left for
   * the device.
   */
  real *lF_Pointers[3] = { aFXX, aFYY, aFZZ };
  int lUnroll_Dp       = (int) (sHD.mPercentageOfIterationsOffloadedToTheDevice * (double) sHD.mUnroll);
  int lUnroll_Hp       = sHD.mUnroll - lUnroll_Dp;
  int lLen_D           = (int) (sHD.mPercentageOfIterationsOffloadedToTheDevice * sCD.mNATMS);
  int lLen_H           = sCD.mNATMS - lLen_D;

  int lIATM_DevBegin, lIATM_HstBegin;
  if (lUnroll_Hp==0) {
    lIATM_DevBegin = 1;
    lIATM_HstBegin = sCD.mNATMS+1;
  } else {
    if (lUnroll_Dp==0) {
      lIATM_HstBegin = 1;
      lIATM_DevBegin = sCD.mNATMS+1;
    } else {
      lIATM_HstBegin = 1;
      lIATM_DevBegin = 1 + lLen_H;
    }
  }

  int lIATM_HstBeginC = lIATM_HstBegin;
  int lIATM_DevBeginC = lIATM_DevBegin;

  /* If the link cell pairs kernel has not run on the device (i.e. when
   * link_cell_pairs_cuda_is_in_valid_context()==0), the initialiser
   * would have already copied the list over. Otherwise, we copy in and
   * out right now.
   */
  if (lUnroll_Hp>0) {
    if (link_cell_pairs_cuda_is_in_valid_context()) {
      link_cell_pairs_cuda_pull_lists(lLen_H);
    }
  }
  if (lUnroll_Dp>0) {
    if (link_cell_pairs_cuda_is_in_valid_context())
      link_cell_pairs_cuda_push_lists(lLen_H+1);
  }

  double lDurationOfK1SpentInHost = 0.0, lDurationOfK1SpentInDevice = 0.0;
  double lDurationOfK2SpentInHost = 0.0, lDurationOfK2SpentInDevice = 0.0;

  for (;;) {
    /* Reduce the unroll factors if necessary to avoid out of bounds errors:
     */
    int lUnroll_D = min(lUnroll_Dp, sCD.mNATMS - lIATM_DevBeginC + 1);
    int lUnroll_H = lUnroll_Hp!=0 ? min(lUnroll_Hp, (lIATM_DevBegin-1) - lIATM_HstBeginC + 1) : 0;
    int lUnroll_H_WithK1 = (int) (sHD.mPercentageOfHostIterationsToBeOverlappedWithK1 * ((double) lUnroll_H));
    int lUnroll_H_WithRD = lUnroll_H - lUnroll_H_WithK1;

    if ((lIATM_HstBeginC>=lIATM_DevBegin || lUnroll_Hp==0) && (lIATM_DevBeginC>sCD.mNATMS || lUnroll_Dp==0))
      break;

    int lLIMIT_MaxInH = 512;
    int lOUT_StripeSizeInItems = 3+3*lLIMIT_MaxInH;
    real *lOUT_HST_FinalRed;
    real *lOUT_DEV;
    real *lOUT_DEV_Red;
    real *lOUT_DEV_FinalRed;

    int lBytesNeeded = (12 + (lOUT_StripeSizeInItems + 12)*lUnroll_D)*sizeof(real);
    if (sHD.mScratchpadSize < lBytesNeeded) {
      printf("%s::%s: Out of global memory: requested %db (lOUT_StripeSizeInItems=%d,"
             " lUnroll_D=%d, lLIMIT_MaxInH=%d), can reserve up to %db\n",
             __FILE__, __FUNCTION__,
             (lOUT_StripeSizeInItems + 16)*lUnroll_D*sizeof(real),
             lOUT_StripeSizeInItems, lUnroll_D, lLIMIT_MaxInH, sHD.mScratchpadSize);
      exit(-1);
    }
    lOUT_DEV          = sHD.mDEV_Scratchpad;
    lOUT_DEV_Red      = lOUT_DEV + lUnroll_D*lOUT_StripeSizeInItems;
    lOUT_DEV_FinalRed = lOUT_DEV_Red + lUnroll_D*12;
    lOUT_HST_FinalRed = sHD.mHST_PinnedScratchpad;

    dim3 lK1_GridDims(1,sHD.mK1_Grid_Dims_Y,1);
    int lNT = CFG_K1_THREADS_PER_BLOCK;

    int lSharedMemorySize = 10*lNT*sizeof(real);// + lNT*4*sizeof(real);
    if (lSharedMemorySize > (16*1024 - 128)) {
      printf("%s::%s: Out of shared memory: requested %db, can reserve up to %db\n",
             __FILE__, __FUNCTION__, lSharedMemorySize, (16*1024 - 128));
      printf("%s::%s: threads/block = %d\n",
             __FILE__, __FUNCTION__, lNT);
      exit(-1);
    }

    /* === K1 begins here
     */
    start_timing_two_body_forces_cuda_k1();
    struct timeval lK1_TV_Begin;
    gettimeofday(&lK1_TV_Begin, NULL);
    cudaEvent_t lCE_K1_DevStart, lCE_K1_DevStop;
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStart));
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStop));

    /* Execute the K1 device portion:
     * TODO: fix the IMCON_ argument!!
     */
    if (lUnroll_D>0) {
      CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStart, 0));

      dim3 lBD = CFG_K1_THREADS_PER_BLOCK;

      two_body_forces_cuda_k1_select<real, CFG_K1_THREADS_PER_BLOCK>
          (lK1_GridDims, lSharedMemorySize, lUnroll_D, lIATM_DevBeginC, lOUT_DEV, lOUT_DEV_Red, lLIMIT_MaxInH);

      CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStop, 0));
    }


    /* Execute the host portion of the workload while K1 is ongoing; this will be the
     * iterations (if any) following those assigned to the device:
     */
    if (lUnroll_H_WithK1 > 0) {
      wrapper_f_two_body_forces_cuda_helper(&lIATM_HstBeginC, &lUnroll_H_WithK1,
                                            &sHD.mKEYFCE, sHD.mMET.mRHO, aSAFE,
                                            &sHD.mIMCON, &sCD.mEWR.mRCUT, &sHD.mVDW.mRVDW, &sCD.mEWR.mALPHA,
                                            &sCD.mEWR.mEPSQ, aSTRS,
                                            aENGMET, aVIRMET, aENGVDW, aVIRVDW, aENGCPE_RL, aVIRCPE_RL);
      lIATM_HstBeginC += lUnroll_H_WithK1;

      struct timeval lK1_TV_Host_End;
      gettimeofday(&lK1_TV_Host_End, NULL);
      lDurationOfK1SpentInHost += secsfromtimeval(lK1_TV_Host_End)-secsfromtimeval(lK1_TV_Begin);
    }

    if (lUnroll_D>0) {
      CUDA_SAFE_CALL(cudaThreadSynchronize());
      CUDA_SAFE_CALL(cudaEventSynchronize(lCE_K1_DevStop));
      float lCE_K1_ElapsedTime;
      CUDA_SAFE_CALL(cudaEventElapsedTime(&lCE_K1_ElapsedTime, lCE_K1_DevStart, lCE_K1_DevStop));
      CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStart));
      CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStop));

      cudaError_t lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);

      struct timeval lK1_TV_Device_End;
      gettimeofday(&lK1_TV_Device_End, NULL);
      lDurationOfK1SpentInDevice += ((double) lCE_K1_ElapsedTime)/1000.0;
    }
    stop_timing_two_body_forces_cuda_k1();



    struct timeval lK2_TV_Begin;
    gettimeofday(&lK2_TV_Begin, NULL);
    cudaEvent_t lCE_K2_DevStart, lCE_K2_DevStop;
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K2_DevStart));
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K2_DevStop));

    start_timing_two_body_forces_cuda_k2();
    if (lUnroll_D>0) {
      CUDA_SAFE_CALL(cudaEventRecord(lCE_K2_DevStart, 0));
      /* K2: update the f{xx, yy, zz} vectors in the device.
       */
      two_body_forces_cuda_k2_bx_switch<real>(lIATM_DevBeginC, lUnroll_D, lOUT_DEV, lOUT_StripeSizeInItems);
      CUDA_SAFE_CALL(cudaEventRecord(lCE_K2_DevStop, 0));

    }

    if (lUnroll_H_WithRD > 0) {
      wrapper_f_two_body_forces_cuda_helper(&lIATM_HstBeginC, &lUnroll_H_WithRD,
                                            &sHD.mKEYFCE, sHD.mMET.mRHO, aSAFE,
                                            &sHD.mIMCON, &sCD.mEWR.mRCUT, &sHD.mVDW.mRVDW, &sCD.mEWR.mALPHA,
                                            &sCD.mEWR.mEPSQ, aSTRS,
                                            aENGMET, aVIRMET, aENGVDW, aVIRVDW, aENGCPE_RL, aVIRCPE_RL);

      lIATM_HstBeginC += lUnroll_H_WithRD;
      struct timeval lK2_TV_Host_End;
      gettimeofday(&lK2_TV_Host_End, NULL);
      lDurationOfK2SpentInHost += secsfromtimeval(lK2_TV_Host_End)-secsfromtimeval(lK2_TV_Begin);
    }


    if (lUnroll_D > 0) {
      CUDA_SAFE_CALL(cudaEventSynchronize(lCE_K2_DevStop));
      float lCE_K2_ElapsedTime;
      CUDA_SAFE_CALL(cudaEventElapsedTime(&lCE_K2_ElapsedTime, lCE_K2_DevStart, lCE_K2_DevStop));
      CUDA_SAFE_CALL(cudaEventDestroy(lCE_K2_DevStart));
      CUDA_SAFE_CALL(cudaEventDestroy(lCE_K2_DevStop));

      CUDA_SAFE_CALL(cudaThreadSynchronize());
      cudaError_t lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);
      lDurationOfK2SpentInDevice += ((double) lCE_K2_ElapsedTime)/1000.0;
    }

    stop_timing_two_body_forces_cuda_k2();

    if (lUnroll_D > 0) {
      start_timing_two_body_forces_cuda_k3();
      /* The length of the y-th dim in lK3_GridDims is the number
       * of individual reductions that must be performed (6 for
       * the stress tensor, 2 for the metalic forces, 2 for the
       * vdw forces, and 2 for the keyfce -- "2" is becuase of
       * the virial and energy).
       */
      dim3 lK3_GridDims(1, 12, 1);
      two_body_forces_cuda_k3<real,512><<<lK3_GridDims, 512, 512*1*sizeof(real)>>>(lK1_GridDims.y, lOUT_DEV_Red, lOUT_DEV_FinalRed);

      CUDA_SAFE_CALL(cudaThreadSynchronize());
      cudaError_t lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);
      stop_timing_two_body_forces_cuda_k3();

      start_timing_two_body_forces_cuda_read();
      CUDA_SAFE_CALL(cudaMemcpy(lOUT_HST_FinalRed, lOUT_DEV_FinalRed, 12*sizeof(real), cudaMemcpyDeviceToHost));
      stop_timing_two_body_forces_cuda_read();

      start_timing_two_body_forces_cuda_finalise();
      aSTRS[0] += lOUT_HST_FinalRed[0]; // stress(1) = stress(1) + strs1
      aSTRS[1] += lOUT_HST_FinalRed[1]; // stress(2) = stress(2) + strs2
      aSTRS[2] += lOUT_HST_FinalRed[2]; // stress(3) = stress(3) + strs3
      aSTRS[3] += lOUT_HST_FinalRed[1]; // stress(4) = stress(4) + strs2
      aSTRS[4] += lOUT_HST_FinalRed[3]; // stress(5) = stress(5) + strs5
      aSTRS[5] += lOUT_HST_FinalRed[4]; // stress(6) = stress(6) + strs6
      aSTRS[6] += lOUT_HST_FinalRed[2]; // stress(7) = stress(7) + strs3
      aSTRS[7] += lOUT_HST_FinalRed[4]; // stress(8) = stress(8) + strs6
      aSTRS[8] += lOUT_HST_FinalRed[5]; // stress(9) = stress(9) + strs9
      if (sHD.mNTPMET > 0) {
        *aENGMET    += lOUT_HST_FinalRed[ 6];
        *aVIRMET    += lOUT_HST_FinalRed[ 7];
      }
      if (sHD.mNTPVDW > 0) {
        *aENGVDW    += lOUT_HST_FinalRed[ 8];
        *aVIRVDW    += lOUT_HST_FinalRed[ 9];
      }
      if (sHD.mKEYFCE != 0 && (sHD.mKEYFCE % 2)==0 && sHD.mKEYFCE<=10) {
        *aENGCPE_RL += lOUT_HST_FinalRed[10];
        *aVIRCPE_RL += lOUT_HST_FinalRed[11];
      }
      stop_timing_two_body_forces_cuda_finalise();

      lIATM_DevBeginC += lUnroll_D;
    }
  }
  if (lUnroll_Dp>0) {

    start_timing_two_body_forces_cuda_k2_reduce();

#define RUN_K3(USLICES)							\
    case USLICES: {							\
      two_body_forces_cuda_k2_reduce<real,CFG_K2_RED_GRID_DIMS_Y,CFG_K2_RED_BLOCK_DIMS_X,USLICES> \
          <<<CFG_K2_RED_GRID_DIMS,CFG_K2_RED_BLOCK_DIMS>>>(0, sHD.mDEV_Spad0); \
      break;\
    }

    switch (sHD.mK2_Unroll_Slices) {
      RUN_K3(1);
      RUN_K3(2);
      RUN_K3(3);
      RUN_K3(4);
      RUN_K3(5);
      RUN_K3(6);
      RUN_K3(7);
      RUN_K3(8);
      RUN_K3(10);
      RUN_K3(12);
      RUN_K3(14);

      default: {
        printf("%s::%s: cannot handle %d numberof k2 slices.\n",
               __FILE__, __FUNCTION__, sHD.mK2_Unroll_Slices);
        exit(-1);
      }
    }
    CUDA_SAFE_CALL(cudaThreadSynchronize());
    stop_timing_two_body_forces_cuda_k2_reduce();

    start_timing_two_body_forces_cuda_read();
    // commit the accummulated updates
    CUDA_SAFE_CALL(cudaMemcpy(sHD.mHST_Spad0, sHD.mDEV_Spad0, sCD.mMXATMS*sizeof(real)*3, cudaMemcpyDeviceToHost));
    for (int lF=0 ; lF<3 ; lF++) {
      for (int lI=0 ; lI<sCD.mMXATMS ; lI++) {
        lF_Pointers[lF][lI] += sHD.mHST_Spad0[lF*sCD.mMXATMS + lI];
      }
    }
    stop_timing_two_body_forces_cuda_read();
  }


  if (sHD.mDynamicRatios) {
    /* Assume that the device and host need t(k1_d) and t(k1_h) units of time
     * per iteration (i-loop). If the total k1 time, k1(tot) is given
     * by k1(tot) = max{ p(k1_d)*I*t(k1_d), (1-p(k1_d))*I*t(k1_h) },
     *
     * Given specimen  t_r(k1_d) and t_r(k1_h), e.g. as observed at runtime
     * for a previously set p(k1_d)), this maximum is minimised when
     * both operands approximate each other, i.e. when their diff tends
     * towards zero; for a p'(k1_d)>0, this is obtained when:
     *   p'(k1_d) = t_r(k1_h)/(t_r(k1_d) + t_r(k1_h))
     */

    double lPI_Dev_K1 = lDurationOfK1SpentInDevice / ((double) lUnroll_Dp/2.0);
    double lPI_Hst_K1 = lDurationOfK1SpentInHost   /
                        ((double) sHD.mPercentageOfHostIterationsToBeOverlappedWithK1*lUnroll_Hp);
    double lPI_Dev_K2 = lDurationOfK2SpentInDevice / ((double) lUnroll_Dp/2.0);
    double lPI_Hst_K2 = lDurationOfK2SpentInHost   /
                        ((double) (lUnroll_Hp - sHD.mPercentageOfHostIterationsToBeOverlappedWithK1*lUnroll_Hp));
    double lPS_Dev_K1 = lPI_Hst_K1 / (lPI_Dev_K1 + lPI_Hst_K1);
    double lPS_Dev_K2 = lPI_Hst_K2 / (lPI_Dev_K2 + lPI_Hst_K2);

    /* Hence, the host will need access to (lPS_Dev_K1 + lPS_Dev_K2)/2 of the
     * iterations; let r = 1 - (lPS_Dev_K1 + lPS_Dev_K2)/2 be the remaining,
     * i.e. host percentage, of that, it needs to be distributed as follows:
     * p(k1,h) = 1-lPS_Dev_K1/2 and (1-lPS_Dev_K2)/2
     */
    double lDevPerc = (lPS_Dev_K1 + lPS_Dev_K2)/2.0;
    //    double lHstK1Perc = (1.0 - lPS_Dev_K1)/(2.0 - lPS_Dev_K1 - lPS_Dev_K2);
    double lHstK1Perc = ((1.0 - lPS_Dev_K1)/(2.0 - lPS_Dev_K1 - lPS_Dev_K2));

    sHD.mPercentageOfIterationsOffloadedToTheDevice = lDevPerc;
    sHD.mPercentageOfHostIterationsToBeOverlappedWithK1 = lHstK1Perc;
  }
}


