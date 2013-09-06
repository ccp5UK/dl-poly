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
#include "cutil.h"

#include "dl_poly_common_cu.cu"

#define CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE (CFG_OVERLAP_WITH_HOST ? 0.8 : 1.0)
#define CFG_K1_DEBUG  0

/* Env var: dlpolycuda_spme_forces_dynamic_ratios
 */
#define CFG_DYNAMIC_RATIOS (CFG_OVERLAP_WITH_HOST)

/* Grid-level block capping.
 *   Env Var: dlpolycuda_spme_forces_grid_dim_y
 */
#define CFG_K1_GRID_DIM_Y                  900

#define CFG_ENABLE_TEXTURE_CACHING         1

extern "C" double spme_forces_cuda_percentage_of_iterations_offloaded_to_the_device() {
  if (dl_poly_cuda_offload_spme_forces()==dl_poly_cuda_fortran_false())
    return(0.0);
  return (CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE);
}

template<typename T_> struct constant_data {
  int  mNATMS, mMXSPL, mMXATMS, mKMAXA, mKMAXB, mKMAXC, mL_CP;
  T_   mRCELL[9];
  T_ mTMP;
  T_  *mBSDX, *mBSDY, *mBSDZ, *mBSPX, *mBSPY, *mBSPZ;
  T_  *mBSPDZYX;
  int *mIXXYYZZ;
  T_  *mQQQ;
  T_  *mCHGE;
  int3 mIXYZDB, mIXYZT, mQQQ_DimLengths;
};

__device__ __constant__ constant_data<real> CONSTANT_DATA;
static                  constant_data<real> sCD;

#if CFG_ENABLE_TEXTURE_CACHING
texture<real_tex, 1, cudaReadModeElementType> TEXTURE_DATA_QQQ;
#endif

template<typename T_> struct host_data {
  int mDynamicRatios;
  int mK1_Grid_Dim_Y;
  int mAreKMAXABCIdentical;

  /* items needed for the herlper's invocation: */
  T_    *mRCELL, *mBSDX, *mBSDY, *mBSDZ, *mBSPX, *mBSPY, *mBSPZ;
  int   *mIXX, *mIYY, *mIZZ;
  T_    *mQQC_DOMAIN;
  int3   mIXYZB;

  double mPercentageOfIterationsOffloadedToTheDevice;
  int    mUnroll;
  T_    *mOUT_HST;

  int   mAreBSPDDataOnline;
};
static                  host_data<real> sHD;

extern "C" void spme_forces_cuda_initialise(
  int *aNATMS, int *aMXSPL, int *aMXATMS, int *aKMAXA, int *aKMAXB, int *aKMAXC,
  int *aL_CP, real *aTMP, real *aRCELL, real *aQQQ, real *aCHGE,
  real *aBSDX, real *aBSDY, real *aBSDZ,
  real *aBSPX, real *aBSPY, real *aBSPZ,
  int *aIXX, int *aIYY, int *aIZZ,
  int *aIXDB, int *aIXT, int *aIYDB, int *aIYT, int *aIZDB, int *aIZT){

  sHD.mDynamicRatios = dl_poly_cuda_getenv("dlpolycuda_spme_forces_dynamic_ratios",
                                           CFG_DYNAMIC_RATIOS);
  sHD.mK1_Grid_Dim_Y = dl_poly_cuda_getenv("dlpolycuda_spme_forces_grid_dim_y",
                                            CFG_K1_GRID_DIM_Y);

  sHD.mUnroll = *aNATMS;

  static int sSetRatios = 1;
  if (sSetRatios) {
    sHD.mPercentageOfIterationsOffloadedToTheDevice = CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE;
    if (sHD.mDynamicRatios) {
      sSetRatios = 0;
    }
  }

  //  sHD.mPercentageOfIterationsOffloadedToTheDevice = CFG_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE;
  sHD.mRCELL = aRCELL;
  sHD.mBSPX  = aBSPX;
  sHD.mBSPY  = aBSPY;
  sHD.mBSPZ  = aBSPZ;
  sHD.mBSDX  = aBSDX;
  sHD.mBSDY  = aBSDY;
  sHD.mBSDZ  = aBSDZ;
  sHD.mIXX   = aIXX;
  sHD.mIYY   = aIYY;
  sHD.mIZZ   = aIZZ;
  sHD.mQQC_DOMAIN = aQQQ;

  sHD.mAreBSPDDataOnline =
      (*spme_container_cuda_bspgen_leave_bspdxyz_data_online())==dl_poly_cuda_fortran_true();

  /* Allocate the host buffer only once as it is large are repetitive costs are too
   * expensive.
   * XXX: We _do_not_ currently ever deallocate the buffer -- this should normally
   *      happen at process destruction time.
   */
  static int sIsFirstInvocation = 1;
  static void *sOUT_HST = NULL;
  if (sIsFirstInvocation) {
    CUDA_SAFE_CALL(cudaHostAlloc(&sOUT_HST, 4*(*aMXATMS)*sizeof(real), 0));
    sIsFirstInvocation = 0;
  }
  sHD.mOUT_HST = (real*)sOUT_HST;

  sCD.mIXYZDB  = make_int3(*aIXDB, *aIYDB, *aIZDB);
  sCD.mIXYZT   = make_int3(*aIXT, *aIYT, *aIZT);
  sCD.mQQQ_DimLengths = sCD.mIXYZT - sCD.mIXYZDB + 1;
  sCD.mNATMS   = *aNATMS;
  sCD.mMXSPL   = *aMXSPL;

  if (sCD.mMXSPL!=8) {
    printf("%s::%s: stub: mxspl(=%d)!=8\n", __FILE__, __FUNCTION__, sCD.mMXSPL);
    exit(-1);
  }

  sCD.mMXATMS  = *aMXATMS;
  sCD.mKMAXA   = *aKMAXA;
  sCD.mKMAXB   = *aKMAXB;
  sCD.mKMAXC   = *aKMAXC;
  sHD.mAreKMAXABCIdentical = (sCD.mKMAXA==sCD.mKMAXB) && (sCD.mKMAXB==sCD.mKMAXC);

  sCD.mTMP     = *aTMP;
  sCD.mL_CP    = *aL_CP;
  for (int lI=0 ; lI<9 ; lI++) {
      sCD.mRCELL[lI] = aRCELL[lI];
  }

  sHD.mIXYZB = sCD.mIXYZDB + sCD.mMXSPL;

  if (!sHD.mAreBSPDDataOnline) {
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mBSDX, sCD.mMXSPL*sCD.mNATMS*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mBSDY, sCD.mMXSPL*sCD.mNATMS*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mBSDZ, sCD.mMXSPL*sCD.mNATMS*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mBSPX, sCD.mMXSPL*sCD.mNATMS*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mBSPY, sCD.mMXSPL*sCD.mNATMS*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mBSPZ, sCD.mMXSPL*sCD.mNATMS*sizeof(real)));
  } else {
    spme_container_cuda_bspgen_grab_bspdxyz_data((void**)&sCD.mBSPX,(void**)&sCD.mBSPY,(void**)&sCD.mBSPZ,
                                                 (void**)&sCD.mBSDX,(void**)&sCD.mBSDY,(void**)&sCD.mBSDZ);
  }
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mCHGE,      sCD.mNATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mIXXYYZZ, 3*sCD.mNATMS*sizeof(int)));
  int  lQQQ_Length  = sCD.mQQQ_DimLengths.x * sCD.mQQQ_DimLengths.y * sCD.mQQQ_DimLengths.z;

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mQQQ, lQQQ_Length*sizeof(real)));
#if CFG_ENABLE_TEXTURE_CACHING
  BIND_TEXTURE_REAL_1D(TEXTURE_DATA_QQQ, sCD.mQQQ, lQQQ_Length);
#endif

  start_timing_spme_forces_write();
  if (!sHD.mAreBSPDDataOnline) {
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mBSDX, aBSDX, sCD.mMXSPL*sCD.mNATMS*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mBSDY, aBSDY, sCD.mMXSPL*sCD.mNATMS*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mBSDZ, aBSDZ, sCD.mMXSPL*sCD.mNATMS*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mBSPX, aBSPX, sCD.mMXSPL*sCD.mNATMS*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mBSPY, aBSPY, sCD.mMXSPL*sCD.mNATMS*sizeof(real), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mBSPZ, aBSPZ, sCD.mMXSPL*sCD.mNATMS*sizeof(real), cudaMemcpyHostToDevice));
  }
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mCHGE, aCHGE, sCD.mNATMS*sizeof(real), cudaMemcpyHostToDevice));

  CUDA_SAFE_CALL(cudaMemcpy(sCD.mIXXYYZZ+0*sCD.mNATMS, aIXX, sCD.mNATMS*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mIXXYYZZ+1*sCD.mNATMS, aIYY, sCD.mNATMS*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mIXXYYZZ+2*sCD.mNATMS, aIZZ, sCD.mNATMS*sizeof(int), cudaMemcpyHostToDevice));

  CUDA_SAFE_CALL(cudaMemcpy(sCD.mQQQ, aQQQ, lQQQ_Length*sizeof(real), cudaMemcpyHostToDevice));

  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA, (void*)&sCD, sizeof(constant_data<real>)));
  stop_timing_spme_forces_write();
}

extern "C" void spme_forces_cuda_finalise() {
  if (!sHD.mAreBSPDDataOnline) {
    CUDA_SAFE_CALL(cudaFree(sCD.mBSDX));
    CUDA_SAFE_CALL(cudaFree(sCD.mBSDY));
    CUDA_SAFE_CALL(cudaFree(sCD.mBSDZ));
    CUDA_SAFE_CALL(cudaFree(sCD.mBSPX));
    CUDA_SAFE_CALL(cudaFree(sCD.mBSPY));
    CUDA_SAFE_CALL(cudaFree(sCD.mBSPZ));
  }
  CUDA_SAFE_CALL(cudaFree(sCD.mCHGE));
  CUDA_SAFE_CALL(cudaFree(sCD.mIXXYYZZ));
  CUDA_SAFE_CALL(cudaFree(sCD.mQQQ));

#if CFG_ENABLE_TEXTURE_CACHING
  CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_QQQ));
#endif
}

__device__ real qqq(int aX, int aY, int aZ) {
  int3 lIdx = make_int3(aX, aY, aZ) - CONSTANT_DATA.mIXYZDB;

#if CFG_ENABLE_TEXTURE_CACHING
  return (fetch_real_1d(TEXTURE_DATA_QQQ,
          lIdx.z*CONSTANT_DATA.mQQQ_DimLengths.x*CONSTANT_DATA.mQQQ_DimLengths.y +
          lIdx.y*CONSTANT_DATA.mQQQ_DimLengths.x + lIdx.x));
#else
  return(*(CONSTANT_DATA.mQQQ +
           lIdx.z*CONSTANT_DATA.mQQQ_DimLengths.x*CONSTANT_DATA.mQQQ_DimLengths.y +
           lIdx.y*CONSTANT_DATA.mQQQ_DimLengths.x + lIdx.x));

#endif
}


template<typename T_, unsigned int MXSPL_>
__global__ void spme_forces_cuda_k1_debug(int aI, int aLast_I, T_ *aOUT) {

  for (int lI=aI + blockIdx.y ; lI<=aLast_I ; lI+=gridDim.y) {
    T_ lTMP = CONSTANT_DATA.mCHGE[lI-1];
    if (fabs(lTMP) == (T_)0) {
      aOUT[4*(lI-aI) + 3] = (T_)0;
      continue;
    }
    T_ lFIX=0, lFIY=0, lFIZ=0;
    T_ lFACX = CONSTANT_DATA.mTMP * (T_)CONSTANT_DATA.mKMAXA;
    T_ lFACY = CONSTANT_DATA.mTMP * (T_)CONSTANT_DATA.mKMAXB;
    T_ lFACZ = CONSTANT_DATA.mTMP * (T_)CONSTANT_DATA.mKMAXC;

    for (int lL=1 ; lL<=MXSPL_ ; lL++) {
      int lLL = CONSTANT_DATA.mIXXYYZZ[2*CONSTANT_DATA.mNATMS + lI -1] - lL + 2;
      T_ lBDXL = lTMP * lFACX * (*F2D_ADDRESS(CONSTANT_DATA.mBSPZ, 1, 1, MXSPL_, lL, lI));
      T_ lBDYL = lTMP * lFACY * (*F2D_ADDRESS(CONSTANT_DATA.mBSPZ, 1, 1, MXSPL_, lL, lI));
      T_ lBDZL = lTMP * lFACZ * (*F2D_ADDRESS(CONSTANT_DATA.mBSDZ, 1, 1, MXSPL_, lL, lI));

      for (int lK=1 ; lK<=MXSPL_ ; lK++) {
        int lKK = CONSTANT_DATA.mIXXYYZZ[1*CONSTANT_DATA.mNATMS + lI -1] - lK + 2;
        T_ lBDXK = lBDXL * (*F2D_ADDRESS(CONSTANT_DATA.mBSPY, 1, 1, MXSPL_, lK, lI));
        T_ lBDYK = lBDYL * (*F2D_ADDRESS(CONSTANT_DATA.mBSDY, 1, 1, MXSPL_, lK, lI));
        T_ lBDZK = lBDZL * (*F2D_ADDRESS(CONSTANT_DATA.mBSPY, 1, 1, MXSPL_, lK, lI));

        for (int lJ=1 ; lJ<=MXSPL_ ; lJ++) {
          int lJJ = CONSTANT_DATA.mIXXYYZZ[0*CONSTANT_DATA.mNATMS + lI -1] - lJ + 2;
          T_ lQSUM = qqq(lJJ, lKK, lLL);
          lFIX += lQSUM * lBDXK * (*F2D_ADDRESS(CONSTANT_DATA.mBSDX, 1, 1, MXSPL_, lJ, lI));
          lFIY += lQSUM * lBDYK * (*F2D_ADDRESS(CONSTANT_DATA.mBSPX, 1, 1, MXSPL_, lJ, lI));
          lFIZ += lQSUM * lBDZK * (*F2D_ADDRESS(CONSTANT_DATA.mBSPX, 1, 1, MXSPL_, lJ, lI));
        }
      }
    }
    T_ lFX = lFIX*CONSTANT_DATA.mRCELL[0] + lFIY*CONSTANT_DATA.mRCELL[1] + lFIZ*CONSTANT_DATA.mRCELL[2];
    T_ lFY = lFIX*CONSTANT_DATA.mRCELL[3] + lFIY*CONSTANT_DATA.mRCELL[4] + lFIZ*CONSTANT_DATA.mRCELL[5];
    T_ lFZ = lFIX*CONSTANT_DATA.mRCELL[6] + lFIY*CONSTANT_DATA.mRCELL[7] + lFIZ*CONSTANT_DATA.mRCELL[8];

    aOUT[4*(lI-aI) + 0] = lFX;
    aOUT[4*(lI-aI) + 1] = lFY;
    aOUT[4*(lI-aI) + 2] = lFZ;
    aOUT[4*(lI-aI) + 3] = (T_)1;
  }
}

template<typename T_, unsigned int MXSPL_, int KAREMAXABCIDENT_>
__global__ void spme_forces_cuda_k1(int aI, int aLast_I, T_ *aOUT) {
  DECLARE_DYNAMIC_SHARED(T_);

  for (int lI=aI+blockIdx.y ; lI<=aLast_I ; lI+=gridDim.y) {
    T_ lTMP = CONSTANT_DATA.mCHGE[lI-1];
    /* Check if we should skip this iteration. We also need to set
     * 'fff(0)' to zero. TODO: zero the vector during init time.
     */
    if (fabs(lTMP) == (T_)0) {
      if (threadIdx.x == 0 && threadIdx.y==0) {
      	aOUT[4*(lI-aI) + 3] = (T_)0;
      }
    } else {
      T_ lFACX = CONSTANT_DATA.mTMP * (T_)CONSTANT_DATA.mKMAXA;
      T_ lFACY = KAREMAXABCIDENT_ ? lFACX : (CONSTANT_DATA.mTMP * (T_)CONSTANT_DATA.mKMAXB);
      T_ lFACZ = KAREMAXABCIDENT_ ? lFACX : (CONSTANT_DATA.mTMP * (T_)CONSTANT_DATA.mKMAXC);

      T_ lBSPZ = *F2D_ADDRESS(CONSTANT_DATA.mBSPZ, 1, 1, MXSPL_, 1+threadIdx.x, lI);
      T_ lBSDZ = *F2D_ADDRESS(CONSTANT_DATA.mBSDZ, 1, 1, MXSPL_, 1+threadIdx.x, lI);

      T_ lBDXL = lTMP * lFACX * lBSPZ;
      T_ lBDYL = lTMP * lFACY * lBSPZ;
      T_ lBDZL = lTMP * lFACZ * lBSDZ;

      T_ lBSPY_K_I = *F2D_ADDRESS(CONSTANT_DATA.mBSPY, 1, 1, MXSPL_, 1+threadIdx.y, lI);
      T_ lBSDY_K_I = *F2D_ADDRESS(CONSTANT_DATA.mBSDY, 1, 1, MXSPL_, 1+threadIdx.y, lI);

      shared[0*MXSPL_*MXSPL_ + threadIdx.x*MXSPL_ + threadIdx.y] = lBDXL * lBSPY_K_I;
      shared[1*MXSPL_*MXSPL_ + threadIdx.x*MXSPL_ + threadIdx.y] = lBDYL * lBSDY_K_I;
      shared[2*MXSPL_*MXSPL_ + threadIdx.x*MXSPL_ + threadIdx.y] = lBDZL * lBSPY_K_I;

      __syncthreads();

      T_ lFIX=0, lFIY=0, lFIZ=0;

      int lL = 1 + threadIdx.y;
      int lLL = CONSTANT_DATA.mIXXYYZZ[2*CONSTANT_DATA.mNATMS + lI -1] - lL + 2;

      int lJ = 1 + threadIdx.x;
      int lJJ = CONSTANT_DATA.mIXXYYZZ[0*CONSTANT_DATA.mNATMS + lI -1] - lJ + 2;
      T_ lBSPX_J_I = *F2D_ADDRESS(CONSTANT_DATA.mBSPX, 1, 1, MXSPL_, lJ, lI);
      T_ lBSDX_J_I = *F2D_ADDRESS(CONSTANT_DATA.mBSDX, 1, 1, MXSPL_, lJ, lI);


#pragma unroll 2
      for (int lK=1 ; lK<=MXSPL_ ; lK++) {
        int lKK = CONSTANT_DATA.mIXXYYZZ[1*CONSTANT_DATA.mNATMS + lI -1] - lK + 2;

        T_ lBDXK = shared[0*MXSPL_*MXSPL_ + (lL-1)*MXSPL_ + lK -1];
        T_ lBDYK = shared[1*MXSPL_*MXSPL_ + (lL-1)*MXSPL_ + lK -1];
        T_ lBDZK = shared[2*MXSPL_*MXSPL_ + (lL-1)*MXSPL_ + lK -1];

        T_ lQSUM  = qqq(lJJ, lKK, lLL);
        T_ lBDXJ  = lQSUM * lBDXK * lBSDX_J_I;
        T_ lBDYJ  = lQSUM * lBDYK * lBSPX_J_I;
        T_ lBDZJ  = lQSUM * lBDZK * lBSPX_J_I;

        lFIX += lBDXJ;
        lFIY += lBDYJ;
        lFIZ += lBDZJ;
      }
      __syncthreads();
      shared[0*MXSPL_*MXSPL_ + threadIdx.x + threadIdx.y*MXSPL_] = lFIX;
      shared[1*MXSPL_*MXSPL_ + threadIdx.x + threadIdx.y*MXSPL_] = lFIY;
      shared[2*MXSPL_*MXSPL_ + threadIdx.x + threadIdx.y*MXSPL_] = lFIZ;
      __syncthreads();

#if !CFG_UNIFIED_ADDRESS_SPACE
      psum2D<T_, MXSPL_, MXSPL_, 3>();
#else
      psum2D_uas<T_, MXSPL_, MXSPL_, 3>(shared);
#endif
      __syncthreads();

      if (MXSPL_>=4) {
        if (threadIdx.x<3 && threadIdx.y==0) {
          lFIX = shared[0*(MXSPL_*MXSPL_)];
          lFIY = shared[1*(MXSPL_*MXSPL_)];
          lFIZ = shared[2*(MXSPL_*MXSPL_)];
          T_ lFXYZ =
              lFIX*CONSTANT_DATA.mRCELL[threadIdx.x*3+0] +
              lFIY*CONSTANT_DATA.mRCELL[threadIdx.x*3+1] +
              lFIZ*CONSTANT_DATA.mRCELL[threadIdx.x*3+2];

          shared[threadIdx.x] = lFXYZ;
        }
        if (threadIdx.x==3) {
          shared[3] = (T_)1;
        }
        if (threadIdx.x<4 && threadIdx.y==0) {
          aOUT[4*(lI-aI) + threadIdx.x] = shared[threadIdx.x];
        }
      } else {

	/* This is a less optimised version of the fragment above. I've kept it
	 * for referencing purposes.
	 */
        if (threadIdx.x==0 && threadIdx.y==0) {
          lFIX = shared[0*(MXSPL_*MXSPL_)];
          lFIY = shared[1*(MXSPL_*MXSPL_)];
          lFIZ = shared[2*(MXSPL_*MXSPL_)];

          T_ lFX = lFIX*CONSTANT_DATA.mRCELL[0] + lFIY*CONSTANT_DATA.mRCELL[1] + lFIZ*CONSTANT_DATA.mRCELL[2];
          T_ lFY = lFIX*CONSTANT_DATA.mRCELL[3] + lFIY*CONSTANT_DATA.mRCELL[4] + lFIZ*CONSTANT_DATA.mRCELL[5];
          T_ lFZ = lFIX*CONSTANT_DATA.mRCELL[6] + lFIY*CONSTANT_DATA.mRCELL[7] + lFIZ*CONSTANT_DATA.mRCELL[8];

          aOUT[4*(lI-aI) + 0] = lFX;
          aOUT[4*(lI-aI) + 1] = lFY;
          aOUT[4*(lI-aI) + 2] = lFZ;
          aOUT[4*(lI-aI) + 3] = (T_)1;
        }
      }
    }
    __syncthreads();
  }
}

extern "C" void spme_forces_cuda_invoke(real *aFFF,
                                        real *aFXX, real *aFYY, real *aFZZ,
                                        real *aFCX, real *aFCY, real *aFCZ) {


  int lI = 1;
  int lUnroll = min(sHD.mUnroll, sCD.mNATMS - lI + 1);
  if (lUnroll<=0)
    return;

  int lUnroll_D = (int) (sHD.mPercentageOfIterationsOffloadedToTheDevice * (double) lUnroll);
  int lUnroll_H = lUnroll - lUnroll_D;
  real *lOUT_DEV, *lOUT_HST;

  float lCE_K1_Hst_ElapsedTime = 0.0f, lCE_K1_Dev_ElapsedTime = 0.0f;
  cudaEvent_t lCE_K1_DevStart, lCE_K1_DevStop;

  if (lUnroll_D>0) {
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStart));
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStop));

    lOUT_HST = sHD.mOUT_HST;
    CUDA_SAFE_CALL(cudaMalloc(&lOUT_DEV, 4*lUnroll_D*sizeof(real)));
    start_timing_spme_forces_k1();

    CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStart, 0));

#if (CFG_K1_DEBUG)
    spme_forces_cuda_k1_debug<real, 8><<<dim3(1,4096,1), dim3(1,1,1)>>>(lI, lI+lUnroll_D-1, lOUT_DEV);
#else
    dim3 lK1_GridDims(1, min(sHD.mK1_Grid_Dim_Y,lUnroll_D), 1);

    if (sHD.mAreKMAXABCIdentical) {
      spme_forces_cuda_k1<real, 8, 1><<<lK1_GridDims,dim3(8,8,1),(3*8*8)*sizeof(real)>>>(lI, lI+lUnroll_D-1, lOUT_DEV);
    } else {
      spme_forces_cuda_k1<real, 8, 0><<<lK1_GridDims,dim3(8,8,1),(3*8*8)*sizeof(real)>>>(lI, lI+lUnroll_D-1, lOUT_DEV);
    }
#endif
    CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStop, 0));
  }

  if (lUnroll_H>0) {
    struct timeval lTV_K1_HstStart, lTV_K1_HstStop;
    gettimeofday(&lTV_K1_HstStart, NULL);

    int lI_H = lI+lUnroll_D;
    wrapper_f_spme_forces_helper(&lI_H, &lUnroll_H, sHD.mRCELL, &sCD.mTMP, aFFF,
                                 &sCD.mIXYZDB.x, &sCD.mIXYZDB.y, &sCD.mIXYZDB.z,
                                 sHD.mIXX, sHD.mIYY, sHD.mIZZ,
                                 sHD.mBSPX, sHD.mBSPY, sHD.mBSPZ,
                                 sHD.mBSDX, sHD.mBSDY, sHD.mBSDZ,
                                 sHD.mQQC_DOMAIN,
                                 &sCD.mIXYZT.x, &sCD.mIXYZT.y, &sCD.mIXYZT.z);

    gettimeofday(&lTV_K1_HstStop, NULL);
    lCE_K1_Hst_ElapsedTime = secsfromtimeval(lTV_K1_HstStop) - secsfromtimeval(lTV_K1_HstStart);
  }

  if (lUnroll_D>0) {
    CUDA_SAFE_CALL(cudaEventSynchronize(lCE_K1_DevStop));
    CUDA_SAFE_CALL(cudaEventElapsedTime(&lCE_K1_Dev_ElapsedTime, lCE_K1_DevStart, lCE_K1_DevStop));
    lCE_K1_Dev_ElapsedTime /= 1000.0f;
    CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStart));
    CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStop));

    CUDA_SAFE_CALL(cudaThreadSynchronize());
    cudaError_t lLastError = cudaGetLastError();
    CUT_CHECK_ERROR(lLastError);
    stop_timing_spme_forces_k1();

    start_timing_spme_forces_read();
    CUDA_SAFE_CALL(cudaMemcpy(lOUT_HST, lOUT_DEV, 4*lUnroll_D*sizeof(real), cudaMemcpyDeviceToHost));
    stop_timing_spme_forces_read();

    start_timing_spme_forces_finalise();

    for (int lU=0 ; lU<lUnroll_D ; lU++) {

      real *lF_HST = &lOUT_HST[lU*4];

      if ((int) lF_HST[3]) {

        aFXX[lI+lU-1] += lF_HST[0];
        aFYY[lI+lU-1] += lF_HST[1];
        aFZZ[lI+lU-1] += lF_HST[2];

        aFFF[0] += (real)1;
        aFFF[1] += lF_HST[0];
        aFFF[2] += lF_HST[1];
        aFFF[3] += lF_HST[2];

        if (sCD.mL_CP) {
          aFCX[lI+lU-1] += lF_HST[0];
          aFCY[lI+lU-1] += lF_HST[1];
          aFCZ[lI+lU-1] += lF_HST[2];
        }
      }
    }

    stop_timing_spme_forces_finalise();

    CUDA_SAFE_CALL(cudaFree(lOUT_DEV));
  }

  if (sHD.mDynamicRatios) {
    double lPI_Dev = lCE_K1_Dev_ElapsedTime / ((double) lUnroll_D);
    double lPI_Hst = lCE_K1_Hst_ElapsedTime / ((double) lUnroll_H);
    double lNewRatio = lPI_Hst / (lPI_Dev + lPI_Hst);

    sHD.mPercentageOfIterationsOffloadedToTheDevice = lNewRatio;
  }
}
