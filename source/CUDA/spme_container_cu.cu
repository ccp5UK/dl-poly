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

#include "dl_poly_common_cu.cu"

template<typename T_> struct constant_data {
  /* Use a union here as more functions may be moved into this file.
   */
  struct {
    struct bspgen {
      int mNLAST, mMXSPL, mNOSPL, mUnroll;
      T_ *mBSPX, *mBSPY, *mBSPZ;
      T_ *mBSDX, *mBSDY, *mBSDZ;
      T_ *mBSPDXYZ;
      T_ *mXXX, *mYYY, *mZZZ;
      T_ *mXXXYYYZZZ;
    } mBSPGEN;
  } mAll;
};

template<typename T_> struct host_data {
  struct {
    struct bspgen {
      int2 mHostRange, mDeviceRange;
      int  mIsHostInvolved, mIsDeviceInvolved;

      double mPercentageOfIterationsOffloadedToTheDevice;
      T_ *mXXX, *mYYY, *mZZZ;
      int mFirstDevIter;
      int mNATMS; // -- spme_forces-specific
    } mBSPGEN;
  } mAll;
};


__device__ __constant__ constant_data<real> CONSTANT_DATA;
#define CONSTANT_DATA_BSPGEN CONSTANT_DATA.mAll.mBSPGEN
static                  constant_data<real> sCD;
static                  host_data<real> sHD;

extern "C" int* spme_container_cuda_bspgen_leave_bspdxyz_data_online() {
  static int lTruth=dl_poly_cuda_fortran_false();
  return (&lTruth);
}

extern "C" int* spme_container_cuda_bspgen_is_under_ewald_spme_forces() {
  static int lTruth=dl_poly_cuda_fortran_false();
  return (&lTruth);
}

extern "C" void spme_container_cuda_bspgen_set_leave_bspdxyz_data_online(int *aTruth) {
  *spme_container_cuda_bspgen_leave_bspdxyz_data_online() = *aTruth;
}

extern "C" void spme_container_cuda_bspgen_set_is_under_ewald_spme_forces(int *aTruth, int *aNATMS) {
  *spme_container_cuda_bspgen_is_under_ewald_spme_forces() = *aTruth;
  sHD.mAll.mBSPGEN.mNATMS=*aNATMS;
}

extern "C" void spme_container_cuda_bspgen_kill_bspdxyz_data() {
  // sanity check:
  if (*spme_container_cuda_bspgen_leave_bspdxyz_data_online()==dl_poly_cuda_fortran_false()) {
    printf("%s::%s:   all data are offline.\n", __FILE__, __FUNCTION__);
    exit(-1);
  }
  CUDA_SAFE_CALL(cudaFree(sCD.mAll.mBSPGEN.mBSPDXYZ));
  // reset to avoid memory leaks:
  *spme_container_cuda_bspgen_leave_bspdxyz_data_online() = dl_poly_cuda_fortran_false();
}

extern "C" void spme_container_cuda_bspgen_grab_bspdxyz_data(void** aBSPX, void ** aBSPY, void **aBSPZ,
							     void** aBSDX, void ** aBSDY, void **aBSDZ) {
  // sanity check:
  if (*spme_container_cuda_bspgen_leave_bspdxyz_data_online()==dl_poly_cuda_fortran_false()) {
    printf("%s::%s:   all data are offline.\n", __FILE__, __FUNCTION__);
    exit(-1);
  }
  *aBSPX = sCD.mAll.mBSPGEN.mBSPX;
  *aBSPY = sCD.mAll.mBSPGEN.mBSPY;
  *aBSPZ = sCD.mAll.mBSPGEN.mBSPZ;
  *aBSDX = sCD.mAll.mBSPGEN.mBSDX;
  *aBSDY = sCD.mAll.mBSPGEN.mBSDY;
  *aBSDZ = sCD.mAll.mBSPGEN.mBSDZ;
}


extern "C" void spme_container_cuda_bspgen_initialise
                (int *aUnroll, int *mNATMS, int *aMXSPL, int *aNOSPL,
                 real *aXXX, real *aYYY, real *aZZZ){

  sCD.mAll.mBSPGEN.mUnroll = *aUnroll;
  sCD.mAll.mBSPGEN.mNLAST  = *mNATMS;
  sCD.mAll.mBSPGEN.mMXSPL  = *aMXSPL;
  sCD.mAll.mBSPGEN.mNOSPL  = *aNOSPL;

  sHD.mAll.mBSPGEN.mPercentageOfIterationsOffloadedToTheDevice = 0.1;
  sHD.mAll.mBSPGEN.mXXX = aXXX;
  sHD.mAll.mBSPGEN.mYYY = aYYY;
  sHD.mAll.mBSPGEN.mZZZ = aZZZ;
  sHD.mAll.mBSPGEN.mFirstDevIter =
      1 + (int) ((1.0 - sHD.mAll.mBSPGEN.mPercentageOfIterationsOffloadedToTheDevice) *
                 (double) sCD.mAll.mBSPGEN.mNLAST);


  CUDA_SAFE_CALL(cudaMalloc(&sCD.mAll.mBSPGEN.mBSPDXYZ,   6*sCD.mAll.mBSPGEN.mMXSPL*sCD.mAll.mBSPGEN.mNLAST*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mAll.mBSPGEN.mXXXYYYZZZ, 3*sCD.mAll.mBSPGEN.mNLAST*sizeof(real)));
  sCD.mAll.mBSPGEN.mBSPX = sCD.mAll.mBSPGEN.mBSPDXYZ + 0*sCD.mAll.mBSPGEN.mMXSPL*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mBSPY = sCD.mAll.mBSPGEN.mBSPDXYZ + 1*sCD.mAll.mBSPGEN.mMXSPL*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mBSPZ = sCD.mAll.mBSPGEN.mBSPDXYZ + 2*sCD.mAll.mBSPGEN.mMXSPL*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mBSDX = sCD.mAll.mBSPGEN.mBSPDXYZ + 3*sCD.mAll.mBSPGEN.mMXSPL*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mBSDY = sCD.mAll.mBSPGEN.mBSPDXYZ + 4*sCD.mAll.mBSPGEN.mMXSPL*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mBSDZ = sCD.mAll.mBSPGEN.mBSPDXYZ + 5*sCD.mAll.mBSPGEN.mMXSPL*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mXXX  = sCD.mAll.mBSPGEN.mXXXYYYZZZ + 0*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mYYY  = sCD.mAll.mBSPGEN.mXXXYYYZZZ + 1*sCD.mAll.mBSPGEN.mNLAST;
  sCD.mAll.mBSPGEN.mZZZ  = sCD.mAll.mBSPGEN.mXXXYYYZZZ + 2*sCD.mAll.mBSPGEN.mNLAST;

  CUDA_SAFE_CALL(cudaMemcpy(sCD.mAll.mBSPGEN.mXXX, aXXX, sCD.mAll.mBSPGEN.mNLAST*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mAll.mBSPGEN.mYYY, aYYY, sCD.mAll.mBSPGEN.mNLAST*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mAll.mBSPGEN.mZZZ, aZZZ, sCD.mAll.mBSPGEN.mNLAST*sizeof(real), cudaMemcpyHostToDevice));

  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA, (void*)&sCD, sizeof(constant_data<real>)));
}

extern "C" void spme_container_cuda_bspgen_finalise() {
  // remove the data if we haven't been instructed otherwise:
  if (*spme_container_cuda_bspgen_leave_bspdxyz_data_online()==dl_poly_cuda_fortran_false()) {
    CUDA_SAFE_CALL(cudaFree(sCD.mAll.mBSPGEN.mBSPDXYZ));
  }
  CUDA_SAFE_CALL(cudaFree(sCD.mAll.mBSPGEN.mXXXYYYZZZ));
}

#define bspgen_bspxyz(X,Y) shared[threadIdx.y*2*NOSPL_ + (X)-1]
#define bspgen_bsdxyz(X,Y) shared[threadIdx.y*2*NOSPL_ + (X)-1+NOSPL_]

template<typename T_, int NOSPL_>
__global__ void spme_container_cuda_bspgen_k1_NOSPLeq8(int aI, int aIEnd) {

  extern __shared__ T_ shared[];

  int lI        = aI + blockIdx.y;

  while (lI<=aIEnd) {

    T_ lXXXYYYZZZ = CONSTANT_DATA_BSPGEN.mXXXYYYZZZ[threadIdx.y*CONSTANT_DATA_BSPGEN.mNLAST + lI - 1];
    T_ lRIXYZ0    = lXXXYYYZZZ - aint(lXXXYYYZZZ);

    if (threadIdx.x==0) {
      bspgen_bspxyz(1, lI) = lRIXYZ0;
      bspgen_bspxyz(2, lI) = ((T_) 1) - bspgen_bspxyz(1, lI);
      bspgen_bsdxyz(1, lI) = (T_)1;
      bspgen_bsdxyz(2, lI) = (T_)-1;
    }
    __syncthreads();

#pragma unroll 5
    for (int lK=3 ; lK<=NOSPL_-1 ; lK++) {

      if ((1+threadIdx.x)==lK) {
        bspgen_bspxyz(lK,lI) = (T_)0;
      }

      T_ lK_R = (T_)lK;
      T_ lKM1_RR = ((T_) 1) / ((T_) (lK-1));

      int lJ = 1+threadIdx.x;

      T_ lBSPXYZ_Ju;
      if (lJ>=2 && lJ<=lK) {
        T_ lJM1_R = (T_) (lJ-1);
        T_ lAAABBBCCC = lRIXYZ0 + lJM1_R;
        lBSPXYZ_Ju = (lAAABBBCCC*bspgen_bspxyz(lJ,lI) + (lK_R - lAAABBBCCC)*bspgen_bspxyz(lJ-1,lI))*lKM1_RR;
      }
      __syncthreads();
      if (lJ>=2 && lJ<=lK) {
        bspgen_bspxyz(lJ,lI) = lBSPXYZ_Ju;
      }

      if ((1+threadIdx.x)==1) {
        bspgen_bspxyz(1,lI) = bspgen_bspxyz(1,lI)*lRIXYZ0*lKM1_RR;
      }
      __syncthreads();
    }

    int lK = NOSPL_;

    if ((1+threadIdx.x)==lK) {
      bspgen_bspxyz(lK,lI) = (T_)0;
    }
    __syncthreads();

    T_ lK_R = (T_)lK;
    T_ lKM1_RR = ((T_) 1) / ((T_) (lK-1));

    int lJ = 1+threadIdx.x;
    T_ lBSPXYZ_J, lBSPXYZ_Jm1;

    if (lJ>=2 && lJ<=NOSPL_) {
      lBSPXYZ_J   = bspgen_bspxyz(lJ,lI);
      lBSPXYZ_Jm1 = bspgen_bspxyz(lJ-1,lI);
    }
    __syncthreads();

    if (lJ>=2 && lJ<=NOSPL_) {
      T_ lJM1_R = (T_) (lJ-1);
      T_ lAAABBBCCC = lRIXYZ0 + lJM1_R;
      bspgen_bsdxyz(lJ,lI) = lBSPXYZ_J - lBSPXYZ_Jm1;
      bspgen_bspxyz(lJ,lI) = (lAAABBBCCC*lBSPXYZ_J + (lK_R - lAAABBBCCC)*lBSPXYZ_Jm1)*lKM1_RR;
    }

    if (threadIdx.x==0) {
      bspgen_bsdxyz(1,lI) = bspgen_bspxyz(1,lI);
      bspgen_bspxyz(1,lI) = bspgen_bspxyz(1,lI)*lRIXYZ0*lKM1_RR;
    }
    __syncthreads();

    T_ *lBSPXYZ = CONSTANT_DATA_BSPGEN.mBSPDXYZ +
                  (0+threadIdx.y)*CONSTANT_DATA_BSPGEN.mMXSPL*CONSTANT_DATA_BSPGEN.mNLAST +
                  (lI-1)*CONSTANT_DATA_BSPGEN.mMXSPL;
    T_ *lBSDXYZ = CONSTANT_DATA_BSPGEN.mBSPDXYZ +
                  (3+threadIdx.y)*CONSTANT_DATA_BSPGEN.mMXSPL*CONSTANT_DATA_BSPGEN.mNLAST +
                  (lI-1)*CONSTANT_DATA_BSPGEN.mMXSPL;

    lBSPXYZ[lJ-1] = bspgen_bspxyz(lJ,lI);
    lBSDXYZ[lJ-1] = bspgen_bsdxyz(lJ,lI);
    lI += gridDim.y;
    __syncthreads();
  }
}


extern "C" void wrapper_f_bspgen_helper
                (int*,int*,int*,real*,real*,real*,
                 real*,real*,real*,real*,real*,real*);


extern "C" void spme_container_cuda_bspgen_invoke
                (real *aBSPX, real *aBSPY, real *aBSPZ,
                 real *aBSDX, real *aBSDY, real *aBSDZ) {
  /* Update: Due to multiple CUDA kernels needing access to the bspgen
   * results and the kernels using dynamic ratios, numerous bugs have
   * popped up here and there. Because bspgen is not particularly
   * expensive, from now on we now have both the host and the device
   * working on the entire dataset separately.
   */

  dim3 lK1_GridDims(1, 64*1024-1, 1);
  int lSharedMemorySize = (3*2*sCD.mAll.mBSPGEN.mNOSPL)*sizeof(real);
  assert(sCD.mAll.mBSPGEN.mNOSPL==8);
  start_timing_bspgen_cuda_k1();
  spme_container_cuda_bspgen_k1_NOSPLeq8<real,8>
      <<<lK1_GridDims, dim3(sCD.mAll.mBSPGEN.mNOSPL,3,1), lSharedMemorySize>>>(1, sCD.mAll.mBSPGEN.mNLAST);

  int lIATMB = 1;
  int lIATME = sCD.mAll.mBSPGEN.mNLAST;
  wrapper_f_bspgen_helper(&lIATMB, &lIATME, &sCD.mAll.mBSPGEN.mNOSPL,
                          sHD.mAll.mBSPGEN.mXXX, sHD.mAll.mBSPGEN.mYYY, sHD.mAll.mBSPGEN.mZZZ,
                          aBSPX, aBSPY, aBSPZ, aBSDX, aBSDY, aBSDZ);


  CUDA_SAFE_CALL(cudaThreadSynchronize());
  cudaError_t lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_bspgen_cuda_k1();
}
