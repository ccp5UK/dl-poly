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
#include <list>
#include <cuda.h>
#include "cutil.h"
#include <cufft.h>

#include "dl_poly_common_cu.cu"

#define CFG_CCCHARGE_K1_BLOCKDIMX_SIZE 64

template<typename T_> struct constant_data_cccharge {
  int3  mBLOCK_XYZ, mKMAXABC, mKMAXABC_Div2;
  complex<T_>  *mQQQ_LOCAL;
  int   mQQQ_LOCAL_NElements;
  complex<real> *mBSCX, *mBSCY, *mBSCZ;
  int  *mINDEX_X, *mINDEX_Y, *mINDEX_Z;
  T_    mRCELL[9];
  T_    mRCPCT2, mRALPH;
  T_   *mSTRS;
};

texture<int4, 1, cudaReadModeElementType> TEXTURE_DATA_BSCX;

__device__ __constant__ constant_data_cccharge<real> CONSTANT_DATA_CCCHARGE;
static                  constant_data_cccharge<real> sCD_CCCHARGE;

template<typename T_> struct host_data_cccharge {
  complex<T_> *mQQQ_LOCAL;
  T_          *mSTRS;
  dim3         mK1_GridDims;
};

static host_data_cccharge<real> sHD_CCCHARGE;

template<typename T_> __host__ __device__ complex<T_> make_complex(T_ aReal, T_ aImag) {
  complex<T_> lC;
  lC.mReal = aReal;
  lC.mImag = aImag;
  return (lC);
}

template<typename T_> __host__ __device__ complex<T_> conjg(const complex<T_>& aC) {
  return(make_complex(aC.mReal, -aC.mImag));
}

template<typename T_> __host__ __device__ complex<T_> operator*(const complex<T_>& aA, const complex<T_>& aB) {
  return(make_complex(aA.mReal*aB.mReal - aA.mImag*aB.mImag, aA.mReal*aB.mImag + aA.mImag*aB.mReal));
}

template<typename T_> __host__ __device__ complex<T_> operator*(const T_& aA, const complex<T_>& aB) {
  return(make_complex(aA*aB.mReal,aA*aB.mImag));
}

__device__ double fetch_double(texture<int2, 1> t, int i)
{
  int2 v = tex1Dfetch(t,i);
  return __hiloint2double(v.y, v.x);
}

template<typename T_> __device__ complex<T_> fetch_complex(texture<int4, 1> aT, int aIndex);

__device__ complex<double> fetch_complex(texture<int4, 1> aT, int aIndex) {
  int4 lV = tex1Dfetch(aT, aIndex);
  return(make_complex(__hiloint2double(lV.y, lV.x), __hiloint2double(lV.w, lV.z)));
}

template<typename T_, unsigned int BXR_, unsigned int BX_>
__global__ void ewald_spme_forces_cuda_cccharge_k1_BB() {

  DECLARE_DYNAMIC_SHARED(T_);

  T_ lTWOPI = (T_)(2.0 * 3.1415926535897932);

  for (int lI=0 ; lI<6 ; lI++)
    shared[lI*BX_ + threadIdx.x] = (T_)0;

  int lL_LOCAL = 1 + blockIdx.y;

  int lL  = CONSTANT_DATA_CCCHARGE.mINDEX_Z[lL_LOCAL-1];
  int lLL = lL-1;

  if (lL > CONSTANT_DATA_CCCHARGE.mKMAXABC_Div2.z) {
    lLL = lL - CONSTANT_DATA_CCCHARGE.mKMAXABC.z -1;
  }

  T_ lTMP_L = lTWOPI * (real)lLL;
  T_ lRKX1  = lTMP_L * CONSTANT_DATA_CCCHARGE.mRCELL[3-1];
  T_ lRKY1  = lTMP_L * CONSTANT_DATA_CCCHARGE.mRCELL[6-1];
  T_ lRKZ1  = lTMP_L * CONSTANT_DATA_CCCHARGE.mRCELL[9-1];
  complex<T_> lBSCZ_L = CONSTANT_DATA_CCCHARGE.mBSCZ[lL-1];
  T_ lBB3 = (lBSCZ_L*conjg(lBSCZ_L)).mReal;

  for (int lK_LOCAL=1 ; lK_LOCAL<=CONSTANT_DATA_CCCHARGE.mBLOCK_XYZ.y ; lK_LOCAL++) {
    int lK  = CONSTANT_DATA_CCCHARGE.mINDEX_Y[lK_LOCAL-1];
    int lKK = lK-1;

    if (lK > CONSTANT_DATA_CCCHARGE.mKMAXABC_Div2.y) {
      lKK = lK - CONSTANT_DATA_CCCHARGE.mKMAXABC.y -1;
    }

    T_ lTMP_K = lTWOPI * (real)lKK;
    T_ lRKX2  = lRKX1 + lTMP_K * CONSTANT_DATA_CCCHARGE.mRCELL[2-1];
    T_ lRKY2  = lRKY1 + lTMP_K * CONSTANT_DATA_CCCHARGE.mRCELL[5-1];
    T_ lRKZ2  = lRKZ1 + lTMP_K * CONSTANT_DATA_CCCHARGE.mRCELL[8-1];
    complex<T_> lBSCY_K = CONSTANT_DATA_CCCHARGE.mBSCY[lK-1];
    T_ lBB2 = lBB3*(lBSCY_K*conjg(lBSCY_K)).mReal;

    for (int lR=0 ; lR<BXR_ ; lR++) {
      int lJ_LOCAL = 1 + threadIdx.x + lR*BX_;

      int lJ  = CONSTANT_DATA_CCCHARGE.mINDEX_X[lJ_LOCAL-1];
      int lJJ = lJ-1;

      if (lJ > CONSTANT_DATA_CCCHARGE.mKMAXABC_Div2.x) {
        lJJ = lJ - CONSTANT_DATA_CCCHARGE.mKMAXABC.x -1;
      }

      T_ lTMP_J = lTWOPI * (real)lJJ;
      T_ lRKX3  = lRKX2 + lTMP_J * CONSTANT_DATA_CCCHARGE.mRCELL[1-1];
      T_ lRKY3  = lRKY2 + lTMP_J * CONSTANT_DATA_CCCHARGE.mRCELL[4-1];
      T_ lRKZ3  = lRKZ2 + lTMP_J * CONSTANT_DATA_CCCHARGE.mRCELL[7-1];
      complex<T_> lBSCX_J = CONSTANT_DATA_CCCHARGE.mBSCX[lJ-1];
      T_ lBB1 = lBB2*(lBSCX_J*conjg(lBSCX_J)).mReal;

      T_ lRKSQ = lRKX3*lRKX3 + lRKY3*lRKY3 + lRKZ3*lRKZ3;
      complex<T_> *lQQQ_LOCAL_JKLp=dev_f3d_address<complex<T_>,0>(CONSTANT_DATA_CCCHARGE.mQQQ_LOCAL, 1,1,1,
                                                                  CONSTANT_DATA_CCCHARGE.mBLOCK_XYZ.x,
                                                                  CONSTANT_DATA_CCCHARGE.mBLOCK_XYZ.y,
                                                                  lJ_LOCAL, lK_LOCAL, lL_LOCAL);

      complex<T_> lVTERM = make_complex((T_)0, (T_)0);
      if (lRKSQ > (T_)1.0e-6 && lRKSQ <=  CONSTANT_DATA_CCCHARGE.mRCPCT2) {

        complex<T_> lQQQ_LOCAL_JKL = *lQQQ_LOCAL_JKLp;
        lVTERM = lBB1*(exp(CONSTANT_DATA_CCCHARGE.mRALPH*lRKSQ)/lRKSQ)*lQQQ_LOCAL_JKL;
        T_ lAKV = ((T_) 2)*((((T_) 1)/lRKSQ) - CONSTANT_DATA_CCCHARGE.mRALPH)*
                   (lVTERM*conjg(lQQQ_LOCAL_JKL)).mReal;
        shared[0*BX_ + threadIdx.x] -= lRKX3*lRKX3*lAKV;
        shared[1*BX_ + threadIdx.x] -= lRKY3*lRKY3*lAKV;
        shared[2*BX_ + threadIdx.x] -= lRKZ3*lRKZ3*lAKV;
        shared[3*BX_ + threadIdx.x] -= lRKX3*lRKY3*lAKV;
        shared[4*BX_ + threadIdx.x] -= lRKX3*lRKZ3*lAKV;
        shared[5*BX_ + threadIdx.x] -= lRKY3*lRKZ3*lAKV;
      }
      *lQQQ_LOCAL_JKLp = lVTERM;

    }
  }
  if (BX_>32)
    __syncthreads();

#if !CFG_UNIFIED_ADDRESS_SPACE
  psum<T_,BX_,6>();
#else
  psum_uas<T_,BX_,6>(shared);
#endif

  if (threadIdx.x==0) {
    for (int lI=0 ; lI<6 ; lI++) {
      CONSTANT_DATA_CCCHARGE.mSTRS[blockIdx.y*6 + lI] = shared[lI*BX_];
    }
  }
}


template<typename T_, unsigned int BX_> __global__ void ewald_spme_forces_cuda_cccharge_k1() {
  DECLARE_DYNAMIC_SHARED(T_);

  T_ lTWOPI = (T_)(2.0 * 3.1415926535897932);

  for (int lI=0 ; lI<6 ; lI++)
    shared[lI*BX_ + threadIdx.x] = (T_)0;

  int lL_LOCAL = 1 + blockIdx.y;

  int lL  = CONSTANT_DATA_CCCHARGE.mINDEX_Z[lL_LOCAL-1];
  int lLL = lL-1;

  if (lL > CONSTANT_DATA_CCCHARGE.mKMAXABC_Div2.z) {
    lLL = lL - CONSTANT_DATA_CCCHARGE.mKMAXABC.z -1;
  }

  T_ lTMP_L = lTWOPI * (real)lLL;
  T_ lRKX1  = lTMP_L * CONSTANT_DATA_CCCHARGE.mRCELL[3-1];
  T_ lRKY1  = lTMP_L * CONSTANT_DATA_CCCHARGE.mRCELL[6-1];
  T_ lRKZ1  = lTMP_L * CONSTANT_DATA_CCCHARGE.mRCELL[9-1];
  complex<T_> lBSCZ_L = CONSTANT_DATA_CCCHARGE.mBSCZ[lL-1];
  T_ lBB3 = (lBSCZ_L*conjg(lBSCZ_L)).mReal;

  for (int lK_LOCAL=1 ; lK_LOCAL<=CONSTANT_DATA_CCCHARGE.mBLOCK_XYZ.y ; lK_LOCAL++) {

    int lK  = CONSTANT_DATA_CCCHARGE.mINDEX_Y[lK_LOCAL-1];
    int lKK = lK-1;

    if (lK > CONSTANT_DATA_CCCHARGE.mKMAXABC_Div2.y) {
      lKK = lK - CONSTANT_DATA_CCCHARGE.mKMAXABC.y -1;
    }

    T_ lTMP_K = lTWOPI * (real)lKK;
    T_ lRKX2  = lRKX1 + lTMP_K * CONSTANT_DATA_CCCHARGE.mRCELL[2-1];
    T_ lRKY2  = lRKY1 + lTMP_K * CONSTANT_DATA_CCCHARGE.mRCELL[5-1];
    T_ lRKZ2  = lRKZ1 + lTMP_K * CONSTANT_DATA_CCCHARGE.mRCELL[8-1];
    complex<T_> lBSCY_K = CONSTANT_DATA_CCCHARGE.mBSCY[lK-1];
    T_ lBB2 = lBB3*(lBSCY_K*conjg(lBSCY_K)).mReal;

    for (int lJ_LOCAL=1+threadIdx.x ; lJ_LOCAL<=CONSTANT_DATA_CCCHARGE.mBLOCK_XYZ.x ; lJ_LOCAL+=BX_) {

      int lJ  = CONSTANT_DATA_CCCHARGE.mINDEX_X[lJ_LOCAL-1];
      int lJJ = lJ-1;

      if (lJ > CONSTANT_DATA_CCCHARGE.mKMAXABC_Div2.x) {
          lJJ = lJ - CONSTANT_DATA_CCCHARGE.mKMAXABC.x -1;
      }

      T_ lTMP_J = lTWOPI * (real)lJJ;
      T_ lRKX3  = lRKX2 + lTMP_J * CONSTANT_DATA_CCCHARGE.mRCELL[1-1];
      T_ lRKY3  = lRKY2 + lTMP_J * CONSTANT_DATA_CCCHARGE.mRCELL[4-1];
      T_ lRKZ3  = lRKZ2 + lTMP_J * CONSTANT_DATA_CCCHARGE.mRCELL[7-1];
      complex<T_> lBSCX_J = CONSTANT_DATA_CCCHARGE.mBSCX[lJ-1];
      T_ lBB1 = lBB2*(lBSCX_J*conjg(lBSCX_J)).mReal;

      T_ lRKSQ = lRKX3*lRKX3 + lRKY3*lRKY3 + lRKZ3*lRKZ3;
      complex<T_> *lQQQ_LOCAL_JKLp=dev_f3d_address<complex<T_>,0>(CONSTANT_DATA_CCCHARGE.mQQQ_LOCAL, 1,1,1,
                                                                  CONSTANT_DATA_CCCHARGE.mBLOCK_XYZ.x,
                                                                  CONSTANT_DATA_CCCHARGE.mBLOCK_XYZ.y,
                                                                  lJ_LOCAL, lK_LOCAL, lL_LOCAL);

      complex<T_> lVTERM = make_complex((T_)0, (T_)0);
      if (lRKSQ > (T_)1.0e-6 && lRKSQ <=  CONSTANT_DATA_CCCHARGE.mRCPCT2) {

          complex<T_> lQQQ_LOCAL_JKL = *lQQQ_LOCAL_JKLp;
          lVTERM = lBB1*(exp(CONSTANT_DATA_CCCHARGE.mRALPH*lRKSQ)/lRKSQ)*lQQQ_LOCAL_JKL;
          T_ lAKV = ((T_) 2)*((((T_) 1)/lRKSQ) - CONSTANT_DATA_CCCHARGE.mRALPH)*
                    (lVTERM*conjg(lQQQ_LOCAL_JKL)).mReal;
          shared[0*BX_ + threadIdx.x] -= lRKX3*lRKX3*lAKV;
          shared[1*BX_ + threadIdx.x] -= lRKY3*lRKY3*lAKV;
          shared[2*BX_ + threadIdx.x] -= lRKZ3*lRKZ3*lAKV;
          shared[3*BX_ + threadIdx.x] -= lRKX3*lRKY3*lAKV;
          shared[4*BX_ + threadIdx.x] -= lRKX3*lRKZ3*lAKV;
          shared[5*BX_ + threadIdx.x] -= lRKY3*lRKZ3*lAKV;
      }
      *lQQQ_LOCAL_JKLp = lVTERM;

    }
  }
  if (BX_>32)
    __syncthreads();

#if !CFG_UNIFIED_ADDRESS_SPACE
  psum<T_,BX_,6>();
#else
  psum_uas<T_,BX_,6>(shared);
#endif


  if (threadIdx.x==0) {
    for (int lI=0 ; lI<6 ; lI++) {
      CONSTANT_DATA_CCCHARGE.mSTRS[blockIdx.y*6 + lI] = shared[lI*BX_];
    }
  }
}

extern "C" void ewald_spme_forces_cuda_cccharge_invoke() {
  start_timing_ewald_spme_forces_cuda_cccharge_k1();
  dim3 lK1_BlockDims = dim3(64,1,1);
  int  lK1_SM        = 6*lK1_BlockDims.x*sizeof(real);

  if ((sCD_CCCHARGE.mBLOCK_XYZ.x % CFG_CCCHARGE_K1_BLOCKDIMX_SIZE)==0) {
    switch (sCD_CCCHARGE.mBLOCK_XYZ.x/CFG_CCCHARGE_K1_BLOCKDIMX_SIZE) {
    case 2: {
      ewald_spme_forces_cuda_cccharge_k1_BB
          <real,2,CFG_CCCHARGE_K1_BLOCKDIMX_SIZE><<<sHD_CCCHARGE.mK1_GridDims,lK1_BlockDims,lK1_SM>>>();
      break;
    }
    case 1: {
      ewald_spme_forces_cuda_cccharge_k1_BB
          <real,1,CFG_CCCHARGE_K1_BLOCKDIMX_SIZE><<<sHD_CCCHARGE.mK1_GridDims,lK1_BlockDims,lK1_SM>>>();
      break;
    }
    default: {
      printf("%s::%s: stub (soft): update the switch case to handle the case sCD_CCCHARGE.mBLOCK_XYZ.x/64=%d\n",
             __FILE__, __FUNCTION__, sCD_CCCHARGE.mBLOCK_XYZ.x/64);
      ewald_spme_forces_cuda_cccharge_k1
          <real,CFG_CCCHARGE_K1_BLOCKDIMX_SIZE><<<sHD_CCCHARGE.mK1_GridDims,lK1_BlockDims,lK1_SM>>>();
    }
    }
  } else {
    ewald_spme_forces_cuda_cccharge_k1
        <real,CFG_CCCHARGE_K1_BLOCKDIMX_SIZE><<<sHD_CCCHARGE.mK1_GridDims,lK1_BlockDims,lK1_SM>>>();
  }


  CUDA_SAFE_CALL(cudaThreadSynchronize());
  cudaError_t lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_ewald_spme_forces_cuda_cccharge_k1();

  start_timing_ewald_spme_forces_cuda_cccharge_read();
  CUDA_SAFE_CALL(cudaMemcpy(sHD_CCCHARGE.mQQQ_LOCAL, sCD_CCCHARGE.mQQQ_LOCAL,
                            sCD_CCCHARGE.mQQQ_LOCAL_NElements*sizeof(complex<real>),
                            cudaMemcpyDeviceToHost));
  real *lAllSTRS = new real[sHD_CCCHARGE.mK1_GridDims.y*6];

  CUDA_SAFE_CALL(cudaMemcpy(lAllSTRS, sCD_CCCHARGE.mSTRS,
                            sHD_CCCHARGE.mK1_GridDims.y*6*sizeof(real), cudaMemcpyDeviceToHost));
  stop_timing_ewald_spme_forces_cuda_cccharge_read();

  start_timing_ewald_spme_forces_cuda_cccharge_finalise();
  for (int lI=0 ; lI<sHD_CCCHARGE.mK1_GridDims.y ; lI++) {
    real *lSTRS = lAllSTRS + lI*6;
    sHD_CCCHARGE.mSTRS[1-1] += lSTRS[0];
    sHD_CCCHARGE.mSTRS[5-1] += lSTRS[1];
    sHD_CCCHARGE.mSTRS[9-1] += lSTRS[2];
    sHD_CCCHARGE.mSTRS[2-1] += lSTRS[3];
    sHD_CCCHARGE.mSTRS[3-1] += lSTRS[4];
    sHD_CCCHARGE.mSTRS[6-1] += lSTRS[5];
  }
  stop_timing_ewald_spme_forces_cuda_cccharge_finalise();
  delete[] lAllSTRS;
}

extern "C" void ewald_spme_forces_cuda_cccharge_initialise
                (int *aBLOCK_X, int *aBLOCK_Y, int *aBLOCK_Z,
                 int *aKMAXA, int *aKMAXB, int *aKMAXC,
                 complex<real> *aQQQ_LOCAL,
                 complex<real> *aBSCX, complex<real> *aBSCY, complex<real> *aBSCZ,
                 int *aINDEX_X, int *aINDEX_Y, int *aINDEX_Z,
                 real *aRCELL, real *aRCPCT2, real *aRALPH, real *aSTRS) {


  sCD_CCCHARGE.mBLOCK_XYZ    = make_int3(*aBLOCK_X, *aBLOCK_Y, *aBLOCK_Z);
  sCD_CCCHARGE.mKMAXABC      = make_int3(*aKMAXA, *aKMAXB, *aKMAXC);
  sCD_CCCHARGE.mKMAXABC_Div2 = make_int3((*aKMAXA)/2, (*aKMAXB)/2, (*aKMAXC)/2);
  sCD_CCCHARGE.mQQQ_LOCAL_NElements =
      sCD_CCCHARGE.mBLOCK_XYZ.x*sCD_CCCHARGE.mBLOCK_XYZ.y*sCD_CCCHARGE.mBLOCK_XYZ.z;

  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mQQQ_LOCAL,
                            sCD_CCCHARGE.mQQQ_LOCAL_NElements*sizeof(complex<real>)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mINDEX_X, sCD_CCCHARGE.mBLOCK_XYZ.x*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mINDEX_Y, sCD_CCCHARGE.mBLOCK_XYZ.y*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mINDEX_Z, sCD_CCCHARGE.mBLOCK_XYZ.z*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mBSCX, sCD_CCCHARGE.mKMAXABC.x*sizeof(complex<real>)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mBSCY, sCD_CCCHARGE.mKMAXABC.y*sizeof(complex<real>)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mBSCZ, sCD_CCCHARGE.mKMAXABC.z*sizeof(complex<real>)));

  for (int lI=0 ; lI<9 ; lI++) {
    sCD_CCCHARGE.mRCELL[lI] = aRCELL[lI];
  }

  sCD_CCCHARGE.mRCPCT2 = *aRCPCT2;
  sCD_CCCHARGE.mRALPH  = *aRALPH;

  sHD_CCCHARGE.mK1_GridDims = dim3(1, sCD_CCCHARGE.mBLOCK_XYZ.z, 1);
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCCHARGE.mSTRS, sHD_CCCHARGE.mK1_GridDims.y*6*sizeof(real)));

  sHD_CCCHARGE.mQQQ_LOCAL = aQQQ_LOCAL;
  sHD_CCCHARGE.mSTRS      = aSTRS;

  start_timing_ewald_spme_forces_cuda_cccharge_write();
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCCHARGE.mQQQ_LOCAL, aQQQ_LOCAL,
                            sCD_CCCHARGE.mQQQ_LOCAL_NElements*sizeof(complex<real>), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCCHARGE.mINDEX_X, aINDEX_X,
                            sCD_CCCHARGE.mBLOCK_XYZ.x*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCCHARGE.mINDEX_Y, aINDEX_Y,
                            sCD_CCCHARGE.mBLOCK_XYZ.y*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCCHARGE.mINDEX_Z, aINDEX_Z,
                            sCD_CCCHARGE.mBLOCK_XYZ.z*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCCHARGE.mBSCX, aBSCX,
                            sCD_CCCHARGE.mKMAXABC.x*sizeof(complex<real>), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCCHARGE.mBSCY, aBSCY,
                            sCD_CCCHARGE.mKMAXABC.y*sizeof(complex<real>), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCCHARGE.mBSCZ, aBSCZ,
                            sCD_CCCHARGE.mKMAXABC.z*sizeof(complex<real>), cudaMemcpyHostToDevice));

  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA_CCCHARGE, (void*)&sCD_CCCHARGE,
                                    sizeof(constant_data_cccharge<real>)));
  stop_timing_ewald_spme_forces_cuda_cccharge_write();
}

extern "C" void ewald_spme_forces_cuda_cccharge_finalise() {
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mQQQ_LOCAL));
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mINDEX_X));
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mINDEX_Y));
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mINDEX_Z));
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mBSCX));
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mBSCY));
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mBSCZ));
  CUDA_SAFE_CALL(cudaFree(sCD_CCCHARGE.mSTRS));
}


/* ======================================
 * Construction of Charge Array (ccarray)
 * ======================================
 */

#define CFG_CCARRAY_DYNAMIC_RATIOS (CFG_OVERLAP_WITH_HOST)

#define CFG_CCARRAY_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE (CFG_OVERLAP_WITH_HOST ? 0.5 : 1.0)

#define CFG_CCARRAY_K1_GRID_DIMS_X   8
#define CFG_CCARRAY_K1_GRID_DIMS_Y  16


extern "C" double ewald_spme_forces_cuda_ccarray_percentage_of_iterations_offloaded_to_the_device() {
  if (dl_poly_cuda_offload_ewald_spme_forces()==dl_poly_cuda_fortran_false())
    return(0.0);
  return (CFG_CCARRAY_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE);
}


template<typename T_> struct constant_data_ccarray {
  int   mNLAST;
  int   mMXATMS;
  int   mMXSPL;
  int3  mBLOCK_XYZ;
  int3  mIXYZB, mIXYZT;
  int  *mIXXYYZZ[3];
  T_   *mBSPXYZ[3];
  int  *mIT;
  T_   *mCHGE;
  T_   *mQQC_LOCAL;
  int   mQQC_LOCAL_NElements;
};

__device__ __constant__ constant_data_ccarray<real> CONSTANT_DATA_CCARRAY;
static                  constant_data_ccarray<real> sCD_CCARRAY;

template<typename T_> struct host_data_ccarray {
  int mDynamicRatios;

  T_   *mQQC_LOCAL;
  dim3  mK1_GridDims;
  double mPercentageOffIterationsOffloadedToTheDevice;
  T_   *mBSPXYZ[3];
  int  *mIXXYYZZ[3];
  int  *mIT;
  T_   *mCHGE;
  int   mAreBSPDDataOnline;
};

static host_data_ccarray<real> sHD_CCARRAY;

template<typename T_> __device__ inline
T_* f3d_address(T_ *aBase, int aXS, int aYS, int aZS, int aXV, int aYV, int aX, int aY, int aZ) {
  return (aBase + __mul24(aXV, __mul24((aZ-aZS),aYV)+(aY-aYS)) +  (aX-aXS));
}


template<typename T_, unsigned int GX_, unsigned int BX_, unsigned int BY_> __device__
void ewald_spme_forces_cuda_ccarray_k1_clear_buffer(){
  T_ *lQQC_LOCAL = CONSTANT_DATA_CCARRAY.mQQC_LOCAL+blockIdx.y*CONSTANT_DATA_CCARRAY.mQQC_LOCAL_NElements;

  int lLBegin = (blockIdx.x==0) ? GX_ : blockIdx.x;
  T_ *lBase = f3d_address(lQQC_LOCAL, 1,1,1,
                          CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.x,
                          CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.y,1,1,lLBegin);
  for (int lL=lLBegin ; lL<=CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.z ; lL += GX_) {

    int lKJ = threadIdx.x + BX_*threadIdx.y;

    while (lKJ<CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.x*CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.y) {
      lBase[lKJ] = (T_)0;
      lKJ += BX_*BY_;
    }
    lBase += GX_*CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.x*CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.y;
  }

  __syncthreads();

}


template<typename T_, int MXSPL_, unsigned int GX_, unsigned int BX_, unsigned int BY_>
__global__ void ewald_spme_forces_cuda_ccarray_k1(int aI, int aLastI) {

  int lI = aI + blockIdx.y;

  T_ *lQQC_LOCAL = f3d_address(CONSTANT_DATA_CCARRAY.mQQC_LOCAL+
                               blockIdx.y*CONSTANT_DATA_CCARRAY.mQQC_LOCAL_NElements,1,1,1,
                               CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.x, CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.y,
                               1 + (1 - CONSTANT_DATA_CCARRAY.mIXYZB.x),
                               1 + (1 - CONSTANT_DATA_CCARRAY.mIXYZB.y),
                               1 + (1 - CONSTANT_DATA_CCARRAY.mIXYZB.z))
                    -CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.x * CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.y
                    -CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.x - 1;

  ewald_spme_forces_cuda_ccarray_k1_clear_buffer<T_,GX_,BX_,BY_>();

  while (lI<=aLastI) {
    int lIT_I = CONSTANT_DATA_CCARRAY.mIT[lI-1];
    if (lIT_I==1) {
      int3 lIXXYYZZ = make_int3(CONSTANT_DATA_CCARRAY.mIXXYYZZ[0][lI-1],
                                CONSTANT_DATA_CCARRAY.mIXXYYZZ[1][lI-1],
                                CONSTANT_DATA_CCARRAY.mIXXYYZZ[2][lI-1]);
      int3 lJJKKLLB = max(CONSTANT_DATA_CCARRAY.mIXYZB, lIXXYYZZ - MXSPL_ + 2);
      int3 lJJKKLLT = min(CONSTANT_DATA_CCARRAY.mIXYZT, lIXXYYZZ + 1);
      T_   lBB3     = CONSTANT_DATA_CCARRAY.mCHGE[lI-1];

      for (int lLL = lJJKKLLB.z ; lLL<=lJJKKLLT.z ; lLL++) {
        if (((lLL - CONSTANT_DATA_CCARRAY.mIXYZB.z +1) % GX_)!=blockIdx.x)
          continue;

        int lKK = lJJKKLLB.y + threadIdx.y;
        while (lKK <= lJJKKLLT.y) {
          int lJJ = lJJKKLLB.x + threadIdx.x;

          while (lJJ<=lJJKKLLT.x) {

            T_ lDET = lBB3*
                      (*F2D_ADDRESS(CONSTANT_DATA_CCARRAY.mBSPXYZ[2], 1, 1, MXSPL_, (lIXXYYZZ.z-lLL+2), lI))*
                      (*F2D_ADDRESS(CONSTANT_DATA_CCARRAY.mBSPXYZ[1], 1, 1, MXSPL_, (lIXXYYZZ.y-lKK+2), lI))*
                      (*F2D_ADDRESS(CONSTANT_DATA_CCARRAY.mBSPXYZ[0], 1, 1, MXSPL_, (lIXXYYZZ.x-lJJ+2), lI));
            *f3d_address(lQQC_LOCAL, 0,0,0,
                         CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.x, CONSTANT_DATA_CCARRAY.mBLOCK_XYZ.y,
                         lJJ, lKK, lLL) += lDET;

            if (MXSPL_!=BX_)
              lJJ += BX_;
            else
              break;
          }

          if (MXSPL_!=BY_)
            lKK += blockDim.y;
          else
            break;
        }
        if (MXSPL_==GX_)
          break;
      }
      __syncthreads();
    }
    lI += gridDim.y;
  }
}

template<typename T_, int PO2_, unsigned int BX_>
__global__ void ewald_spme_forces_cuda_ccarray_k2(int aI, dim3 aK1_GridDims) {

  if (aK1_GridDims.y==1)
    return;
  int lI = aI + BX_*blockIdx.y + threadIdx.x - 1;
  while (lI<CONSTANT_DATA_CCARRAY.mQQC_LOCAL_NElements) {
    T_ lV = (CONSTANT_DATA_CCARRAY.mQQC_LOCAL+0*CONSTANT_DATA_CCARRAY.mQQC_LOCAL_NElements)[lI];
    if (PO2_) {
#pragma unroll 2
      for (int lJ=1 ; lJ<aK1_GridDims.y ; lJ++) {
        lV += (CONSTANT_DATA_CCARRAY.mQQC_LOCAL+lJ*CONSTANT_DATA_CCARRAY.mQQC_LOCAL_NElements)[lI];
      }
    } else {
      for (int lJ=1 ; lJ<aK1_GridDims.y ; lJ++) {
        lV += (CONSTANT_DATA_CCARRAY.mQQC_LOCAL+lJ*CONSTANT_DATA_CCARRAY.mQQC_LOCAL_NElements)[lI];
      }
    }
    (CONSTANT_DATA_CCARRAY.mQQC_LOCAL+0*CONSTANT_DATA_CCARRAY.mQQC_LOCAL_NElements)[lI] = lV;
    lI += BX_*gridDim.y;
  }
}

extern "C" void ewald_spme_forces_cuda_ccarray_invoke() {
  cudaError_t lLastError;

  int lUnroll_D = (int) (sHD_CCARRAY.mPercentageOffIterationsOffloadedToTheDevice * (double) sCD_CCARRAY.mNLAST);
  int lUnroll_H = sCD_CCARRAY.mNLAST-lUnroll_D;
  int lI_H=lUnroll_D+1;
  int lI_D=1;
  real *lQQC_LOCALL_H = NULL;

  float lCE_K1_Hst_ElapsedTime = 0.0f, lCE_K1_Dev_ElapsedTime = 0.0f;
  cudaEvent_t lCE_K1_DevStart, lCE_K1_DevStop;

  if (lUnroll_D==0)
    lQQC_LOCALL_H = sHD_CCARRAY.mQQC_LOCAL;
  else
    lQQC_LOCALL_H = new real[sCD_CCARRAY.mQQC_LOCAL_NElements];

  if (lUnroll_D>0) {
    start_timing_ewald_spme_forces_cuda_ccarray_k1();
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStart));
    CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStop));

#define INVOKE_K1(GX)\
    ewald_spme_forces_cuda_ccarray_k1<real,8,GX,8,8><<<sHD_CCARRAY.mK1_GridDims, dim3(8,8,1)>>>(lI_D, lI_D+lUnroll_D-1)

    CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStart, 0));

    switch (sHD_CCARRAY.mK1_GridDims.x) {
    case  2: { INVOKE_K1(2); break; }
    case  4: { INVOKE_K1(4); break; }
    case  8: { INVOKE_K1(8); break; }
    case  16: { INVOKE_K1(16); break; }
    case  32: { INVOKE_K1(32); break; }
    default: {
      printf("%s::%s: stub: sHD_CCARRAY.mK1_GridDims.x=%d.\n",
             __FILE__, __FUNCTION__,sHD_CCARRAY.mK1_GridDims.x);
      exit(-1); }
    }

    CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStop, 0));
  }

  if (lUnroll_H>0) {
    for (int lI=0 ; lI<sCD_CCARRAY.mQQC_LOCAL_NElements ; lI++) {
      lQQC_LOCALL_H[lI] = (real)0;
    }

    struct timeval lTV_K1_HstStart, lTV_K1_HstStop;
    gettimeofday(&lTV_K1_HstStart, NULL);

    wrapper_f_ewald_spme_forces_ccarray_helper(&lI_H, &lUnroll_H,
				      &sCD_CCARRAY.mBLOCK_XYZ.x,&sCD_CCARRAY.mBLOCK_XYZ.y,&sCD_CCARRAY.mBLOCK_XYZ.z,
				      &sCD_CCARRAY.mIXYZB.x,&sCD_CCARRAY.mIXYZB.y,&sCD_CCARRAY.mIXYZB.z,
				      &sCD_CCARRAY.mIXYZT.x,&sCD_CCARRAY.mIXYZT.y,&sCD_CCARRAY.mIXYZT.z,
				      sHD_CCARRAY.mIXXYYZZ[0], sHD_CCARRAY.mIXXYYZZ[1], sHD_CCARRAY.mIXXYYZZ[2], sHD_CCARRAY.mIT,
				      sHD_CCARRAY.mBSPXYZ[0], sHD_CCARRAY.mBSPXYZ[1], sHD_CCARRAY.mBSPXYZ[2],
				      lQQC_LOCALL_H);

    gettimeofday(&lTV_K1_HstStop, NULL);
    lCE_K1_Hst_ElapsedTime = secsfromtimeval(lTV_K1_HstStop) - secsfromtimeval(lTV_K1_HstStart);
  }

  if (lUnroll_D>0) {
    CUDA_SAFE_CALL(cudaEventSynchronize(lCE_K1_DevStop));
    CUDA_SAFE_CALL(cudaEventElapsedTime(&lCE_K1_Dev_ElapsedTime, lCE_K1_DevStart, lCE_K1_DevStop));
    lCE_K1_Dev_ElapsedTime /= 1000.0f;
//malysaght310512
    CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStart));
    CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStop));
//end_malysaght310512

    CUDA_SAFE_CALL(cudaThreadSynchronize());
    lLastError = cudaGetLastError();
    CUT_CHECK_ERROR(lLastError);
    stop_timing_ewald_spme_forces_cuda_ccarray_k1();

    start_timing_ewald_spme_forces_cuda_ccarray_k2();

    int lK2_NOfThreadBlocks = 128;
    int lK2_NOfBlockThreads = 64;

#define INVOKE_K2(NT)\
    do {\
      if ((sHD_CCARRAY.mK1_GridDims.y % 2)==0)\
        ewald_spme_forces_cuda_ccarray_k2<real,1,NT>                  \
            <<<dim3(1,lK2_NOfThreadBlocks,1), dim3((NT),1,1)>>>(1,sHD_CCARRAY.mK1_GridDims); \
      else\
        ewald_spme_forces_cuda_ccarray_k2<real,0,NT>                  \
            <<<dim3(1,lK2_NOfThreadBlocks,1), dim3((NT),1,1)>>>(1,sHD_CCARRAY.mK1_GridDims); \
    } while (0)

    switch (lK2_NOfBlockThreads) {
    case 64: { INVOKE_K2(64); break; }
    case 32: { INVOKE_K2(32); break; }
    case 16: { INVOKE_K2(16); break; }
    case  8: { INVOKE_K2(8); break; }
    default: {
      printf("%s::%s: stub: sHD_CCARRAY.mK1_GridDims.y=%d.\n",
	     __FILE__, __FUNCTION__,sHD_CCARRAY.mK1_GridDims.y);
      exit(-1); }
    };

    CUDA_SAFE_CALL(cudaThreadSynchronize());
    lLastError = cudaGetLastError();
    CUT_CHECK_ERROR(lLastError);
    stop_timing_ewald_spme_forces_cuda_ccarray_k2();
  }


  if (lUnroll_D>0) {
    start_timing_ewald_spme_forces_cuda_ccarray_read();
    CUDA_SAFE_CALL(cudaMemcpy(sHD_CCARRAY.mQQC_LOCAL, sCD_CCARRAY.mQQC_LOCAL,
                              sCD_CCARRAY.mQQC_LOCAL_NElements*sizeof(real),
                              cudaMemcpyDeviceToHost));
    stop_timing_ewald_spme_forces_cuda_ccarray_read();
  }

  if (lUnroll_D!=0 && lUnroll_H!=0) {
    /* Final reduction with those iterations that the host handled.
     */
    wrapper_f_ewald_spme_forces_ccarray_final_reduction
        (&sCD_CCARRAY.mBLOCK_XYZ.x,
         &sCD_CCARRAY.mBLOCK_XYZ.y,
         &sCD_CCARRAY.mBLOCK_XYZ.z,
         sHD_CCARRAY.mQQC_LOCAL, lQQC_LOCALL_H);
//malysaght110612: note memory leak here -- fix
//    delete[] lQQC_LOCALL_H;
  }
  delete[] lQQC_LOCALL_H;

  if (sHD_CCARRAY.mDynamicRatios) {
    double lPI_Dev = lCE_K1_Dev_ElapsedTime / ((double) lUnroll_D);
    double lPI_Hst = lCE_K1_Hst_ElapsedTime / ((double) lUnroll_H);
    double lNewRatio = lPI_Hst / (lPI_Dev + lPI_Hst);

    sHD_CCARRAY.mPercentageOffIterationsOffloadedToTheDevice = lNewRatio;
  }
}

extern "C" void ewald_spme_forces_cuda_ccarray_initialise
                (int *aNLAST, int *aMXSPL, int *aMXATMS,
                 int *aIXB, int *aIYB, int *aIZB, int *aIXT, int *aIYT, int *aIZT,
                 int *aIXX, int *aIYY, int *aIZZ,
                 real *aBSPX, real *aBSPY, real *aBSPZ,
                 real *aBSDX, real *aBSDY, real *aBSDZ,
                 int *aIT, real *aCHGE,
                 int *aBLOCK_X, int *aBLOCK_Y, int *aBLOCK_Z, real *aQQC_LOCAL) {


  sHD_CCARRAY.mDynamicRatios =
      dl_poly_cuda_getenv("dlpolycuda_ewald_spme_forces_ccarray_dynamic_ratios",
                          CFG_CCARRAY_DYNAMIC_RATIOS);
  static int sSetRatios = 1;
  if (sSetRatios) {
    sHD_CCARRAY.mPercentageOffIterationsOffloadedToTheDevice =
        CFG_CCARRAY_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE;
    if (sHD_CCARRAY.mDynamicRatios) {
      sSetRatios = 0;
    }
  }

  int  *lIntPointers_HST[3]  = { aIXX,  aIYY,  aIZZ  };
  real *lRealPointers_HST[3] = { aBSPX, aBSPY, aBSPZ };

  sCD_CCARRAY.mNLAST     = *aNLAST;
  sCD_CCARRAY.mMXSPL     = *aMXSPL;
  sCD_CCARRAY.mMXATMS    = *aMXATMS;
  sCD_CCARRAY.mBLOCK_XYZ = make_int3(*aBLOCK_X, *aBLOCK_Y, *aBLOCK_Z);
  sCD_CCARRAY.mIXYZB     = make_int3(*aIXB, *aIYB, *aIZB);
  sCD_CCARRAY.mIXYZT     = make_int3(*aIXT, *aIYT, *aIZT);
  sCD_CCARRAY.mQQC_LOCAL_NElements =
      sCD_CCARRAY.mBLOCK_XYZ.x*sCD_CCARRAY.mBLOCK_XYZ.y*sCD_CCARRAY.mBLOCK_XYZ.z;

  int lK1_Grid_Dims_X = dl_poly_cuda_getenv("dlpolycuda_ewald_spme_forces_ccarray_k1_grid_dims_x",
                                            CFG_CCARRAY_K1_GRID_DIMS_X);
                                            //CFG_CCARRAY_K1_GRID_DIMS_Y); //Above line is updated by Maria
  int lK1_Grid_Dims_Y = dl_poly_cuda_getenv("dlpolycuda_ewald_spme_forces_ccarray_k1_grid_dims_y",
                                            CFG_CCARRAY_K1_GRID_DIMS_Y);

  sHD_CCARRAY.mK1_GridDims = dim3(lK1_Grid_Dims_X,lK1_Grid_Dims_Y,1);
  sHD_CCARRAY.mQQC_LOCAL = aQQC_LOCAL;
  sHD_CCARRAY.mPercentageOffIterationsOffloadedToTheDevice =
      CFG_CCARRAY_PERCENTAGE_OF_ITERATIONS_OFFLOADED_TO_THE_DEVICE;

  for (int lI=0 ; lI<3 ; lI++) {
    sHD_CCARRAY.mIXXYYZZ[lI] = lIntPointers_HST[lI];
    sHD_CCARRAY.mBSPXYZ[lI]  = lRealPointers_HST[lI];
  }

  sHD_CCARRAY.mIT   = aIT;
  sHD_CCARRAY.mCHGE = aCHGE;
  sHD_CCARRAY.mAreBSPDDataOnline =
      (*spme_container_cuda_bspgen_leave_bspdxyz_data_online())==dl_poly_cuda_fortran_true();


  /* When the device's work is overlapped with that of the host, as some iterations will
   * not be executed by the device, the relevant data need no copying over to the device.
   * Note that the device will be working on the range 1:lUnroll_D-1, i.e. it is the
   * first lUnroll_D that we need to make space for and copy over.
   */
  int lUnroll_D = (int) (sHD_CCARRAY.mPercentageOffIterationsOffloadedToTheDevice *
                         (double) sCD_CCARRAY.mNLAST);

  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCARRAY.mQQC_LOCAL,
                            sCD_CCARRAY.mQQC_LOCAL_NElements*sHD_CCARRAY.mK1_GridDims.y*sizeof(real)));


  start_timing_ewald_spme_forces_cuda_ccarray_write();
  for (int lI=0 ; lI<3 ; lI++) {
    CUDA_SAFE_CALL(cudaMalloc(&sCD_CCARRAY.mIXXYYZZ[lI],
                              lUnroll_D*sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpy(sCD_CCARRAY.mIXXYYZZ[lI], lIntPointers_HST[lI],
                              lUnroll_D*sizeof(int), cudaMemcpyHostToDevice));
    if (!sHD_CCARRAY.mAreBSPDDataOnline) {
      CUDA_SAFE_CALL(cudaMalloc(&sCD_CCARRAY.mBSPXYZ[lI],
                                sCD_CCARRAY.mMXSPL*lUnroll_D*sizeof(real)));
      CUDA_SAFE_CALL(cudaMemcpy(sCD_CCARRAY.mBSPXYZ[lI], lRealPointers_HST[lI],
                                sCD_CCARRAY.mMXSPL*lUnroll_D*sizeof(real),  cudaMemcpyHostToDevice));
    } else {
      void *lDummy;
      spme_container_cuda_bspgen_grab_bspdxyz_data((void**)&sCD_CCARRAY.mBSPXYZ[0],
                                                   (void**)&sCD_CCARRAY.mBSPXYZ[1],
                                                   (void**)&sCD_CCARRAY.mBSPXYZ[2],
                                                   &lDummy, &lDummy, &lDummy);
    }
  }
  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCARRAY.mIT,
                            lUnroll_D*sizeof(int)));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCARRAY.mIT, aIT,
                            lUnroll_D*sizeof(int), cudaMemcpyHostToDevice));

  CUDA_SAFE_CALL(cudaMalloc(&sCD_CCARRAY.mCHGE,
                            lUnroll_D*sizeof(real)));
  CUDA_SAFE_CALL(cudaMemcpy(sCD_CCARRAY.mCHGE, aCHGE,
                            lUnroll_D*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA_CCARRAY, (void*)&sCD_CCARRAY,
                                    sizeof(constant_data_ccarray<real>)));
  stop_timing_ewald_spme_forces_cuda_ccarray_write();
}

extern "C" void ewald_spme_forces_cuda_ccarray_finalise() {
  for (int lI=0 ; lI<3 ; lI++) {
    CUDA_SAFE_CALL(cudaFree(sCD_CCARRAY.mIXXYYZZ[lI]));
    if (!sHD_CCARRAY.mAreBSPDDataOnline) {
      CUDA_SAFE_CALL(cudaFree(sCD_CCARRAY.mBSPXYZ[lI]));
    }
  }
  CUDA_SAFE_CALL(cudaFree(sCD_CCARRAY.mQQC_LOCAL));
  CUDA_SAFE_CALL(cudaFree(sCD_CCARRAY.mIT));
  CUDA_SAFE_CALL(cudaFree(sCD_CCARRAY.mCHGE));
}




