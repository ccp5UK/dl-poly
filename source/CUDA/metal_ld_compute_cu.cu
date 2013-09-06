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

#define CFG_K2_THREADBLOCKS                        64
#define CFG_K2_THREADS_PER_BLOCK                   64
#define CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY    1
#define CFG_K1_MXMET_MAX_VALUE                      4
#define CFG_K1_TEXTURE_CACHE_DMET_1234_K0_1         1

template<typename T_> struct constant_data {
  int mNATMS, mKEYPOT, mMXLIST, mMXMET, mMXGRID, mMXATMS;
#if (CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
  int mLTPMET[CFG_K1_MXMET_MAX_VALUE];
  int mLSTMET[CFG_K1_MXMET_MAX_VALUE];
  T_  mDMET_1234_K0_1[CFG_K1_MXMET_MAX_VALUE][4];
#else
  int *mLTPMET, *mLSTMET;
#endif

  int *mLTYPE, *mLIST;
  T_ *mVMET, *mDMET, *mRHO;
  T_ *mXXX, *mYYY, *mZZZ;
  T_ mCELL[9], mCELL_REC[9], mCELL_INVERTED[9];
  T_ *mRHO_JATM_ACCU;
  T_ *mRHO_IATM_ACCU;

};

__device__ __constant__ constant_data<real> CONSTANT_DATA;
static constant_data<real>                  sCD;

template<typename T_> struct host_data {
  int mIsListAlreadyOnline, mFreeList;
  int mNTPMET, mMXMET, mMXATMS, mMXGRID, mMXATDM, mMXLIST;
  T_ *mRHO;
};

static host_data<real> sHD;

extern "C" void metal_ld_compute_cuda_initialise
                (int *aIsListOnline, int *aMXATMS, int *aNATMS,
                 int *aMXGRID, int *aNTPMET, int *aMXMET, int *aMXATDM, int *aMXLIST,
                 real *aXXX, real *aYYY, real *aZZZ,
                 int *aLIST, int *aLTYPE, int *aLTPMET, int *aLSTMET,
                 real *aVMET, real *aDMET, real *aCELL, real *aRHO)
{

  sHD.mIsListAlreadyOnline = *aIsListOnline;
  sHD.mMXATMS = sCD.mMXATMS = *aMXATMS;
  sCD.mNATMS  = *aNATMS;
  sHD.mMXGRID = sCD.mMXGRID = *aMXGRID;
  sHD.mMXMET  = sCD.mMXMET  = *aMXMET;
  sHD.mNTPMET = *aNTPMET;
  sHD.mMXLIST = sCD.mMXLIST = *aMXLIST;
  sHD.mMXATDM = *aMXATDM;
  sHD.mRHO    = aRHO;

  for (int lI=0 ; lI<9 ; lI++) {
    real lCELL        = aCELL[lI];
    sCD.mCELL[lI]     = lCELL;
    sCD.mCELL_REC[lI] = ((real) 1) / lCELL;
  }
  real lDummy;
  wrapper_f_invert(sCD.mCELL, sCD.mCELL_INVERTED, &lDummy);

  if (CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY && sHD.mMXMET > CFG_K1_MXMET_MAX_VALUE) {
    printf("%s::%s: can only handle mxmet<%d when CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY"
	   " is enabled; found %d. Try setting CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY"
           " to zero in %s\n",
	   __FILE__, __FUNCTION__, CFG_K1_MXMET_MAX_VALUE, sHD.mMXMET, __FILE__);
    exit(-1);
  }

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mXXX, sHD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mYYY, sHD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mZZZ, sHD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLTYPE, sHD.mMXATMS*sizeof(int)));

#if (!CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLSTMET, sHD.mMXMET*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLTPMET, sHD.mMXMET*sizeof(int)));
#endif

#if (CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
  /* The 'collect fst' part of the kernel, for every secondary atoms ietration,
   * loads four items (dmet(1:4,k0,1). If mxmet is small (e.g. 4 in the case of
   * TEST12), or generally the {lst,ltp}met structures fit in cmem, those
   * 4*mxmet possible accesses against gmem can be potentially serviced from
   * cmem.
   */


  for (int lI=1 ; lI<=sHD.mMXMET ; lI++) {
    if (aLSTMET[lI-1]==0)
      continue;
    int lLSTMET_I = aLSTMET[lI-1];
    for (int lJ=1 ; lJ<=4 ; lJ++) {
      sCD.mDMET_1234_K0_1[lLSTMET_I-1][lJ-1] =
        *F3D_ADDRESS(aDMET, 1,1,1, sCD.mMXGRID, sCD.mMXMET, lJ, lLSTMET_I, 1);
      if (lJ==3) {
        sCD.mDMET_1234_K0_1[lLSTMET_I-1][lJ-1] *=sCD.mDMET_1234_K0_1[lLSTMET_I-1][lJ-1];
      }
      if (lJ==4) {
        sCD.mDMET_1234_K0_1[lLSTMET_I-1][lJ-1]
          = ((real) 1) /  sCD.mDMET_1234_K0_1[lLSTMET_I-1][lJ-1];
      }
    }
  }
#endif

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mVMET, sHD.mMXGRID*sHD.mMXMET*2*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDMET, sHD.mMXGRID*sHD.mMXMET*2*sizeof(real)));

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mRHO_JATM_ACCU, sCD.mNATMS*(sHD.mMXLIST)*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mRHO_IATM_ACCU, sCD.mNATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mRHO, CFG_K2_THREADBLOCKS*sCD.mMXATMS*sizeof(real)));

  static int sHasKeypotBeenComputed = 0;
  if (!sHasKeypotBeenComputed) {
    int lKEYPOT = 0;
    for (int lL=1 ; lL<=sHD.mNTPMET ; lL++) {
      lKEYPOT = aLTPMET[lL-2];
      if (lL>1) {
        if (lKEYPOT != aLTPMET[(lL-2)]) {
          int lError = 92;
          wrapper_f_error(&lError);
        }
      }
    }
    sCD.mKEYPOT = lKEYPOT;
    sHasKeypotBeenComputed = 1;

    if (sCD.mKEYPOT==0) {
      // TODO -- easy to implement as code is already there but *has not been tested*
      printf("%s::%s: stub: keypot==0.\n", __FILE__, __FUNCTION__);
      exit(-1);
    }
  }

  start_timing_metal_ld_compute_cuda_write();

  if (link_cell_pairs_cuda_is_in_valid_context()) {
    sHD.mFreeList = 0;
    sCD.mLIST = (int*) link_cell_pairs_cuda_get_list();
  } else {
    CUDA_SAFE_CALL(cudaMalloc((void**)&sCD.mLIST, (512+sHD.mMXATDM*(1+2+sHD.mMXLIST))*sizeof(int)));
    sHD.mFreeList = 1;
  }

  CUDA_SAFE_CALL(cudaMemcpy(sCD.mXXX, aXXX, sHD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mYYY, aYYY, sHD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mZZZ, aZZZ, sHD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLTYPE, aLTYPE, sHD.mMXATMS*sizeof(int), cudaMemcpyHostToDevice));

#if (CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
  for (int lI=0 ; lI<sHD.mMXMET ; lI++) {
    sCD.mLSTMET[lI] = aLSTMET[lI];
    sCD.mLTPMET[lI] = aLTPMET[lI];
  }
#else
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLSTMET, aLSTMET, sHD.mMXMET*sizeof(int),
			    cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLTPMET, aLTPMET, sHD.mMXMET*sizeof(int),
			    cudaMemcpyHostToDevice));
#endif

  CUDA_SAFE_CALL(cudaMemcpy(sCD.mVMET, aVMET, sHD.mMXGRID*sHD.mMXMET*2*sizeof(real),
			    cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mDMET, aDMET, sHD.mMXGRID*sHD.mMXMET*2*sizeof(real),
			    cudaMemcpyHostToDevice));

  if (link_cell_pairs_cuda_is_in_valid_context()) {
    link_cell_pairs_cuda_push_lists(1);
    sHD.mIsListAlreadyOnline = 1;
  } else {
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mLIST, aLIST, sHD.mMXATDM*(1+2+sHD.mMXLIST)*sizeof(int),
                              cudaMemcpyHostToDevice));
    sHD.mIsListAlreadyOnline = 1;
  }

  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA, (void*)&sCD, sizeof(constant_data<real>)));
  stop_timing_metal_ld_compute_cuda_write();
}

extern "C" void metal_ld_compute_cuda_finalise()
{
  CUDA_SAFE_CALL(cudaFree(sCD.mXXX));
  CUDA_SAFE_CALL(cudaFree(sCD.mYYY));
  CUDA_SAFE_CALL(cudaFree(sCD.mZZZ));
  CUDA_SAFE_CALL(cudaFree(sCD.mVMET));
  CUDA_SAFE_CALL(cudaFree(sCD.mDMET));
#if (!CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
  CUDA_SAFE_CALL(cudaFree(sCD.mLTPMET));
  CUDA_SAFE_CALL(cudaFree(sCD.mLSTMET));
#endif
  CUDA_SAFE_CALL(cudaFree(sCD.mLTYPE));
  CUDA_SAFE_CALL(cudaFree(sCD.mRHO_JATM_ACCU));
  CUDA_SAFE_CALL(cudaFree(sCD.mRHO_IATM_ACCU));
  if (sHD.mFreeList) {
    CUDA_SAFE_CALL(cudaFree(sCD.mLIST));
  }
  CUDA_SAFE_CALL(cudaFree(sCD.mRHO));
}

template<typename T_>
__device__ void obtain_iatm_specific(int aIATM, T_& aXXX_I, T_& aYYY_I, T_& aZZZ_I,
                                     int& aLIMIT, int& aAI)
{

  if (threadIdx.x==0) {
    aXXX_I = CONSTANT_DATA.mXXX[aIATM-1];
    aYYY_I = CONSTANT_DATA.mYYY[aIATM-1];
    aZZZ_I = CONSTANT_DATA.mZZZ[aIATM-1];
  }
  __syncthreads();
  aLIMIT = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0+2, aIATM);
  aAI    = CONSTANT_DATA.mLTYPE[aIATM-1];
}

template<typename T_, int IMCON_>
__device__ void obtain_jatm_rsqdf(int aJATM, T_ aXXX_I, T_ aYYY_I, T_ aZZZ_I, T_& aRSQDF)
{

  T_ lXXX_J = CONSTANT_DATA.mXXX[aJATM-1];
  T_ lYYY_J = CONSTANT_DATA.mYYY[aJATM-1];
  T_ lZZZ_J = CONSTANT_DATA.mZZZ[aJATM-1];

  T_ aXDF      = aXXX_I - lXXX_J;
  T_ aYDF      = aYYY_I - lYYY_J;
  T_ aZDF      = aZZZ_I - lZZZ_J;

  if (IMCON_==2) { // rectangular (slab) boundary conditions:
    aXDF = msub(aXDF, CONSTANT_DATA.mCELL[0], anint(CONSTANT_DATA.mCELL_REC[0]*aXDF));
    aYDF = msub(aYDF, CONSTANT_DATA.mCELL[4], anint(CONSTANT_DATA.mCELL_REC[4]*aYDF));
    aZDF = msub(aZDF, CONSTANT_DATA.mCELL[8], anint(CONSTANT_DATA.mCELL_REC[8]*aZDF));
  }

  if (IMCON_==3) { // parallelepiped boundary conditions
    T_ lXSS = CONSTANT_DATA.mCELL_INVERTED[0]*aXDF + CONSTANT_DATA.mCELL_INVERTED[3]*aYDF +
              CONSTANT_DATA.mCELL_INVERTED[6]*aZDF;
    T_ lYSS = CONSTANT_DATA.mCELL_INVERTED[1]*aXDF + CONSTANT_DATA.mCELL_INVERTED[4]*aYDF +
              CONSTANT_DATA.mCELL_INVERTED[7]*aZDF;
    T_ lZSS = CONSTANT_DATA.mCELL_INVERTED[2]*aXDF + CONSTANT_DATA.mCELL_INVERTED[5]*aYDF +
              CONSTANT_DATA.mCELL_INVERTED[8]*aZDF;

    lXSS = lXSS - anint(lXSS);
    lYSS = lYSS - anint(lYSS);
    lZSS = lZSS - anint(lZSS);

    aXDF = CONSTANT_DATA.mCELL[0]*lXSS + CONSTANT_DATA.mCELL[3]*lYSS + CONSTANT_DATA.mCELL[6]*lZSS;
    aYDF = CONSTANT_DATA.mCELL[1]*lXSS + CONSTANT_DATA.mCELL[4]*lYSS + CONSTANT_DATA.mCELL[7]*lZSS;
    aZDF = CONSTANT_DATA.mCELL[2]*lXSS + CONSTANT_DATA.mCELL[5]*lYSS + CONSTANT_DATA.mCELL[8]*lZSS;
  }
  aRSQDF = add(mul(aXDF,aXDF), add(mul(aYDF,aYDF), mul(aZDF,aZDF)));
}

template<typename T_> __device__ T_* vmet(int aX, int aY, int aZ)
{
  return(dev_f3d_address<T_,0>(CONSTANT_DATA.mVMET, 1,1,1,
                               CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX, aY, aZ));
}

template<typename T_> __device__ T_* dmet(int aX, int aY, int aZ)
{
  return(dev_f3d_address<T_,0>(CONSTANT_DATA.mDMET, 1,1,1,
                               CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX, aY, aZ));
}

template<typename T_> __device__ T_ metal_forces_calc_gamma(T_ aPPx, T_ aT1, T_ aT2, int lL)
{
  T_ lSelect_T1T2 = (aPPx < (T_)0) ? aT1 : aT2;
  T_ lPlusMinus1  = (aPPx < (T_)0) ? ((T_) 1) : ((T_) -1);
  T_ lGamma = madd(lSelect_T1T2, ((T_)0.5), mul(add(aT2,-aT1),add(aPPx, lPlusMinus1)));
//mlysaght
  if (lL == 5) {
     lGamma = aT2;
  }
  return (lGamma);
}

template<typename T_>
__device__ T_ metal_forces_calc_gamma(T_ *aBase, int aX, int aY, int aZ, T_ aPPx)
{
//mlysaght
  T_ lGVK0 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX-1, aY, aZ);
  T_ lGVK1 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX, aY, aZ);
  T_ lGVK2 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX+1, aY, aZ);
//  T_ lGVK0 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX+3, aY, aZ);
//  T_ lGVK1 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX+4, aY, aZ);
//  T_ lGVK2 = *dev_f3d_address<T_,0>(aBase, 1,1,1, CONSTANT_DATA.mMXGRID, CONSTANT_DATA.mMXMET, aX+5, aY, aZ);
  T_ lT1   = madd(lGVK1, aPPx, (lGVK1 - lGVK0));
  T_ lT2   = madd(lGVK1, aPPx, (lGVK2 - lGVK1));
  return(metal_forces_calc_gamma(aPPx, lT1, lT2, aX));
}


template<typename T_>
__device__ T_ metal_forces_calc_gamma_vmet(int aX, int aY, int aZ, T_ aPPx)
{
  return(metal_forces_calc_gamma(CONSTANT_DATA.mVMET, aX, aY, aZ, aPPx));
}

template<typename T_>
__device__ T_ metal_forces_calc_gamma_dmet(int aX, int aY, int aZ, T_ aPPx)
{
  return(metal_forces_calc_gamma(CONSTANT_DATA.mDMET, aX, aY, aZ, aPPx));
}

template<typename T_,unsigned int BX_>
__global__ void metal_ld_compute_cuda_k0(int aI)
{
  for (int lIATM=aI+blockIdx.y*BX_+threadIdx.x;
       lIATM<=CONSTANT_DATA.mNATMS ; lIATM+=gridDim.y*BX_) {
    CONSTANT_DATA.mRHO_IATM_ACCU[lIATM-1] = (T_)0;
  }
}

template<typename T_, unsigned int GY_, unsigned int BX_, int IMCON_>
__global__ void metal_ld_compute_cuda_k1(int aI)
{
  extern __shared__ T_ shared[];

  // +1 reg for the gridDim.y increment (instead of a constant)
  for (int lIATM=aI+blockIdx.y ; lIATM<=CONSTANT_DATA.mNATMS ; lIATM+=GY_) {

#define lXXX_I shared[BX_]
#define lYYY_I shared[BX_+1]
#define lZZZ_I shared[BX_+2]

    int lLIMIT, lAI;
    shared[threadIdx.x] = (T_)0;

    obtain_iatm_specific(lIATM, lXXX_I, lYYY_I, lZZZ_I, lLIMIT, lAI);

    // -1 reg with this
    if (lLIMIT<=0) {
      if (threadIdx.x==0) {
        CONSTANT_DATA.mRHO_IATM_ACCU[lIATM-1] = (T_)0;
      }
      __syncthreads();
      continue;
    }

    for (int lJ=1+threadIdx.x ; lJ<=lLIMIT ; lJ+=BX_) {
      int lJATM = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), lJ+2, lIATM);
      int lAJ;

      lAJ = CONSTANT_DATA.mLTYPE[lJATM - 1];

      int lKEY = max(lAI,lAJ) * (max(lAI,lAJ)-1) / 2 + min(lAI,lAJ);
      int lK0  = CONSTANT_DATA.mLSTMET[lKEY - 1];

#if (!CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
      T_  lDMET_3_K0_1 = *dmet<T_>(3,lK0,1);
#endif
      int lLTPMET_K0   = CONSTANT_DATA.mLTPMET[lK0-1];

      // +1 extra register for this var
      T_  lRHO_JATM_Update = (T_)0;


      if (lLTPMET_K0>0) {
        T_ lRSQ ;
        obtain_jatm_rsqdf<T_,IMCON_>(lJATM, lXXX_I, lYYY_I, lZZZ_I, lRSQ);

#if (CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
        T_  lDMET_3_K0_1_Pow2 = CONSTANT_DATA.mDMET_1234_K0_1[lK0-1][3-1];
        if (lRSQ<=lDMET_3_K0_1_Pow2) {
#else
        if (lRSQ<=(lDMET_3_K0_1*lDMET_3_K0_1)) {
#endif

#if (CFG_LSTLTPMET_FETCH_FROM_CONSTANT_MEMORY)
          T_  lDMET_2_K0_1 = CONSTANT_DATA.mDMET_1234_K0_1[lK0-1][2-1];
//mlysaght
          T_  lDMET_1_K0_1 = CONSTANT_DATA.mDMET_1234_K0_1[lK0-1][1-1];
          T_  lRDR         = CONSTANT_DATA.mDMET_1234_K0_1[lK0-1][4-1];
#else
//mlysaght
          T_  lDMET_1_K0_1 = *dmet<T_>(1,lK0,1);
          T_  lDMET_2_K0_1 = *dmet<T_>(2,lK0,1);
          T_  lDMET_4_K0_1 = *dmet<T_>(4,lK0,1);
          T_ lRDR         = ((T_) 1) / lDMET_4_K0_1;
#endif
          T_ lRRR         = sqrt(lRSQ) - lDMET_2_K0_1;

//          int lL          = nint(lRRR*lRDR);
//mlysaght
          int lL          = min(nint(lRRR*lRDR),nint(lDMET_1_K0_1-1));
          if (lL < 5) {
//             aSAFE = .false.;
             lL = 6;
          }
//end mlysaght
          T_ lPPP         = mul(lRRR,lRDR) - (T_)lL;
          T_ lDENSITY     = metal_forces_calc_gamma_dmet(lL, lK0, 1, lPPP);

          int4 lXZ        = lAI>lAJ ? make_int4(1,2,2,2) : make_int4(2,2,1,2);

          shared[threadIdx.x] += lDENSITY*(*dmet<T_>(lXZ.x,lK0,lXZ.y));

          if (lJATM<=CONSTANT_DATA.mNATMS) {
            // +4 registers for this
            //	  *F2D_ADDRESS(CONSTANT_DATA.mRHO_JATM_ACCU, 1, 1, CONSTANT_DATA.mMXLIST, lJ, lIATM) =
            lRHO_JATM_Update = lDENSITY*(*dmet<T_>(lXZ.z,lK0,lXZ.w));
          }
        }
      }
      *F2D_ADDRESS(CONSTANT_DATA.mRHO_JATM_ACCU, 1, 1, CONSTANT_DATA.mMXLIST, lJ, lIATM) = lRHO_JATM_Update;
    }
    __syncthreads();
#if !CFG_UNIFIED_ADDRESS_SPACE
    psum<T_,BX_,1>();
#else
    psum_uas<T_,BX_,1>(shared);
#endif
    if (threadIdx.x==0) {
        // +2 registers for this
        CONSTANT_DATA.mRHO_IATM_ACCU[lIATM-1] = shared[0];
    }
    __syncthreads();
  }
}

template<typename T_, unsigned int BX_>
__global__ void metal_ld_compute_cuda_k2(int aI)
{
  T_ *lRHO = CONSTANT_DATA.mRHO + blockIdx.y*CONSTANT_DATA.mMXATMS;

  /* Cleanup the 'rho' array; this is what the host does as well as necessary
   * for the += operator.
   */
  for (int lI=threadIdx.x ; lI<CONSTANT_DATA.mMXATMS ; lI+=BX_) {
    lRHO[lI] = (T_)0;
  }
  __syncthreads();

  int2 lIATM_Bounds = make_int2(1,CONSTANT_DATA.mNATMS+1);

  if (blockIdx.y!=0)
    lIATM_Bounds.x = blockIdx.y * (CONSTANT_DATA.mNATMS / gridDim.y);
  if ((blockIdx.y+1)!=gridDim.y)
    lIATM_Bounds.y = (blockIdx.y + 1)*(CONSTANT_DATA.mNATMS / gridDim.y);

  for (int lIATM=lIATM_Bounds.x ; lIATM<lIATM_Bounds.y ; lIATM++) {
    int lLIMIT = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0+2, lIATM);

    if (lLIMIT<=0)
      continue;

    for (int lJ=1+threadIdx.x ; lJ<=lLIMIT ; lJ+=BX_) {
      int lJATM = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), lJ+2, lIATM);
      if (lJATM<=CONSTANT_DATA.mNATMS) {
        lRHO[lJATM-1] +=
          *F2D_ADDRESS(CONSTANT_DATA.mRHO_JATM_ACCU, 1, 1, CONSTANT_DATA.mMXLIST, lJ, lIATM);
      }
    }
    if (threadIdx.x==0) {
      lRHO[lIATM-1] += CONSTANT_DATA.mRHO_IATM_ACCU[lIATM-1];
    }
    __syncthreads();
  }
}

template<typename T_, unsigned int GX_, unsigned int BX_, int NS_>
__global__ void metal_ld_compute_cuda_k3b(int aI)
{
  for (int lI=aI+threadIdx.x+blockIdx.x*BX_ ; lI<=CONSTANT_DATA.mNATMS ; lI+=GX_*BX_) {
    T_ lV = CONSTANT_DATA.mRHO[lI-1];

    for (int lU=1 ; lU<NS_ ; lU++) {
      lV += (CONSTANT_DATA.mRHO + lU*CONSTANT_DATA.mMXATMS)[lI-1];
    }

    CONSTANT_DATA.mRHO[lI-1] = lV;
  }
}

extern "C" void metal_ld_compute_cuda_invoke()
{
  cudaError_t lLastError;
  start_timing_metal_ld_compute_cuda_k0();
  metal_ld_compute_cuda_k0<real,64><<<dim3(1,8*30,1), 64>>>(1);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_metal_ld_compute_cuda_k0();

  start_timing_metal_ld_compute_cuda_k1();
  metal_ld_compute_cuda_k1<real,3000,64,3><<<dim3(1,3000,1), 64, (64+3)*sizeof(real)>>>(1);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_metal_ld_compute_cuda_k1();

  start_timing_metal_ld_compute_cuda_k2();
  metal_ld_compute_cuda_k2<real, CFG_K2_THREADS_PER_BLOCK><<<dim3(1,CFG_K2_THREADBLOCKS,1), CFG_K2_THREADS_PER_BLOCK>>>(1);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_metal_ld_compute_cuda_k2();

  start_timing_metal_ld_compute_cuda_k3();
  metal_ld_compute_cuda_k3b<real,30*8,64,CFG_K2_THREADBLOCKS><<<30*8,64>>>(1);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_metal_ld_compute_cuda_k3();

  start_timing_metal_ld_compute_cuda_read();
  CUDA_SAFE_CALL(cudaMemcpy(sHD.mRHO, sCD.mRHO, sCD.mMXATMS*sizeof(real), cudaMemcpyDeviceToHost));
  stop_timing_metal_ld_compute_cuda_read();
}
