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

/* Enables host/device self load-balancing at runtime.
 *  Env var: dlpolycuda_linkcellpairs_dynamic_ratios
 */
#define CFG_DYNAMIC_RATIOS              (CFG_OVERLAP_WITH_HOST)

/* These are the threadblock dimensions of kernel k1. As usual, ensure
 * you don't run out of rmem space.
 *  Env var: dlpolycuda_linkcellpairs_k1_block_dims_{x,y,z}
 */
#define CFG_K1_BLOCK_DIMS_X             1
#define CFG_K1_BLOCK_DIMS_Y             1

#if (CFG_COMPUTE_MAJOR==2)
#  define CFG_K1_BLOCK_DIMS_Z            64
#else
#  define CFG_K1_BLOCK_DIMS_Z            32
#endif

/* Allows link_cell_pairs_cuda_pull_lists to transfer the list from the
 * device to host in a packed form. This optimisation is practical *only*
 * in the context of multiple GPUs, where contention over the data link
 * may occur.
 *  Env var: dlpolycuda_linkcellpairs_sparse_data_transfers
 */
#define CFG_SPARSE_DATA_TRANSFERS       1

/* Certain procedures declared here "make sense" if and only if (1) offloading
 * is enabled, and (2) the link_cell_pairs_cuda_initialise_ procedure has been
 * previously invoked (because it configures sCD). Ignoring this leads to
 * myriads of bugs.
 */

int* link_cell_pairs_cuda_get_pointer_is_in_valid_context() {
  static int sIsInValidcontext = 0;
  return(&sIsInValidcontext);
}

int link_cell_pairs_cuda_is_in_valid_context() {
  return(*link_cell_pairs_cuda_get_pointer_is_in_valid_context());
}

void link_cell_pairs_cuda_set_is_in_valid_context(int aIs) {
  //printf ("Calling link cell pairs\n");
	int *lIsInValidcontext = link_cell_pairs_cuda_get_pointer_is_in_valid_context();
	if (dl_poly_cuda_offload_link_cell_pairs()==dl_poly_cuda_fortran_false()) {
    printf("%s::%s: can change context when acceleration has been disabled.\n",
	   __FILE__, __FUNCTION__);
    exit(-1);
  }
  if ((*lIsInValidcontext) == aIs) {
    printf("%s::%s: can only toggle this item.\n", __FILE__, __FUNCTION__);
    exit(-1);
  }
  *lIsInValidcontext = aIs;
}
#define ASSERT_ENABLED(A)\
  if (!(*link_cell_pairs_cuda_get_pointer_is_in_valid_context())) {\
    printf("%s::%s: bogus call -- method requires link cell pairs'"\
           " accelerator to be in a valid context.\n",             \
           __FILE__, __FUNCTION__);                                \
    exit(-1);                                                      \
  }

//malysaght240112 : modification of constant_data structure
template<typename T_> struct constant_data {
  T_ *mXXX, *mYYY, *mZZZ;

  int *mNIR;
  int *mLFRZN;
  int *mAT_LIST, *mLCT_START;
  int *mLIST;
  int *mLTG;
  int *mLEXATM;

  int3 *mNIXYZ;

  T_ mRCSQ;

  int mNATMS, mMXATMS, mNCELLS, mMXLIST, mNSBCLL;
  int mNLP, mNLP3,mNLP3b, mNLX, mNLY, mMXATDM;
  int mMEGFRZ;
  int mNLX_Plus_2xNLP, mNLY_Plus_2xNLP;
  int mMXEXCL;
  int3 mNLXYZ0E, mNLXYZ1S;
};
//end_malysaght240112




__device__ __constant__ constant_data<real> CONSTANT_DATA;
static                  constant_data<real> sCD;

texture<real_tex, 1, cudaReadModeElementType> TEXTURE_DATA_XXX, TEXTURE_DATA_YYY, TEXTURE_DATA_ZZZ;

template<typename T_> struct host_data {
  int mDynamicRatios;
  dim3 mK1_BlockDims;
  int mSparseDataTransfers;

  int *mNIR;
  int *mAT_LIST, *mLCT_START, *mLIST;
  int *mDEV_ScratchBuffer;
  size_t  mDEV_ScratchBuffer_Size;

  double mPercentageOffloadedToTheDevice;
  double mPercentageOffloadedToTheDevice_NextValue;

  /* This is the host's "nl" values, possibly different than what is found
   * in the CD structure.
   */
  int3 mNLXYZ0E, mNLXYZ1S;
  int *mNIX, *mNIY, *mNIZ;
  int  mIsLBOOKTrue;
};
static host_data<real> sHD;

extern "C" void* link_cell_pairs_cuda_get_list() {
  return(sCD.mLIST);
}

extern "C" void link_cell_pairs_cuda_initialise(
    int *aNATMS, int *aMXATMS, int *aNCELLS, int *aMXLIST, int *aMXATDM, int *aNSBCLL,
    int *aNLP, int *aNLX, int *aNLY,
    int *aNLX0E, int *aNLY0E, int *aNLZ0E,
    int *aNLX1S, int *aNLY1S, int *aNLZ1S,
    real *aXXX, real *aYYY, real *aZZZ, int *aLFRZN, int *aNIR,
    int *aAT_LIST, int *aLCT_START, int *aLIST, int *aNLP3, int *aNLP3b,
    int *aNIX, int *aNIY, int *aNIZ, real *aRCSQ,
    int *aLBOOK, int *aMXEXCL, int *aLEXATM, int *aLTG, int *aMEGFRZ) {



  sHD.mDynamicRatios = dl_poly_cuda_getenv("dlpolycuda_linkcellpairs_dynamic_ratios",
                                           CFG_DYNAMIC_RATIOS);

  sHD.mK1_BlockDims  = dim3(dl_poly_cuda_getenv("dlpolycuda_linkcellpairs_k1_block_dims_x",
                                                CFG_K1_BLOCK_DIMS_X),
                            dl_poly_cuda_getenv("dlpolycuda_linkcellpairs_k1_block_dims_y",
                                                CFG_K1_BLOCK_DIMS_Y),
                            dl_poly_cuda_getenv("dlpolycuda_linkcellpairs_k1_block_dims_z",
                                                CFG_K1_BLOCK_DIMS_Z));

  sHD.mSparseDataTransfers = dl_poly_cuda_getenv("dlpolycuda_linkcellpairs_sparse_data_transfers",
                                                 CFG_SPARSE_DATA_TRANSFERS);


  static int sSetRatios = 1;
  if (sSetRatios) {
    sHD.mPercentageOffloadedToTheDevice = sHD.mPercentageOffloadedToTheDevice_NextValue = CFG_OVERLAP_WITH_HOST ? 0.95 : 1.0;
    if (sHD.mDynamicRatios) {
      sSetRatios = 0;
    }
  }
  /* When dynamic ratios is enabled, the _invoke function will compute the
   * new ratio and store it in the ..._NextValue
   */
  if (sHD.mDynamicRatios) {
    sHD.mPercentageOffloadedToTheDevice = sHD.mPercentageOffloadedToTheDevice_NextValue;

    if (sHD.mPercentageOffloadedToTheDevice==0.0 || sHD.mPercentageOffloadedToTheDevice==1.0) {
      sHD.mDynamicRatios = 0;
      printf("%s::%s: warning: disabled host/device load-balancing\n", __FILE__, __FUNCTION__);
    }
  }

  sHD.mIsLBOOKTrue = *aLBOOK == dl_poly_cuda_fortran_true();

  sHD.mDEV_ScratchBuffer_Size = 32*1024*1024;
  CUDA_SAFE_CALL(cudaMalloc(&sHD.mDEV_ScratchBuffer, sHD.mDEV_ScratchBuffer_Size));

  sHD.mAT_LIST   = aAT_LIST;
  sHD.mLCT_START = aLCT_START;
  sHD.mLIST      = aLIST;
  sHD.mNIX       = aNIX;
  sHD.mNIY       = aNIY;
  sHD.mNIZ       = aNIZ;
//malysaght300412
  sHD.mNIR 	 = aNIR;
  sCD.mNIR	 = aNIR;
//end_malysaght300412
  sCD.mRCSQ      = *aRCSQ;
  sCD.mNATMS     = *aNATMS;
  sCD.mMXATMS    = *aMXATMS;
  sCD.mMXATDM    = *aMXATDM;
  sCD.mNCELLS    = *aNCELLS;
  sCD.mMXLIST    = *aMXLIST;
  sCD.mNSBCLL    = *aNSBCLL;
  sCD.mNLP       = *aNLP;
  sCD.mNLP3      = *aNLP3;
  sCD.mNLP3b     = *aNLP3b;
  sCD.mNLX       = *aNLX;
  sCD.mNLY       = *aNLY;
  sCD.mNLX_Plus_2xNLP = sCD.mNLX+2*sCD.mNLP;
  sCD.mNLY_Plus_2xNLP = sCD.mNLY+2*sCD.mNLP;
  sCD.mMXEXCL = *aMXEXCL;
  sCD.mMEGFRZ = *aMEGFRZ;

  /* Load distribution; recall from the FORTRAN code the outmost loop:
   *    Do iz=nlz0e+1,nlz1s-1 ...
   * To cover the range correctly, the device part needs to start where
   * the host part finishes _minus_ 1
   */

  sCD.mNLXYZ0E   = make_int3(*aNLX0E, *aNLY0E, *aNLZ0E);
  sCD.mNLXYZ1S   = make_int3(*aNLX1S, *aNLY1S, *aNLZ1S);


  CUDA_SAFE_CALL(cudaMalloc(&sCD.mXXX,       sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mYYY,       sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mZZZ,       sCD.mMXATMS*sizeof(real)));
  BIND_TEXTURE_REAL_1D(TEXTURE_DATA_XXX, sCD.mXXX, sCD.mMXATMS);
  BIND_TEXTURE_REAL_1D(TEXTURE_DATA_YYY, sCD.mYYY, sCD.mMXATMS);
  BIND_TEXTURE_REAL_1D(TEXTURE_DATA_ZZZ, sCD.mZZZ, sCD.mMXATMS);

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mNIR,       2*sCD.mNLP3*sizeof(int)));  
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLFRZN,     sCD.mMXATMS*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mAT_LIST,   sCD.mMXATMS*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLCT_START, (sCD.mNCELLS+1+1)*sizeof(int)));
//malysaght240112
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLIST, (512 + (1+2+sCD.mMXLIST)*sCD.mMXATDM)*sizeof(int)));

  /* ni{x,y,z} constants; the second half of this array holds the
   * negated values. */
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mNIXYZ, 2*sCD.mNLP3b*sizeof(int3)));


  if (sHD.mIsLBOOKTrue) {
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mLTG, sCD.mMXATMS*sizeof(int)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mLEXATM, (1+sCD.mMXEXCL)*sCD.mMXATDM*sizeof(int)));
  }


  start_timing_link_cell_pairs_cuda_write();
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mXXX,       aXXX,
                            sCD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mYYY,       aYYY,
                            sCD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mZZZ,       aZZZ,
                            sCD.mMXATMS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLFRZN,     aLFRZN,
                            sCD.mMXATMS*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mAT_LIST,   aAT_LIST,
                            sCD.mMXATMS*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLCT_START, aLCT_START,
                            (sCD.mNCELLS+1+1)*sizeof(int), cudaMemcpyHostToDevice));

  if (sHD.mIsLBOOKTrue) {
    // TODO: this can be collected from the two body force's CUDA copy.
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mLTG, aLTG,
                              sCD.mMXATMS*sizeof(int), cudaMemcpyHostToDevice));

    CUDA_SAFE_CALL(cudaMemcpy(sCD.mLEXATM, aLEXATM,
			      (1+sCD.mMXEXCL)*sCD.mMXATDM*sizeof(int), cudaMemcpyHostToDevice));
  }

//  int lNLP3 = *aNLP3;
//  for (int lI=0 ; lI<lNLP3 ; lI++) {
//   sCD.mNIXYZ[lI]          = make_int3( aNIX[lI],  aNIY[lI],  aNIZ[lI]);
//   sCD.mNIXYZ[NLP3_MAX+lI] = make_int3(-aNIX[lI], -aNIY[lI], -aNIZ[lI]);
//  }


//malysaght240112
  // Create a temporary NIXYZ array - the second half stores the negated values
  int3 *lNIXYZ = (int3*) malloc(2*sCD.mNLP3b*sizeof(int3));
  for (int lI=0 ; lI < sCD.mNLP3b ; lI++) {
    lNIXYZ[lI]           = make_int3( aNIX[lI],  aNIY[lI],  aNIZ[lI]);
    lNIXYZ[sCD.mNLP3b+lI] = make_int3(-aNIX[lI], -aNIY[lI], -aNIZ[lI]);
  }
  cudaMemcpy(sCD.mNIXYZ, lNIXYZ, 2*sCD.mNLP3b*sizeof(int3), cudaMemcpyHostToDevice);
  free(lNIXYZ);
//end_malysaght240112


  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA, (void*)&sCD, sizeof(constant_data<real>)));
  stop_timing_link_cell_pairs_cuda_write();

  link_cell_pairs_cuda_set_is_in_valid_context(1);
}

extern "C" void link_cell_pairs_cuda_finalise() {
  link_cell_pairs_cuda_set_is_in_valid_context(0);
  CUDA_SAFE_CALL(cudaFree(sHD.mDEV_ScratchBuffer));
  CUDA_SAFE_CALL(cudaFree(sCD.mXXX));
  CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_XXX));
  CUDA_SAFE_CALL(cudaFree(sCD.mYYY));
  CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_YYY));
  CUDA_SAFE_CALL(cudaFree(sCD.mZZZ));
  CUDA_SAFE_CALL(cudaUnbindTexture(&TEXTURE_DATA_ZZZ));

  CUDA_SAFE_CALL(cudaFree(sCD.mNIR));
  CUDA_SAFE_CALL(cudaFree(sCD.mLFRZN));
  CUDA_SAFE_CALL(cudaFree(sCD.mAT_LIST));
  CUDA_SAFE_CALL(cudaFree(sCD.mLCT_START));
  CUDA_SAFE_CALL(cudaFree(sCD.mLIST));
//malysaght240112
  CUDA_SAFE_CALL(cudaFree(sCD.mNIXYZ));
//end_malysaght240112

  if (sHD.mIsLBOOKTrue) {
    CUDA_SAFE_CALL(cudaFree(sCD.mLTG));
    CUDA_SAFE_CALL(cudaFree(sCD.mLEXATM));
  }
}

/* @return Either a valid item that can be added to the list, or -1.
 */
template<typename T_>  __device__ int jjresult(int aJJ, int aMEGFRZ, int aLFRZN_I,
                                               T_ aXXX_I, T_ aYYY_I, T_ aZZZ_I) {

  int lJ       = CONSTANT_DATA.mAT_LIST[aJJ-1];
  int lLFRZN_J = CONSTANT_DATA.mLFRZN[lJ-1];

  T_ lXXX_J    = fetch_real_1d(TEXTURE_DATA_XXX, lJ-1);//CONSTANT_DATA.mXXX[lJ-1];
  T_ lYYY_J    = fetch_real_1d(TEXTURE_DATA_YYY, lJ-1);//CONSTANT_DATA.mYYY[lJ-1];
  T_ lZZZ_J    = fetch_real_1d(TEXTURE_DATA_ZZZ, lJ-1);//CONSTANT_DATA.mZZZ[lJ-1];

    T_ lRSQ = (lXXX_J-aXXX_I)*(lXXX_J-aXXX_I) +
              (lYYY_J-aYYY_I)*(lYYY_J-aYYY_I) +
              (lZZZ_J-aZZZ_I)*(lZZZ_J-aZZZ_I);
    if (lRSQ <= CONSTANT_DATA.mRCSQ)
      return (lJ);
  return (-1);
}

template<typename T_> __global__ void link_cell_pairs_cuda_k0() {
  /* Intialise the length fields of the list to zero:
   */
  int lIdx = 1 + blockIdx.x*blockDim.x + threadIdx.x;
  while (lIdx<=CONSTANT_DATA.mMXATDM) {
    *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0, lIdx) = (T_)0;
//malysaght240112
    *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 1, lIdx) = (T_)0;
    *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 2, lIdx) = (T_)0;
//end_malysaght240112
    lIdx += blockDim.x*gridDim.x;
  }
}


template<typename T_, unsigned int BX_, unsigned int BY_, unsigned int BZ_, int IPASS_>
__device__ void link_cell_pairs_cuda_k1_0(int3 aDims, int aIATMBegin) {

  __shared__ int shared[BZ_];

  // declamp the 1D grid to aDims:

  int3 lBlockIdx;
  lBlockIdx.x = blockIdx.x % aDims.x;
  lBlockIdx.y = (blockIdx.x % (aDims.x*aDims.y)) / aDims.x;
  lBlockIdx.z = blockIdx.x / (aDims.x*aDims.y);

  int lMEGFRZ = CONSTANT_DATA.mMEGFRZ;
  int lNLP3   = CONSTANT_DATA.mNLP3;
  int lNLP3b   = CONSTANT_DATA.mNLP3b;

  // compute i{x,y,z}, then i{x,y,z}{1,2}, from the blockIdx vector:
  int3 lIXYZ  = CONSTANT_DATA.mNLXYZ0E + 1 + lBlockIdx;

  int3 lIXYZ1 = lIXYZ - CONSTANT_DATA.mNLXYZ0E;
  int3 lIXYZ2 = lIXYZ - CONSTANT_DATA.mNLXYZ1S;


  if (!(IPASS_==1 || (lIXYZ1.x >= 1 && lIXYZ1.x <= CONSTANT_DATA.mNLP)   ||
                     (lIXYZ2.x <= -1 && lIXYZ2.x >= -CONSTANT_DATA.mNLP) ||
                     (lIXYZ1.y >= 1 && lIXYZ1.y <= CONSTANT_DATA.mNLP)   ||
                     (lIXYZ2.y <= -1 && lIXYZ2.y >= -CONSTANT_DATA.mNLP) ||
                     (lIXYZ1.z >= 1 && lIXYZ1.z <= CONSTANT_DATA.mNLP)   ||
                     (lIXYZ2.z <= -1 && lIXYZ2.z >= -CONSTANT_DATA.mNLP)))
     return;

  int  lIC        = 1 + lIXYZ.x + CONSTANT_DATA.mNLX_Plus_2xNLP*(lIXYZ.y + CONSTANT_DATA.mNLY_Plus_2xNLP*lIXYZ.z);
  int2 lII_Bounds = make_int2(CONSTANT_DATA.mLCT_START[lIC-1], CONSTANT_DATA.mLCT_START[lIC+1-1]-1);

  for (int lII = lII_Bounds.x + threadIdx.z ; lII<=lII_Bounds.y ; lII += BZ_) {
    int lI = CONSTANT_DATA.mAT_LIST[lII-1];

    if (lI>= aIATMBegin && lI <= CONSTANT_DATA.mNATMS) {
      if (threadIdx.x==0) {
//malysaght240112
        shared[threadIdx.z] = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0+2, lI);
//end_malysaght240112
      }

      if (BX_>32) __syncthreads();

      int aLFRZN_I = CONSTANT_DATA.mLFRZN[lI - 1];

      T_  aXXX_I   = fetch_real_1d(TEXTURE_DATA_XXX, lI-1);//CONSTANT_DATA.mXXX[lI - 1];
      T_  aYYY_I   = fetch_real_1d(TEXTURE_DATA_YYY, lI-1);//CONSTANT_DATA.mYYY[lI - 1];
      T_  aZZZ_I   = fetch_real_1d(TEXTURE_DATA_ZZZ, lI-1);//CONSTANT_DATA.mZZZ[lI - 1];

      int lOff = (IPASS_-1)*lNLP3b;

      for (int lKK = IPASS_ + threadIdx.y ; lKK <= CONSTANT_DATA.mNSBCLL ; lKK += BY_) {
        int3 lNIXYZ = CONSTANT_DATA.mNIXYZ[lOff+lKK-1];
        int3 lJXYZ = lIXYZ + lNIXYZ;

        if (IPASS_==1 || (lJXYZ.x <= CONSTANT_DATA.mNLXYZ0E.x) || (lJXYZ.x >= CONSTANT_DATA.mNLXYZ1S.x) ||
            (lJXYZ.y <= CONSTANT_DATA.mNLXYZ0E.y) || (lJXYZ.y >= CONSTANT_DATA.mNLXYZ1S.y) ||
            (lJXYZ.z <= CONSTANT_DATA.mNLXYZ0E.z) || (lJXYZ.z >= CONSTANT_DATA.mNLXYZ1S.z)) {


             int lJC = 1 + lJXYZ.x + CONSTANT_DATA.mNLX_Plus_2xNLP * (lJXYZ.y + CONSTANT_DATA.mNLY_Plus_2xNLP*lJXYZ.z);

             int lJ_START = (lJC != lIC) ? CONSTANT_DATA.mLCT_START[lJC-1] : (lII+1);
             int lJ_END   = CONSTANT_DATA.mLCT_START[lJC+1-1]-1;

             for (int lJJ = lJ_START + threadIdx.x ; lJJ<=lJ_END ; lJJ += BX_) {
                 int lJ     = jjresult(lJJ, lMEGFRZ, aLFRZN_I, aXXX_I, aYYY_I, aZZZ_I);

                 if (lJ != -1) {
                    int lSlot = atomicAdd(&shared[threadIdx.z], 1);
///malysaght240112
                    *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 1+2+lSlot, lI) = lJ;
///end_malysaght240112
                 }
             }
         } 
      }

      if (BX_>32) __syncthreads();

      if (threadIdx.x==0) {
//malysaght240112
          *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0+2, lI) = shared[threadIdx.z];
//end_malysaght240112
      }
    }
  }
}



template<typename T_,  unsigned int BX_, unsigned int BY_, unsigned int BZ_>
__global__ void link_cell_pairs_cuda_k1(int3 aDims, int aIATMBegin) {

  link_cell_pairs_cuda_k1_0<T_, BX_, BY_, BZ_, 1>(aDims, aIATMBegin);

  __threadfence();

  link_cell_pairs_cuda_k1_0<T_, BX_, BY_, BZ_, 2>(aDims, aIATMBegin);

}


extern "C" void link_cell_pairs_cuda_pull_lists(int aLastIteration) {
  ASSERT_ENABLED();
  int lHostATMs = 0; // default, if no offloading
  if (sHD.mPercentageOffloadedToTheDevice!=1.0) {
      lHostATMs = (int)((1.0 - sHD.mPercentageOffloadedToTheDevice)*sCD.mNATMS);
  }
  if (lHostATMs >= aLastIteration)
    return;

  if (lHostATMs==0) {
    start_timing_link_cell_pairs_cuda_read();
//malysaght240112
    CUDA_SAFE_CALL(cudaMemcpy(F2D_ADDRESS(sHD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, 1),
                              F2D_ADDRESS(sCD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, 1),
                              aLastIteration*(1+2+sCD.mMXLIST)*sizeof(int),
                              cudaMemcpyDeviceToHost));
//end_malysaght240112
    stop_timing_link_cell_pairs_cuda_read();
  } else {
    if (sHD.mSparseDataTransfers) {
      link_cell_pairs_sparseListTransfer_any(lHostATMs+1, aLastIteration - lHostATMs);
    } else {
      start_timing_link_cell_pairs_cuda_read();
//malysaght240112
      CUDA_SAFE_CALL(cudaMemcpy(F2D_ADDRESS(sHD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, lHostATMs+1),
                                F2D_ADDRESS(sCD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, lHostATMs+1),
                                (aLastIteration - lHostATMs)*(1+2+sCD.mMXLIST)*sizeof(int),
                                cudaMemcpyDeviceToHost));
//end_malysaght240112
      stop_timing_link_cell_pairs_cuda_read();
    }
  }
}

extern "C" void link_cell_pairs_cuda_push_lists(int aFirstIteration) {
  ASSERT_ENABLED();
  int lHostATMs = 0; // default, if no offloading

  if (sHD.mPercentageOffloadedToTheDevice!=1.0) {
      lHostATMs = (int)((1.0 - sHD.mPercentageOffloadedToTheDevice)*sCD.mNATMS);
  }

  if (lHostATMs <= aFirstIteration)
    return;

  start_timing_link_cell_pairs_cuda_write();

//malysaght240112
  CUDA_SAFE_CALL(cudaMemcpy(F2D_ADDRESS(sCD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, aFirstIteration),
                            F2D_ADDRESS(sHD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, aFirstIteration),
                            (lHostATMs - aFirstIteration + 1)*(1+2+sCD.mMXLIST)*sizeof(int),
                            cudaMemcpyHostToDevice));
//end_malysaght240112
  stop_timing_link_cell_pairs_cuda_write();
}


__global__ void link_cell_pairs_sparse_list_transfer_reorder(int aIATM_Begin, int aN, int *aOUT_Lengths) {

  for (int lI=blockIdx.x*blockDim.x + threadIdx.x ; lI<aN ; lI += gridDim.x*blockDim.x) {
//malysaght240112
    aOUT_Lengths[lI] = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0+2, aIATM_Begin+lI);
//end_malysaght240112
  }

}


/* Double-buffered inclusive prefix scan for host-based coordination.
 */
template<unsigned int GX_, unsigned int BX_, int W_>
__global__ void link_cell_pairs_sparse_list_transfer_gxscan(int aN, int *aOUT_Lengths, int aS, int *aOUT_LengthsTmp) {

  for (int lI=aS + blockIdx.x*BX_+threadIdx.x ; lI<aN ; lI+=GX_*BX_) {
    if (!W_) {
      aOUT_LengthsTmp[lI] = aOUT_Lengths[lI-aS];
    }
    if (W_) {
      aOUT_Lengths[lI] += aOUT_LengthsTmp[lI];
    }
  }
}

template<unsigned int BX_>
__global__ void link_cell_pairs_sparse_list_transfer_gxscan_nonstop(int aN, int *aOUT_Lengths) {

  __shared__ int shared[BX_];

  shared[threadIdx.x] = aOUT_Lengths[threadIdx.x];
  __syncthreads();

  for (int lI=aN/2, lS=1 ; lI>0 ; lI/=2, lS*=2) {
    int lItem;
    if (threadIdx.x >= lS) {
      lItem = shared[threadIdx.x] + shared[threadIdx.x - lS];
    }
    __syncthreads();

    if (threadIdx.x >= lS) {
      shared[threadIdx.x] = lItem;
    }
    __syncthreads();
  }
  aOUT_Lengths[threadIdx.x] = shared[threadIdx.x];
}

__global__ void link_cell_pairs_sparse_list_transfer_pack(int aIATM_Begin, int aN, int *aLength, int *aOUT_PackedLists) {
  for (int lI=blockIdx.y ; lI<aN ; lI+=gridDim.y) {
    int lPB_Offset = lI==0 ? 0 : aLength[lI-1];
    int lLength = aLength[lI] - lPB_Offset;

    for (int lJ=threadIdx.x ; lJ<lLength ; lJ+=blockDim.x) {

      aOUT_PackedLists[lPB_Offset + lJ] =
          *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 1+2+lJ, aIATM_Begin+lI); //malysaght240112
    }
  }
}

extern "C" void link_cell_pairs_sparseListTransfer(int aIATM_Begin, int aN);

extern "C" void link_cell_pairs_sparseListTransfer_any(int aIATM_Begin, int aN) {
  int lN = aN;
  int lOffset = 0;
  while (lN > 0) {
    int lK = 1 << 30;
    while (lK > lN) {
      lK/=2;
    }
    link_cell_pairs_sparseListTransfer(aIATM_Begin + lOffset, lK);
    lOffset += lK;
    lN -= lK;

    if (lN<=2048) { // contiguous strategy:
//malysaght240112
      CUDA_SAFE_CALL(cudaMemcpy(F2D_ADDRESS(sHD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0+2, aIATM_Begin + lOffset),
                                F2D_ADDRESS(sCD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0+2, aIATM_Begin + lOffset),
                                lN*(1+sCD.mMXLIST+2)*sizeof(int), cudaMemcpyDeviceToHost));
//end_malysaght240112
      break;
    }
  }
}

extern "C" void link_cell_pairs_sparseListTransfer(int aIATM_Begin, int aN) {
  start_timing_link_cell_pairs_sparse_list_transfer();

  assert(aN>0);
  cudaError_t lLastError;
  /* Do an exclusive scan of the lengths of the sub-list's lists and store
   * it in lDEV_Lengths. lNElements will holds the total size, i.e. the value of
   * the last summation.
   */
  int *lDEV_Lengths, lNElements=0;
  static int *lHST_Lengths;
  static int sHasHostLengthsBeenAllocated=0;
  static int sN = 0;
  if (sN < aN) {
    if (sHasHostLengthsBeenAllocated) {
      CUDA_SAFE_CALL(cudaFreeHost(lHST_Lengths));
    }
    CUDA_SAFE_CALL(cudaHostAlloc(&lHST_Lengths, aN*sizeof(int), 0));
    sN = aN;
    sHasHostLengthsBeenAllocated = 1;
  }
  CUDA_SAFE_CALL(cudaMalloc(&lDEV_Lengths, aN*sizeof(int)));


  start_timing_link_cell_pairs_sparse_list_transfer_reorder();
  link_cell_pairs_sparse_list_transfer_reorder<<<90,128>>>(aIATM_Begin, aN, lDEV_Lengths);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_link_cell_pairs_sparse_list_transfer_reorder();

  int *lDEV_LengthsTmp;
  CUDA_SAFE_CALL(cudaMalloc(&lDEV_LengthsTmp, aN*sizeof(int)));

  start_timing_link_cell_pairs_sparse_list_transfer_gxscan();

  for (int lI=aN/2, lS=1 ; lI>0 ; lI/=2, lS*=2) {
    if (aN<=512) {
      link_cell_pairs_sparse_list_transfer_gxscan_nonstop<512><<<1,512>>>(aN, lDEV_Lengths);

      CUDA_SAFE_CALL(cudaThreadSynchronize());
      lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);
      break;
    }

    link_cell_pairs_sparse_list_transfer_gxscan<30,512,0><<<30,512>>>(aN, lDEV_Lengths, lS, lDEV_LengthsTmp);

    CUDA_SAFE_CALL(cudaThreadSynchronize());
    lLastError = cudaGetLastError();
    CUT_CHECK_ERROR(lLastError);

    link_cell_pairs_sparse_list_transfer_gxscan<30,512,1><<<30,512>>>(aN, lDEV_Lengths, lS, lDEV_LengthsTmp);

    CUDA_SAFE_CALL(cudaThreadSynchronize());
    lLastError = cudaGetLastError();
    CUT_CHECK_ERROR(lLastError);
  }

  stop_timing_link_cell_pairs_sparse_list_transfer_gxscan();

  CUDA_SAFE_CALL(cudaFree(lDEV_LengthsTmp));

  CUDA_SAFE_CALL(cudaMemcpy(&lNElements, lDEV_Lengths+aN-1, sizeof(int), cudaMemcpyDeviceToHost));

  /* l{DEV,HST}_PackedList will hold the lists in a compressed format (dev/hst copies)
   * while lDEV_Lengths is the means for localy "decompressing" the list. Do the
   * packing:
   */
  int *lDEV_PackedLists;
  static int *lHST_PackedLists;
  CUDA_SAFE_CALL(cudaMalloc(&lDEV_PackedLists, lNElements*sizeof(int)));

  start_timing_link_cell_pairs_sparse_list_transfer_pack();
  link_cell_pairs_sparse_list_transfer_pack<<<dim3(1,900,1),32>>>(aIATM_Begin, aN, lDEV_Lengths, lDEV_PackedLists);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_link_cell_pairs_sparse_list_transfer_pack();

  start_timing_link_cell_pairs_sparse_list_transfer_read();
  start_timing_link_cell_pairs_cuda_read();

  static int sHasPackedListsBeenAllocated = 0;
  static int sNElements= 0 ;
  if (sNElements < lNElements) {
    if (sHasPackedListsBeenAllocated)
      CUDA_SAFE_CALL(cudaFreeHost(lHST_PackedLists));
    CUDA_SAFE_CALL(cudaHostAlloc(&lHST_PackedLists, lNElements*sizeof(int), 0));
    sNElements = lNElements;
    sHasPackedListsBeenAllocated = 1;
  }

  CUDA_SAFE_CALL(cudaMemcpy(lHST_PackedLists, lDEV_PackedLists, lNElements*sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaFree(lDEV_PackedLists));

  CUDA_SAFE_CALL(cudaMemcpy(lHST_Lengths, lDEV_Lengths, aN*sizeof(int), cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaFree(lDEV_Lengths));

  stop_timing_link_cell_pairs_cuda_read();
  stop_timing_link_cell_pairs_sparse_list_transfer_read();

  /* Unpack:
   */
  start_timing_link_cell_pairs_sparse_list_transfer_unpack();


#pragma omp parallel for private(lI)
  for (int lI=0 ; lI<aN ; lI++) {
    int lOffset = lI==0 ? 0 : lHST_Lengths[lI-1];
    int lLength = lHST_Lengths[lI] - lOffset;

//malysaght240112
    *F2D_ADDRESS(sHD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0+2, aIATM_Begin+lI) = lLength;
//end_malysaght240112

    if (lLength>0) {
//malysaght240112
      __builtin_memcpy(F2D_ADDRESS(sHD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 1+2, aIATM_Begin+lI),
	     lHST_PackedLists+lOffset,
	     lLength*sizeof(int));
    }
  }

  stop_timing_link_cell_pairs_sparse_list_transfer_unpack();
  stop_timing_link_cell_pairs_sparse_list_transfer();
}

extern "C" void link_cell_pairs_push_the_lists_back_if_necessary(int aIATM_DevBegin, int aHostATMs) {
  ASSERT_ENABLED();
  /* If the various forces are executed exclusively on the host,
   * move the lists to host memory.
   */
  if (!dl_poly_cuda_offload_tbforces() || !dl_poly_cuda_offload_link_cell_pairs_re() ) {
    start_timing_link_cell_pairs_cuda_read();
//malysaght240112
    CUDA_SAFE_CALL(cudaMemcpy(F2D_ADDRESS(sHD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, aIATM_DevBegin),
                              F2D_ADDRESS(sCD.mLIST, 0, 1, (sCD.mMXLIST+1+2), 0, aIATM_DevBegin),
                              (sCD.mNATMS-aHostATMs)*(1+2+sCD.mMXLIST)*sizeof(int),
                              cudaMemcpyDeviceToHost));
//end_malysaght240112

    stop_timing_link_cell_pairs_cuda_read();
  }
}

template<typename T_, unsigned int BX_, unsigned int BY_>
inline void link_cell_pairs_cuda_k1_switch_bz(int aBlocks, int3 aK1_GridDims, int aIATM_DevBegin) {

#define CASE_BZ(BZ) \
  case BZ: { link_cell_pairs_cuda_k1<real, BX_, BY_, BZ>		\
            <<<dim3(aBlocks,1,1), dim3(BX_,BY_,BZ)>>>(aK1_GridDims, aIATM_DevBegin); break; }

  switch (sHD.mK1_BlockDims.z) {
    CASE_BZ(512);
    CASE_BZ(256);
    CASE_BZ(128);
    CASE_BZ(64);
    CASE_BZ(32);
    CASE_BZ(16);
    CASE_BZ(8);
    CASE_BZ(4);
    CASE_BZ(2);
    CASE_BZ(1);
    default: {
      printf("%s::%s: (e) unsupported sHD.mK1_BlockDims.z==%d\n",
             __FILE__, __FUNCTION__, sHD.mK1_BlockDims.z);
      exit(-1);
    }
  }
#undef CASE_BZ
}

template<typename T_, unsigned int BX_>
inline void link_cell_pairs_cuda_k1_switch_byz(int aBlocks, int3 aK1_GridDims, int aIATM_DevBegin) {

#define CASE_BY(BY) \
  case BY: { link_cell_pairs_cuda_k1_switch_bz<real, BX_, BY>(aBlocks, aK1_GridDims, aIATM_DevBegin); break; }

  switch (sHD.mK1_BlockDims.y) {
    CASE_BY(4);
    CASE_BY(2);
    CASE_BY(1);
    default: {
      printf("%s::%s: (e) unsupported sHD.mK1_BlockDims.y==%d\n",
             __FILE__, __FUNCTION__, sHD.mK1_BlockDims.y);
      exit(-1);
    }
  }
#undef CASE_BY
}


template<typename T_> inline void link_cell_pairs_cuda_k1_switch_bxyz(int aBlocks, int3 aK1_GridDims, int aIATM_DevBegin) {

#define CASE_BX(BX) \
  case BX: { link_cell_pairs_cuda_k1_switch_byz<real, BX>(aBlocks, aK1_GridDims, aIATM_DevBegin); break; }

  switch (sHD.mK1_BlockDims.x) {
    CASE_BX(4);
    CASE_BX(2);
    CASE_BX(1);
    default: {
      printf("%s::%s: (e) unsupported sHD.mK1_BlockDims.x==%d\n",
             __FILE__, __FUNCTION__, sHD.mK1_BlockDims.x);
      exit(-1);
    }
  }
#undef CASE_BX
}

extern "C" void link_cell_pairs_cuda_invoke() {
  ASSERT_ENABLED();

  float lCE_K1_Hst_ElapsedTime = 0.0f, lCE_K1_Dev_ElapsedTime = 0.0f;

  if (sHD.mPercentageOffloadedToTheDevice==0.0) {
    int lFirstATM = 1;
    wrapper_f_link_cell_pairs_helper(&lFirstATM, &sCD.mNATMS, sHD.mNIX, sHD.mNIY, sHD.mNIZ,
                                     &sCD.mNLX, &sCD.mNLY,
                                     &sCD.mNLXYZ0E.x, &sCD.mNLXYZ0E.y, &sCD.mNLXYZ0E.z,
                                     &sCD.mNLXYZ1S.x, &sCD.mNLXYZ1S.y, &sCD.mNLXYZ1S.z,
                                     sHD.mLCT_START, sHD.mAT_LIST,
                                     &sCD.mNCELLS, &sCD.mNLP, &sCD.mNLP3, &sCD.mNSBCLL,
                                     &sCD.mMEGFRZ, &sCD.mRCSQ);
  } else {
      int lHostATMs = (int)((1.0 - sHD.mPercentageOffloadedToTheDevice)*sCD.mNATMS);
    int lIATM_HstBegin, lIATM_HstEnd, lIATM_DevBegin;
    if (sHD.mPercentageOffloadedToTheDevice==1.0) {
      lIATM_DevBegin = 1;
      lIATM_HstBegin=-1;
      lIATM_HstEnd=-1;
    } else {
      lIATM_HstBegin = 1;
      lIATM_HstEnd   = lIATM_HstBegin + lHostATMs - 1;
      lIATM_DevBegin = lIATM_HstEnd+1;
    }

    cudaEvent_t lCE_K1_DevStart, lCE_K1_DevStop;

    if (lIATM_DevBegin<=sCD.mNATMS) {
      /* This will intialise the "length" fields of the list structure to 0 so
       * that if we never touch them, the won't contain garbage.
       */
      start_timing_link_cell_pairs_cuda_k0();
      link_cell_pairs_cuda_k0<real><<<30,512>>>();
      CUDA_SAFE_CALL(cudaThreadSynchronize());
      cudaError_t lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);
      stop_timing_link_cell_pairs_cuda_k0();

      start_timing_link_cell_pairs_cuda_k1();

      CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStart));
      CUDA_SAFE_CALL(cudaEventCreate(&lCE_K1_DevStop));

      int3 lK1_GridDims = (sCD.mNLXYZ1S-1) - (sCD.mNLXYZ0E+1) +1;
      int lBlocks = lK1_GridDims.x * lK1_GridDims.y * lK1_GridDims.z;

      if (lBlocks > CFG_GRID_DIMENSION_MAX_SIZE) {
        printf("%s::%s: cannot support %d blocks (max=%d)\n",
               __FILE__, __FUNCTION__, lBlocks, CFG_GRID_DIMENSION_MAX_SIZE);
        exit(-1);
      }
      uint3 lTPB = sHD.mK1_BlockDims;//make_uint3(2,1,32);

      CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStart, 0));

      link_cell_pairs_cuda_k1_switch_bxyz<real>(lBlocks, lK1_GridDims, lIATM_DevBegin);

      CUDA_SAFE_CALL(cudaEventRecord(lCE_K1_DevStop, 0));
    }

    if (sHD.mPercentageOffloadedToTheDevice!=1.0) {
      struct timeval lTV_K1_HstStart, lTV_K1_HstStop;
      gettimeofday(&lTV_K1_HstStart, NULL);

      wrapper_f_link_cell_pairs_helper(&lIATM_HstBegin, &lIATM_HstEnd, sHD.mNIX, sHD.mNIY, sHD.mNIZ,
                                       &sCD.mNLX, &sCD.mNLY,
                                       &sCD.mNLXYZ0E.x, &sCD.mNLXYZ0E.y, &sCD.mNLXYZ0E.z,
                                       &sCD.mNLXYZ1S.x, &sCD.mNLXYZ1S.y, &sCD.mNLXYZ1S.z,
                                       sHD.mLCT_START, sHD.mAT_LIST,
                                       &sCD.mNCELLS, &sCD.mNLP, &sCD.mNLP3, &sCD.mNSBCLL,
                                       &sCD.mMEGFRZ, &sCD.mRCSQ);

      gettimeofday(&lTV_K1_HstStop, NULL);
      lCE_K1_Hst_ElapsedTime = secsfromtimeval(lTV_K1_HstStop) - secsfromtimeval(lTV_K1_HstStart);
    }

    if (lIATM_DevBegin<=sCD.mNATMS) {
      CUDA_SAFE_CALL(cudaEventSynchronize(lCE_K1_DevStop));
      CUDA_SAFE_CALL(cudaEventElapsedTime(&lCE_K1_Dev_ElapsedTime, lCE_K1_DevStart, lCE_K1_DevStop));
      lCE_K1_Dev_ElapsedTime /= 1000.0f;
      CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStart));
      CUDA_SAFE_CALL(cudaEventDestroy(lCE_K1_DevStop));

      CUDA_SAFE_CALL(cudaThreadSynchronize());
      cudaError_t lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);
      stop_timing_link_cell_pairs_cuda_k1();
    }

    /* If the various forces are executed exclusively on the host,
     * move the lists to host memory. However, if exclusions removal
     * must take place (i.e. when lbook==.true.) postpone this process
     * for link_cell_pairs_cuda_invoke_remove_exclusions_(..) invocation
     * time.
     */
//malysaght240112
    if (!sHD.mIsLBOOKTrue || !dl_poly_cuda_offload_link_cell_pairs_re()) {
      link_cell_pairs_push_the_lists_back_if_necessary(lIATM_DevBegin, lHostATMs);
    }

    if (sHD.mDynamicRatios) {
      double lPI_Dev = lCE_K1_Dev_ElapsedTime / ((double) (sCD.mNATMS - lIATM_DevBegin + 1));
      double lPI_Hst = lCE_K1_Hst_ElapsedTime / ((double) (lIATM_HstEnd - lIATM_HstBegin + 1));
      double lNewRatio = lPI_Hst / (lPI_Dev + lPI_Hst);

      sHD.mPercentageOffloadedToTheDevice_NextValue = lNewRatio;
    }
  }
}

__device__ void link_cell_pairs_cuda_invoke_remove_exclusions_search(int aJJ, int aII) {
  extern __shared__ int shared[];

}

template<unsigned int GY_, unsigned int BX_> __global__
void link_cell_pairs_cuda_invoke_remove_exclusions(int aI) {
  /* Structure:
   *   0     - the list counter
   *   1-512 - the list itself
   */
  extern __shared__ int shared[];


  for (int lI=aI+blockIdx.y ; lI<=CONSTANT_DATA.mNATMS ; lI+= GY_) {
    int lII = *F2D_ADDRESS(CONSTANT_DATA.mLEXATM, 0, 1, (CONSTANT_DATA.mMXEXCL+1), 0, lI);
    if (lII<=0)
      continue;

    /* Cache the excluded ones:
     */
    for (int lU=1+threadIdx.x ; lU<=lII ; lU+=BX_) {
      shared[1+512+lU-1] = *F2D_ADDRESS(CONSTANT_DATA.mLEXATM, 0, 1, (CONSTANT_DATA.mMXEXCL+1), lU, lI);
    }

    if (threadIdx.x==0) {
      shared[0] = 0; // reset the counter
    }

    if (BX_>32) __syncthreads();
//malysaght240112
    int lLL=*F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 2, lI);
    // note : Filling in list(-1,i)
    *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 1, lI)=lLL;


//end_malysaght240112

    for (int lKK=1+threadIdx.x ; lKK<=lLL ; lKK+=BX_) {
//malysaght240112
      int lJATM = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), lKK+2, lI);
//end_malysaght240112
      /* This acecss here is expensive; what happens is because the items in
       * lexatm use global identifiers, and "jatm" is a local identifier, a
       * translation is necessary (ltg=local-to-global).
       */
      int lJJ = CONSTANT_DATA.mLTG[lJATM-1];

      int lMatched = 0;

      /* Linear search: takes advantage of smem broadcast performance but
       * threads have to wait for others (of the same warp) to complete
       * the search.
       * Check if we do indeed have to enter the loop.
       */
      if (shared[1+512+lII-1]>=lJJ) {
	      for (int lU=1 ; lU<=lII ; lU++) {
	        int lLEXATM = shared[1+512+lU-1];

	  /* It triggers better predicaion logic to include the second test
	   * here than leaving it outside. The second OR test signals early
	   * exit from the loop (that is, as the lexatm's contents are
	   * already sorted, we know there won't be a match with the rest
	   * if this item is greater than jj).
	   */
	        if (lLEXATM>=lJJ) {
	          if (lLEXATM==lJJ)
	            lMatched = 1;
	          break;
          }
        }
      }

      if (!lMatched) {
      	shared[atomicAdd(&shared[0], 1)+1] = lJATM;
      }
    }

    if (BX_>32) __syncthreads();
    /* If the list length has changed, write the new list back:
     */
    if (shared[0]<lLL) {
      for (int lU=threadIdx.x ; lU<=shared[0] ; lU+=BX_) {
//malysaght240112
        *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), lU+2, lI)= shared[lU];
//end_malysaght240112
      }
    }
  }
}


extern "C" void link_cell_pairs_cuda_invoke_remove_exclusions() {

  ASSERT_ENABLED();

  if (sHD.mPercentageOffloadedToTheDevice==0.0) {
    int lFirstATM = 1;
    wrapper_f_link_cell_pairs_remove_exclusions_helper(&lFirstATM, &sCD.mNATMS);
  } else {
      int lHostATMs = (int)((1.0 - sHD.mPercentageOffloadedToTheDevice)*sCD.mNATMS);
    int lIATM_HstBegin, lIATM_HstEnd, lIATM_DevBegin;

    if (sHD.mPercentageOffloadedToTheDevice==1.0) {
      lIATM_DevBegin = 1;
      lIATM_HstBegin=-1;
      lIATM_HstEnd=-1;
    } else {
      lIATM_HstBegin = 1;
      lIATM_HstEnd   = lIATM_HstBegin + lHostATMs - 1;
      lIATM_DevBegin = lIATM_HstEnd+1;
    }

    if (sHD.mPercentageOffloadedToTheDevice!=0.0) {
      link_cell_pairs_cuda_invoke_remove_exclusions
          <900,64><<<dim3(1,900,1),dim3(64,1,1),(1+512+sCD.mMXEXCL)*sizeof(int)>>>(lIATM_DevBegin);
    }

    if (sHD.mPercentageOffloadedToTheDevice!=1.0) {
      wrapper_f_link_cell_pairs_remove_exclusions_helper(&lIATM_HstBegin, &lIATM_HstEnd);
    }

    if (sHD.mPercentageOffloadedToTheDevice!=0.0) {
      CUDA_SAFE_CALL(cudaThreadSynchronize());
      cudaError_t lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);
    }

    link_cell_pairs_push_the_lists_back_if_necessary(lIATM_DevBegin, lHostATMs);
  }
}

template<unsigned int BX_>
__global__ void link_cell_pairs_cuda_invoke_find_max_list_length(int aI, int aIEnd, int *aOut) {

  extern __shared__ int shared[];

  int lMax=0;
  for (int lI=aI+threadIdx.x ; lI<=aIEnd ; lI+=BX_) {
    int lNextLIMIT = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0+2, lI);
    lMax = max(lMax, lNextLIMIT);
  }
  shared[threadIdx.x] = lMax;
  __syncthreads();
  pmax<int,BX_,1>();
  if (threadIdx.x==0)
    *aOut = shared[0];
}

template<unsigned int BX_>
__global__ void link_cell_pairs_cuda_invoke_sort_atoms(int aI, int aIEnd) {

  extern __shared__ int shared[];

  for (int lI=aI+blockIdx.y ; lI<=aIEnd ; lI+=gridDim.y) {
    int lLIMIT = *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), 0+2, lI);
    if (lLIMIT<=16) // mealingless for half-warp work (or less)
      continue;

    int lNUM = next_pow2_ge(lLIMIT);

    for (int lU=1+threadIdx.x ; lU<=lNUM ; lU+=BX_) {
      if (lU<=lLIMIT)
        shared[lU-1] =
            *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), lU+2, lI);
      else
        shared[lU-1] = INT_MAX;
    }
    __syncthreads();

    unsigned int lTID = threadIdx.x;

    /* Bitonic sort based on
     * http://developer.download.nvidia.com/compute/cuda/sdk/Projects/bitonic.tar.gz
     * It has been modified to work with less threads (but still a pow of 2).
     */

    for (unsigned int lK=2 ; lK<=lNUM ; lK*=2) {
      for (unsigned int lJ=lK/2 ; lJ>0 ; lJ/=2) {
        if (threadIdx.x<lNUM) {
          unsigned int lIXJ = lTID ^ lJ;
          if (lIXJ > lTID) {
            if ((lTID & lK) == 0) {
              if (shared[lTID] > shared[lIXJ]) {
                int lTmp = shared[lTID];
                shared[lTID] = shared[lIXJ];
                shared[lIXJ] = lTmp;
              }
            } else {
              if (shared[lTID] < shared[lIXJ]) {
                int lTmp = shared[lTID];
                shared[lTID] = shared[lIXJ];
                shared[lIXJ] = lTmp;
              }
            }
          }
        }
        __syncthreads();
      }
    }

    __syncthreads();
    for (int lU=1+threadIdx.x ; lU<=lLIMIT ; lU+=BX_) {
      *F2D_ADDRESS(CONSTANT_DATA.mLIST, 0, 1, (CONSTANT_DATA.mMXLIST+1+2), lU+2, lI) = shared[lU-1];
    }
    __syncthreads();
  }
}

/* This kernel was developed to test what effect the sorting of secondary
 * atoms (in ascending order) has for two_body_forces_cuda_{k1,k2}. Sorting
 * them only yields ~2GBytes of both load/store gmem bandwidth; however, the
 * operation is too expecnsive (~ 50% of two_body_forces_cuda_k1's cost) to
 * be of any practical value.
 *
 * The kernel can be invoked right after link_cell_pairs_cuda_invoke.
 */
extern "C" void link_cell_pairs_cuda_invoke_sort_atoms() {
  ASSERT_ENABLED();
  if (sHD.mPercentageOffloadedToTheDevice!=0.0) {
      int lHostATMs = (int)((1.0 - sHD.mPercentageOffloadedToTheDevice)*sCD.mNATMS);
    int lIATM_HstBegin, lIATM_HstEnd, lIATM_DevBegin;
    if (sHD.mPercentageOffloadedToTheDevice==1.0) {
      lIATM_DevBegin = 1;
      lIATM_HstBegin=-1;
      lIATM_HstEnd=-1;
    } else {
      lIATM_HstBegin = 1;
      lIATM_HstEnd   = lIATM_HstBegin + lHostATMs - 1;
      lIATM_DevBegin = lIATM_HstEnd+1;
    }

    if (sHD.mPercentageOffloadedToTheDevice!=0.0) {
      start_timing_link_cell_pairs_cuda_sort_atoms();
      int *lDEV_Max;
      CUDA_SAFE_CALL(cudaMalloc(&lDEV_Max, sizeof(int)));

      link_cell_pairs_cuda_invoke_find_max_list_length
          <512><<<1,512,512*sizeof(int)>>>(lIATM_DevBegin, sCD.mNATMS, lDEV_Max);
      CUDA_SAFE_CALL(cudaThreadSynchronize());
      cudaError_t lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);

      int lHST_Max;
      CUDA_SAFE_CALL(cudaMemcpy(&lHST_Max, lDEV_Max, sizeof(int), cudaMemcpyDeviceToHost));
      if (lHST_Max<=0)
        return;

      lHST_Max = next_pow2_ge(lHST_Max);

      switch (lHST_Max) {
      case 512: {
          link_cell_pairs_cuda_invoke_sort_atoms
              <512><<<dim3(1,3000,1),512,512*sizeof(int)>>>(lIATM_DevBegin, sCD.mNATMS);
          break;
      }
      case 256: {
          link_cell_pairs_cuda_invoke_sort_atoms
              <256><<<dim3(1,3000,1),256,256*sizeof(int)>>>(lIATM_DevBegin, sCD.mNATMS);
          break;
      }
      case 128: {
          link_cell_pairs_cuda_invoke_sort_atoms
              <128><<<dim3(1,3000,1),128,128*sizeof(int)>>>(lIATM_DevBegin, sCD.mNATMS);
          break;
      }
      case  64: {
          link_cell_pairs_cuda_invoke_sort_atoms
              < 64><<<dim3(1,3000,1), 64, 64*sizeof(int)>>>(lIATM_DevBegin, sCD.mNATMS);
          break;
      }
      default : {
          link_cell_pairs_cuda_invoke_sort_atoms
              < 32><<<dim3(1,3000,1), 32, 32*sizeof(int)>>>(lIATM_DevBegin, sCD.mNATMS);
          break;
      }
      }

      CUDA_SAFE_CALL(cudaThreadSynchronize());
      lLastError = cudaGetLastError();
      CUT_CHECK_ERROR(lLastError);

      CUDA_SAFE_CALL(cudaFree(lDEV_Max));
      stop_timing_link_cell_pairs_cuda_sort_atoms();
    }
  }
}

