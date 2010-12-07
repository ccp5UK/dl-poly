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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <list>
#include <cuda.h>
#include <cutil.h>

/* Some notes:
 * (1) As access to both local & shared position data is required by
 *     this MPI process, we need to copy mxatms (instead of natms)-
 *     proportionate data to the device implying some imbalance when
 *     writing to the devices.
 */

#define CFG_ACCURRATE 1
#include "dl_poly_common_cu.cu"

/**
 * CUDA acceleration will be invoked if and only if the 'ntcons' value
 * is greater or equal to this threshold.
 *   dlpolycuda_constraints_shake_ntcons_threshold
 */
#define CFG_LEAST_NTCONS_THRESHOLD 0

extern "C" int dl_poly_cuda_constraints_shake_ntcons_enough_to_offload(int *aNTCONS) {
  int lTh = dl_poly_cuda_getenv("dlpolycuda_constraints_shake_ntcons_threshold",
				CFG_LEAST_NTCONS_THRESHOLD);
  return ((*aNTCONS)>=lTh ? dl_poly_cuda_fortran_true() : dl_poly_cuda_fortran_false());
}

template<typename T_> class constraints_shake_pack_bufs {
private:
  host_pinned_buffer<int> mIL_List;
  host_pinned_buffer<T_> mXXX_Packed, mYYY_Packed, mZZZ_Packed;
  int mAllocatedLength, mNumberOfLegitItems;
public:
  __host__ inline constraints_shake_pack_bufs() {
    mAllocatedLength    = 0;
    mNumberOfLegitItems = 0;
  }

  /**
   * Modifies the buffers so that they can host aNewLength number of
   * items.
   */
  inline void setLength(int aNewLength) {
    if (mAllocatedLength < aNewLength) {
      mIL_List.realloc(0, aNewLength);
      mXXX_Packed.realloc(0, aNewLength);
      mYYY_Packed.realloc(0, aNewLength);
      mZZZ_Packed.realloc(0, aNewLength);
      mAllocatedLength = aNewLength;
    }
    mNumberOfLegitItems = aNewLength;
  }
  inline int* list() { return (mIL_List.pointer()); }
  inline T_* xxx() { return (mXXX_Packed.pointer()); }
  inline T_* yyy() { return (mYYY_Packed.pointer()); }
  inline T_* zzz() { return (mZZZ_Packed.pointer()); }
  inline int length() { return(mAllocatedLength); }
  inline int numberOfItems() { return(mNumberOfLegitItems); }
  inline int isEmpty() { return(numberOfItems()<=0); }
  __host__ inline ~constraints_shake_pack_bufs(){}
};


template<typename T_> struct constant_data {
  int mK1_GridDims_Y, mK2_GridDims_Y;
  int mNTCONS, mMXCONS, mMXATMS, mNATMS;
  T_  mCELL[9], mRCELL[9], mICELL[9];
  T_ mTSTEP2;
  int *mLSTOPT, *mLISTOT, *mLFRZN, *mLISTCON;
  T_ *mXXX, *mYYY, *mZZZ;
  T_ *mXXT, *mYYT, *mZZT;
  T_ *mDXX, *mDYY, *mDZZ;
  T_ *mDXT, *mDYT, *mDZT;
  T_ *mWEIGHT, *mDT2, *mPRMCON_K;
  T_ *mAMTI, *mAMTJ;
  T_ *mSTRCON, *mSTRCON_Final;
  T_ *mESIG;

  int2 *mIJSo;
  int  *mIJKPos;

  // "packed transfer" of local positions from device to host
  int *mIL_SendList;
  T_ *mXXX_PackedSend, *mYYY_PackedSend, *mZZZ_PackedSend;

  // "packed transfer" of shared positions from host to device
  int *mIL_RecvList;
  T_ *mXXX_PackedRecv, *mYYY_PackedRecv, *mZZZ_PackedRecv;

  // A set of multiplication factors which are set in the initialiser
  // whose values depend on whether the lfv or vv algorithm is used
  T_ mGamma_mul, mGammIJ_mul;
};

__device__ __constant__ constant_data<real> CONSTANT_DATA;
static constant_data<real>                  sCD;

template<typename T_> struct host_data {
  int mIsFirstIteration;
  int mHaveIterationsFinished;
  int mShouldCollectStrcon;
  int mIterationId;
  dim3 mK1_GridDims, mK2_GridDims;
  int mIMCON, mNLAST;
  T_ *mXXT, *mYYT, *mZZT, *mXXX, *mYYY, *mZZZ;
  T_ *mSTRCON, *mSTRCON_Final;
  int *mLSTOPT, *mLFRZN, *mLISTOT;
  T_ mTOLNCE;
  T_ *mDXT, *mDYT, *mDZT, *mDT2;

  int2 *mIJSo;
  int  *mIJKPos;

  constraints_shake_pack_bufs<T_> *mPackedSend, *mPackedRecv;

  int mHasRecvListBeenConstructed;
};

static struct host_data<real> sHD;


/* @param aIList The list of local position indices that the update_shared_units.f90
 * will need access in order to send them over to other MPI processes. The structure
 * is set by this function; if more than zero indices have been gathered, then it is
 * the caller's responsibility to release the memory (mIndices member).
 */
template<typename T_> void constraints_shake_cuda_compile_send_list
(int *aLSI, int *aLSA, int *aLISHP, int *aLASHP, int *aMOP, int *aMXBUFF, int *aNLAST) {

  std::list<int> lList;
  int lJ0=0;

  /* These are the loops found in update_shared_units.f90 in order for
   * the buffers to be filled in with local positions.
   */
  for (int lK=1 ; lK<=26 ; lK++) {
    if (aMOP[lK-1]==0) {
      for (int lJ=lJ0+1 ; lJ<=aLASHP[lK-1] ; lJ++) {
        int lM = wrapper_f_local_index(&aLISHP[lJ-1], aNLAST, aLSI, aLSA);
        if (lM<=sCD.mNATMS) {
          lList.push_back(lM);
        }
      }
      lJ0 = aLASHP[lK-1];
    }
  }

  /* Now that the length is known, convert the list into an array:
   */
  sHD.mPackedSend->setLength(lList.size());
  int lLength = sHD.mPackedSend->numberOfItems();
  if (lLength>0) {
    int *lIndicesArray = sHD.mPackedSend->list();
    int lJ=-1;
    for (std::list<int>::iterator lIter=lList.begin(); lIter!=lList.end(); ++lIter) {
      lIndicesArray[++lJ] = *lIter;
    }

    CUDA_SAFE_CALL(cudaMalloc(&sCD.mXXX_PackedSend, lLength*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mYYY_PackedSend, lLength*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mZZZ_PackedSend, lLength*sizeof(real)));
    CUDA_SAFE_CALL(cudaMalloc(&sCD.mIL_SendList,    lLength*sizeof(int)));

    // copy the send list to the device:
    CUDA_SAFE_CALL(cudaMemcpy(sCD.mIL_SendList, sHD.mPackedSend->list(),
                              lLength*sizeof(int), cudaMemcpyHostToDevice));
  }
}

template<typename T_> void constraints_shake_cuda_compile_receive_list(int *aIndices, int aN) {
  /* The list exists in the linked list form that was constructed gradually
   * as constraints_shake_cuda_add_to_receive_list_ was invoked from
   * update_shared_units.f90. We now only need to construct the indices list
   * by flattening it out (sHD.mIL_Recv)
   */

  sHD.mPackedRecv->setLength(aN);
  if (aN<=0)
    return;

  memcpy(sHD.mPackedRecv->list(), aIndices, aN*sizeof(int));

  /* (1) allocate the buffers; this we would normally be doing from within
   * the initialiser.
   */
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mXXX_PackedRecv, aN*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mYYY_PackedRecv, aN*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mZZZ_PackedRecv, aN*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mIL_RecvList, aN*sizeof(int)));

  /* (2) copy the list over:
   */
  start_timing_constraints_shake_cuda_write();
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mIL_RecvList, sHD.mPackedRecv->list(), aN*sizeof(int), cudaMemcpyHostToDevice));
  stop_timing_constraints_shake_cuda_write();

  // (3) update the constant data:
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA, (void*)&sCD, sizeof(constant_data<real>)));
}

extern "C" void constraints_shake_cuda_commit_receive_list(int *aIndices, int *aN) {
  if (sHD.mHasRecvListBeenConstructed) {
    printf("%s::%s: bogus call; the list has already been commited\n");
    exit(-1);
  }
  sHD.mHasRecvListBeenConstructed = 1;
  constraints_shake_cuda_compile_receive_list<real>(aIndices, *aN);
}

template<typename T_, unsigned int BX_> __global__
void constraints_shake_cuda_gather_hs_scatter_dv_k1(int aI, int aN) {
  for (int lI=aI+blockIdx.x*BX_+threadIdx.x ; lI<aN ; lI+=gridDim.x*BX_) {
    int lIndex = CONSTANT_DATA.mIL_RecvList[lI]-1;
    CONSTANT_DATA.mXXX[lIndex] = CONSTANT_DATA.mXXX_PackedRecv[lI];
    CONSTANT_DATA.mYYY[lIndex] = CONSTANT_DATA.mYYY_PackedRecv[lI];
    CONSTANT_DATA.mZZZ[lIndex] = CONSTANT_DATA.mZZZ_PackedRecv[lI];
  }
}


template<typename T_, unsigned int BX_> __global__
void constraints_shake_cuda_gather_dv_scatter_hs_k1(int aI, int aN) {
  for (int lI=aI+blockIdx.x*BX_+threadIdx.x ; lI<aN ; lI+=gridDim.x*BX_) {
    int lIndex = CONSTANT_DATA.mIL_SendList[lI]-1;
    CONSTANT_DATA.mXXX_PackedSend[lI] = CONSTANT_DATA.mXXX[lIndex];
    CONSTANT_DATA.mYYY_PackedSend[lI] = CONSTANT_DATA.mYYY[lIndex];
    CONSTANT_DATA.mZZZ_PackedSend[lI] = CONSTANT_DATA.mZZZ[lIndex];
  }
}

extern "C" void  constraints_shake_cuda_gather_hs_scatter_dv() {
  start_timing_constraints_shake_cuda_gather_hs_scatter_dv();
  assert(sHD.mHasRecvListBeenConstructed==1);

  if (sHD.mPackedRecv->isEmpty())
    return;

  int lLength = sHD.mPackedRecv->numberOfItems();

  // (1) pack (offline):
  for (int lI=0 ; lI<lLength ; lI++) {
    int lIndex = sHD.mPackedRecv->list()[lI]-1;
    sHD.mPackedRecv->xxx()[lI] = sHD.mXXX[lIndex];
    sHD.mPackedRecv->yyy()[lI] = sHD.mYYY[lIndex];
    sHD.mPackedRecv->zzz()[lI] = sHD.mZZZ[lIndex];
  }
  // (2) copy to the device:
  start_timing_constraints_shake_cuda_write();
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mXXX_PackedRecv, sHD.mPackedRecv->xxx(), lLength*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mYYY_PackedRecv, sHD.mPackedRecv->yyy(), lLength*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mZZZ_PackedRecv, sHD.mPackedRecv->zzz(), lLength*sizeof(real), cudaMemcpyHostToDevice));
  stop_timing_constraints_shake_cuda_write();

  // (3) unpack (online):
  constraints_shake_cuda_gather_hs_scatter_dv_k1<real, 64><<<300,64>>>(0, lLength);
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  cudaError_t lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_constraints_shake_cuda_gather_hs_scatter_dv();
}

/* These positions that are going to appear in the send buffers
 */
extern "C" void constraints_shake_cuda_gather_dv_scatter_hs() {
  start_timing_constraints_shake_cuda_gather_dv_scatter_hs();

  if (sHD.mPackedSend->numberOfItems()<=0)
    return;

  // (1) pack (online):
  constraints_shake_cuda_gather_dv_scatter_hs_k1<real, 64><<<300,64>>>(0, sHD.mPackedSend->numberOfItems());
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  cudaError_t lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);

  // (2) copy back to host:
  start_timing_constraints_shake_cuda_read();
  CUDA_SAFE_CALL(cudaMemcpy(sHD.mPackedSend->xxx(), sCD.mXXX_PackedSend,
			    sHD.mPackedSend->numberOfItems()*sizeof(real), cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(sHD.mPackedSend->yyy(), sCD.mYYY_PackedSend,
			    sHD.mPackedSend->numberOfItems()*sizeof(real), cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(sHD.mPackedSend->zzz(), sCD.mZZZ_PackedSend,
			    sHD.mPackedSend->numberOfItems()*sizeof(real), cudaMemcpyDeviceToHost));
  stop_timing_constraints_shake_cuda_read();

  // (3) unpack (offline):
  for (int lI=0 ; lI<sHD.mPackedSend->numberOfItems() ; lI++) {
    int lIndex = sHD.mPackedSend->list()[lI]-1;
    sHD.mXXX[lIndex] = sHD.mPackedSend->xxx()[lI];
    sHD.mYYY[lIndex] = sHD.mPackedSend->yyy()[lI];
    sHD.mZZZ[lIndex] = sHD.mPackedSend->zzz()[lI];
  }

  stop_timing_constraints_shake_cuda_gather_dv_scatter_hs();
}

/* @return FORTRAN true if the receive indices list has been constructed; otherwise
 * FORTRAN false. If true, then we can safe append indices to the list.
 */
extern "C" int constraints_shake_cuda_has_recv_list_been_constructed() {
  return (sHD.mHasRecvListBeenConstructed ? dl_poly_cuda_fortran_true() : dl_poly_cuda_fortran_false());
}

template<typename T_> class simplelist {
  public:
  int mPrevSlot;
  T_ *mList;

  inline void init(int aMaxSize) {
    mList     = aMaxSize>0 ? (T_*) malloc(aMaxSize*sizeof(T_)) : NULL;
    mPrevSlot = -1;
  }

  inline void push_back(T_ aItem) {
    mList[++mPrevSlot] = aItem;
  }
  inline int size() { return (mPrevSlot+1); }
  ~simplelist() {
    if (mList!=NULL) free(mList);
  }
};



template<typename T_> void
constraints_shake_cuda_initialise_install_red_struct(int *aLSTOPT, int *aLFRZN) {

  start_timing_constraints_shake_cuda_install_red_struct();

  /* (i,j) groupping-for-reductions logic: group k-indices by atom identifier so that
   * the updates as well as their corrections can be done in parallel.
   *
   * In constraints_shake.f90, the constraints forces loop updates the {xx,yy,zz}t
   * array; it is possible that an arbitrary {xx,yy,zz}t(i) might have been updated
   * multiply -- i.e. an accummulation might have happened.
   *
   * Later, in the corrections loop, each element of the xxx,yyy & zzz arrays gets
   * updated as many times as the number of accummulations during the constraint
   * forces loops.
   */

  simplelist<int> *lOcc = (simplelist<int>*) malloc(sCD.mNATMS*sizeof(simplelist<int>));
  for (int lI=0 ; lI<sCD.mNATMS ; lI++) {
    /* The maximum number of constraints is encoded in the 'listot' array and thus
     * we can use it for the lists' initialisation.
     */
    lOcc[lI].init(sHD.mLISTOT[lI]);
  }

  /* During ccforces' K1, the i/j position updates are stored in adjacent positions, i.e.
   * if k in 1..ntcons, then they will be in 2*k and 2*k+1. We scan the lstopt structure
   * sort out these indices by atom id (i & j). The total number of updates is obviously
   * ntcons*2. Once sorted, we can perform the position update reductions without rescanning
   * the lstopt structure.
   *
   * We have two options: (1) use the original "k" indices (i.e. as we iterate the
   * k loop's elements), or (2) change them. As ccforces' k1 loop needs to store
   * data in the sCD.m{XX,YY,ZZ}T buffers, if we take the first option, then the successive
   * writes (i.e. k-wise) will require a striding, as each thread needs to store
   * 2 elements (at most, for i & j). What we can do, instead, is to re-arrange the
   * store indices so that threads can store the i,js in a coallesced fashion.
   */

  for (int lK=0 ; lK<sCD.mNTCONS ; lK++) {
    int2 lIJ = ((int2*) aLSTOPT)[lK];
    int lLFRZN_I = aLFRZN[lIJ.x-1], lLFRZN_J = aLFRZN[lIJ.y-1];

    if (lIJ.x>0 && lIJ.y>0 && (lIJ.x<=sCD.mNATMS || lIJ.y<=sCD.mNATMS) && (lLFRZN_I==0 || lLFRZN_J==0)) {

      if (lLFRZN_I==0 && lIJ.x<=sCD.mNATMS) {
        lOcc[lIJ.x-1].push_back(lK);
      }

      if (lLFRZN_J==0 && lIJ.y<=sCD.mNATMS) {
        lOcc[lIJ.y-1].push_back(sCD.mNTCONS + lK);
      }
    }
  }

  /* We now need to 'flatten out' the lOcc structure so that it is transportable to the
   * device. We will create two structures: lCounts (to be writen in sCD.mIJSo) and
   * lIJKpos (to be writen in sCD.mIJKPos), such that: for each atom i (in 1..natms),
   * if lCounts[i-1].x>0, then this means that there are that many position updates for
   * atom i, whose indices in the updates structure (see prev. paragraph) start at
   * lCounts[i-1].y; that is, lCounts describes count & offset, with lCounts[0].y==0.
   */
  int2 *lCounts;
  int *lIJKPos;
  lCounts = (int2*) dl_poly_cuda_get_buffer_constraints_shake_ijso();
  lIJKPos = (int*)  dl_poly_cuda_get_buffer_constraints_shake_ijkpos();

  sHD.mIJSo = lCounts;
  sHD.mIJKPos = lIJKPos;

  int lTotalLen = 0;

  int *lIJKPosRef = lIJKPos;

  for (int lI=0 ; lI<sCD.mNATMS ; lI++) {
    int lSize = lOcc[lI].size();
    lCounts[lI] = make_int2(lSize, lTotalLen);

    for (int lJ=0 ; lJ<lSize ; lJ++, lIJKPosRef++) {
      *lIJKPosRef = lOcc[lI].mList[lJ];
    }

    lTotalLen += lSize;
    lOcc[lI].~simplelist(); // need to call the destructor explicitly
  }
  assert(lTotalLen<=(2*sCD.mNTCONS)); // this must always hold
  free(lOcc);

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mIJSo, sCD.mNATMS*sizeof(int2)));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mIJSo, lCounts, sCD.mNATMS*sizeof(int2), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mIJKPos, lTotalLen*sizeof(int)));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mIJKPos, lIJKPos, lTotalLen*sizeof(int), cudaMemcpyHostToDevice));
  stop_timing_constraints_shake_cuda_install_red_struct();
}


extern "C" void constraints_shake_cuda_initialise
(int *aNTCONS, int *aMXCONS, int *aMXATMS, int *aNATMS, int *aIMCON,
 int *aLSI, int *aLSA, int *aLISHP, int *aLASHP, int *aMOP, int *aMXBUFF, int *aNLAST,
 int *aLSTOPT, int *aLFRZN, int *aLISTCON, int *aLISTOT, real *aPRMCON, real *aWEIGHT,
 real *aDXX, real *aDYY, real *aDZZ, real *aDXT, real *aDYT, real *aDZT,
 real *aDT2, real *aTSTEP2, real *aTOLNCE, real *aCELL,
 real *aXXX, real *aYYY, real *aZZZ,
 real *aXXT, real *aYYT, real *aZZT, real *aSTRCON, bool aIS_VV) {

  start_timing_constraints_shake_cuda_initialise();

  if (aIS_VV){
    sCD.mGamma_mul  = (real) 1.0;
    sCD.mGammIJ_mul = (real) 0.5;
  }else{
    sCD.mGamma_mul  = (real) 0.5;
    sCD.mGammIJ_mul = (real) 1.0;
  }

  sCD.mNTCONS = *aNTCONS;
  sCD.mMXCONS = *aMXCONS;
  sCD.mMXATMS = *aMXATMS;
  sCD.mNATMS  = *aNATMS;
  sCD.mTSTEP2 = *aTSTEP2;
  sHD.mIMCON  = *aIMCON;
  sHD.mK1_GridDims = dim3(1,300,1);
  sHD.mK2_GridDims = dim3(1,300,1);
  sCD.mK1_GridDims_Y = sHD.mK1_GridDims.y;
  sCD.mK2_GridDims_Y = sHD.mK2_GridDims.y;
  sHD.mIsFirstIteration = 1;
  sHD.mLSTOPT = aLSTOPT;
  sHD.mLFRZN = aLFRZN;
  sHD.mLISTOT = aLISTOT;
  sHD.mNLAST = *aNLAST;
  sHD.mXXX = aXXX;
  sHD.mYYY = aYYY;
  sHD.mZZZ = aZZZ;
  sHD.mXXT = aXXT;
  sHD.mYYT = aYYT;
  sHD.mZZT = aZZT;
  sHD.mDXT = aDXT;
  sHD.mDYT = aDYT;
  sHD.mDZT = aDZT;
  sHD.mDT2 = aDT2;
  sHD.mTOLNCE = *aTOLNCE;
  sHD.mHaveIterationsFinished = 0;
  sHD.mSTRCON = aSTRCON;
  sHD.mIterationId=0;
  sHD.mShouldCollectStrcon = 0;
  sHD.mHasRecvListBeenConstructed = 0;

  sHD.mPackedSend = new constraints_shake_pack_bufs<real>();
  sHD.mPackedRecv = new constraints_shake_pack_bufs<real>();

  for (int lI=0 ; lI<9 ; lI++) {
    real lCELL        = aCELL[lI];
    sCD.mCELL[lI]     = lCELL;
    sCD.mRCELL[lI] = ((real) 1) / lCELL;
  }
  real lDummy;
  wrapper_f_invert(sCD.mCELL, sCD.mICELL, &lDummy);

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mXXX,   sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mYYY,   sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mZZZ,   sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDXT,   sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDYT,   sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDZT,   sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDT2,   sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDXX,   sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDYY,   sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mDZZ,   sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mXXT,  2*sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mYYT,  2*sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mZZT,  2*sCD.mNTCONS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mPRMCON_K, sCD.mNTCONS*sizeof(real)));

  real *lPRMCON_K = (real*)dl_poly_cuda_get_buffer_constraints_shake_prmcon_k();
  for (int lK=1 ; lK<=sCD.mNTCONS ; lK++) {
    lPRMCON_K[lK-1] = aPRMCON[(*F2D_ADDRESS(aLISTCON, 0,1, 3, 0,lK))-1];
  }

  CUDA_SAFE_CALL(cudaMalloc(&sCD.mWEIGHT, sCD.mMXATMS*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLSTOPT, 2*sCD.mNTCONS*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLFRZN,  sCD.mMXATMS*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLISTOT, sCD.mMXATMS*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mLISTCON,3*sCD.mNTCONS*sizeof(int)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mESIG, sHD.mK1_GridDims.y*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mSTRCON, sHD.mK2_GridDims.y*6*sizeof(real)));
  CUDA_SAFE_CALL(cudaMalloc(&sCD.mSTRCON_Final, 6*sizeof(real)));

  CUDA_SAFE_CALL(cudaHostAlloc(&sHD.mSTRCON_Final, 6*sizeof(real), 0));

  start_timing_constraints_shake_cuda_write();
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mXXX,     aXXX,     sHD.mNLAST*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mYYY,     aYYY,     sHD.mNLAST*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mZZZ,     aZZZ,     sHD.mNLAST*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mDXX,     aDXX,     sCD.mNTCONS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mDYY,     aDYY,     sCD.mNTCONS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mDZZ,     aDZZ,     sCD.mNTCONS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mPRMCON_K,  lPRMCON_K,  sCD.mNTCONS*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mWEIGHT,  aWEIGHT,  sHD.mNLAST*sizeof(real), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLSTOPT,  aLSTOPT,  2*sCD.mNTCONS*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLFRZN,   aLFRZN,   sHD.mNLAST*sizeof(int), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(sCD.mLISTOT,  aLISTOT,  sHD.mNLAST*sizeof(int), cudaMemcpyHostToDevice));

  /* Construct the "send-list"; these are the identifiers of the shared positions
   * that processes around us will need (if any).
   *
   * TODO: Pin down the buffers sent later over to the device and overlap with
   *       the host code below.
   */
  constraints_shake_cuda_compile_send_list<real>(aLSI, aLSA, aLISHP, aLASHP, aMOP, aMXBUFF, aNLAST);
  constraints_shake_cuda_initialise_install_red_struct<real>(aLSTOPT, aLFRZN);


  CUDA_SAFE_CALL(cudaMemcpyToSymbol(CONSTANT_DATA, (void*)&sCD, sizeof(constant_data<real>)));
  stop_timing_constraints_shake_cuda_write();

  stop_timing_constraints_shake_cuda_initialise();
}

extern "C" void constraints_shake_cuda_finalise() {
  CUDA_SAFE_CALL(cudaFree(sCD.mXXX));
  CUDA_SAFE_CALL(cudaFree(sCD.mYYY));
  CUDA_SAFE_CALL(cudaFree(sCD.mZZZ));
  CUDA_SAFE_CALL(cudaFree(sCD.mDXT));
  CUDA_SAFE_CALL(cudaFree(sCD.mDYT));
  CUDA_SAFE_CALL(cudaFree(sCD.mDZT));
  CUDA_SAFE_CALL(cudaFree(sCD.mDT2));
  CUDA_SAFE_CALL(cudaFree(sCD.mDXX));
  CUDA_SAFE_CALL(cudaFree(sCD.mDYY));
  CUDA_SAFE_CALL(cudaFree(sCD.mDZZ));
  CUDA_SAFE_CALL(cudaFree(sCD.mXXT));
  CUDA_SAFE_CALL(cudaFree(sCD.mYYT));
  CUDA_SAFE_CALL(cudaFree(sCD.mZZT));
  CUDA_SAFE_CALL(cudaFree(sCD.mPRMCON_K));
  CUDA_SAFE_CALL(cudaFree(sCD.mWEIGHT));
  CUDA_SAFE_CALL(cudaFree(sCD.mLSTOPT));
  CUDA_SAFE_CALL(cudaFree(sCD.mLFRZN));
  CUDA_SAFE_CALL(cudaFree(sCD.mLISTOT));
  CUDA_SAFE_CALL(cudaFree(sCD.mLISTCON));
  CUDA_SAFE_CALL(cudaFree(sCD.mESIG));
  CUDA_SAFE_CALL(cudaFree(sCD.mSTRCON));
  CUDA_SAFE_CALL(cudaFree(sCD.mSTRCON_Final));
  CUDA_SAFE_CALL(cudaFree(sCD.mIJKPos));
  CUDA_SAFE_CALL(cudaFree(sCD.mIJSo));

  if (!sHD.mPackedSend->isEmpty()) {
    CUDA_SAFE_CALL(cudaFree(sCD.mXXX_PackedSend));
    CUDA_SAFE_CALL(cudaFree(sCD.mYYY_PackedSend));
    CUDA_SAFE_CALL(cudaFree(sCD.mZZZ_PackedSend));
    CUDA_SAFE_CALL(cudaFree(sCD.mIL_SendList));
  }

  if (!sHD.mPackedRecv->isEmpty()) {
    CUDA_SAFE_CALL(cudaFree(sCD.mXXX_PackedRecv));
    CUDA_SAFE_CALL(cudaFree(sCD.mYYY_PackedRecv));
    CUDA_SAFE_CALL(cudaFree(sCD.mZZZ_PackedRecv));
    CUDA_SAFE_CALL(cudaFree(sCD.mIL_RecvList));
    constraints_shake_cuda_get_recvlisttmp()->clear();
  }

  CUDA_SAFE_CALL(cudaFreeHost(sHD.mSTRCON_Final));
}

template<typename T_, unsigned int BX_, int IMCON_, int ISFINALRED_>
__global__ void constraints_shake_cuda_th() {

  DECLARE_DYNAMIC_SHARED(T_);
  shared[threadIdx.x] = (T_)0;

  if (!ISFINALRED_) {
    for (int lK=1+blockIdx.y*BX_+threadIdx.x ; lK<=CONSTANT_DATA.mNTCONS ; lK+=gridDim.y*BX_) {
      int lI = *F2D_ADDRESS(CONSTANT_DATA.mLSTOPT, 1,1, 2, 1,lK);
      int lJ = *F2D_ADDRESS(CONSTANT_DATA.mLSTOPT, 1,1, 2, 2,lK);

      T_ lDXT = (T_)0, lDYT=(T_)0, lDZT=(T_)0;

      if ((lI>0 && lJ>0) && (lI<=CONSTANT_DATA.mNATMS || lJ<=CONSTANT_DATA.mNATMS)) {
        lDXT = CONSTANT_DATA.mXXX[lI-1] - CONSTANT_DATA.mXXX[lJ-1];
        lDYT = CONSTANT_DATA.mYYY[lI-1] - CONSTANT_DATA.mYYY[lJ-1];
        lDZT = CONSTANT_DATA.mZZZ[lI-1] - CONSTANT_DATA.mZZZ[lJ-1];
      }
      IMAGES(CONSTANT_DATA, T_, IMCON_, lDXT, lDYT, lDZT);

      CONSTANT_DATA.mDXT[lK-1] = lDXT;
      CONSTANT_DATA.mDYT[lK-1] = lDYT;
      CONSTANT_DATA.mDZT[lK-1] = lDZT;

      T_  lDT2     = (T_)0;
      int lLFRZN_I = CONSTANT_DATA.mLFRZN[lI-1];
      int lLFRZN_J = CONSTANT_DATA.mLFRZN[lJ-1];

      if ((lI>0 && lJ>0) && (lI<=CONSTANT_DATA.mNATMS || lJ<=CONSTANT_DATA.mNATMS) &&
          lLFRZN_I*lLFRZN_J==0) {
        lDT2 = addp3(lDXT,lDXT,lDYT,lDYT,lDZT,lDZT);
        T_ lPRMCON = CONSTANT_DATA.mPRMCON_K[lK-1];
        T_ lESIG1  = ((T_) 0.5)*fabs(msub(lDT2, lPRMCON, lPRMCON)) / lPRMCON;
        shared[threadIdx.x] = max(shared[threadIdx.x], lESIG1);
        CONSTANT_DATA.mDT2[lK-1] = lDT2;
      }

    }
  } else {
    for (int lI=threadIdx.x ; lI<CONSTANT_DATA.mK1_GridDims_Y ; lI+=BX_) {
      shared[threadIdx.x] = max(shared[threadIdx.x], CONSTANT_DATA.mESIG[lI]);
    }
  }

  __syncthreads();
  pmax<real, BX_, 1>();
  if (threadIdx.x==0) {
    CONSTANT_DATA.mESIG[blockIdx.y] = shared[0];
  }
}

template<typename T_, unsigned int BX_> inline void constraints_shake_cuda_invoke_th_switch_imcon() {
  cudaError_t lLastError;

#define INVOKE_IMCON(IMCON) \
  do {									\
    constraints_shake_cuda_th<real, BX_, IMCON, 0><<<sHD.mK1_GridDims, BX_, BX_*sizeof(real)>>>(); \
    CUDA_SAFE_CALL(cudaThreadSynchronize());                            \
    lLastError = cudaGetLastError();                                    \
    CUT_CHECK_ERROR(lLastError);                                        \
    constraints_shake_cuda_th<real, BX_, IMCON, 1><<<dim3(1,1,1), BX_, BX_*sizeof(real)>>>(); \
    CUDA_SAFE_CALL(cudaThreadSynchronize());                            \
    lLastError = cudaGetLastError();                                    \
    CUT_CHECK_ERROR(lLastError);                                        \
  } while (0)

  switch (sHD.mIMCON) {
    case 2: { INVOKE_IMCON(2); break; }
    case 3: { INVOKE_IMCON(3); break; }
    default: { printf("%s::%s: stub reached.\n", __FILE__, __FUNCTION__); exit(-1); }
  }

#undef INVOKE_IMCON
}

extern "C" void constraints_shake_cuda_invoke_th(real *aESIG) {
  start_timing_constraints_shake_cuda_k1_th();

  constraints_shake_cuda_invoke_th_switch_imcon<real, 64>();

  CUDA_SAFE_CALL(cudaMemcpy(aESIG, sCD.mESIG, sizeof(real), cudaMemcpyDeviceToHost));
  sHD.mHaveIterationsFinished = (*aESIG) < sHD.mTOLNCE;
  stop_timing_constraints_shake_cuda_k1_th();
}


template<typename T_, unsigned int BX_, int ISFIRSTITER_, int HAVEITERSFINISHED_>
__global__ void constraints_shake_cuda_invoke_bh() {
  DECLARE_DYNAMIC_SHARED(T_);

  T_ lGamma_mul  = CONSTANT_DATA.mGamma_mul;
  T_ lGammIJ_mul = CONSTANT_DATA.mGammIJ_mul;

  if (!HAVEITERSFINISHED_) {
    for (int lI=0 ; lI<6 ; lI++) {
      shared[lI*BX_ + threadIdx.x] = (T_)0;
    }

    for (int lK=1+blockIdx.y*BX_+threadIdx.x ; lK<=CONSTANT_DATA.mNTCONS ; lK+=gridDim.y*BX_) {
      int lI = *F2D_ADDRESS(CONSTANT_DATA.mLSTOPT, 1,1, 2, 1,lK);
      int lJ = *F2D_ADDRESS(CONSTANT_DATA.mLSTOPT, 1,1, 2, 2,lK);
      int lLFRZN_I = CONSTANT_DATA.mLFRZN[lI-1];
      int lLFRZN_J = CONSTANT_DATA.mLFRZN[lJ-1];

      if ((lI>0 && lJ>0) && (lI<=CONSTANT_DATA.mNATMS || lJ<=CONSTANT_DATA.mNATMS) &&
          lLFRZN_I*lLFRZN_J==0) {
        T_ lAMTI = CONSTANT_DATA.mTSTEP2 / CONSTANT_DATA.mWEIGHT[lI-1];
        T_ lAMTJ = CONSTANT_DATA.mTSTEP2 / CONSTANT_DATA.mWEIGHT[lJ-1];
        if (lLFRZN_I!=0)
          lAMTI = (T_)0;
        if (lLFRZN_J!=0)
          lAMTJ = (T_)0;

        T_ lDXX = CONSTANT_DATA.mDXX[lK-1];
        T_ lDYY = CONSTANT_DATA.mDYY[lK-1];
        T_ lDZZ = CONSTANT_DATA.mDZZ[lK-1];
        T_ lDXT = CONSTANT_DATA.mDXT[lK-1];
        T_ lDYT = CONSTANT_DATA.mDYT[lK-1];
        T_ lDZT = CONSTANT_DATA.mDZT[lK-1];
        T_ lPRMCON = CONSTANT_DATA.mPRMCON_K[lK-1];
        T_ lDT2 = addp3(lDXT,lDXT,lDYT,lDYT,lDZT,lDZT);

#define lGAMMA shared[6*BX_+threadIdx.x]	
        lGAMMA = lGamma_mul * msub(lDT2,lPRMCON,lPRMCON)/
                 ((lAMTI+lAMTJ) * addp3(lDXX, lDXT, lDYY, lDYT, lDZZ, lDZT));
	
        if (lI <= CONSTANT_DATA.mNATMS) {
#define UPDSTRCON(I,A1,A2) shared[(I)*BX_ + threadIdx.x] = m3add(shared[(I)*BX_ + threadIdx.x], lGAMMA, lD##A1, lD##A2)
          UPDSTRCON(0, XX, XX); UPDSTRCON(1, XX, YY); UPDSTRCON(2, XX, ZZ);
          UPDSTRCON(3, YY, YY); UPDSTRCON(4, YY, ZZ);
          UPDSTRCON(5, ZZ, ZZ);
#undef UPDSTRCON

          if (lLFRZN_I==0) {
            T_ lGAMMI = -lGammIJ_mul * lGAMMA * lAMTI;
            CONSTANT_DATA.mXXT[lK-1] = lDXX*lGAMMI;
            CONSTANT_DATA.mYYT[lK-1] = lDYY*lGAMMI;
            CONSTANT_DATA.mZZT[lK-1] = lDZZ*lGAMMI;
          }
        }
        if (lJ <= CONSTANT_DATA.mNATMS && lLFRZN_J==0) {
          T_ lGAMMJ = lGammIJ_mul * lGAMMA * lAMTJ;
          CONSTANT_DATA.mXXT[CONSTANT_DATA.mNTCONS+lK-1] = lDXX*lGAMMJ;
          CONSTANT_DATA.mYYT[CONSTANT_DATA.mNTCONS+lK-1] = lDYY*lGAMMJ;
          CONSTANT_DATA.mZZT[CONSTANT_DATA.mNTCONS+lK-1] = lDZZ*lGAMMJ;
        }
      }
    }
#undef lGAMMA
    __syncthreads();

#if !CFG_UNIFIED_ADDRESS_SPACE
    psum<real, BX_, 6>();
#else
    psum_uas<real, BX_, 6>(shared);
#endif

  } else {

    T_ lSum = (T_)0;
    for (int lI=threadIdx.x ; lI<CONSTANT_DATA.mK2_GridDims_Y ; lI+=BX_) {
      lSum += CONSTANT_DATA.mSTRCON[blockIdx.x*CONSTANT_DATA.mK2_GridDims_Y + lI];
    }
    shared[threadIdx.x] = lSum;
    __syncthreads();
#if !CFG_UNIFIED_ADDRESS_SPACE
    psum<real, BX_, 1>();
#else
    psum_uas<real, BX_, 1>(shared);
#endif
  }

  if (threadIdx.x==0) {
    if (!HAVEITERSFINISHED_) {
      for (int lI=0 ; lI<6 ; lI++) {
        if (ISFIRSTITER_) {
          CONSTANT_DATA.mSTRCON[lI*CONSTANT_DATA.mK2_GridDims_Y+blockIdx.y] = shared[lI*BX_];
        } else {
          CONSTANT_DATA.mSTRCON[lI*CONSTANT_DATA.mK2_GridDims_Y+blockIdx.y] += shared[lI*BX_];
        }
      }
    } else {
      CONSTANT_DATA.mSTRCON_Final[blockIdx.x] = shared[0];
    }
  }
}


extern "C" void constraints_shake_cuda_invoke_bh() {
  start_timing_constraints_shake_cuda_k1_bh();
  cudaError_t lLastError;

  if (sHD.mIsFirstIteration) {
    constraints_shake_cuda_invoke_bh<real, 64, 1, 0><<<sHD.mK2_GridDims, 64, 64*7*sizeof(real)>>>();
    sHD.mIsFirstIteration = 0;
  } else {
    constraints_shake_cuda_invoke_bh<real, 64, 0, 0><<<sHD.mK2_GridDims, 64, 64*7*sizeof(real)>>>();
  }
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_constraints_shake_cuda_k1_bh();

  sHD.mShouldCollectStrcon = 1;
}

extern "C" void constraints_shake_cuda_invoke_strcon_extraction() {
  if (!sHD.mShouldCollectStrcon) {
    return;
  }

  constraints_shake_cuda_invoke_bh<real, 64, 0, 1><<<dim3(6,1,1), 64, 64*sizeof(real)>>>();
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  cudaError_t lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);

  CUDA_SAFE_CALL(cudaMemcpy(sHD.mSTRCON_Final, sCD.mSTRCON_Final, 6*sizeof(real), cudaMemcpyDeviceToHost));

  sHD.mSTRCON[1-1] -= sHD.mSTRCON_Final[0];
  sHD.mSTRCON[2-1] -= sHD.mSTRCON_Final[1];
  sHD.mSTRCON[3-1] -= sHD.mSTRCON_Final[2];
  sHD.mSTRCON[5-1] -= sHD.mSTRCON_Final[3];
  sHD.mSTRCON[6-1] -= sHD.mSTRCON_Final[4];
  sHD.mSTRCON[9-1] -= sHD.mSTRCON_Final[5];
}

/* Apply corrections and update the position vectors sCD.m{XXX,YYY,ZZZ}
 * from sCD.m{XX,YY,ZZ}T using sCD.mIJ{So, KPos} as a guide.
 */
template<typename T_, unsigned int BX_>
__global__ void constraints_shake_cuda_correct_positions_k1() {

  for (int lI=blockIdx.y*BX_ + threadIdx.x ; lI<CONSTANT_DATA.mNATMS ; lI+=gridDim.y*BX_) {
    int2 lB = CONSTANT_DATA.mIJSo[lI];
    if (lB.x>0) {
      T_ lDL = (((T_) 1) / ((real) CONSTANT_DATA.mLISTOT[lI])) * ((T_) lB.x);
      T_ lX=(T_)0, lY=(T_)0, lZ=(T_)0;
      for (int lQ=lB.y ; lQ<(lB.y+lB.x) ; lQ++) {
        int lIndex = CONSTANT_DATA.mIJKPos[lQ];
        lX += CONSTANT_DATA.mXXT[lIndex];
        lY += CONSTANT_DATA.mYYT[lIndex];
        lZ += CONSTANT_DATA.mZZT[lIndex];
      }
      CONSTANT_DATA.mXXX[lI] = madd(CONSTANT_DATA.mXXX[lI], lX, lDL);
      CONSTANT_DATA.mYYY[lI] = madd(CONSTANT_DATA.mYYY[lI], lY, lDL);
      CONSTANT_DATA.mZZZ[lI] = madd(CONSTANT_DATA.mZZZ[lI], lZ, lDL);
    }
  }
}

extern "C" void constraints_shake_cuda_invoke_correct_positions() {
  start_timing_constraints_shake_cuda_invoke_correct_positions();
  constraints_shake_cuda_correct_positions_k1<real, 64><<<dim3(1,300,1), 64>>>();
  CUDA_SAFE_CALL(cudaThreadSynchronize());
  cudaError_t lLastError = cudaGetLastError();
  CUT_CHECK_ERROR(lLastError);
  stop_timing_constraints_shake_cuda_invoke_correct_positions();
}

extern "C" void constraints_shake_cuda_copy_local_positions_to_host() {
  if (!sHD.mShouldCollectStrcon) { // if not, positions haven't changed
    return;
  }
  CUDA_SAFE_CALL(cudaMemcpy(sHD.mXXX, sCD.mXXX, sCD.mNATMS*sizeof(real), cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(sHD.mYYY, sCD.mYYY, sCD.mNATMS*sizeof(real), cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(sHD.mZZZ, sCD.mZZZ, sCD.mNATMS*sizeof(real), cudaMemcpyDeviceToHost));
}

