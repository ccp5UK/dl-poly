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

#include "dl_poly_cu.h"

extern "C" void start_writing(int);
extern "C" void stop_writing();

/* We use this macro where 24-bit IA ops are sufficient in order
 * to handle cases transparently.
 */
#if CFG_USE_24BIT_IAOPS
#  define IMULXX(A,B)  __mul24(A,B)
#  define IUMULXX(A,B) __umul24(A,B)
#else
#  define IMULXX(A,B)  (A)*(B)
#  define IUMULXX(A,B) (A)*(B)
#endif

template<typename T_> void bindTexture1D(texture<T_,1,cudaReadModeElementType> *aT, T_ *aGmemPointer, int aLen) {
  cudaChannelFormatDesc lCD = cudaCreateChannelDesc<T_>();
  size_t lOffset = 0;
  CUDA_SAFE_CALL(cudaBindTexture(&lOffset, aT, aGmemPointer, &lCD, ((size_t) (aLen*sizeof(T_)))));
}


#define BIND_TEXTURE_REAL_1D(TREF,GMPOINTER,LEN)			\
  do {                                              \
      cudaChannelFormatDesc lCD_;                   \
      lCD_ = cudaCreateChannelDesc<real_tex>();                         \
      size_t lOffset_ = 0;                                              \
      CUDA_SAFE_CALL(cudaBindTexture(&lOffset_, &TREF, GMPOINTER, &lCD_, \
                                     ((size_t) (LEN))*sizeof(real)));   \
  } while (0)

inline real __device__ fetch_real_1d(texture<real_tex> aTexture, int aIndex) {
  real_tex lD = tex1Dfetch(aTexture, aIndex);
#if (CFG_DOUBLE_PRECISION)
  return (__hiloint2double(lD.y, lD.x));
#else
  return (lD);
#endif
}

template<typename T_, unsigned int D_> __device__ void psum_NFS() {
    extern __shared__ T_ shared[];

    for (unsigned int lS=blockDim.x/2 ; lS>0 ; lS>>=1) {
       if (threadIdx.x < lS) {

         for (unsigned int lI=0U ; lI<D_ ; lI++) {
           shared[(lI*blockDim.x) + threadIdx.x] += shared[(lI*blockDim.x) + threadIdx.x + lS];
         }
       }
       __syncthreads();
    }
}


/* reduction(+) of D_ distinct sets and B_ being a power of 2. Results wil appear
 * in shared[lI*blockDim.x] for lI=0..D_-1. An extra __syncthreads() call after
 * the function if necessary if shared memory ius either updated or read by trheads
 * outside the first warp are to be executed. This implementation is based on the
 * NVIDIA slides on fast reductions.
 *  NOTE: Due to the complete unrolling, this kernel may increase register
 *        pressure.
 */
template<typename T_, unsigned int B_, unsigned int D_> __device__ void psum() {
    extern __shared__ T_ shared[];

#define REDATA(I)\
    if (B_ >= (I)) {\
       if (threadIdx.x < ((I)/2)) {\
          for (unsigned int lI=0U ; lI<D_ ; lI++) {\
             shared[(lI*B_) + threadIdx.x] += shared[(lI*B_) + threadIdx.x + ((I)/2)];\
          }\
          __syncthreads();\
       }\
    }

    REDATA(512); REDATA(256); REDATA(128);
#undef REDATA

    if (threadIdx.x < 32) {
#define REDATB(I)\
       if (B_ >= (I)) {\
          for (unsigned int lI=0U ; lI<D_ ; lI++) {\
             shared[(lI*B_) + threadIdx.x] += shared[(lI*B_) + threadIdx.x + ((I)/2)];\
          }\
       }
       REDATB(64); REDATB(32); REDATB(16); REDATB(8); REDATB(4); REDATB(2);
    }
#undef REDATB
}

template<typename T_, unsigned int BX_, unsigned int D_, int S_>
__device__ inline void psum_uas0(volatile T_* smem) {
  if (BX_ >= S_) {
    if (threadIdx.x < (S_/2)) {
      for (int lI=0 ; lI<D_ ; lI++) {
	smem[lI*BX_ + threadIdx.x] += smem[lI*BX_ + threadIdx.x + (S_/2)];
      }
    }
    if ((S_/2)>32)
      __syncthreads();
  }
}


template<typename T_, unsigned int BX_, unsigned int D_>
__device__ inline void psum_uas(volatile T_ *smem) {
  psum_uas0<T_, BX_, D_, 512>(smem);
  psum_uas0<T_, BX_, D_, 256>(smem);
  psum_uas0<T_, BX_, D_, 128>(smem);
  psum_uas0<T_, BX_, D_,  64>(smem);
  psum_uas0<T_, BX_, D_,  32>(smem);
  psum_uas0<T_, BX_, D_,  16>(smem);
  psum_uas0<T_, BX_, D_,   8>(smem);
  psum_uas0<T_, BX_, D_,   4>(smem);
  psum_uas0<T_, BX_, D_,   2>(smem);
}

template<typename T_, unsigned int BX_, unsigned int BY_, unsigned int D_, int S_>
__device__ inline void psum2D_uas0(volatile T_* smem, int tid) {
  if ((BX_*BY_) >= S_) {
    if (tid < (S_/2)) for (int lI=0 ; lI<D_ ; lI++) smem[lI*BX_*BY_ + tid] += smem[lI*BX_*BY_ + tid + (S_/2)];
    if ((S_/2)>32)    __syncthreads();
  }
}

template<typename T_, unsigned int BX_, unsigned int BY_, unsigned int D_>
__device__ inline void psum2D_uas(volatile T_* smem) {
  int lTI = threadIdx.y*BX_ + threadIdx.x;
  psum2D_uas0<T_, BX_, BY_, D_, 512>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_, 256>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_, 128>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_,  64>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_,  32>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_,  16>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_,   8>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_,   4>(smem, lTI);
  psum2D_uas0<T_, BX_, BY_, D_,   2>(smem, lTI);
}

template<typename T_, unsigned int BX_, unsigned int BY_, unsigned int D_> __device__ void psum2D() {
  extern __shared__ T_ shared[];

  // clamp the thread indentifier:
  int lTI = threadIdx.x + threadIdx.y*BX_;
#define REDATA(I)             \
  if ((BX_*BY_) >= (I)) {			\
    if (lTI < ((I)/2)) {      \
      for (unsigned int lI=0U ; lI<D_ ; lI++) {\
        shared[(lI*(BX_*BY_)) + lTI] += shared[(lI*(BX_*BY_)) + lTI + ((I)/2)]; \
      }\
      __syncthreads();\
    }\
  }

  REDATA(512); REDATA(256); REDATA(128);
#undef REDATA

  if (lTI < 32) {
#define REDATB(I)\
    if ((BX_*BY_) >= (I)) {			\
      for (unsigned int lI=0U ; lI<D_ ; lI++) {\
        shared[(lI*(BX_*BY_)) + lTI] += shared[(lI*(BX_*BY_)) + lTI + ((I)/2)]; \
      }\
    }
    REDATB(64); REDATB(32); REDATB(16); REDATB(8); REDATB(4); REDATB(2);
  }
#undef REDATB
}


template<typename T_, unsigned int BX_, unsigned int BY_, unsigned int BZ_, unsigned int D_> __device__ void psum3D() {
  extern __shared__ T_ shared[];
  // clamp the thread indentifier:
  int lTI = threadIdx.x + threadIdx.y*BX_ + threadIdx.z*BX_*BY_;
#define REDATA(I)\
  if ((BX_*BY_*BZ_) >= (I)) {                       \
    if (lTI < ((I)/2)) {\
      for (unsigned int lI=0U ; lI<D_ ; lI++) {\
        shared[(lI*(BX_*BY_*BZ_)) + lTI] += shared[(lI*(BX_*BY_*BZ_)) + lTI + ((I)/2)]; \
      }\
      __syncthreads();\
    }\
  }

  REDATA(512); REDATA(256); REDATA(128);
#undef REDATA

  if (lTI < 32) {
#define REDATB(I)\
    if ((BX_*BY_*BZ_) >= (I)) {                     \
      for (unsigned int lI=0U ; lI<D_ ; lI++) {\
        shared[(lI*(BX_*BY_*BZ_)) + lTI] += shared[(lI*(BX_*BY_*BZ_)) + lTI + ((I)/2)]; \
      }\
    }
    REDATB(64); REDATB(32); REDATB(16); REDATB(8); REDATB(4); REDATB(2);
  }
#undef REDATB
}



/* psum variant where each D_ loop ishandled by threads from the y-th block dim. It
 * is required that blockDim.y==D_
 */
template<typename T_, unsigned int B_, unsigned int D_> __device__ void psumye() {
    extern __shared__ T_ shared[];
#undef REDATA
#define REDATA(I)\
    if (B_ >= (I)) {\
       if (threadIdx.x < ((I)/2) && threadIdx.y<D_) {\
          shared[(threadIdx.y*B_) + threadIdx.x] += shared[(threadIdx.y*B_) + threadIdx.x + ((I)/2)];\
          __syncthreads();\
       }\
    }

    REDATA(512); REDATA(256); REDATA(128);


    if (threadIdx.x < 32 && threadIdx.y<D_) {
#undef REDATB
#define REDATB(I)\
       if (B_ >= (I)) {\
          shared[(threadIdx.y*B_) + threadIdx.x] += shared[(threadIdx.y*B_) + threadIdx.x + ((I)/2)];\
       }
       REDATB(64); REDATB(32); REDATB(16); REDATB(8); REDATB(4); REDATB(2);
    }
}

template<typename T_, unsigned int B_, unsigned int D_, unsigned int S_> __device__ void psum_full_k1() {
  extern __shared__ T_ shared[];
  if (S_<B_) {
    if ((B_/S_)<D_) {
      unsigned int lTI = threadIdx.x % S_;
      unsigned int lDO = threadIdx.x / S_;

      for (unsigned int lD=0 ; lD<D_ ; lD+=B_/S_) {
        shared[(lD+lDO)*B_ + lTI] += shared[(lD+lDO)*B_ + lTI + S_];
      }
      __syncthreads();
    }
  }
}

template<typename T_, unsigned int B_, unsigned int D_, unsigned int SO_, unsigned int S_>
__device__ void psum_full_group_mode0(unsigned int lDO, unsigned int lTI, unsigned int lATI) {
  extern __shared__ T_ shared[];
  if (S_<B_ && (B_/S_)>=D_) {
    if ((B_/S_)==D_) {
      shared[lTI] += shared[lTI + S_];
    } else {
      if ((B_/S_)>=D_ && lATI < S_) {
        shared[lTI] += shared[lTI + S_];
      }
    }
    if (S_>32)
      __syncthreads();
  }
}


template<typename T_, unsigned int B_, unsigned int D_, unsigned int S_> __device__ void psum_full_group_mode() {
  extern __shared__ T_ shared[];

  unsigned int lDO = threadIdx.x / (B_/D_);
  unsigned int lATI = threadIdx.x % (B_/D_);
  unsigned int lTI = (lDO)*B_ + (threadIdx.x % (B_/D_));

  psum_full_group_mode0<T_,B_,D_,S_,256>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,128>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,64>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,32>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,16>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,8>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,4>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,2>(lDO,lTI,lATI);
  psum_full_group_mode0<T_,B_,D_,S_,1>(lDO,lTI,lATI);
}

template<typename T_, unsigned int B_, unsigned int D_> __device__ void psum_full() {
  extern __shared__ T_ shared[];
  if (B_==1) {
    return;
  }

  psum_full_k1<T_, B_, D_, 256>();
  psum_full_k1<T_, B_, D_, 128>();
  psum_full_k1<T_, B_, D_,  64>();
  psum_full_k1<T_, B_, D_,  32>();
  psum_full_k1<T_, B_, D_,  16>();
  psum_full_k1<T_, B_, D_,   8>();
  psum_full_k1<T_, B_, D_,   4>();
  psum_full_k1<T_, B_, D_,   2>();
  psum_full_k1<T_, B_, D_,   1>();
  psum_full_group_mode<T_, B_, D_,  256>();
}



template<typename T_, unsigned int B_, unsigned int D_> __device__ void pmax() {
    extern __shared__ T_ shared[];

#define MXREDATA(I)\
    if (B_ >= (I)) {\
       if (threadIdx.x < ((I)/2)) {\
          for (unsigned int lI=0U ; lI<D_ ; lI++) {\
             shared[(lI*B_) + threadIdx.x] = max(shared[(lI*B_) + threadIdx.x], shared[(lI*B_) + threadIdx.x + ((I)/2)]);\
          }\
          __syncthreads();\
       }\
    }

    MXREDATA(512); MXREDATA(256); MXREDATA(128);


    if (threadIdx.x < 32) {
#define MXREDATB(I)\
       if (B_ >= (I)) {\
          for (unsigned int lI=0U ; lI<D_ ; lI++) {\
             shared[(lI*B_) + threadIdx.x] = max(shared[(lI*B_) + threadIdx.x], shared[(lI*B_) + threadIdx.x + ((I)/2)]);\
          }\
       }
       MXREDATB(64); MXREDATB(32); MXREDATB(16); MXREDATB(8); MXREDATB(4); MXREDATB(2);
    }
#undef MXREDATA
#undef MXREDATB
}

/* FORTRAN AINT(A, KIND) implementation. From the manual:
 *   If |A| < 1, the result is zero.
 *   If |A| >= 1, the result has a value equal to the integer whose magnitude is the largest integer that
 *   does not exceed the magnitude of A and whose sign is the same as the sign of A.
 */
template<typename T_> inline __device__ T_ aint(T_ aX) {
  T_ lABS_X = fabs(aX);
  return (lABS_X < ((T_)1) ? ((T_)0) : trunc(aX));//trunc(lABS_X));
  //return(trunc(aX));
}

/* FORTRAN ANINT(A, KIND) implementation. Fromm the manual:
 *   If A > 0, ANINT(A) = AINT(A + 0.5)
 *   If A <= 0, ANINT(A) = AINT(A - 0.5)
 * Note:
 *   The addition and subtraction of 0.5 are done in round-to-zero mode.
 */
template<typename T_> inline __device__ T_ anint(T_ aX);
template<> inline __device__ double anint(double aX) {
  return(aint(__dadd_rz(aX, (aX > 0.0) ? 0.5 : -0.5)));
}
template<> inline __device__ float anint(float aX) {
  return(aint(__fadd_rz(aX, (aX > 0.0f) ? 0.5f : -0.5f)));
}


/* FORTRAN NINT(A, KIND)
 */
template<typename T_> inline __device__ int nint(T_ aArg);
template<> inline __device__ int nint(double aArg) {
  return((int) trunc((aArg>0.0) ? (aArg+0.5) : (aArg-0.5)));
}
template<> inline __device__ int nint(float aArg) {
  return((int) aint(aArg));
}

#define F3D_ADDRESS(ABASE, AXS, AYS, AZS, AXV, AYV, AX, AY, AZ)\
  ((ABASE) + ((AZ)-(AZS))*(AXV)*(AYV) + ((AY)-(AYS))*(AXV) + ((AX)-(AXS)))

/* Calculates the address for a FORTRAN 3D array. If USE24_ equals 1, then the
 * intrinsic __mul24 is used for the offset calculation prior to scaling by the
 * size of T_. The FORTRAN array is expected to be of shape (aXS:aXV,aYS:aYV,aZS:).
 */
template<typename T_, int USE24_> inline __device__
T_* dev_f3d_address(T_ *aBase, int aXS, int aYS, int aZS, int aXV, int aYV, int aX, int aY, int aZ) {

  if (USE24_)
    return (aBase + __mul24(aXV, __mul24((aZ-aZS),aYV)+(aY-aYS)) +  (aX-aXS));
  else
    return (aBase + (aZ-aZS)*aXV*aYV + (aY-aYS)*aXV + (aX-aXS));
}

#define F2D_ADDRESS(BASE, ROWSTART, COLUMNSTART, ROWS, ROW, COLUMN)\
  ((BASE) + ((COLUMN) - (COLUMNSTART))*(ROWS) + ((ROW)-(ROWSTART)))

template<typename T_, int USE24_> inline __device__
T_* dev_f2d_address(T_ *aBase, int aXS, int aYS, int aXV, int aX, int aY) {
  if (USE24_)
    return (aBase + __mul24(aY-aYS,aXV) + (aX-aXS));
  else
    return (F2D_ADDRESS(aBase, aXS, aYS, aXV, aX, aY));
}

/* Find an integer Q >= aV, such that I is a power of 2
 */
static unsigned int __host__ __device__ next_pow2_ge(unsigned int aV) {
  unsigned int lQ = 1;
  while (lQ < aV)
        lQ *= 2;
  return (lQ);
}

static __host__ __device__ int3 operator-(int3 aA, int3 aB) {
  return(make_int3(aA.x - aB.x, aA.y - aB.y, aA.z - aB.z));
}


static __host__ __device__ int3 operator+(int3 aA, int3 aB) {
  return(make_int3(aA.x + aB.x, aA.y + aB.y, aA.z + aB.z));
}


static __host__ __device__ int3 operator+(int3 aA, uint3 aB) {
  return(make_int3(aA.x + aB.x, aA.y + aB.y, aA.z + aB.z));
}

static __host__ __device__ int3 operator+(int3 aA, dim3 aB) {
  return(make_int3(aA.x + aB.x, aA.y + aB.y, aA.z + aB.z));
}

static __host__ __device__ int3 operator+(int3 aA, int aS) {
  return(make_int3(aA.x + aS, aA.y + aS, aA.z + aS));
}

static __host__ __device__ int3 operator/(int3 aA, int aS) {
  return(make_int3(aA.x / aS, aA.y / aS, aA.z / aS));
}


static __host__ __device__ int3 operator-(int3 aA, int aS) {
  return(make_int3(aA.x - aS, aA.y - aS, aA.z - aS));
}

static __host__ __device__ int3 max(int3 aA, int3 aB) {
  return (make_int3(max(aA.x, aB.x), max(aA.y, aB.y), max(aA.z, aB.z)));
}

static __host__ __device__ int3 min(int3 aA, int3 aB) {
  return (make_int3(min(aA.x, aB.x), min(aA.y, aB.y), min(aA.z, aB.z)));
}

#if 0
/* returns the truth of aA.i >= aB.i for all i=x,y,z
 */
static __host__ __device__ bool operator>=(const int3& aA, const int3& aB) {
  return (aA.x >= aB.x && aA.y >= aB.y && aA.z >= aB.z);
}

static __host__ __device__ bool operator<=(const int3& aA, const int3& aB) {
  return(aA.x <= aB.x && aA.y <= aB.y && aA.z <= aB.z);
}

static __host__ __device__ bool operator>(const int3& aA, const int3& aB) {
  return (aA.x > aB.x && aA.y > aB.y && aA.z > aB.z);
}

static __host__ __device__ bool operator<(const int3& aA, const int3& aB) {
  return(aA.x < aB.x && aA.y < aB.y && aA.z < aB.z);
}
#endif

static dim3 convert(int3 aV) { return(dim3(aV.x, aV.y, aV.z)); }
static __host__ __device__ int3 convert(uint3 aV) { return(make_int3(aV.x, aV.y, aV.z)); }


/* When CFG_ACCURRATE==1, the following functions will prevent nvcc from using
 * the non-fused multiply add instructions.
 */
#if (!CFG_DISABLE_ACCURATE_VERSIONS && CFG_ACCURATE)
// r = x + (y*z)
template<typename T_> __device__ T_ madd(T_ aX, T_ aY, T_ aZ);
template<> __device__ float  madd(float aX,  float aY,  float aZ) { return (__fadd_rn(aX, __fmul_rn(aY,aZ))); }
template<> __device__ double madd(double aX, double aY, double  aZ) { return (__dadd_rn(aX, __dmul_rn(aY,aZ))); }

// r = x - (y*z)
template<typename T_> __device__ T_ msub(T_ aX, T_ aY, T_ aZ);
template<> __device__ float  msub(float aX,  float aY,  float aZ) { return (aX - __fmul_rn(aY,aZ)); }
template<> __device__ double msub(double aX, double aY, double  aZ) { return (aX - __dmul_rn(aY,aZ)); }

// r = x*y
template<typename T_> __device__ T_ mul(T_ aX, T_ aY);
template<> __device__ float  mul(float aX, float aY) { return(__fmul_rn(aX, aY)); }
template<> __device__ double  mul(double aX, double aY) { return(__dmul_rn(aX, aY)); }

template<typename T_> __device__ T_ add(T_ aX, T_ aY);
template<> __device__ float  add(float aX, float aY) { return(__fadd_rn(aX, aY)); }
template<> __device__ double  add(double aX, double aY) { return(__dadd_rn(aX, aY)); }
#else
template<typename T_> __device__ T_ madd(T_ aX, T_ aY, T_ aZ) { return(aX + aY*aZ); }
template<typename T_> __device__ T_ msub(T_ aX, T_ aY, T_ aZ) { return(aX - aY*aZ); }
template<typename T_> __device__ T_ mul(T_ aX, T_ aY) { return(aX*aY); }
template<typename T_> __device__ T_ add(T_ aX, T_ aY) { return(aX+aY); }
#endif

// sum of three products: computes aI1*aI2 + aI3*aI4 + aI5*aI6
template<typename T_> __device__ T_ addp3(T_ aI1, T_ aI2, T_ aI3, T_ aI4, T_ aI5, T_ aI6) {
  return (add(mul(aI1,aI2), add(mul(aI3,aI4), mul(aI5,aI6))));
}


// aX - aMO1*aMO2*aMO3
template<typename T_> __device__ T_ subm3(T_ aX, T_ aMO1, T_ aMO2, T_ aMO3) {
  return (msub(aX, mul(aMO1, aMO2), aMO3));
}

// Multiply three operands (aY,aZ & aW) and add to another (aX): aX + aY*aZ*aW
template<typename T_> __device__ T_ m3add(T_ aX, T_ aY, T_ aZ, T_ aW) {
  return (madd(aX, aY, mul(aZ, aW)));
}

/**
 * As __constant__ data cannot be passed as args to CUDA functions, and the
 * CUDA "images" needs them, this macro is provided for inlining. Insofar,
 * it's needed by both two_body_forces_cuda & constraints_shake_vv_cuda.
 * @param CONSTANT_DATA A data structure with the mCELL, mRCELL (reciprocal) and
 *                      mICELL (inverted) as accessible members.
 * @param TYPE_T The floating point data type.
 * @param IMCON  The FORTRAN images(..)'s 'imcon' argument.
 */
#define IMAGES(CONSTANT_DATA, TYPE_T, IMCON, XDF, YDF, ZDF)		\
  do {									\
    if (IMCON==2) { /* rectangular (slab) boundary conditions: */	\
      XDF = msub(XDF, CONSTANT_DATA.mCELL[0], anint(CONSTANT_DATA.mRCELL[0]*XDF)); \
      YDF = msub(YDF, CONSTANT_DATA.mCELL[4], anint(CONSTANT_DATA.mRCELL[4]*YDF)); \
      ZDF = msub(ZDF, CONSTANT_DATA.mCELL[8], anint(CONSTANT_DATA.mRCELL[8]*ZDF)); \
    }									\
    if (IMCON_==3) {							\
      TYPE_T lXSS = addp3(CONSTANT_DATA.mICELL[0], XDF, CONSTANT_DATA.mICELL[3], YDF, CONSTANT_DATA.mICELL[6], ZDF); \
      TYPE_T lYSS = addp3(CONSTANT_DATA.mICELL[1], XDF, CONSTANT_DATA.mICELL[4], YDF, CONSTANT_DATA.mICELL[7], ZDF); \
      TYPE_T lZSS = addp3(CONSTANT_DATA.mICELL[2], XDF, CONSTANT_DATA.mICELL[5], YDF, CONSTANT_DATA.mICELL[8], ZDF); \
									\
      lXSS = lXSS - anint(lXSS);					\
      lYSS = lYSS - anint(lYSS);					\
      lZSS = lZSS - anint(lZSS);					\
									\
      XDF = addp3(CONSTANT_DATA.mCELL[0], lXSS, CONSTANT_DATA.mCELL[3], lYSS, CONSTANT_DATA.mCELL[6], lZSS); \
      YDF = addp3(CONSTANT_DATA.mCELL[1], lXSS, CONSTANT_DATA.mCELL[4], lYSS, CONSTANT_DATA.mCELL[7], lZSS); \
      ZDF = addp3(CONSTANT_DATA.mCELL[2], lXSS, CONSTANT_DATA.mCELL[5], lYSS, CONSTANT_DATA.mCELL[8], lZSS); \
    }									\
  } while (0)


