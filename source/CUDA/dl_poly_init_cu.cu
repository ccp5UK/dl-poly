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
#include <sys/types.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <sched.h>
#include <cuda.h>
#include <cutil.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "dl_poly_common_cu.cu"

#if (PARALLEL_BUILD==1)
#  ifndef MPICH_IGNORE_CXX_SEEK
#    define MPICH_IGNORE_CXX_SEEK
#  endif
# include <mpi.h>
#endif

static int sTRUE, sFALSE;
static int sIsCudaCapable, sDeviceBoundTo;

extern "C" __host__ int dl_poly_cuda_process() {
#if PARALLEL_BUILD
  int lRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &lRank);
  return (lRank);
#else
  return (0);
#endif
}

extern "C" __host__ int dl_poly_cuda_fortran_true() { return(sTRUE); }
extern "C" __host__ int dl_poly_cuda_fortran_false() { return(sFALSE); }

extern "C" __host__ void dl_poly_cuda_allocate_pinned_default(
  int *aN, int *aUnitSize, void **aAddrOfPointer) {
  CUDA_SAFE_CALL(cudaHostAlloc(aAddrOfPointer, (*aN)*(*aUnitSize), 0));
}



// Perform some checks at initialisation to decide what's offloaded
// to the device (to trap unimplemented functionality):

static bool offload_link_cell_pairs   = true;
static bool offload_tbforces          = true;
static bool offload_constraints_shake = true;
static bool offload_metal_ld_compute  = true;
static bool offload_ewald_spme_forces = true;
static bool offload_spme_forces       = true;
static bool offload_bspgen            = true;
static bool offload_ewald_spme_forces_ccarray  = true;
static bool offload_ewald_spme_forces_cccharge = true;

extern "C" void dl_poly_cuda_offload_set(bool offload_link_cell_pairs_f90,
                                         bool offload_tbforces_f90,                  
                                         bool offload_constraints_shake_f90,         
                                         bool offload_metal_ld_compute_f90,          
                                         bool offload_ewald_spme_forces_f90,         
                                         bool offload_spme_forces_f90,  
                                         bool offload_bspgen_f90,
                                         bool offload_ewald_spme_forces_ccarray_f90, 
                                         bool offload_ewald_spme_forces_cccharge_f90)
{
  // The offload_cuda array comes from the f90 subroutine
  // dl_poly_cuda_check_offload_conditions()
  offload_link_cell_pairs            = offload_link_cell_pairs_f90;
  offload_tbforces                   = offload_tbforces_f90;
  offload_constraints_shake          = offload_constraints_shake_f90;
  offload_metal_ld_compute           = offload_metal_ld_compute_f90;
  offload_ewald_spme_forces          = offload_ewald_spme_forces_f90;
  offload_spme_forces                = offload_spme_forces_f90;
  offload_bspgen                     = offload_bspgen_f90;
  offload_ewald_spme_forces_ccarray  = offload_ewald_spme_forces_ccarray_f90;
  offload_ewald_spme_forces_cccharge = offload_ewald_spme_forces_cccharge_f90;

}

extern "C" __host__ int dl_poly_cuda_is_cuda_capable() {
   return (sIsCudaCapable ? sTRUE : sFALSE);
}

extern "C" __host__ int dl_poly_cuda_offload_tbforces() {
  if (offload_tbforces)
    return (sTRUE);
  return (sFALSE);
}

extern "C" __host__ int dl_poly_cuda_offload_link_cell_pairs() {
  if (offload_link_cell_pairs)
    return (sTRUE);
  return (sFALSE);
}

extern "C" __host__ int dl_poly_cuda_offload_constraints_shake() {
  if (offload_constraints_shake)
    return (sTRUE);
  return (sFALSE);
}

/* If .false., so will return the functions below.
 */
extern "C" __host__ int dl_poly_cuda_offload_metal_ld_compute() {
  if (offload_metal_ld_compute)
    return (sTRUE);
  return (sFALSE);
}

extern "C" __host__ int dl_poly_cuda_offload_ewald_spme_forces() {
  if (offload_ewald_spme_forces)
    return (sTRUE);
  return (sFALSE);
}

extern "C" __host__ int dl_poly_cuda_offload_ewald_spme_forces_ccarray() {
  if (dl_poly_cuda_offload_ewald_spme_forces() == sFALSE)
    return (sFALSE);
  if (offload_ewald_spme_forces_ccarray)
    return (sTRUE);
  return (sFALSE);
}


extern "C" __host__ int dl_poly_cuda_offload_ewald_spme_forces_cccharge() {
  if (dl_poly_cuda_offload_ewald_spme_forces() == sFALSE)
    return(sFALSE);
  if (offload_ewald_spme_forces_cccharge)
    return (sTRUE);
  return (sFALSE);
}

extern "C" __host__ int dl_poly_cuda_offload_spme_forces() {
  if (offload_spme_forces)
    return (sTRUE);
  return (sFALSE);
}

extern "C" __host__ int dl_poly_cuda_offload_bspgen() {
  if (offload_bspgen)
    return (sTRUE);
  return (sFALSE);
}


template<> int dl_poly_cuda_getenv(const char *aName, int aDefault) {
  char *lVStr = getenv(aName);
  if (lVStr!=NULL) {
    return (atoi(lVStr));
  } else {
    return (aDefault);
  }
}


/* This data structure will be holding buffers of pinned memory that have
 * simulation-wide lifetime. We do that as allocating pinned memory _is_
 * expensive.
 */
template<typename T_> struct host_data {
  int mMXATMS, mMXCONS;

  // constraints_shake (lfv or vv) pinned buffers:
  struct {
    int2 *mIJSo;
    int  *mIJKPos;
    T_   *mLISTOT_R;
    T_   *mPRMCON_K;
  } mConstraintsShakeBuffers;
};
static struct host_data<real> sHD;

/* Getter methgods for those buffers:
 */
extern "C" void* dl_poly_cuda_get_buffer_constraints_shake_ijso() { return(sHD.mConstraintsShakeBuffers.mIJSo); }
extern "C" void* dl_poly_cuda_get_buffer_constraints_shake_ijkpos() { return(sHD.mConstraintsShakeBuffers.mIJKPos); }
extern "C" void* dl_poly_cuda_get_buffer_constraints_shake_prmcon_k() { return(sHD.mConstraintsShakeBuffers.mPRMCON_K); }

extern "C" __host__ void dl_poly_cuda_initialise2(int *aMXATMS, int *aMXCONS){
  sHD.mMXATMS = *aMXATMS;
  sHD.mMXCONS = *aMXCONS;

  CUDA_SAFE_CALL(cudaHostAlloc(&sHD.mConstraintsShakeBuffers.mIJSo,     sHD.mMXATMS*sizeof(int2),   0));
  CUDA_SAFE_CALL(cudaHostAlloc(&sHD.mConstraintsShakeBuffers.mIJKPos,   2*sHD.mMXCONS*sizeof(int),  0));
  CUDA_SAFE_CALL(cudaHostAlloc(&sHD.mConstraintsShakeBuffers.mLISTOT_R, sHD.mMXATMS*sizeof(real),   0));
  CUDA_SAFE_CALL(cudaHostAlloc(&sHD.mConstraintsShakeBuffers.mPRMCON_K, sHD.mMXCONS*sizeof(real),   0));
}


#define CHARS_MAX 256

extern "C" __host__ void dl_poly_cuda_initialise1(int *aTRUE, int *aFALSE){

#if PARALLEL_BUILD
  int lRank, lSize, i;
  char lNodeName[CHARS_MAX];
  int lNodeNameLength;

  MPI_Comm_rank(MPI_COMM_WORLD, &lRank);
  MPI_Comm_size(MPI_COMM_WORLD, &lSize);

  if (lRank == 0) {
      printf("DL_POLY_4 CUDA port build Info: \n");
      printf("CFG_OVERLAP_WITH_HOST = %d\n", CFG_OVERLAP_WITH_HOST);
      printf("Built for compute capability = %d.%d\n\n", CFG_COMPUTE_MAJOR, CFG_COMPUTE_MINOR);
  }

  //Get the name of the node that this process is running on
  MPI_Get_processor_name(lNodeName, &lNodeNameLength);

  //Gather the lNodeNames from all MPI tasks
  char *lNodeNameRbuf;
  lNodeNameRbuf = (char*) malloc(lSize * CHARS_MAX * sizeof(char));

  MPI_Allgather(lNodeName, CHARS_MAX, MPI_CHAR, lNodeNameRbuf, CHARS_MAX, MPI_CHAR, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  //Count the total number of MPI processes running on this node
  //by comparing the node names in lNodeNameRbuf with the current node name
  int lNumRanksThisNode = 0;
  int *lRanksThisNode;

  //lRanksThisNode is a list of the global ranks running on this node
  lRanksThisNode = (int*) malloc(lSize * sizeof(int));

  for(i=0; i<lSize; i++)
  {
      if(strncmp(lNodeName, (lNodeNameRbuf + i * CHARS_MAX),
                 CHARS_MAX) == 0)
      {
          lRanksThisNode[lNumRanksThisNode] = i;
          lNumRanksThisNode++;
      }
  }

  //Create a communicator consisting of the ranks running on this node.
  MPI_Group lWorldGroup;
  MPI_Group lThisNodeGroup;
  MPI_Comm  lThisNodeComm;
  int lRankThisNode=0, lSizeThisNode=0;

  MPI_Comm_group(MPI_COMM_WORLD, &lWorldGroup);
  MPI_Group_incl(lWorldGroup, lNumRanksThisNode, lRanksThisNode, &lThisNodeGroup);
  MPI_Comm_create(MPI_COMM_WORLD, lThisNodeGroup, &lThisNodeComm);

  MPI_Comm_rank(lThisNodeComm, &lRankThisNode);
  MPI_Comm_size(lThisNodeComm, &lSizeThisNode);

  //Attach all MPI processes on this node to the available GPUs in round-robin fashion
  int lNumDevicesThisNode = 0;

  CUDA_SAFE_CALL(cudaGetDeviceCount(&lNumDevicesThisNode));

  if (lNumDevicesThisNode == 0 && lRankThisNode == 0) {
      printf("***ERROR: no CUDA-capable devices were found on node %s.\n", lNodeName);
      exit(-1);
  }

  sIsCudaCapable = 1;
  sDeviceBoundTo = lRankThisNode % lNumDevicesThisNode;
  CUDA_SAFE_CALL(cudaSetDevice(sDeviceBoundTo));

  //The rest of the routine is for printing information to the user.
  //Gather info about what CPU cores this MPI process is affined to
  int lNumCores;
  int lCore;
  int lOpenMPThreads;
  char lAffinityInfo[CHARS_MAX] = {};
  char lBuffer[CHARS_MAX] = {};
  cpu_set_t lMask;

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp single
      {
      lOpenMPThreads = omp_get_num_threads();
      }
    }
#else
  lOpenMPThreads = 1;
#endif

  sched_getaffinity(getpid(), sizeof(lMask), &lMask);
  lNumCores = sysconf(_SC_NPROCESSORS_CONF);

  for(lCore=0; lCore < lNumCores; lCore++)
  {
      sprintf(lBuffer, " %d", lCore);
      if (CPU_ISSET(lCore, &lMask))
          strcat(lAffinityInfo, lBuffer);
  }

  //Print process-device binding info to stdout
  for(i=0; i < lSize; i++)
  {
      if(i == lRank)
      {
          printf("Bound MPI process %d (pid=%d; affined to CPU(s) %s; %d OpenMP thread(s)) to device %d@%s\n",
                 lRank, (int) getpid(), lAffinityInfo, lOpenMPThreads, sDeviceBoundTo, lNodeName);
      }
      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
  }

  //Print warning message if number of MPI processes on a given node
  //is greater than the number of available devices on that node
  if(lSizeThisNode > lNumDevicesThisNode && lRankThisNode == 0)
  {
      printf("\n** WARNING: The number of MPI processes (%d) on node %s is greater than the number of devices (%d) on that node\n   Ideally, the number of MPI processes on a given node should match the number of available devices on that node.\n\n",
             lSizeThisNode, lNodeName, lNumDevicesThisNode);
  }

  free(lNodeNameRbuf);
  free(lRanksThisNode);

#else  //Non-parallel build
  printf("DL_POLY_4 CUDA port build Info: \n");
  printf("CFG_OVERLAP_WITH_HOST = %d\n", CFG_OVERLAP_WITH_HOST);
  printf("Built for compute capability = %d.%d\n\n", CFG_COMPUTE_MAJOR, CFG_COMPUTE_MINOR);

  sIsCudaCapable = 1;
  sDeviceBoundTo = 0;
  CUDA_SAFE_CALL(cudaSetDevice(sDeviceBoundTo));
  printf("Bound process to device %d.\n", sDeviceBoundTo);
#endif //PARALLEL_BUILD

  // grab some FORTRAN90 constants:
  sTRUE  = *aTRUE;
  sFALSE = *aFALSE;
}

extern "C" __host__ void dl_poly_cuda_finalise(){
  CUDA_SAFE_CALL(cudaFreeHost(sHD.mConstraintsShakeBuffers.mIJSo));
  CUDA_SAFE_CALL(cudaFreeHost(sHD.mConstraintsShakeBuffers.mIJKPos));
  CUDA_SAFE_CALL(cudaFreeHost(sHD.mConstraintsShakeBuffers.mLISTOT_R));
  CUDA_SAFE_CALL(cudaFreeHost(sHD.mConstraintsShakeBuffers.mPRMCON_K));
}



