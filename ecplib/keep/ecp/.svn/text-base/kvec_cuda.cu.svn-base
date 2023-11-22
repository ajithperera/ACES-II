/* kvec_cuda.cu: interface between C++ and CUDA for kvecs
 *
 * Written by Tom Grimes, 9-April 2010
 *
 * This call is blocking.
 */

// need k_struct
#include "ecpints.h"

#include "cuda.h"
#include "cutil_inline.h"

// kernel prototypes
__global__ void kvec_kernel
               (const double *CApows, const int dim_1, const int dim_2, const int dim_3,
                const k_struct *kvecs, const int nkvecs,
                const double *uniq_exp, 
                double *kvec_vals);
__global__ void kvp_kernel
               (const double *CApows, const int dim_1, const int dim_2, const int dim_3,
                const kvpair_struct *kvecs, const int nkvpairs,
                const double *uniq_exp,
                double *kvp_vals, double *screen);

// call this
void kvec_cuda(double *h_CApows, int nat, int maxdim,
               k_struct *h_kvecs, int nkvecs,
               kvpair_struct *h_kvpairs, int nkvpairs,
               double *h_exps, int nexp,
               double *kvec_vals,
               double *kvp_vals,
               double *scr_vals)
{
   // get device memory
   double *d_CApows, *d_exps, *d_kv_vals, *d_kvp_vals, *d_scr_vals;
   k_struct *d_kvecs;
   kvpair_struct *d_kvpairs;
   cutilSafeCall(cudaMalloc((void**)&d_CApows, 3*nat*nat*maxdim*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&d_exps, nexp*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&d_kvecs, nkvecs*sizeof(k_struct)));
   cutilSafeCall(cudaMalloc((void**)&d_kvpairs, nkvpairs*sizeof(kvpair_struct)));

   cutilSafeCall(cudaMalloc((void**)&d_kv_vals, 4*nkvecs*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&d_kvp_vals, 4*nkvpairs*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&d_scr_vals, nkvpairs*sizeof(double)));

   // move input data to device
   cutilSafeCall(cudaMemcpy(d_CApows, h_CApows, 3*nat*nat*maxdim*sizeof(double), cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(d_exps, h_exps, nexp*sizeof(double), cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(d_kvecs, h_kvecs, nkvecs*sizeof(k_struct), cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(d_kvpairs, h_kvpairs, nkvpairs*sizeof(kvpair_struct), cudaMemcpyHostToDevice));

   // get block sizes
   const int nthread = 10;
   const int kv_blocks = nkvecs/nthread + (nkvecs%nthread==0 ? 0 : 1),
             kvp_blocks = nkvpairs/nthread + (nkvpairs%nthread==0 ? 0 : 1);

   // make the calls
   kvec_kernel<<<kv_blocks, nthread>>>(d_CApows, maxdim, 3*maxdim, nat*3*maxdim, d_kvecs, nkvecs, d_exps,
                                       d_kv_vals);
   cutilCheckMsg("kernel launch failure");

   kvp_kernel<<<kvp_blocks, nthread>>>(d_CApows, maxdim, 3*maxdim, nat*3*maxdim, d_kvpairs, nkvpairs, d_exps,
                                       d_kvp_vals, d_scr_vals);
   cutilCheckMsg("kernel launch failure");

   // copy information back to host
   cutilSafeCall(cudaMemcpy(kvec_vals, d_kv_vals, 4*nkvecs*sizeof(double), cudaMemcpyDeviceToHost));
   cutilSafeCall(cudaMemcpy(kvp_vals, d_kvp_vals, 4*nkvpairs*sizeof(double), cudaMemcpyDeviceToHost));
   cutilSafeCall(cudaMemcpy(scr_vals, d_scr_vals, nkvpairs*sizeof(double), cudaMemcpyDeviceToHost));

   // free device memory
   cutilSafeCall(cudaFree(d_CApows));
   cutilSafeCall(cudaFree(d_exps));
   cutilSafeCall(cudaFree(d_kvecs));
   cutilSafeCall(cudaFree(d_kvpairs));
   cutilSafeCall(cudaFree(d_kv_vals));
   cutilSafeCall(cudaFree(d_kvp_vals));
   cutilSafeCall(cudaFree(d_scr_vals));

   return;
};

