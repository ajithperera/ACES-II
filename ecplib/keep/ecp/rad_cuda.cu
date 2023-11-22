/* rad_cuda.cu: interface between C++ and CUDA for radial integrals
 *
 * Written by Tom Grimes, 9-April 2010
 *
 * This call is blocking.
 */

#include "ecpints.h"

#include "cuda.h"
#include "cutil_inline.h"

// kernel prototypes
__global__ void local_radial_kernel
                             (const rad_int *rad_ints, const int *rad_list, 
                              const int lower, const int upper,
                              const double *kvp_vals, double *output);

__global__ void nonlocal_radial_kernel
                                (const rad_int *rad_ints, const int *rad_list,
                                 const int lower, const int upper,
                                 const double *kvec_vals, double *output);

// this function should be ELSEWHERE
void cuda_load_kvecs (double *h_kv_vals, int nkvecs, double *h_kvp_vals, int nkvpairs)
{
   cutilSafeCall(cudaMemcpy(ecp_data.d_kvec_vals, h_kv_vals, 4*nkvecs*sizeof(double), cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(ecp_data.d_kvp_vals, h_kvp_vals, 4*nkvpairs*sizeof(double), cudaMemcpyHostToDevice));
   return;
};

// call this
void radial_cuda(int lower, int upper, double *h_output)
{
   // get lower and upper limits for local and nonlocal
   int *local_ptr = &(ecp_data.c_local_rad.front()),
       *nonlocal_ptr = &(ecp_data.c_nonlocal_rad.front()),
       l_lower = lower_bound(local_ptr,local_ptr+ecp_data.c_local_rad.size(),lower)-local_ptr,
       l_upper = lower_bound(local_ptr,local_ptr+ecp_data.c_local_rad.size(),upper)-local_ptr,
       nl_lower = lower_bound(nonlocal_ptr,nonlocal_ptr+ecp_data.c_nonlocal_rad.size(),lower)-nonlocal_ptr,
       nl_upper = lower_bound(nonlocal_ptr,nonlocal_ptr+ecp_data.c_nonlocal_rad.size(),upper)-nonlocal_ptr;

   // get dimensions
   int nlocal = l_upper-l_lower,
       nnlocal = nl_upper-nl_lower,
       nlthreadx = 4,
       nlthready = 4,
       nlthread = nlthreadx*nlthready,
       nnlthreadx = 4,
       nnlthready = 4,
       nnlthread = nnlthreadx*nnlthready,
       nlblocks = nlocal/nlthread + (nlocal%nlthread==0? 0 : 1),
       nnlblocks = nnlocal/nnlthread + (nnlocal%nnlthread==0? 0 : 1),
       bdim1l, bdim2l, bdim1n, bdim2n;
   if (nlblocks>0)
   {
      bdim1l = (int) sqrt(nlblocks);
      bdim2l = nlblocks/bdim1l + (nlblocks%bdim1l==0? 0 : 1);
   };
   if (nnlblocks>0)
   {
      bdim1n = (int) sqrt(nnlblocks);
      bdim2n = nnlblocks/bdim1n + (nnlblocks%bdim1n==0? 0 : 1);
   };
   dim3 lblocks(bdim1l,bdim2l),
        nblocks(bdim1n,bdim2n),
        lthreads(nlthreadx,nlthready),
        nthreads(nnlthreadx,nnlthready);

   // make the calls
   if (l_upper>l_lower)
   {
      local_radial_kernel<<<lblocks,lthreads>>>(ecp_data.d_radints, ecp_data.d_local_rad,
                                               l_lower, l_upper,
                                               ecp_data.d_kvp_vals, ecp_data.d_radvals);
      cutilCheckMsg("kernel launch failure");
   };
   if (nl_upper>nl_lower)
   {
      nonlocal_radial_kernel<<<nblocks,nthreads>>>(ecp_data.d_radints, ecp_data.d_nonlocal_rad,
                                                  nl_lower, nl_upper,
                                                  ecp_data.d_kvec_vals, ecp_data.d_radvals);
      cutilCheckMsg("kernel launch failure");
   };

   // copy information back to host
   cutilSafeCall(cudaMemcpy(h_output+lower, ecp_data.d_radvals+lower, (upper-lower)*sizeof(double),
                            cudaMemcpyDeviceToHost));

   return;
};

