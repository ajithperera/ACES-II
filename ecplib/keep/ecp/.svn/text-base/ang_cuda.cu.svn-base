/* ang_cuda.cu: interface between C++ and CUDA for angular integrals
 *
 * Written by Tom Grimes, 27-April 2010
 *
 * This call is blocking.
 */

#include "ecpints.h"

#include "cuda.h"
#include "cutil_inline.h"

// kernel prototypes
__global__ void local_angular_kernel
                             (const ang_struct *ang_ints, const int *ang_list,
                              const int lower, const int upper,
                              const double *kvp_vals,
                              const double *Ylm, const int max_l,
                              const double *local_tab, const int max_IJK,
                              double *output);

__global__ void nonlocal_angular_kernel
                             (const ang_struct *ang_ints, const int *ang_list,
                              const int lower, const int upper,
                              const double *kvec_vals,
                              const double *Ylm, const int max_l,
                              const double *nonlocal_tab, const int max_IJK,
                              double *output);

// call this
void angular_cuda(int lower, int upper, double *h_output)
{
   // get lower and upper limits for local and nonlocal angular ints
   int *local_ptr = &(ecp_data.c_local_ang.front()),
       *nonlocal_ptr = &(ecp_data.c_nonlocal_ang.front()),
       l_lower = lower_bound(local_ptr,local_ptr+ecp_data.c_local_ang.size(),lower)-local_ptr,
       l_upper = lower_bound(local_ptr,local_ptr+ecp_data.c_local_ang.size(),upper)-local_ptr,
       nl_lower = lower_bound(nonlocal_ptr,nonlocal_ptr+ecp_data.c_nonlocal_ang.size(),lower)-nonlocal_ptr,
       nl_upper = lower_bound(nonlocal_ptr,nonlocal_ptr+ecp_data.c_nonlocal_ang.size(),upper)-nonlocal_ptr;

   // get dimensions
   const int nlocal = l_upper-l_lower,
             nnonlocal = nl_upper-nl_lower,
             nlthreadx = 4,
             nlthready = 4,
             nlthreads = nlthreadx*nlthready,
             nnlthreadx = 16,
             nnlthready = 2,
             nlblocks = nlocal/nlthreads+(nlocal%nlthreads==0? 0 : 1),
             nnlblocks = nnonlocal/nnlthreadx+(nnonlocal%nnlthreadx==0? 0 : 1);
   int bdim1l = 0, bdim2l = 0, bdim1n = 0, bdim2n = 0;
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
   if (nlocal>0)
   {
      local_angular_kernel<<<lblocks,lthreads>>>(ecp_data.d_angints, ecp_data.d_local_ang, l_lower, l_upper,
                                                 ecp_data.d_kvp_vals,
                                                 ecp_data.d_Ylm, ecp_data.Ylm_maxl,
                                                 ecp_data.d_l_ang, ecp_data.max_IJK,
                                                 ecp_data.d_angvals);
      cutilCheckMsg("kernel launch failure");
   };
   if (nnonlocal>0)
   {
      nonlocal_angular_kernel<<<nblocks,nthreads>>>(ecp_data.d_angints, ecp_data.d_nonlocal_ang, nl_lower, nl_upper,
                                                    ecp_data.d_kvec_vals,
                                                    ecp_data.d_Ylm, ecp_data.Ylm_maxl,
                                                    ecp_data.d_nl_ang, ecp_data.max_nl_IJK,
                                                    ecp_data.d_angvals);
      cutilCheckMsg("kernel launch failure");
   };

   // copy information back to host
   cutilSafeCall(cudaMemcpy(h_output+lower, ecp_data.d_angvals+lower, (upper-lower)*sizeof(double),
                            cudaMemcpyDeviceToHost));

   return;
};

