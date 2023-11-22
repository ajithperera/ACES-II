/* ecp_data_cuda_free.cu: frees the CUDA pointers in ecp_data
 *
 * Written by Tom Grimes, 13-April 2010
 *
 */

#include "ecpints.h"
#include "cuda.h"

void ecp_data_cuda_free(ecp_data_type *obj)
{
   cudaFree(obj->d_radvals);
   cudaFree(obj->d_kvec_vals);
   cudaFree(obj->d_kvp_vals);
   cudaFree(obj->d_radints);
   cudaFree(obj->d_Ylm);
   cudaFree(obj->d_l_ang);
   cudaFree(obj->d_nl_ang);
   cudaFree(obj->d_local_rad);
   cudaFree(obj->d_local_ang);
   cudaFree(obj->d_nonlocal_rad);
   cudaFree(obj->d_nonlocal_ang);

   return;
};

