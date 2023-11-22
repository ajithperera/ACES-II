/* ecp_init_cuda.cu: initialize storage on device
 *
 * Written by Tom Grimes, 13-April 2010
 *
 */

#include "ecpints.h"
#include "cuda.h"
#include "cutil_inline.h"

void ecp_init_cuda()
{
   // allocate device memory for radial ints
   int nradints=ecp_data.radints.size(),
       nkvecs=ecp_data.kvecs.size(),
       nkvp=ecp_data.kvpairs.size();
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_radints, nradints*sizeof(rad_int)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_kvec_vals, 4*nkvecs*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_kvp_vals, 4*nkvp*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_radvals, nradints*sizeof(double)));

   // transfer the radial int list
   cutilSafeCall(cudaMemcpy(ecp_data.d_radints, &(ecp_data.radints.front()),
                            nradints*sizeof(rad_int), cudaMemcpyHostToDevice));

   // get local/nonlocal ordering of radial ints
   ecp_data.c_local_rad.clear();
   ecp_data.c_nonlocal_rad.clear();
   for (int i=0; i<ecp_data.radints.size(); i++)
   {
      if (ecp_data.radints[i].lambdaP==-1)
      {
         ecp_data.c_local_rad.push_back(i);
      } else {
         ecp_data.c_nonlocal_rad.push_back(i);
      };
   };

   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_local_rad,ecp_data.c_local_rad.size()*sizeof(int)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_nonlocal_rad,ecp_data.c_nonlocal_rad.size()*sizeof(int)));
   cutilSafeCall(cudaMemcpy(ecp_data.d_local_rad, &(ecp_data.c_local_rad.front()),
                            ecp_data.c_local_rad.size()*sizeof(int), cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(ecp_data.d_nonlocal_rad, &(ecp_data.c_nonlocal_rad.front()),
                            ecp_data.c_nonlocal_rad.size()*sizeof(int), cudaMemcpyHostToDevice));

   // transfer angular integral data
   int nangints = ecp_data.angints.size(),
       max_l = ecp_data.Ylm_maxl,
       nYlm = ((max_l+1)*(max_l+2)+1)*(max_l+1)*(max_l+1)*(max_l+1),
       max_IJK = ecp_data.max_IJK,
       nloctab = ((max_l+1)*(max_l+2)+1)*(max_IJK+1)*(max_IJK+1)*(max_IJK+1),
       max_nl_IJK = ecp_data.max_nl_IJK,
       nnltab = ((max_l+1)*(max_l+2)+1)*((max_l+1)*(max_l+2)+1)*(max_nl_IJK+1)*(max_nl_IJK+1)*(max_nl_IJK+1);
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_angvals, nangints*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_angints, nangints*sizeof(ang_struct)));

   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_Ylm, nYlm*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_l_ang, nloctab*sizeof(double)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_nl_ang, nnltab*sizeof(double)));
   cutilSafeCall(cudaMemcpy(ecp_data.d_angints, &(ecp_data.angints.front()),
                 nangints*sizeof(ang_struct),cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(ecp_data.d_Ylm, ecp_data.Ylm,
                 nYlm*sizeof(double),cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(ecp_data.d_l_ang, ecp_data.local_ang,
                 nloctab*sizeof(double),cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(ecp_data.d_nl_ang, ecp_data.nonlocal_ang,
                 nnltab*sizeof(double),cudaMemcpyHostToDevice));

   // separate local and nonlocal 
   ecp_data.c_local_ang.clear();
   ecp_data.c_nonlocal_ang.clear();
   for (int i=0; i<ecp_data.angints.size(); i++)
   {
      if (ecp_data.angints[i].ang2_lambda==-1)
      {
         ecp_data.c_local_ang.push_back(i);
      } else {
         ecp_data.c_nonlocal_ang.push_back(i);
      };
   };

   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_local_ang, ecp_data.c_local_ang.size()*sizeof(int)));
   cutilSafeCall(cudaMalloc((void**)&ecp_data.d_nonlocal_ang, ecp_data.c_nonlocal_ang.size()*sizeof(int)));
   cutilSafeCall(cudaMemcpy(ecp_data.d_local_ang, &(ecp_data.c_local_ang.front()),
                 ecp_data.c_local_ang.size()*sizeof(int),cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(ecp_data.d_nonlocal_ang, &(ecp_data.c_nonlocal_ang.front()),
                 ecp_data.c_nonlocal_ang.size()*sizeof(int),cudaMemcpyHostToDevice));

   return;
};

