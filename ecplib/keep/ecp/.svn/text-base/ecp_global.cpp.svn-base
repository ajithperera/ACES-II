/* ecp_global.cpp: implements the global data class
 *
 * Written by Tom Grimes, 13-April 2010
 *
 */

#include "ecpints.h"
#ifdef USING_CUDA
void ecp_data_cuda_free(ecp_data_type*);
#endif

// explicit constructor
ecp_data_type::ecp_data_type()
{
   maxshell = 0;

   Ylm = NULL;
   Ylm_maxl = 0;

   local_ang = NULL;
   max_IJK = 0;

   nonlocal_ang = NULL;
   max_nl_IJK = 0;

   angvals = NULL; 
   radvals = NULL;
   screen = NULL;
   CApows = NULL;
   kvec_vals = NULL;
   kvp_vals = NULL;
   norm_klen = 0;
   norm_kvplen = 0;
   norm_anglen = 0;
   norm_radlen = 0;

   // cuda data items
   d_radvals = NULL; 
   d_kvec_vals = NULL;
   d_kvp_vals = NULL;
   d_radints = NULL;

   return;
};

// explicit destructor
ecp_data_type::~ecp_data_type()
{
   delete [] Ylm;
   delete [] local_ang;
   delete [] nonlocal_ang;
   delete [] angvals;
   delete [] radvals;
   delete [] screen;
   delete [] CApows;
   delete [] kvec_vals;
   delete [] kvp_vals;

#ifdef USING_CUDA
   ecp_data_cuda_free(this);
#endif

   return;
};

