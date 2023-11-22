/* kvec_kernel.cu: implements the CUDA kernel for computing k-vectors
 *
 * Written by Tom Grimes, 9-April 2010
 *
 */

// need to get k_struct
#include "ecpints.h"

#include "stdio.h"

// compute k-vectors
__global__ void kvec_kernel
               (const double *CApows, const int dim_1, const int dim_2, const int dim_3,
                const k_struct *kvecs, const int nkvecs,
                const double *uniq_exp,
                double *kvec_vals)
{
   int i=blockIdx.x*blockDim.x + threadIdx.x,
       xyz;
   double klen, vec[3];

   // limit to number of data items
   if (i>=nkvecs) return;

   if (kvecs[i].a!=kvecs[i].c)
   {
      klen=0.;
      #pragma unroll 3
      for (xyz=0; xyz<3; xyz++)
      {
         vec[xyz]=-2*uniq_exp[kvecs[i].alpha]
                    *CApows[dim_3*kvecs[i].c+dim_2*kvecs[i].a+dim_1*xyz+1];
         klen+=vec[xyz]*vec[xyz];
      };
      klen=sqrt(klen);
      kvec_vals[4*i+3]=klen;
      if (klen>0)
      {
         kvec_vals[4*i+0]=vec[0]/klen;
         kvec_vals[4*i+1]=vec[1]/klen;
         kvec_vals[4*i+2]=vec[2]/klen;
      } else {
         kvec_vals[4*i+0]=0.;
         kvec_vals[4*i+1]=0.;
         kvec_vals[4*i+2]=0.;
         kvec_vals[4*i+3]=0.;
      };
   } else {
      kvec_vals[4*i+0]=0.;
      kvec_vals[4*i+1]=0.;
      kvec_vals[4*i+2]=0.;
      kvec_vals[4*i+3]=0.;
   };
   return;
};

// compute k-vector pairs and screening values
__global__ void kvp_kernel
               (const double *CApows, const int dim_1, const int dim_2, const int dim_3,
                const kvpair_struct *kvecs, const int nkvpairs,
                const double *uniq_exp,
                double *kvp_vals, double *screen)
{
   int i=blockIdx.x*blockDim.x + threadIdx.x,
       xyz;
   double klen, arg, vec[3];
   const kvpair_struct *k=kvecs+i;

   // limit to number of data items
   if (i>=nkvpairs) return;

   if ((k->a!=k->c)||(k->b!=k->c))
   {
      klen=0;
      #pragma unroll 3
      for (xyz=0; xyz<3; xyz++)
      {
         vec[xyz]=-2*uniq_exp[k->alpha_a]*CApows[dim_3*k->c+dim_2*k->a+dim_1*xyz+1]
                  -2*uniq_exp[k->alpha_b]*CApows[dim_3*k->c+dim_2*k->b+dim_1*xyz+1];
         klen+=vec[xyz]*vec[xyz];
      };
      klen=sqrt(klen);
      kvp_vals[4*i+3]=klen;
      if (klen>0)
      {
         kvp_vals[4*i+0]=vec[0]/klen;
         kvp_vals[4*i+1]=vec[1]/klen;
         kvp_vals[4*i+2]=vec[2]/klen;
      } else {
         kvp_vals[4*i+0]=0;
         kvp_vals[4*i+1]=0;
         kvp_vals[4*i+2]=0;
         kvp_vals[4*i+3]=0;
      };

      // screening term
      arg=uniq_exp[k->alpha_a]*
           ( CApows[dim_3*k->c+dim_2*k->a+dim_1*0+2]
            +CApows[dim_3*k->c+dim_2*k->a+dim_1*1+2]
            +CApows[dim_3*k->c+dim_2*k->a+dim_1*2+2])
           +uniq_exp[k->alpha_b]*
           ( CApows[dim_3*k->c+dim_2*k->b+dim_1*0+2]
            +CApows[dim_3*k->c+dim_2*k->b+dim_1*1+2]
            +CApows[dim_3*k->c+dim_2*k->b+dim_1*2+2]);
      screen[i]=exp(-arg);

   } else {
      kvp_vals[4*i+0]=0.;
      kvp_vals[4*i+1]=0.;
      kvp_vals[4*i+2]=0.;
      kvp_vals[4*i+3]=0.;
      screen[i]=1.;
   };

   return;
};

