/* rad_kernel.cu: implements the CUDA kernel for computing k-vectors
 *
 * Written by Tom Grimes, 9-April 2010
 *
 */

#include "ecpints.h"

// define infinities
#define CUDA_INF __int_as_float(0x7f800000)
#define CUDA_NINF __int_as_float(0xff800000)

// forward decl for the confluent hypergeometric function
__device__ double cuda_chgf (const int N, const int l, const double z);

// compute local integrals
__global__ void local_radial_kernel
                             (const rad_int *rad_ints, const int *rad_list, 
                              const int lower, const int upper,
                              const double *kvp_vals, double *output)
{
   const int idx = lower + threadIdx.x + threadIdx.y*blockDim.x
                         +(blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
   // skip if out of limits
   if (idx>=upper) return;

   // get real index into total rad_ints list
   const rad_int *rad = rad_ints+rad_list[idx];

   double retval = cuda_chgf(rad->N, rad->lambda,
                             kvp_vals[4*rad->ktot+3]*kvp_vals[4*rad->ktot+3]/(4*rad->alpha))
                   *pow(kvp_vals[4*rad->ktot+3],rad->lambda);

   if ((retval==CUDA_INF)||(retval==CUDA_NINF)||isnan(retval)||(retval<0)) retval=0;
   output[rad_list[idx]]=retval;

   return;
};

// nonlocal integrals
__global__ void nonlocal_radial_kernel
                                (const rad_int *rad_ints, const int *rad_list,
                                 const int lower, const int upper,
                                 const double *kvec_vals, double *output) 
{
   const double sqrt_pi = 1.77245385091;
   const int max_nl_terms = 150;

   const int idx = lower + threadIdx.x + threadIdx.y*blockDim.x
                         +(blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
   // skip if out of limits
   if (idx>=upper) return;

   // get real index into total rad_ints list
   const rad_int *rad = rad_ints+rad_list[idx];

#define N rad->N
#define l rad->lambda
#define lp rad->lambdaP
#define K(x) kvec_vals[4*rad->k##x+3]
#define Z(x) K(x)*K(x)/(4*rad->alpha)

   // get CHGFs
   double Q0 = cuda_chgf(N+l,lp,Z(b))*pow(K(b),lp)*sqrt_pi*pow(rad->alpha,-0.5*(N+l+lp+1))/(1<<(lp+2))
               *tgamma(0.5*(N+l+lp+1))/tgamma(lp+1.5),
          Q2 = cuda_chgf(N+l+2,lp,Z(b))*pow(K(b),lp)*sqrt_pi*pow(rad->alpha,-0.5*(N+l+2+lp+1))/(1<<(lp+2))
               *tgamma(0.5*(N+l+2+lp+1))/tgamma(lp+1.5);

   double thisterm = pow(K(a),l);

   // compute double factorial
   for (int i=2*l+1; i>1; i-=2)
      thisterm/=i;

   double retval=thisterm*Q0;
   thisterm*=K(a)*K(a)/(2*(2*l+1+2));
   retval+=thisterm*Q2;

   double term = 1.;
   for (int j=2; ((j<max_nl_terms)&&(term>ECP_TOL)&&(retval>0)); j++)
   {
#define n  (N + l + 2*j)
      thisterm*=K(a)*K(a)/(2*j*(2*l+1+2*j));
      double Q4 = ((K(b)*K(b)+2*rad->alpha*(2*n-5))*Q2+(n+lp-3)*(lp-n+4)*Q0)/(4*rad->alpha*rad->alpha);
      term=thisterm*Q4;
      retval+=term;
      Q0=Q2;
      Q2=Q4;
   };

   if ((retval==CUDA_INF)||(retval==CUDA_NINF)||isnan(retval)||(retval<0)) retval=0;
   output[rad_list[idx]]=retval;

#undef n
#undef Z
#undef K
#undef lp
#undef l
#undef N

   return;
};

// compute the confluent hypergeometric function
__device__ double cuda_chgf (const int N, const int l, const double z)
{
   const int maxterms = 200;
   double a = 0.5*(N+l+1),
          b = l+1.5;

   // decide which series
   int series = 1;
   const int max_z[9] = {31,28,25,23,22,20,19,18,15};
   if (z>=(N>8? 15 : max_z[N])) series = 3;
   int trunc = -1;
   if (((N+l)%2==0)&&(N>=(l+2)))
   {
      trunc=(N-l-2)/2;
      series=2;
   }

   double retval = 1,
          thisterm = 1;
   int i;
   if (series==1)
   {
      // normal series
      for (i=0; (thisterm>ECP_TOL)&&(i<maxterms); i++)
      {
         thisterm*=(a+i)*z/((b+i)*(i+1));
         retval+=thisterm;
      };

   } else if (series==2) {
      // truncating series
      for (i=0; i<trunc; i++)
      {
         thisterm*=(b-a+i)*-z/((b+i)*(i+1));
         retval+=thisterm;
      };
      retval*=exp(z);

   } else if (series==3) {
      // asymptotic series for large z (see bounds above)
      // loop is in order so that only the minimum is reached, not beyond
      int sign = ((b-a)*(1-a)>0? 1 : -1);
      double lastterm = 2*sign;
      retval=0;
      for (i=0; ((sign*thisterm<sign*lastterm)&&(thisterm>ECP_TOL)); i++)
      {
         retval+=thisterm;
         lastterm=thisterm;
         thisterm*=(b-a+i)*(1-a+i)/(z*(i+1));
      };
      retval *= tgamma(b)/tgamma(a)*pow(z,a-b)*exp(z);
   };

   return retval;
};

