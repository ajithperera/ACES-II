/* ang_kernel.cu: CUDA computation of angular integrals
 *
 * Written by Tom Grimes, 27-April 2010
 *
 */

#include "ecpints.h"

#include "cuda.h"
#include "cutil_inline.h"

// kernel prototypes
__global__ void local_angular_kernel
                             (const ang_struct *ang_ints, const int *ang_list,
                              const int lower, const int upper,
                              const double *kvp_vals,
                              const double *Ylm_, const int max_l,
                              const double *local_tab_, const int max_IJK,
                              double *output)
{
   // get index and pointer to integral
   int list_idx = lower+threadIdx.x+threadIdx.y*blockDim.x +
                  (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
   if (list_idx>=upper) return;

   // load ang struct and k vector
   const ang_struct ang = ang_ints[ang_list[list_idx]];
   short int lambda = ang.ang1_lambda,
             I = ang.ang1_I,
             J = ang.ang1_J,
             K = ang.ang1_K;
   int kptr = 4*ang.ktot;
   float k0 = kvp_vals[kptr],
         k1 = kvp_vals[kptr+1],
         k2 = kvp_vals[kptr+2];

   // get pointer into Ylm
   double (*Ylm)[max_l+1][max_l+1][max_l+1] = (double(*)[max_l+1][max_l+1][max_l+1]) Ylm_;
   double (*l_angtab)[max_IJK+1][max_IJK+1][max_IJK+1] = 
                                  (double(*)[max_IJK+1][max_IJK+1][max_IJK+1]) local_tab_;

#define kpow(x,n) pow(k##x,n)

   double angterm = 0;
   for (short int mu=-lambda; mu<=lambda; mu++)
   {
      double eval = 0.;

      // loop over xyz exponents for Y(;lambda,mu)
      for (int nx=0; nx<=lambda; nx++)
      for (int ny=0; (nx+ny)<=lambda; ny++)
      for (int nz=lambda-nx-ny; nz>=0; nz-=2)
      {
         // evaluate the monomial for Y(Omega_k; lambda,mu)
         eval+=Ylm[lambda*(lambda+1)+lambda+mu][nx][ny][nz]*kpow(0,nx)*kpow(1,ny)*kpow(2,nz);
      };
      angterm+=eval*l_angtab[lambda*(lambda+1)+lambda+mu][I][J][K];
   };
   output[ang_list[list_idx]]=angterm;

#undef kpow

   return;
};

// forward decl
__device__ double nonlocal_angular_term (const short int &lambda, const short int &l, const short int &m,
                                         const short int &I, const short int &J, const short int &K,
                                         const double *kvec,
                                         const double *Ylm, const int &max_l,
                                         const double *nonlocal_tab, const int &max_nl_IJK);

__global__ void nonlocal_angular_kernel
                             (const ang_struct *ang_ints, const int *ang_list,
                              const int lower, const int upper,
                              const double *kvec_vals,
                              const double *Ylm, const int max_l,
                              const double *nonlocal_tab, const int max_nl_IJK,
                              double *output)
{
   // get index and pointer to integral
   int idx = lower+threadIdx.x+(blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x;
   if (idx>=upper) return;

   const ang_struct *ang = ang_ints+ang_list[idx];

   // do individual terms in parallel
   __shared__ double work[32];
   short int lambda, I, J, K,
             l = ang->l,
             tid = 2*threadIdx.x+threadIdx.y;
   typedef struct {
      double val[3];
   } double3;
   double3 kvec;
   if (threadIdx.y==0)
   {
      lambda=ang->ang1_lambda;
      I=ang->ang1_I;
      J=ang->ang1_J;
      K=ang->ang1_K;
      kvec=*((double3*)(kvec_vals+4*ang->ka));
   } else {
      lambda=ang->ang2_lambda;
      I=ang->ang2_I;
      J=ang->ang2_J;
      K=ang->ang2_K;
      kvec=*((double3*)(kvec_vals+4*ang->kb));
   };

   double intval=0.;

   for (short int m=-ang->l; m<=ang->l; m++)
   {
      work[tid]=nonlocal_angular_term(lambda, l, m, I, J, K, kvec.val,
                                      Ylm, max_l, nonlocal_tab, max_nl_IJK);
      __threadfence_block();
      if (threadIdx.y==0) intval+=work[tid]*work[tid+1];
   };

   if (threadIdx.y==0) output[ang_list[idx]]=intval;

   return;
};

__device__ double nonlocal_angular_term (const short int &lambda, const short int &l, const short int &m,
                                         const short int &I, const short int &J, const short int &K,
                                         const double *kvec,
                                         const double *Ylm_, const int &max_l,
                                         const double *nonlocal_tab, const int &max_nl_IJK)
{
   // get pointer into Ylm and nonlocal_tab
   double (*Ylm)[max_l+1][max_l+1][max_l+1] = (double(*)[max_l+1][max_l+1][max_l+1]) Ylm_;
   double (*nl_angtab)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1] =
       (double(*)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1]) nonlocal_tab;

   int l_idx=l*(l+1)+l+m;
   double angterm = 0.;

   // double loop: two angular c_factors to expand:
   // sum[mu=-lambda..lambda](Y(Omega_k;lambda,mu)*int((1/4pi)*x^a*y^b*z*c*Y(;lambda,mu)*Y(;l,m))
   for (int mu=-lambda; mu<=lambda; mu++)
   {
      int lambda_idx=lambda*(lambda+1)+lambda+mu;
      double eval = 0.;

      // loop over xyz exponents for Y(;lambda,mu)
      for (int nx=0; nx<=lambda; nx++)
      for (int ny=0; (nx+ny)<=lambda; ny++)
      for (int nz=lambda-nx-ny; nz>=0; nz-=2)
      {
         // evaluate the monomial for Y(Omega_k; lambda,mu)
         eval+=Ylm[lambda_idx][nx][ny][nz]*pow(kvec[0],nx)*pow(kvec[1],ny)*pow(kvec[2],nz);
      };
      angterm+=eval*nl_angtab[lambda_idx][l_idx][I][J][K];
   };

   return angterm;
};

/**************************************************************************/
// Try calculating Ylm constants
/*
// helper functions
__device__ short int c_fact(short int x)
{
   if (x<2) return 1;
   switch(x)
   {
   case 2: return 2;
   case 3: return 6;
   case 4: return 24;
   case 5: return 120;
   case 6: return 720;
   case 7: return 5040;
   };
};

__device__ short int c_fact2(short int x)
{
   if (x<3) return 1;
   switch(x)
   {
   case 3: return 3;
   case 5: return 15;
   case 7: return 105;
   case 9: return 945;
   case 11: return 10395;
   };
};

__device__ short int c_nCk(short int n, short int k)
{
   return c_fact(n)/(c_fact(k)*c_fact(n-k));
};

// this takes k and p, internal counters: cf. ecp_spec.cpp
// note: if m>=0 use cos_fact[(abs_m-p)%4] (otherwise add 3 to get sine series)
// constants are 4pi, sqrt(2)
__device__ double Ylm_coef(short int l, short int &abs_m, short int &k, short int &p)
{
   double retval = sqrt(c_fact(l-abs_m)*(2*l+1)/(12.5663706144*c_fact(l+abs_m)))
                   *c_nCk(l,k)*c_nCk(2*(l-k),l)*c_fact(l-2*k)/((1<<l)*c_fact(l-2*k-abs_m))
                   *c_nCk(abs_m,p);
   if (abs_m!=0) retval*=1.41421356237;
   if (k%2==1) retval=-retval;
   return retval;
};

// integral of xyz hat monomial over solid angle
// constant is 4pi
__device__ double monomial(short int I, short int J, short int K)
{
   if ((I%2==1)||(J%2==1)||(K%2==1)) return 0;
   return 12.5663706144*c_fact2(I-1)*c_fact2(J-1)*c_fact2(K-1)/c_fact2(I+J+K+1);
};

__device__ double nonlocal_angular_term_alternate
                                        (const short int &lambda, const short int &l, const short int &m,
                                         const short int &I, const short int &J, const short int &K,
                                         const double *kvec)
{
   short int abs_m = (m>=0? m : -m);

   // main loop
   double angterm = 0;
   for (short int mu=0; mu<=lambda; mu++)
   {
      double intval[2] = {0,0},
             Ylmval[2] = {0,0};

      for (short int k=0; k<=(lambda-mu)/2; k++)
      for (short int p=0; p<=mu; p++)
      {
         // total Ylm coefficient, including sin/cos factor
         double coef = Ylm_coef(lambda, mu, k, p);
         if ((mu-p)%4>1) coef=-coef;

         Ylmval[(mu-p)%2] += coef*pow(kvec[0],p)*pow(kvec[1],mu-p)*pow(kvec[2],lambda-2*k-mu);

         // compute Ylm (not lambda mu) terms
         for (short int kk=0; kk<=(l-abs_m)/2; k++)
         for (short int pp=(m>=0? abs_m%2 : 1-abs_m%2); pp<=abs_m; pp+=2)
         {
            double coef2 = Ylm_coef(l,abs_m,kk,pp)*((abs_m-pp)%4>1? -1 : 1);
            intval[(mu-p)%2] += coef*coef2*monomial(I+p+pp,J+mu-p+abs_m-pp,K+lambda-2*k-mu+l-2*kk-abs_m);
         };
      };

      angterm += Ylmval[0]*intval[0];
      if (mu>0) angterm += Ylmval[1]*intval[1];

   };

   return angterm;
};

__global__ void NEW_local_angular_kernel
                             (const ang_struct *ang_ints, const int *ang_list,
                              const int lower, const int upper,
                              const double *kvp_vals,
                              double *output)
{
   // get index and pointer to integral
   int list_idx = lower+threadIdx.x+threadIdx.y*blockDim.x +
                  (blockIdx.x+blockIdx.y*gridDim.x)*blockDim.x*blockDim.y;
   if (list_idx>=upper) return;

   // load ang struct and k vector
   ang_struct ang = ang_ints[ang_list[list_idx]];
#define lambda ang.ang1_lambda
#define I ang.ang1_I
#define J ang.ang1_J
#define K ang.ang1_K

   // attempt to coalesce memory access
   typedef struct {
      double val[3];
   } double3;
   double3 kvec = *((double3*)(kvp_vals+4*ang.ktot));

   // main loops: positive and negative mu have been compactified
   double angterm = 0;
   for (short int mu=0; mu<=lambda; mu++)
   {
      double intval[2] = {0,0},
             Ylmval[2] = {0,0};

      for (short int k=0; k<=(lambda-mu)/2; k++)
      for (short int p=0; p<=mu; p++)
      {
         // total Ylm coefficient, including sin/cos factor
         double coef = Ylm_coef(lambda, mu, k, p);
         if ((mu-p)%4>1) coef=-coef;

         Ylmval[(mu-p)%2] += coef*pow(kvec.val[0],p)*pow(kvec.val[1],mu-p)*pow(kvec.val[2],lambda-2*k-mu);
         intval[(mu-p)%2] += coef*monomial(I+p,J+mu-p,K+lambda-2*k-mu);

      };

      angterm += Ylmval[0]*intval[0];
      if (mu>0) angterm += Ylmval[1]*intval[1];

   };

   output[ang_list[list_idx]]=angterm;

   return;
};
*/


