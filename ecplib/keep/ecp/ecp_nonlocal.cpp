// ecp_nonlocal.cpp
//
// Written by Tom Grimes, 17-Dec 2009
//
// This file implements nonlocal ECP integrals
//

#include "ecpints.h"

/* uses precalculated coefficients */
double ecp_nonlocal_angular
(int lambda, int l, int m, int I, int J, int K, double kvec[3])
{
   static const double pi = acos(-1.);

   // connect to Y(;lambda,mu) coefficient tables
   int max_l=ecp_data.Ylm_maxl,
       max_nl_IJK=ecp_data.max_nl_IJK;
   double (*Ylm)[max_l+1][max_l+1][max_l+1] = (double(*)[max_l+1][max_l+1][max_l+1]) ecp_data.Ylm;
   double (*nl_angtab)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1] =
       (double(*)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1]) ecp_data.nonlocal_ang;

   int l_idx=l*(l+1)+l+m;
   double angterm = 0.;

   // optimized method of computing norm factors
   double kvpows[3][lambda+1];
   kvpows[0][0]=1.;
   kvpows[1][0]=1.;
   kvpows[2][0]=1.;
   for (int i=1; i<=lambda; i++)
      for (int xyz=0; xyz<3; xyz++)
         kvpows[xyz][i]=kvpows[xyz][i-1]*kvec[xyz];

   // double loop: two angular factors to expand:
   // sum[mu=-lambda..lambda](Y(Omega_k;lambda,mu)*int((1/4pi)*x^a*y^b*z*c*Y(;lambda,mu)*Y(;l,m))
   for (int mu=-lambda; mu<=lambda; mu++)
   {
      int lambda_idx=lambda*(lambda+1)+lambda+mu;
      int abs_mu=abs(mu);
      double eval = 0.;

      // loop over xyz exponents for Y(;lambda,mu)
      for (int nx=0; nx<=lambda; nx++)
      for (int ny=0; (nx+ny)<=lambda; ny++)
      for (int nz=lambda-nx-ny; nz>=0; nz-=2)
      {
         // evaluate the monomial for Y(Omega_k; lambda,mu)
         eval+=Ylm[lambda_idx][nx][ny][nz]*kvpows[0][nx]*kvpows[1][ny]*kvpows[2][nz];
      };
      angterm+=eval*nl_angtab[lambda_idx][l_idx][I][J][K];
   };
   return angterm;
};

