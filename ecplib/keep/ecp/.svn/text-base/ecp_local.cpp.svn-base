// ecp_local.cpp
//
// Written by Tom Grimes, 17-Dec 2009
//
// This file implements the local part of ECP integrals
//

#include "ecpints.h"

double ecp_local_angular
(int lambda, int I, int J, int K, double kvec[3])
{
   // connect to coefficient tables
   int max_l=ecp_data.Ylm_maxl,
       max_IJK=ecp_data.max_IJK;
   double (*Ylm)[max_l+1][max_l+1][max_l+1] = (double(*)[max_l+1][max_l+1][max_l+1]) ecp_data.Ylm;
   double (*l_angtab)[max_IJK+1][max_IJK+1][max_IJK+1] = (double(*)[max_IJK+1][max_IJK+1][max_IJK+1]) ecp_data.local_ang;

   double angterm = 0.;

   // optimized method of computing norm factors
   double kvpows[3][lambda+1];
   kvpows[0][0]=1.;
   kvpows[1][0]=1.;
   kvpows[2][0]=1.;
   for (int i=1; i<=lambda; i++)
      for (int xyz=0; xyz<3; xyz++)
         kvpows[xyz][i]=kvpows[xyz][i-1]*kvec[xyz];

   for (int mu=-lambda; mu<=lambda; mu++)
   {
      int idx=lambda*(lambda+1)+lambda+mu;
      int abs_mu=abs(mu);
      double eval = 0.;

      // loop over xyz exponents for Y(;lambda,mu)
      for (int nx=0; nx<=lambda; nx++)
      for (int ny=0; (nx+ny)<=lambda; ny++)
      for (int nz=lambda-nx-ny; nz>=0; nz-=2)
      {
         // evaluate the monomial for Y(Omega_k; lambda,mu)
         eval+=Ylm[idx][nx][ny][nz]*kvpows[0][nx]*kvpows[1][ny]*kvpows[2][nz];
      };
      angterm+=eval*l_angtab[idx][I][J][K];
   };
   
   return angterm;
};
   
