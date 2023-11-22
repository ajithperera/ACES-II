/* ecp_drivers.cpp
 *
 * Written by Tom Grimes, 18-March 2010
 *
 * This file contains the TBB classes for computing the integrals in parallel
 *
 * Beware of the static global ecp specification!
 *
 * Parallelism lives here.
 *
 */

#include "ecpints.h"

// compute the angular integrals
void ang_ints_driver::operator() ( const tbb::blocked_range<int> &range ) const
{
   for (int i=range.begin(); i!=range.end(); ++i)
   {
      ang_struct *omega = list+i;
      if (omega->ang2_lambda==-1)
      {
         // local
         output[i]=ecp_local_angular(omega->ang1_lambda,omega->ang1_I,omega->ang1_J,omega->ang1_K,
                                      kvp_vals[omega->ktot]);
      } else {
         // nonlocal
         output[i]=0.;
         for (int m=-omega->l; m<=omega->l; m++)
         {
            double term1=ecp_nonlocal_angular(omega->ang1_lambda,omega->l,m,
                                              omega->ang1_I,omega->ang1_J,omega->ang1_K,kvec_vals[omega->ka]),
                   term2=ecp_nonlocal_angular(omega->ang2_lambda,omega->l,m,
                                              omega->ang2_I,omega->ang2_J,omega->ang2_K,kvec_vals[omega->kb]);
            output[i]+=term1*term2;
         };
      };
   };
   return;
};

// compute the radial integrals
void rad_ints_driver::operator() ( const tbb::blocked_range<int> &range ) const
{
   // maximum number of terms in nonlocal series
   const int max_nl = 150;

   for (int idx=range.begin(); idx!=range.end(); ++idx)
   {
      rad_int *rad = list+idx;

      // check for local or nonlocal
      if (rad->lambdaP==-1)
      {
         // exact local solution
         double k = kvp_vals[rad->ktot][3],
                z = k*k/(4*rad->alpha);
         output[idx]=intpwr_mod(k,rad->lambda)*chgf(rad->N,rad->lambda,z);

      } else {
         // exact nonlocal solution
         static const double sqrt_pi = sqrt(acos(-1.));
         int N = rad->N,
             l = rad->lambda,
             lp = rad->lambdaP;
         double ka = kvec_vals[rad->ka][3],
                ka2 = ka*ka/2,
                kb = kvec_vals[rad->kb][3],
                kb2 = kb*kb,
                zb = kb2/(4*rad->alpha),
                a2 = 4*rad->alpha*rad->alpha,
                Q0 = sqrt_pi*pow(kb,lp)*pow(rad->alpha,-0.5*(N+l+lp+1))/(1<<(lp+2))
                     *tgamma(0.5*(N+l+lp+1))/tgamma(lp+1.5)
                     *chgf(N+l,lp,zb),
                Q2 = sqrt_pi*pow(kb,lp)*pow(rad->alpha,-0.5*(N+l+2+lp+1))/(1<<(lp+2))
                     *tgamma(0.5*(N+l+2+lp+1))/tgamma(lp+1.5)
                     *chgf(N+l+2,lp,zb),
                thisterm = pow(ka,l)/fact2(2*l+1);

         output[idx]=thisterm*Q0;
         thisterm*=ka2/(2*l+1+2);
         output[idx]+=thisterm*Q2;

         double term = 1.;
         int j;
         for (j=2; ((j<max_nl)&&(term>ECP_TOL)); j++)
         {
            int n = N + l + 2*j;
            thisterm*=ka2/(j*(2*l+1+2*j));
            double Q4 = ((kb2+2*rad->alpha*(2*n-5))*Q2+(n+lp-3)*(lp-n+4)*Q0)/a2;
            term=thisterm*Q4;
            output[idx]+=term;
            Q0=Q2;
            Q2=Q4;
         };
      };

      if (isnan(output[idx])||isinf(output[idx])) output[idx]=0.;

   };
   return;
};

// combine the results
void combine_expr_driver::operator() ( const tbb::blocked_range<int> &range ) const
{
   // recast distance structure
   double (*CApows)[nat][3][maxdim] = (double(*)[nat][3][maxdim]) CApows_;

   // get indices for our range
   int begin = blocks[range.begin()],
       end = blocks[range.end()];

   for (int i=begin; i<end; i++)
   {
      int a = prims[i].atid_a,
          b = prims[i].atid_b, 
          c = prims[i].atid_c;

      double val=prims[i].prefact*screen[prims[i].scr];
      val*=CApows[c][a][0][prims[i].a]*CApows[c][a][1][prims[i].b]*CApows[c][a][2][prims[i].c];
      val*=CApows[c][b][0][prims[i].d]*CApows[c][b][1][prims[i].e]*CApows[c][b][2][prims[i].f];

      // unified for local and nonlocal
      // note the radial part may be constant and included in prefact
      if (prims[i].omega!=-1) val*=angvals[prims[i].omega];
      if (prims[i].Q!=-1) val*=radvals[prims[i].Q];

      output[prims[i].idx]+=val;
   };
   return;
};

// combine the results of a gradient calculation
// calculate <di/dX(A)|V(C)|j> and translate to CSD structure
void combine_grad_driver::operator() ( const tbb::blocked_range<int> &range ) const
{
   // recast data
   double (*CApows)[nat][3][maxdim] = (double(*)[nat][3][maxdim]) CApows_;
   double (*gradout)[nat][norbp] = (double(*)[nat][norbp]) output;

   // get indices for our range
   int begin = blocks[range.begin()],
       end = blocks[range.end()];

   for (int i=begin; i<end; i++)
   {
      int a = prims[i].atid_a,
          b = prims[i].atid_b, 
          c = prims[i].atid_c;

      double val=prims[i].prefact*screen[prims[i].scr];
      val*=CApows[c][a][0][prims[i].a]*CApows[c][a][1][prims[i].b]*CApows[c][a][2][prims[i].c];
      val*=CApows[c][b][0][prims[i].d]*CApows[c][b][1][prims[i].e]*CApows[c][b][2][prims[i].f];

      // unified for local and nonlocal
      if (prims[i].omega!=-1) val*=angvals[prims[i].omega];
      if (prims[i].Q!=-1) val*=radvals[prims[i].Q];

      // calculate indices
      int n = prims[i].idx/dim1,
          orb_i = (prims[i].idx-n*dim1)/dim2,
          orb_j = (prims[i].idx-n*dim1-orb_i*dim2)/dim3,
          xyz = prims[i].idx%dim3,
          pack_ij=orb_i*(orb_i+1)/2+orb_j,
          pack_ji=orb_j*(orb_j+1)/2+orb_i;

      if (orb_i>=orb_j)
      {
         gradout[xyz][n][pack_ij]+=val;
         gradout[xyz+3][n][pack_ij]-=val;
      };
      if (orb_j>=orb_i)
      {
         gradout[xyz+3][n][pack_ji]-=val;
      };
   };
   return;
};

