/* ecp_expr_opt.cpp: optimize level 1 expressions
 *
 * Written by Tom Grimes, 25-March 2010
 *
 */

#include "ecpints.h"
#include "ecp_opt.h"

// simplify comparison
#define LEXCMP(member) if (a.##member!=b.##member) return (a.##member<b.##member);
bool operator<(const k_struct &a, const k_struct &b)
{
   LEXCMP(alpha)
   LEXCMP(a)
   LEXCMP(c)

   return false;
};
bool operator<(const kvpair_struct &a, const kvpair_struct &b)
{
   LEXCMP(a)
   LEXCMP(b)
   LEXCMP(c)
   LEXCMP(alpha_a)
   LEXCMP(alpha_b)

   return false;
};
bool operator<(const rad_int &a, const rad_int &b)
{
   LEXCMP(N)
   LEXCMP(lambda)
   LEXCMP(lambdaP)

   if (a.lambdaP==-1)
   {
      LEXCMP(ktot)
   } else {
      LEXCMP(ka)
      LEXCMP(kb)
   };

   if (fabs(a.alpha-b.alpha)>ECP_TOL) return (a.alpha<b.alpha);

   return false;
};
bool operator<(const ang_struct &a, const ang_struct &b)
{
   LEXCMP(l)
   LEXCMP(ka)
   LEXCMP(kb)
   LEXCMP(ktot)
   LEXCMP(ang1_I)
   LEXCMP(ang1_J)
   LEXCMP(ang1_K)
   LEXCMP(ang1_lambda)
   LEXCMP(ang2_I)
   LEXCMP(ang2_J)
   LEXCMP(ang2_K)
   LEXCMP(ang2_lambda)

   return false;
};

// simplify sorting
bool primints_predicate(const dryrun_lvl2 &a, const dryrun_lvl2 &b)
{
   return (a.idx<b.idx);
};

// actual optimization routine prototype
void opt_list(vector<dryrun_lvl1> &primlist, vector<dryrun_lvl2> &outprims,
              int &cvt_rad, int &const_ang, int &const_rad, int &cache_k, int &cache_kvp,
              int &cache_ang, int &cache_rad);

void ecp_expr_opt(vector<dryrun_lvl1> &primints1, vector<dryrun_lvl1> &primgrad1, FILE *outf)
{
   // declare vars
   int cvt_rad, const_ang, const_rad, cache_k, cache_kvp, cache_ang, cache_rad;

   // do the normal expressions
   opt_list(primints1, ecp_data.primints2,
            cvt_rad, const_ang, const_rad, cache_k, cache_kvp, cache_ang, cache_rad);

   fprintf(outf,"  The normal primitive expressions generate:\n");
   fprintf(outf,"      %i k-vectors\n",ecp_data.kvecs.size());
   fprintf(outf,"      %i k-vector pairs\n",ecp_data.kvpairs.size());
   fprintf(outf,"      %i angular integrals\n",ecp_data.angints.size());
   fprintf(outf,"      %i radial integrals\n",ecp_data.radints.size());
   fprintf(outf,"  (%i angular integrals determine to be const)\n",const_ang);
   fprintf(outf,"  (%i nonlocal radial integrals converted to local)\n",cvt_rad);
   fprintf(outf,"  (%i local radial integrals determined to be const)\n",const_rad);
   fprintf(outf,"\n");

   // save the normal integral extents
   ecp_data.norm_klen = ecp_data.kvecs.size();
   ecp_data.norm_kvplen = ecp_data.kvpairs.size();
   ecp_data.norm_anglen = ecp_data.angints.size();
   ecp_data.norm_radlen = ecp_data.radints.size();

   // now do the gradient integrals
   opt_list(primgrad1, ecp_data.primgrad2,
            cvt_rad, const_ang, const_rad, cache_k, cache_kvp, cache_ang, cache_rad);

   fprintf(outf,"  The gradient primitive expressions generate an additional:\n");
   fprintf(outf,"      %i k-vectors\n",ecp_data.kvecs.size()-ecp_data.norm_klen);
   fprintf(outf,"      %i k-vector pairs\n",ecp_data.kvpairs.size()-ecp_data.norm_kvplen);
   fprintf(outf,"      %i angular integrals\n",ecp_data.angints.size()-ecp_data.norm_anglen);
   fprintf(outf,"      %i radial integrals\n",ecp_data.radints.size()-ecp_data.norm_radlen);
   fprintf(outf,"  (%i angular integrals determine to be const)\n",const_ang);
   fprintf(outf,"  (%i nonlocal radial integrals converted to local)\n",cvt_rad);
   fprintf(outf,"  (%i local radial integrals determined to be const)\n",const_rad);
   fprintf(outf,"\n");

   fprintf(outf,"  The gradient primitive expressions cache:\n");
   fprintf(outf,"      %i k-vectors\n",cache_k);
   fprintf(outf,"      %i k-vector pairs\n",cache_kvp);
   fprintf(outf,"      %i angular integrals\n",cache_ang);
   fprintf(outf,"      %i radial integrals\n",cache_rad);
   fprintf(outf,"  (from the normal integrals)\n\n");

   // sort the vector for blocking and create the block list
   fprintf(outf,"  Creating blocking lists... ");
   sort(ecp_data.primints2.begin(),ecp_data.primints2.end(),primints_predicate);
   int lastidx=-1;
   for (int i=0; i<ecp_data.primints2.size(); i++)
   {
      if (ecp_data.primints2[i].idx!=lastidx)
      {
         lastidx=ecp_data.primints2[i].idx;
         ecp_data.blocks.push_back(i);
      };
   };
   ecp_data.blocks.push_back(ecp_data.primints2.size());

   sort(ecp_data.primgrad2.begin(),ecp_data.primgrad2.end(),primints_predicate);
   lastidx=-1;
   for (int i=0; i<ecp_data.primgrad2.size(); i++)
   {
      if (ecp_data.primgrad2[i].idx!=lastidx)
      {
         lastidx=ecp_data.primgrad2[i].idx;
         ecp_data.gradblocks.push_back(i);
      };
   };
   ecp_data.gradblocks.push_back(ecp_data.primgrad2.size());

   fprintf(outf,"done.\n\n");
   fflush(outf);

   return;
};

/* optimizing code */
// note that primlist is consumed
void opt_list(vector<dryrun_lvl1> &primlist, vector<dryrun_lvl2> &outprims,
              int &cvt_rad, int &const_ang, int &const_rad, int &cache_k, int &cache_kvp,
              int &cache_ang, int &cache_rad)
{
   // get uniquifying tool
   uniquify<k_struct> uniq_kvec(ecp_data.kvecs);
   uniquify<kvpair_struct> uniq_kvp(ecp_data.kvpairs);
   uniquify<ang_struct> uniq_ang(ecp_data.angints);
   uniquify<rad_int> uniq_rad(ecp_data.radints);

   // clear variables
   int dummy;
   cvt_rad=0;
   const_ang=0;
   const_rad=0;
   cache_k=0;
   cache_kvp=0;
   cache_ang=0;
   cache_rad=0;

   // reserve space for expressions
   outprims.reserve(primlist.size());

   const double pi=acos(-1.);

   while (!primlist.empty())
   {
      // get the top element
      dryrun_lvl1 prim = primlist.back();
      primlist.pop_back();

      // check the k-vectors
      k_struct k_a = {prim.alpha_a,prim.atid_a,prim.atid_ecp},
               k_b = {prim.alpha_b,prim.atid_b,prim.atid_ecp};
      int idx_a = uniq_kvec.index(k_a,cache_k),
          idx_b = uniq_kvec.index(k_b,cache_k);

      // and pairwise k-vectors ("total" k-vectors for local angular terms)
      kvpair_struct kvp_1 = {prim.alpha_a,prim.alpha_b,prim.atid_a,prim.atid_b,prim.atid_ecp},
                    kvp_2 = {prim.alpha_b,prim.alpha_a,prim.atid_b,prim.atid_a,prim.atid_ecp},
                    kvp = (kvp_1<kvp_2 ? kvp_1 : kvp_2);
      int ktot = uniq_kvp.index(kvp,cache_kvp);

      // check the angular integrals
      ang_struct ang_1 = { prim.l, idx_a, idx_b, ktot,
                           prim.ang1_I, prim.ang1_J, prim.ang1_K, prim.ang1_lambda,
                           prim.ang2_I, prim.ang2_J, prim.ang2_K, prim.ang2_lambda },
                 ang_2 = { prim.l, idx_b, idx_a, ktot,
                           prim.ang2_I, prim.ang2_J, prim.ang2_K, prim.ang2_lambda,
                           prim.ang1_I, prim.ang1_J, prim.ang1_K, prim.ang1_lambda },
                 ang = ((prim.l==-1)||(ang_1<ang_2) ? ang_1 : ang_2);
      int ang_idx = -1;

      // strength reduction on lambda==0 (isotropic)
      if ((ang.ang2_lambda==-1)&&(ang.ang1_lambda==0))
      {
         // compute constant value for local angular
         // monomial should never evaluate to zero, but check just in case
         if ((ang.ang1_I%2==1)||(ang.ang1_J%2==1)||(ang.ang1_K%2==1)) continue;

         prim.prefact*=fact2(ang.ang1_I-1)*fact2(ang.ang1_J-1)*fact2(ang.ang1_K-1)
                      /fact2(ang.ang1_I+ang.ang1_J+ang.ang1_K+1);

         const_ang++;

/*
      } else if ((prim.l==0)&&(ang.ang1_lambda==0)&&(ang.ang2_lambda==0)) {
         // constant value for nonlocal angular
         // monomials should never evaluate to zero, but check just in case
         if ((ang.ang1_I%2==1)||(ang.ang1_J%2==1)||(ang.ang1_K%2==1)) continue;
         if ((ang.ang2_I%2==1)||(ang.ang2_J%2==1)||(ang.ang2_K%2==1)) continue;

         prim.prefact*=fact2(ang.ang1_I-1)*fact2(ang.ang1_J-1)*fact2(ang.ang1_K-1)
                      *fact2(ang.ang2_I-1)*fact2(ang.ang2_J-1)*fact2(ang.ang2_K-1)
                      /(4*pi*fact2(ang.ang1_I+ang.ang1_J+ang.ang1_K+1)
                            *fact2(ang.ang2_I+ang.ang2_J+ang.ang2_K+1));

         const_ang++;
*/

      } else if ((ang.ang1_lambda==0)&&(ang.ang2_lambda==0)) {
         // dirty! steal integrals from global ecp_data
         int max_l=ecp_data.Ylm_maxl,
              max_nl_IJK=ecp_data.max_nl_IJK;
         double (*nl_angtab)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1] =
              (double(*)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1]) ecp_data.nonlocal_ang;

         double angterm = 0;
         for (int m=-prim.l; m<=prim.l; m++)
         {
            angterm+=nl_angtab[0][prim.l*(prim.l+1)+prim.l+m][ang.ang1_I][ang.ang1_J][ang.ang1_K]
                    *nl_angtab[0][prim.l*(prim.l+1)+prim.l+m][ang.ang2_I][ang.ang2_J][ang.ang2_K]
                    /(4*pi);
         };
         prim.prefact*=angterm;
         const_ang++;

      } else {
         // default case: compute value as needed
         ang_idx=uniq_ang.index(ang,cache_ang);

      };

      // construct the radial factors
      double rad_alpha=ecp_data.uniq_exp[prim.alpha_a]+ecp_data.uniq_exp[prim.alpha_b]+prim.zeta;
      rad_int Q_1 = {prim.Q_lambda, prim.Q_lambdaP, prim.Q_N, idx_a, idx_b, ktot, rad_alpha},
              Q_2 = {prim.Q_lambdaP, prim.Q_lambda, prim.Q_N, idx_b, idx_a, ktot, rad_alpha},
              Q = ((prim.Q_lambdaP==-1)||(Q_1<Q_2) ? Q_1 : Q_2);

      // convert nonlocal radial to local, if possible
      // due to possible order switching, need to go back to ecp_data struct
      if (Q.lambdaP!=-1)
      {
         if ((ecp_data.kvecs[Q.ka].a==prim.atid_ecp)&&(Q.lambda==0))
         {
            // first Bessel is always unity: translate
            Q.lambda=Q.lambdaP;
            Q.lambdaP=-1;
            int temp = Q.ka;
            Q.ka=Q.kb;
            Q.kb=temp;
            cvt_rad++;
         } else if ((ecp_data.kvecs[Q.kb].a==prim.atid_ecp)&&(Q.lambdaP==0)) {
            // second Bessel is always unity: translate
            Q.lambdaP=-1;
            cvt_rad++;
         };
      };

      // check for constant radial int and get index if not
      // only applies to local and converted nonlocal (=local)
      int Q_idx = -1;
      if ((Q.lambdaP==-1)&&(prim.atid_a==prim.atid_ecp)&&(prim.atid_b==prim.atid_ecp))
      {
         if (prim.Q_N%2==0)
         {
            prim.prefact*=sqrt(pi/2)*fact2(Q.N-1)/pow(2*rad_alpha,0.5*(Q.N+1));
         } else {
            prim.prefact*=fact2(Q.N-1)/pow(2*rad_alpha,0.5*(Q.N+1));
         };
         const_rad++;

      } else {
         // general case: get index
         Q_idx = uniq_rad.index(Q,cache_rad);
         
         // if local, set normalization
         if (Q.lambdaP==-1)
            prim.prefact*=sqrt(pi)*pow(Q.alpha,-0.5*(Q.N+Q.lambda+1))/(1<<(Q.lambda+2))
                          *tgamma(0.5*(Q.N+Q.lambda+1))/tgamma(Q.lambda+1.5);

      };

      // generate the element
      dryrun_lvl2 newprim = {prim.idx,
                             ecp_data.uniq_exp[prim.alpha_a], ecp_data.uniq_exp[prim.alpha_b],
                             prim.atid_a, prim.atid_b, prim.atid_ecp,
                             prim.prefact, ktot,
                             prim.a, prim.b, prim.c, prim.d, prim.e, prim.f,
                             ang_idx, Q_idx};
      outprims.push_back(newprim);
   };

   return;
};


