/* ecp_prepare.cpp: generates and optimizes the primitive integral expressions
 *
 * Written by Tom Grimes, 25-March 2010
 *
 */

#include "ecpints.h"
#include "ecp_opt.h"

// generate expressions in parallel
class gener_prim_expr {
   // private data
   int *exp_ind, *primofs, *aoofs, *shelltype, *atid;
   double *cc;
   int nunique, norb, nat, norbp;

   public:
      // output
      vector<dryrun_lvl1> primints1, primgrad1;
      int zero_pre, zero_ang, zero_rad;

      // new constructor
      gener_prim_expr(int *exp_ind_,double *cc_,int *primofs_,int *aoofs_,int *shelltype_,
                      int *atid_, int nunique_,int norb_):
                      exp_ind(exp_ind_),cc(cc_),primofs(primofs_),aoofs(aoofs_),
                      shelltype(shelltype_),atid(atid_),
                      nunique(nunique_),norb(norb_),norbp(norb_*(norb_+1)/2),nat(atid_[nunique_-1]),
                      zero_pre(0),zero_ang(0),zero_rad(0) {};

      // split constructor
      gener_prim_expr( gener_prim_expr &x, tbb::split ):
                      exp_ind(x.exp_ind),cc(x.cc),primofs(x.primofs),aoofs(x.aoofs),
                      shelltype(x.shelltype),atid(x.atid),
                      nunique(x.nunique),norb(x.norb),norbp(x.norbp),nat(x.nat),
                      zero_pre(0),zero_ang(0),zero_rad(0) {};

      // join method
      void join( gener_prim_expr &x )
      {
         // reserve space in the vectors
         primints1.reserve(primints1.size()+x.primints1.size());
         primgrad1.reserve(primgrad1.size()+x.primgrad1.size());

         // concat the vectors
         primints1.insert(primints1.end(), x.primints1.begin(), x.primints1.end());
         primgrad1.insert(primgrad1.end(), x.primgrad1.begin(), x.primgrad1.end());

         // update the counters
         zero_pre+=x.zero_pre;
         zero_ang+=x.zero_ang;
         zero_rad+=x.zero_rad;
      };

      // do the actual work
      void operator() ( tbb::blocked_range2d<int,int> &range );
};

void ecp_prepare(int maxshell,
               double *alpha, double *cc, int *primofs, int *aoofs, int *shelltype,
               int *atid, int nunique, int norb, FILE* outf)
{
   // get a timer for this procedure
   timer preptime;
   preptime.punch("Total in prepare");

   const double pi=acos(-1.);

   /* uniquify exponents */
   // exp_ind vars maps alpha->uniq_exp
   const int nprimtot=primofs[nunique]-1;
   int exp_ind[nprimtot];

   fprintf(outf,"  Uniquifying basis exponents...");
   preptime.punch("Uniq exponents");
   for (int i=0; i<nprimtot; i++)
   {
      double test=alpha[i];

      int idx=0;
      for (; (idx<ecp_data.uniq_exp.size())&&(fabs(ecp_data.uniq_exp[idx]-test)>ECP_TOL); idx++);
      exp_ind[i]=idx;
      if (exp_ind[i]==ecp_data.uniq_exp.size()) ecp_data.uniq_exp.push_back(test);
   };
   fprintf(outf," done.\n\n");
   preptime.punch("Uniq exponents");

   // get all the unique ecp centers
   int dummy;
   uniquify<int> ecpatid(ecp_data.ecpcenters);
   for (vector<ecpdat>::iterator ecp=ecp_data.info.begin(); ecp!= ecp_data.info.end();
        ecp++) ecpatid.index(ecp->atid,dummy);

   /* perform dry run */
   fprintf(outf,"  Constructing primitive integral list...");
   preptime.punch("Gen. primitives");

   gener_prim_expr expr(exp_ind,cc,primofs,aoofs,shelltype,atid,nunique,norb);
   tbb::parallel_reduce( tbb::blocked_range2d<int,int>(0,nunique,0,nunique), expr);

   preptime.punch("Gen. primitives");
   fprintf(outf," done.\n\n");

   fprintf(outf,"  %i normal primitive expressions\n",expr.primints1.size());
   fprintf(outf,"  %i gradient primitive expressions\n",expr.primgrad1.size());
   fprintf(outf,"      %i removed by prefactor screening\n",expr.zero_pre);
   fprintf(outf,"      %i removed by angular symmetry\n",expr.zero_ang);
   fprintf(outf,"      %i removed by radial integral screening\n\n",expr.zero_rad);

   // optimize expressions
   preptime.punch("Strength reduction");
   ecp_expr_opt(expr.primints1, expr.primgrad1, outf);
   preptime.punch("Strength reduction");

   preptime.print("Timing inside prepare:",outf);
   fprintf(outf,"\n");
   fflush(outf);

   return;
};

// parallel reduction operator implementation
void gener_prim_expr::operator() ( tbb::blocked_range2d<int,int> &range )
{
   // constants and recasting
   const double pi = acos(-1.);
   int max_l=ecp_data.Ylm_maxl,
       max_IJK=ecp_data.max_IJK,
       max_nl_IJK=ecp_data.max_nl_IJK;
   double (*l_angtab)[max_IJK+1][max_IJK+1][max_IJK+1] = (double(*)[max_IJK+1][max_IJK+1][max_IJK+1]) ecp_data.local_ang;
   double (*nl_angtab)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1] =
       (double(*)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1]) ecp_data.nonlocal_ang;

   for (int i=range.rows().begin(); i!=range.rows().end(); i++)
   for (int j=range.cols().begin(); j!=range.cols().end(); j++)
   {
      // get i-params
      int orb_i=aoofs[i]-2,
          shella=shelltype[i]-1;

      // enumerate shells in function i
      for (int na=shella; na>=0; na--)
      for (int la=shella-na; la>=0; la--)
      {

         orb_i++;

         // cycle through the normal shell and its derivatives
         for (int dx=0; dx<3; dx++)
         for (int dir=-1; dir<2; dir++)
         {

         // only get one normal shell
         if ((dir==0)&&(dx>0)) continue;

         // get the new powers in a
         int pwrs_a[3] = {na,la,shella-na-la};
         pwrs_a[dx]+=dir;

         // skip unphysical -1 component
         if (pwrs_a[dx]<0) continue;

         // get j-params
         int orb_j=aoofs[j]-2,
             shellb=shelltype[j]-1;

         // enumerate shells
         for (int nb=shellb; (nb>=0)&&(((dir==0)&&(orb_i>=orb_j))||(dir!=0)); nb--)
         for (int lb=shellb-nb; lb>=0; lb--)
         {
            int mb=shellb-nb-lb,
                pwrs_b[3]={nb,lb,mb};
            orb_j++;
            if ((dir==0)&&(orb_j>orb_i)) break;

            // loop over ECPs
            for (vector<ecpdat>::iterator ecp=ecp_data.info.begin();
                 ecp!=ecp_data.info.end(); ecp++)
            {
               // don't compute known-zero gradients
               if ((dir!=0)&&(atid[i]==atid[j])&&((atid[i]-1)==ecp->atid)) continue;

               // get output idx
               // gradient element matrix <di/dr|V|j> computed first
               // dimensions [nat][norb(i)][norb(j)][3]
               int idx = ecp->atid*norbp+orb_i*(orb_i+1)/2+orb_j;
               if (dir!=0) idx = ecp->atid*norb*norb*3+orb_i*norb*3+orb_j*3+dx;

               // calculate mock CA and CB
               double eq[3] = {0.,0.,0.},
                      ne[3] = {1.,1.,1.};
               double *CA = ne,
                      *CB = ne;
               if ((atid[i]-1)==ecp->atid) CA=eq;
               if ((atid[j]-1)==ecp->atid) CB=eq;

               // loop over primitive gaussians
               for (int prima=primofs[i]-1; prima<(primofs[i+1]-1); prima++)
               for (int primb=primofs[j]-1; primb<(primofs[j+1]-1); primb++)
               {
                  // prefactor loop
                  for (int a=0; a<=pwrs_a[0]; a++)
                  for (int b=0; b<=pwrs_a[1]; b++)
                  for (int c=0; c<=pwrs_a[2]; c++)
                  for (int d=0; d<=nb; d++)
                  for (int e=0; e<=lb; e++)
                  for (int f=0; f<=mb; f++)
                  {
                     double prefactor=4*pi*cc[prima]*cc[primb]*ecp->coef;
                     if (dir==-1) prefactor*=-(pwrs_a[dx]+1);
                     if (dir==1) prefactor*=2*ecp_data.uniq_exp[exp_ind[prima]];

                     int idxa[3] = {a,b,c},
                         idxb[3] = {d,e,f};
                     for (int xyz=0; xyz<3; xyz++)
                     {
                        prefactor*=nCk(pwrs_a[xyz],idxa[xyz])*intpwr_mod(CA[xyz],pwrs_a[xyz]-idxa[xyz])
                                  *nCk(pwrs_b[xyz],idxb[xyz])*intpwr_mod(CB[xyz],pwrs_b[xyz]-idxb[xyz]);
                     };
                     if (prefactor==0.)
                     {
                        zero_pre++;
                        continue;
                     };

                     // generate local or nonlocal term
                     if (ecp->l==-1)
                     {
                        // local term
                        int lambda_max=a+b+c+d+e+f,
                            lambda_min=lambda_max%2;
                        for (int lambda=lambda_min; lambda<=lambda_max; lambda+=2)
                        {
                           // check for zero angular terms
                           bool zero=true;
                           for (int mu=-lambda; (mu<=lambda)&&zero; mu++)
                              zero=(l_angtab[lambda*(lambda+1)+lambda+mu][a+d][b+e][c+f]==0.);
                           if (zero)
                           {
                              zero_ang++;
                              continue;
                           };

                           // check for radial zero
                           if ((atid[i]==atid[j])&&((atid[i]-1)==ecp->atid)&&(lambda>0))
                           {
                              zero_rad++;
                              continue;
                           };

                           // add the primitive
                           dryrun_lvl1 prim =
                                {idx,
                                 -1, ecp->zeta, ecp->atid,
                                 prefactor, pwrs_a[0]-a, pwrs_a[1]-b, pwrs_a[2]-c, nb-d, lb-e, mb-f,
                                 atid[i]-1, exp_ind[prima],
                                 atid[j]-1, exp_ind[primb],
                                 a+d, b+e, c+f, lambda,
                                 -1, -1, -1, -1,
                                 lambda, -1, a+b+c+d+e+f+ecp->r};
                           if (dir==0) primints1.push_back(prim);
                           else primgrad1.push_back(prim);

                        };
                     } else {
                        // nonlocal term

                        prefactor*=4*pi;

                        int lambda_max=ecp->l+a+b+c,
                            lambda_min=ecp->l-a-b-c;
                        if (lambda_min<0) lambda_min=0;
                        if ((lambda_max-lambda_min)%2==1) lambda_min+=1;

                        int lambdaP_max=ecp->l+d+e+f,
                            lambdaP_min=ecp->l-d-e-f;
                        if (lambdaP_min<0) lambdaP_min=0;
                        if ((lambdaP_max-lambdaP_min)%2==1) lambdaP_min+=1;

                        for (int lambda=lambda_min; lambda<=lambda_max; lambda+=2)
                        for (int lambdaP=lambdaP_min; lambdaP<=lambdaP_max; lambdaP+=2)
                        {
                           // check for zero angular terms
                           bool zero=true;
                           for (int m=-ecp->l; (m<=ecp->l)&&zero; m++)
                           {
                              int zero1=true;
                              for (int mu=-lambda; (mu<=lambda)&&zero1; mu++)
                                 zero1=(nl_angtab[lambda*(lambda+1)+lambda+mu][ecp->l*(ecp->l+1)+ecp->l+m][a][b][c]==0.);
                              int zero2=true;
                              for (int mu=-lambdaP; (mu<=lambdaP)&&zero2; mu++)
                                 zero2=(nl_angtab[lambdaP*(lambdaP+1)+lambdaP+mu][ecp->l*(ecp->l+1)+ecp->l+m][d][e][f]==0.);
                              zero=zero1||zero2;
                           };
                           if (zero)
                           {
                              zero_ang++;
                              continue;
                           };

                           // check radial term
                           if (((atid[i]-1)==ecp->atid)&&(lambda>0)||
                               ((atid[j]-1)==ecp->atid)&&(lambdaP>0))
                           {
                              zero_rad++;
                              continue;
                           };

                           dryrun_lvl1 prim =
                                 {idx,
                                  ecp->l, ecp->zeta, ecp->atid,
                                  prefactor, pwrs_a[0]-a, pwrs_a[1]-b, pwrs_a[2]-c, nb-d, lb-e, mb-f,
                                  atid[i]-1, exp_ind[prima],
                                  atid[j]-1, exp_ind[primb],
                                  a, b, c, lambda,
                                  d, e, f, lambdaP,
                                  lambda, lambdaP, a+b+c+d+e+f+ecp->r};
                           if (dir==0) primints1.push_back(prim);
                           else primgrad1.push_back(prim);

                        };
                     };
                  }; // prefactor loop
               }; // primitive loop
            }; // ECP loop
         }; // b-shell loop
         }; // derivatives loop
      }; // a-shell loop
   }; // parallel loop

   return;
}; // end of method
