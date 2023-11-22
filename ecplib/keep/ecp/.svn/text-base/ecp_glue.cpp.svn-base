/* ecp_glue.cpp
 *
 * Written by Tom Grimes, 25-Jan 2010
 *
 * This file is the partner of add_ecp.F:
 * Given the basis set as exposed by add_ecp, generate the integral
 * batches to be handled by the actual ECP routines.
 *
 * Beware of the static global ecp specification!

 * Parallelism lives here.
 *
 */

// turn on this define to report balancing
//define BALANCE_VERBOSE

#include "ecpints.h"

// prototype: eventually move to header
void compute_kvecs(double *coords, int &nat, int maxdim, int nkvecs, int nkvp);
void balance_radial(int lower, int upper);
void balance_angular(int lower, int upper);

extern "C" {

// compute normal integrals
void ecp_glue_(double *coords, double *output, int &norb, int &nat)
{ 
   ecp_data.time.punch("Normal ints");

   // useful constants
   const int norbp=norb*(norb+1)/2,
             maxdim=ecp_data.maxshell+1,
             nkvecs=ecp_data.kvecs.size(),
             nkvp=ecp_data.kvpairs.size();

   // connect data
   compute_kvecs(coords, nat, maxdim, nkvecs, nkvp);

   // get angular and radial integrals
   ecp_data.time.punch("Normal ints:angular");

#ifndef USING_CUDA
   tbb::parallel_for( tbb::blocked_range<int>(0,ecp_data.norm_anglen),
                      ang_ints_driver(&(ecp_data.angints.front()),ecp_data.kvec_vals,
                                      ecp_data.kvp_vals,ecp_data.angvals) );

#else
   balance_angular(0,ecp_data.norm_anglen);

#endif

/* DEBUG */
/*
   tbb::parallel_for( tbb::blocked_range<int>(0,ecp_data.norm_anglen),
                      ang_ints_driver(&(ecp_data.angints.front()),ecp_data.kvec_vals,
                                      ecp_data.kvp_vals,ecp_data.angvals) );
   double angvals[ecp_data.norm_anglen];
   angular_cuda(0, ecp_data.norm_anglen, angvals);
   for (int i=0; i<ecp_data.norm_anglen; i++)
      if ((fabs(ecp_data.angvals[i]-angvals[i])/ecp_data.angvals[i])>1.e-4)
      printf("%f %f\n",angvals[i],ecp_data.angvals[i]);
   exit(1);
*/

   ecp_data.time.punch("Normal ints:angular");

   ecp_data.time.punch("Normal ints:radial");


#ifndef USING_CUDA
   tbb::parallel_for( tbb::blocked_range<int>(0,ecp_data.norm_radlen),
                      rad_ints_driver(&(ecp_data.radints.front()),&(ecp_data.kvecs.front()),
                                      ecp_data.kvec_vals,ecp_data.kvp_vals,
                                      ecp_data.radvals) );

#else
   balance_radial(0,ecp_data.norm_radlen);

#endif

/* DEBUG SECTION */
/*
   tbb::parallel_for( tbb::blocked_range<int>(0,ecp_data.norm_radlen),
                      rad_ints_driver(&(ecp_data.radints.front()),&(ecp_data.kvecs.front()),
                                      ecp_data.kvec_vals,ecp_data.kvp_vals,
                                      ecp_data.radvals) );

   double radvals[ecp_data.norm_radlen];
   radial_cuda( 0, ecp_data.norm_radlen, radvals);
   for (int i=0; i<ecp_data.norm_radlen; i++)
   {
      if ((fabs(radvals[i]-ecp_data.radvals[i])/ecp_data.radvals[i]>1.e-5))
      {
         printf("%e %e %e\n",radvals[i],ecp_data.radvals[i],
                fabs(radvals[i]-ecp_data.radvals[i])/ecp_data.radvals[i]);
      };
   };
   exit(1);
*/

   ecp_data.time.punch("Normal ints:radial");

   // finish computation in primints2
   ecp_data.time.punch("Normal ints:combine");
   tbb::parallel_for( tbb::blocked_range<int>(0,ecp_data.blocks.size()-1),
                      combine_expr_driver(&(ecp_data.primints2.front()),&(ecp_data.blocks.front()),
                                          ecp_data.CApows,ecp_data.angvals,ecp_data.radvals,ecp_data.screen,
                                          output,nat,maxdim) );
   ecp_data.time.punch("Normal ints:combine");

   ecp_data.time.punch("Normal ints");
   return;
};

// add derivatives directly to passed structure
// output is (norbp,nat,6)<->[6][nat][norbp]
// The last two indices are packed and the order of the first index
// is {i(x,y,z) V_c(x,y,z)} on center c in {nat}
void ecp_glue_deriv_(double *coords, double *output, int &norb, int &nat)
{
   ecp_data.time.punch("Grad ints");

   // useful constants
   const int norbp=norb*(norb+1)/2,
             maxdim=ecp_data.maxshell+1,
             nkvecs=ecp_data.kvecs.size(),
             nkvp=ecp_data.kvpairs.size();

   // ensure data is connected
   int ang_lower = ecp_data.norm_anglen,
       rad_lower = ecp_data.norm_radlen;
   if (ecp_data.CApows==NULL)
   {
      compute_kvecs(coords,nat,maxdim,nkvecs,nkvp);
      ang_lower=0;
      rad_lower=0;
   };

   // compute the remaining angular and radial integrals
   ecp_data.time.punch("Grad ints:angular");

#ifndef USING_CUDA
   tbb::parallel_for( tbb::blocked_range<int>(ang_lower,ecp_data.angints.size()),
                      ang_ints_driver(&(ecp_data.angints.front()),ecp_data.kvec_vals,
                                      ecp_data.kvp_vals,ecp_data.angvals) );

#else
   balance_angular(ang_lower,ecp_data.angints.size());

#endif

   ecp_data.time.punch("Grad ints:angular");

   ecp_data.time.punch("Grad ints:radial");

#ifndef USING_CUDA
   tbb::parallel_for( tbb::blocked_range<int>(rad_lower,ecp_data.radints.size()),
                      rad_ints_driver(&(ecp_data.radints.front()),&(ecp_data.kvecs.front()),
                                      ecp_data.kvec_vals,ecp_data.kvp_vals,
                                      ecp_data.radvals) );

#else
   balance_radial(rad_lower,ecp_data.radints.size());

#endif

   ecp_data.time.punch("Grad ints:radial");

   // finish combining expressions
   ecp_data.time.punch("Grad ints:combine");
   tbb::parallel_for( tbb::blocked_range<int>(0,ecp_data.gradblocks.size()-1),
                      combine_grad_driver(&(ecp_data.primgrad2.front()),&(ecp_data.gradblocks.front()),
                                          ecp_data.CApows,ecp_data.angvals,ecp_data.radvals,ecp_data.screen,
                                          (double*)output,nat,norb,maxdim) );
   ecp_data.time.punch("Grad ints:combine");

   // consume the integral lists as a safeguard against
   // changing geometry before recomputing gradients w/o normal integrals first
   delete [] ecp_data.CApows;
   delete [] ecp_data.kvec_vals;
   delete [] ecp_data.kvp_vals;
   delete [] ecp_data.angvals;
   delete [] ecp_data.radvals;
   ecp_data.CApows=NULL;
   ecp_data.kvec_vals=NULL;
   ecp_data.kvp_vals=NULL;
   ecp_data.angvals=NULL;
   ecp_data.radvals=NULL;

   ecp_data.time.punch("Grad ints");

   return;
};

}; // extern C

// this procedure makes R and the k-vectors available
void compute_kvecs(double *coords, int &nat, int maxdim, int nkvecs, int nkvp)
{
   ecp_data.time.punch("K-vectors");

   delete [] ecp_data.CApows;
   delete [] ecp_data.kvec_vals;
   delete [] ecp_data.kvp_vals;
   delete [] ecp_data.angvals;
   delete [] ecp_data.radvals;
   delete [] ecp_data.screen;

   // angular and radial integrals are not actually computed here
   // but they need to be allocated
   ecp_data.kvec_vals = new double[nkvecs][4];
   ecp_data.kvp_vals = new double[nkvp][4];
   ecp_data.angvals = new double[ecp_data.angints.size()];
   ecp_data.radvals = new double[ecp_data.radints.size()];
   ecp_data.screen = new double[nkvp];
   ecp_data.CApows = new double[nat*nat*3*maxdim];

   double (*CApows)[nat][3][maxdim] = (double(*)[nat][3][maxdim]) ecp_data.CApows;
   double (*kvec_vals)[4] = ecp_data.kvec_vals;
   double (*kvp_vals)[4] = ecp_data.kvp_vals;
   double *angvals = ecp_data.angvals;
   double *radvals = ecp_data.radvals;

   // first compute CA^x
   double (*R)[3] = (double(*)[3]) coords;
   for (int a=0; a<nat; a++)
   for (int c=0; c<nat; c++)
   for (int xyz=0; xyz<3; xyz++)
   {
      CApows[c][a][xyz][0]=1.;
      CApows[c][a][xyz][1]=R[c][xyz]-R[a][xyz];
      for (int p=2; p<maxdim; p++)
      {
         CApows[c][a][xyz][p]=CApows[c][a][xyz][p-1]*CApows[c][a][xyz][1];
      };
   };

   // compute single k-vectors
   for (int i=0; i<nkvecs; i++)
   {
      k_struct k=ecp_data.kvecs[i];
      if (k.a!=k.c)
      {
         double klen=0;
         for (int xyz=0; xyz<3; xyz++)
         {
            kvec_vals[i][xyz]=-2*ecp_data.uniq_exp[k.alpha]*CApows[k.c][k.a][xyz][1];
            klen+=kvec_vals[i][xyz]*kvec_vals[i][xyz];
         };
         klen=sqrt(klen);
         kvec_vals[i][3]=klen;
         if (klen>0)
         {
            kvec_vals[i][0]/=klen;
            kvec_vals[i][1]/=klen;
            kvec_vals[i][2]/=klen;
         };
      } else {
         kvec_vals[i][0]=0.;
         kvec_vals[i][1]=0.;
         kvec_vals[i][2]=0.;
         kvec_vals[i][3]=0.;
      };
   };

   // compute pairwise k-vectors
   // also, screening term exp(-alpha_a|CA|^2-alpha_b|CB|^2)
   for (int i=0; i<nkvp; i++)
   {
      kvpair_struct k=ecp_data.kvpairs[i];

      if ((k.a!=k.c)||(k.b!=k.c))
      {
         double klen=0;
         for (int xyz=0; xyz<3; xyz++)
         {
            kvp_vals[i][xyz]=-2*ecp_data.uniq_exp[k.alpha_a]*CApows[k.c][k.a][xyz][1]
                             -2*ecp_data.uniq_exp[k.alpha_b]*CApows[k.c][k.b][xyz][1];
            klen+=kvp_vals[i][xyz]*kvp_vals[i][xyz];
         };
         klen=sqrt(klen);
         kvp_vals[i][3]=klen;
         if (klen>0)
         {
            kvp_vals[i][0]/=klen;
            kvp_vals[i][1]/=klen;
            kvp_vals[i][2]/=klen;
         };

         // screening term
         double arg=ecp_data.uniq_exp[k.alpha_a]*
                    (CApows[k.c][k.a][0][2]+CApows[k.c][k.a][1][2]+CApows[k.c][k.a][2][2])
                   +ecp_data.uniq_exp[k.alpha_b]*
                    (CApows[k.c][k.b][0][2]+CApows[k.c][k.b][1][2]+CApows[k.c][k.b][2][2]);
         ecp_data.screen[i]=exp(-arg);

      } else {
         kvp_vals[i][0]=0.;
         kvp_vals[i][1]=0.;
         kvp_vals[i][2]=0.;
         kvp_vals[i][3]=0.;
         ecp_data.screen[i]=1.;
      };
   };

   // if using CUDA, copy over the k-vec data
#ifdef USING_CUDA
   cuda_load_kvecs((double*) kvec_vals, nkvecs,
                   (double*) kvp_vals, nkvp);
#endif

   ecp_data.time.punch("K-vectors");

   return;
};

// balance radial integral load between CUDA and CPU
// run tbb threads from another thread
// (CUDA must use the main host thread to keep device context)
class cpu_thread_rad {
   int lower, upper;
public:
   double *t;
   cpu_thread_rad(int lower_, int upper_, double &t_):
                 lower(lower_), upper(upper_), t(&t_) {};
   void operator()()
   {
      tbb::tick_count t0 = tbb::tick_count::now();
      tbb::parallel_for( tbb::blocked_range<int>(lower,upper),
                      rad_ints_driver(&(ecp_data.radints.front()),&(ecp_data.kvecs.front()),
                                      ecp_data.kvec_vals,ecp_data.kvp_vals,
                                      ecp_data.radvals) );
      tbb::tick_count t1 = tbb::tick_count::now();
      *t=(t1-t0).seconds();
      return;
   };
};

#ifndef USING_CUDA
void balance_radial(int lower, int upper) {};

#else
void balance_radial(int lower, int upper)
{
   const int nkvecs=ecp_data.kvecs.size(),
             nkvp=ecp_data.kvpairs.size(),
             nradint = upper-lower;

   // debug!!
   const int nstep = 1;
   static int counter = nstep;

   // balancing params
   static double p_cuda = 0.5;
   int ncuda = p_cuda*nradint;

   // start TBB in its own thread
   double cpu_t;
   cpu_thread_rad  cpu_tobj(lower+ncuda, upper, cpu_t);
   tbb::tbb_thread t_cpu(cpu_tobj);

   // run CUDA
   tbb::tick_count cuda_t0 = tbb::tick_count::now();
   if (ncuda>0)
   {
      radial_cuda(lower, lower+ncuda, ecp_data.radvals);
   };
   double cuda_t=(tbb::tick_count::now()-cuda_t0).seconds();

   // wait for TBB
   t_cpu.join();

   // re-balance for next call (nc = num cuda, tc = time cuda, etc.)
   // new cuda = (3*old + measured)/4
   double nctp = ncuda*cpu_t,
          nptc = (nradint-ncuda)*cuda_t;
   if (((nptc+nctp)>0)&&(ncuda>0)) p_cuda = (nctp/(nptc+nctp)+3*p_cuda)/4;

#ifdef BALANCE_VERBOSE
   if (counter==0)
   {
      printf("Balancing params (radial): CUDA %e s, CPU %e s, new p_cuda %f\n",
           cuda_t,cpu_t,p_cuda);
      counter=nstep;
   } else {
      counter--;
   };
#endif

   return;
};

#endif

// balance angular integral load between CUDA and CPU
// run tbb threads from another thread
// (CUDA must use the main host thread to keep device context)
class cpu_thread_ang {
   int lower, upper;
public:
   double *t;
   cpu_thread_ang(int lower_, int upper_, double &t_):
                 lower(lower_), upper(upper_), t(&t_) {};
   void operator()()
   {
      tbb::tick_count t0 = tbb::tick_count::now();
      tbb::parallel_for( tbb::blocked_range<int>(lower,upper),
                      ang_ints_driver(&(ecp_data.angints.front()),ecp_data.kvec_vals,
                                      ecp_data.kvp_vals,ecp_data.angvals) );
      tbb::tick_count t1 = tbb::tick_count::now();
      *t=(t1-t0).seconds();
      return;
   };
};

#ifndef USING_CUDA
void balance_angular(int lower, int upper) {};

#else
void balance_angular(int lower, int upper)
{
   int nangint = upper-lower;

   // debug!!
   const int nstep = 1;
   static int counter = nstep;

   // balancing params
   static double p_cuda = 0.5;
   int ncuda = p_cuda*nangint;

   // start TBB in its own thread
   double cpu_t;
   cpu_thread_ang  cpu_tobj(lower+ncuda, upper, cpu_t);
   tbb::tbb_thread t_cpu(cpu_tobj);

   // run CUDA
   tbb::tick_count cuda_t0 = tbb::tick_count::now();
   if (ncuda>0)
   {
      angular_cuda(lower, lower+ncuda, ecp_data.angvals);
   };
   double cuda_t=(tbb::tick_count::now()-cuda_t0).seconds();

   // wait for TBB
   t_cpu.join();

   // re-balance for next call (nc = num cuda, tc = time cuda, etc.)
   // new cuda = (3*old + measured)/4
   double nctp = ncuda*cpu_t,
          nptc = (nangint-ncuda)*cuda_t;
   if (((nptc+nctp)>0)&&(ncuda>0)) p_cuda = (nctp/(nptc+nctp)+3*p_cuda)/4;

#ifdef BALANCE_VERBOSE
   if (counter==0)
   {
      printf("Balancing params (angular): CUDA %e s, CPU %e s, new p_cuda %f\n",
           cuda_t,cpu_t,p_cuda);
      counter=nstep;
   } else {
      counter--;
   };
#endif

   return;
};

#endif

