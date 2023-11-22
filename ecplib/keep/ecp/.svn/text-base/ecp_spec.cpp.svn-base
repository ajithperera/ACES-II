/* ecp_spec.cpp: initializes the ECP specification
 *
 * Written by Tom Grimes, 27-Jan 2010
 *
 */

// ecpdat is r,l,start,end,atid

#include "ecpints.h"

// global ECP specification
/*
class ecp_data_type {
public:
   vector<ecpdat> info;
   int maxshell;

   double *Ylm;
   int Ylm_maxl;

   double *local_ang;
   int max_IJK;

   double *nonlocal_ang;
   int max_nl_IJK;

   vector<dryrun_lvl2> primints2;
   vector<dryrun_lvl2> primgrad2;

   vector<double> uniq_exp;
   vector<k_struct> kvecs;
   vector<ang_struct> angints;
   vector<kvpair_struct> kvpairs;
   vector<rad_int> radints;

   vector<int> blocks;
   vector<int> gradblocks;
   vector<int> ecpcenters;

   double *angvals, *radvals, *screen, *CApows, (*kvec_vals)[4], (*kvp_vals)[4];
   int norm_klen, norm_kvplen, norm_anglen, norm_radlen;

   timer time;

};
*/

ecp_data_type ecp_data;

extern "C" {

void ecp_spec_(int &maxshell, const char *bas_str_, const char *libname_, int *nel_rm,
               double *alpha, double *cc, int *primofs, int *aoofs, int *shelltype,
               int *atid, int &nunique, int &norb,
               unsigned int baslen, unsigned int liblen)
{
   // begin timing
   timer ecptime;
   ecptime.punch("Total init");

   // reset data structures
   ecp_data.info.clear();
   ecp_data.primints2.clear();
   ecp_data.primgrad2.clear();
   ecp_data.uniq_exp.clear();

   ecp_data.kvecs.clear();
   ecp_data.norm_klen = 0;

   ecp_data.angints.clear();
   ecp_data.norm_anglen = 0;

   ecp_data.kvpairs.clear();
   ecp_data.norm_kvplen = 0;

   ecp_data.radints.clear();
   ecp_data.norm_radlen = 0;

   ecp_data.blocks.clear();
   ecp_data.gradblocks.clear();
   ecp_data.ecpcenters.clear();

   delete [] ecp_data.Ylm;
   delete [] ecp_data.local_ang;
   delete [] ecp_data.nonlocal_ang;

   delete [] ecp_data.angvals;
   delete [] ecp_data.radvals;
   delete [] ecp_data.screen;
   delete [] ecp_data.CApows;
   delete [] ecp_data.kvec_vals;
   delete [] ecp_data.kvp_vals;

   ecp_data.Ylm_maxl=0;
   ecp_data.maxshell=(maxshell+1>2 ? maxshell+1 : 2);

   double *angvals, *radvals, *CApows, (*kvec_vals)[4], (*kvp_vals)[4];

   // initialize output data connection
   // and write header
   const char *ang[] = {"local","s","p","d","f","g","h","i"};
   FILE *outf = fdopen(output_fd,"a");
   if (outf == NULL) { printf("Error connecting to output file!\n"); exit(1); };
   fprintf(outf,"  *** ECP module initialization ***\n\n");
   fprintf(outf,"  ECP data read:\n\n");
   fprintf(outf,"  Atom  Num  Basis       L      N   Coefficient      Exponent      \n");
   fprintf(outf,"  -----------------------------------------------------------------\n");

   // load up known ECPs
   ecptime.punch("Load ecp");
   map<string, map<string, vector<ecpdat> > > ecp_lib;
   map<string, map<string, int> > ecp_core;
   load_ecps(ecp_lib, ecp_core, string(libname_,liblen));

   // parse the basis spec
   // note that lowercase is used exclusively
   vector<string> tokens;
   int atnum = 0;
   string bas_str(bas_str_,baslen);
   transform(bas_str.begin(),bas_str.end(),bas_str.begin(),static_cast<int(*)(int)>(tolower));
   Tokenize(bas_str,tokens);
   for (vector<string>::iterator iter=tokens.begin(); iter!=tokens.end(); iter++)
   {
      // loop vars
      int num;
      string atom, bset;

      // decompose token of form [#]At:Basis
      vector<string> tok2;
      Tokenize(*iter,tok2,":");
      istringstream isstr;
      isstr.str(tok2[0]);

      if (!(isstr >> num >> atom))
      {
         isstr.clear();
         num=1;
         isstr >> atom;
      };
      bset=tok2[1];

      // add the ECPs
      map<string,map<string,vector<ecpdat> > >::iterator map_basis = ecp_lib.find(bset);
      if (map_basis!=ecp_lib.end())
      {
         map<string,vector<ecpdat> >::iterator map_atom = map_basis->second.find(atom);
         if (map_atom!=map_basis->second.end())
         {
            for (int cpy=0; cpy<num; cpy++)
            {
               // record the core electrons
               nel_rm[atnum]=ecp_core[bset][atom];

               // copy basis
               for (vector<ecpdat>::iterator ecpiter=map_atom->second.begin();
                    ecpiter!=map_atom->second.end(); ecpiter++)
               {
                  // copy data from library, updating atom numbers
                  ecpdat addecp = *(ecpiter);
                  addecp.atid=atnum;
                  ecp_data.info.push_back(addecp);

                  // print data
                  fprintf(outf,"  %-4s  %3i  %-10s  %-5s  %1i  % 10.8e  % 10.8e\n",
                                atom.c_str(), atnum+1, bset.c_str(), ang[addecp.l+1],
                                addecp.r, addecp.coef, addecp.zeta);
               };
               atnum++;
            };
         } else {
            // atom not found in basis
            atnum+=num;

            fprintf(outf,"  %s%i-%s%i -- not in set %s --\n",atom.c_str(),atnum-num+1,atom.c_str(),atnum,bset.c_str());

         };
      } else {
         // ECP basis not found
         atnum+=num;

         fprintf(outf,"  %s%i-%s%i -- basis %s not found --\n",atom.c_str(),atnum-num+1,atom.c_str(),atnum,bset.c_str());
      };
   };
   fprintf(outf,"\n");
   ecptime.punch("Load ecp");

   // get table dimensions
   // keep in mind that maxshell is in normal convention (s=0)
   ecptime.punch("Precompute const");
   int max_l = 0,
       necp=ecp_data.info.size();
   for (vector<ecpdat>::iterator iter=ecp_data.info.begin(); iter!=ecp_data.info.end(); iter++)
      max_l=(iter->l>max_l ? iter->l : max_l);
   max_l+=2*maxshell+1;   // maximum sph harm (Ylm) needed

   // build the Ylm coefficients
   ecp_data.Ylm = new double[((max_l+1)*(max_l+2)+1)*(max_l+1)*(max_l+1)*(max_l+1)];
   memset((void*)ecp_data.Ylm,0,((max_l+1)*(max_l+2)+1)*(max_l+1)*(max_l+1)*(max_l+1)*sizeof(double));
   ecp_data.Ylm_maxl=max_l;
   double (*Ylm)[max_l+1][max_l+1][max_l+1] = (double(*)[max_l+1][max_l+1][max_l+1]) ecp_data.Ylm;
   const double pi = acos(-1.0),
                cosfact[4]={1.,0.,-1.,0.};

   for (int l=0; l<=max_l; l++)
   {
      // compute the coefficients
      for (int m=-l; m<=l; m++)
      {
         // get first index
         int idx=l*(l+1)+l+m;

         int abs_m=abs(m);
         double Nlm=sqrt(fact(l-abs_m)*(2*l+1)/(4*pi*fact(l+abs_m)));
         if (m!=0) Nlm*=sqrt(2.);
         for (int k=0; k<=(l-abs_m)/2; k++)
         {
            double gnorm=Nlm*nCk(l,k)*nCk(2*(l-k),l)*fact(l-2*k)/((1<<l)*fact(l-2*k-abs_m));
            if (k%2==1) gnorm=-gnorm;

            for (int p=0; p<=abs_m; p++)
            {
               double coef=gnorm*nCk(abs_m,p);

               if (m>=0) { coef*=cosfact[(abs_m-p)%4]; }
               else { coef*=cosfact[(abs_m-p+3)%4]; };

               Ylm[idx][p][abs_m-p][l-2*k-abs_m]+=coef;
            };
         };
      };
   };

   // build the monomial integral table: int(xyz (hat))
   int mon_len=2*max_l+1;
   double monomials[mon_len][mon_len][mon_len];
   for (int i=0; i<mon_len; i++)
   for (int j=0; j<mon_len; j++)
   for (int k=0; k<mon_len; k++)
   {
      if ((i%2==1)||(j%2==1)||(k%2==1))
      {
         monomials[i][j][k]=0.;
      } else {
         monomials[i][j][k]=4*pi*fact2(i-1)*fact2(j-1)*fact2(k-1)/fact2(i+j+k+1);
      };
   };
   
   // build local angular integrals
   // address as [lambda(lambda+1)+lambda+mu][I][J][K]
   int max_IJK=2*(maxshell+1);
   ecp_data.local_ang = new double[((max_l+1)*(max_l+2)+1)*(max_IJK+1)*(max_IJK+1)*(max_IJK+1)];
   memset((void*)ecp_data.local_ang,0,((max_l+1)*(max_l+2)+1)*(max_IJK+1)*(max_IJK+1)*(max_IJK+1)*sizeof(double));
   ecp_data.max_IJK=max_IJK;
   double (*l_angtab)[max_IJK+1][max_IJK+1][max_IJK+1] = (double(*)[max_IJK+1][max_IJK+1][max_IJK+1]) ecp_data.local_ang;

   for (int lambda=0; lambda<=max_l; lambda++)
   {
      // compute the terms
      for (int mu=-lambda; mu<=lambda; mu++)
      {
         // get first index
         int idx=lambda*(lambda+1)+lambda+mu;

         for (int nx=0; nx<=lambda; nx++)
         for (int ny=0; (nx+ny)<=lambda; ny++)
         for (int nz=lambda-nx-ny; nz>=0; nz-=2)
         for (int I=0; I<=max_IJK; I++)
         for (int J=0; J<=max_IJK; J++)
         for (int K=0; K<=max_IJK; K++)
         {
            l_angtab[idx][I][J][K]+=Ylm[idx][nx][ny][nz]*monomials[nx+I][ny+J][nz+K];
         };
      };
   };

   // build nonlocal angular integrals
   // address as [lambda(lambda+1)+lambda+mu][l(l+1)+l+m][I][J][K]
   int max_nl_IJK=maxshell+1;
   ecp_data.nonlocal_ang = 
          new double[((max_l+1)*(max_l+2)+1)*((max_l+1)*(max_l+2)+1)*(max_nl_IJK+1)*(max_nl_IJK+1)*(max_nl_IJK+1)];
   memset((void*)ecp_data.nonlocal_ang,0,
          ((max_l+1)*(max_l+2)+1)*((max_l+1)*(max_l+2)+1)*(max_nl_IJK+1)*(max_nl_IJK+1)*(max_nl_IJK+1)*sizeof(double));
   ecp_data.max_nl_IJK=max_nl_IJK;

   double (*nl_angtab)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1] = 
       (double(*)[(max_l+1)*(max_l+2)+1][max_nl_IJK+1][max_nl_IJK+1][max_nl_IJK+1]) ecp_data.nonlocal_ang;
   for (int lambda=0; lambda<=max_l; lambda++)
   for (int l=0; l<=max_l; l++)
   {
      for (int mu=-lambda; mu<=lambda; mu++)
      for (int m=-l; m<=l; m++)
      {
         // calculate indices
         int lambda_idx=lambda*(lambda+1)+lambda+mu,
             l_idx=l*(l+1)+l+m;

         for (int I=0; I<=max_nl_IJK; I++)
         for (int J=0; J<=max_nl_IJK; J++)
         for (int K=0; K<=max_nl_IJK; K++)
         for (int nx=0; nx<=lambda; nx++)
         for (int ny=0; (nx+ny)<=lambda; ny++)
         for (int nz=lambda-nx-ny; nz>=0; nz-=2)
         for (int nxp=0; nxp<=l; nxp++)
         for (int nyp=0; (nxp+nyp)<=l; nyp++)
         for (int nzp=l-nxp-nyp; nzp>=0; nzp-=2)
         {
            nl_angtab[lambda_idx][l_idx][I][J][K]+=Ylm[lambda_idx][nx][ny][nz]*Ylm[l_idx][nxp][nyp][nzp]
                                                   *monomials[nx+nxp+I][ny+nyp+J][nz+nzp+K]; 
         };
      };
   };
   ecptime.punch("Precompute const");

   // prepare the integral expressions
   // flush the stream because prepare also writes to it
   fflush(outf);
   ecptime.punch("Calling prepare");
   ecp_prepare(maxshell,alpha,cc,primofs,aoofs,shelltype,atid,nunique,norb,outf);
   ecptime.punch("Calling prepare");

   // if using CUDA, allocate memory and transfer radial integral list
#ifdef USING_CUDA
   ecp_init_cuda();
#endif

   // print out our level of timing info
   ecptime.print("Overall ECP initialization timing:",outf);
   fprintf(outf,"\n");

   // write closing message and flush the stream
   fprintf(outf,"  *** ECP module finished initializing ***\n\n");
   fflush(outf);

   return;
};

}; // extern C

