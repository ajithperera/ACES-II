
#ifndef _SCFFLAGS_COM_
#define _SCFFLAGS_COM_

c This file contains flags that are set in the SCF namelist.  See the file
c initscfflags.F for a description of each of them.

c The following flags are not read in but keep track of some useful
c information:
c
c     scflastiter   : true if this is the final iteration of the SCF
c                     procedure (the final fock matrix is being
c                     calculated)
c     scfksiter     : true if the fock matrix should be calculated using
c                     a Kohn-Sham procedure this iteration
c     scfksfirsiter : true if this is the first Kohn-Sham iteration
c     scfkslastiter : true if this is the last Kohn-Sham iteration

      integer scfprint,scfmaxit,scftdlst,scfrppord,
     &    scfrppbeg,scfdamptyp,
     &    scfprintevec,scfprintevec_i,scfprinteval,scfprinteval_i
      M_REAL
     &    scfconvtol,scfstadamp,scfalpha1,
     &    scfbeta1,scfdavdamp
      logical scfaofil,scfdamp,scfks,scfrpp,scflastiter,scfksiter,
     &    scfksfirsiter,scfkslastiter

      common /scfflags/  scfprint,scfmaxit,scftdlst,
     &    scfrppord,scfrppbeg,scfdamptyp,
     &    scfprintevec,scfprintevec_i,scfprinteval,scfprinteval_i
      common /scfflagsd/ scfconvtol,scfdavdamp,scfstadamp,
     &    scfalpha1,scfbeta1
      common /scfflagsl/ scfaofil,scfdamp,scfrpp,scfks,
     &    scflastiter,scfksiter,scfksfirsiter,scfkslastiter

      save /scfflags/
      save /scfflagsd/
      save /scfflagsl/

c The following parameters are used in the namelist.

      integer scfeig_none,scfeig_occ,scfeig_ov,scfeig_occv,scfeig_all
      parameter (scfeig_none     =1)
      parameter (scfeig_occ      =2)
      parameter (scfeig_ov       =3)
      parameter (scfeig_occv     =4)
      parameter (scfeig_all      =5)

      integer scfdamp_none,scfdamp_davidson,scfdamp_static
      parameter (scfdamp_none    =1)
      parameter (scfdamp_davidson=2)
      parameter (scfdamp_static  =3)

#endif /* _SCFFLAGS_COM_ */

