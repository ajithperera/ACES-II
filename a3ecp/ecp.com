#ifndef _ECP_COM_
#define _ECP_COM_
C
C This file contain all the ECP variables that need to be known
C across multiple files.
C
      common /ECP_INT_VARS/Zlm(Lmnpwr), Lmnval(3,84),
     &                     Istart(0:Maxang),Iend(0:Maxang),
     &                     Ideg(0:Maxang),Lmf(Maxangpwr),
     &                     Lml(Maxangpwr),
     &                     Lmx(Lmnpwr),Lmy(Lmnpwr),Lmz(Lmnpwr),
     &                     Pi,Fpi,Sqpi2,Sqrt_Fpi,R_intcutoff

      common/ECP_POT_VARS/clp(Mxecpprim),zlp(Mxecpprim),
     &                    nlp(Mxecpprim),kfirst(Maxang,Max_centers),
     &                    klast(Maxang,Max_centers),
     &                    llmax(Max_centers)

      common /pseud / nelecp(Max_centers),ipseux(Max_centers),ipseud 

      common /nshel / expnt(Max_prims),contr(Max_prims,Max_prims),
     &                numcon(Max_prims),katom(max_shells),
     &                ktype(max_shells),
     &                kprim(max_shells),kbfn(max_shells),
     &                kmini(max_shells),
     &                kmaxi(max_shells),nprims(max_shells),
     &                ndegen(max_shells),
     &                nshell,nbf

      Common /Qstore/Alpha,Beta,Xval
     
      Common /RadAng_sums/Rad_Sum(Maxang,Maxang), 
     &                    Ang_sum(Maxang,Maxang)
   
      Common /Fints/Fijk(0:4*Maxang,0:4*Maxang,0:4*Maxang)

      common /factorials/Fact(0:2*Maxang),Fac2(-1:4*Maxang),
     &                   Faco(0:2*Maxang),
     &                   Bcoefs(0:2*Maxang,0:2*Maxang),
     &                   Fprod(2*Maxang, 2*maxang)
  
#endif /* _ECP_COM_ */

