#include "flags.h"
      SUBROUTINE GET_PERTURB_RSPNS_DENS(DOO,DVV,WORK,MAXCOR,IRREPX,IUHF,
     +                                  LENVV,LENOO)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      Integer T2ln_aa,T2ln_bb,T2ln_ab,T2ln
      Integer T1ln_aa,T1ln_bb,T1ln
C
      DIMENSION WORK(MAXCOR),LENVV(2),LENOO(2)
      DIMENSION DOO(LENOO(1)+IUHF*LENOO(2))
      DIMENSION DVV(LENVV(1)+IUHF*LENVV(2))
C
#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"
C
      DATA IONE, ONE /1, 1.0D0/    

      Listr1 = 490
      Listr2 = 444 

      T1ln_aa = Irpdpd(Irrepx,9)
      If (Iuhf .Ne. 0) T1ln_bb =Irpdpd(Irrepx,10)
      T2ln_aa = Idsymsz(Irrepx,1,3)
      If (Iuhf .Ne. 0) T2ln_bb = Idsymsz(Irrepx,2,4)
      T2ln_ab = Idsymsz(Irrepx,13,14)

      T1ln = T1ln_aa + Iuhf*T1ln_bb
      T2ln = T2ln_aa + Iuhf*T2ln_bb + T2ln_ab

      I000 = Ione
      I010 = I000 + T1ln
      Iend = I010 + T2ln 

      If (Iend .Gt. Maxcor) Call Insmem(get_perturb_rspn_dens,Iend,
     +                                  Maxcor)
      Maxcor = Maxcor - Iend

      Call Getlst(Work(I000),1,1,1,1,Listr1)
      call checksum("T1aa",Work(I000),T1ln_aa)

      Ioff = T1ln_aa 
      If (Iuhf .Ne. 0) Call Getlst(Work(I000+Ioff),1,1,1,2,Listr1)
      call checksum("T1aa",Work(I000+Ioff),T1ln_bb)

      Call Getall(Work(I010),T2ln_aa,Irrepx,Listr2)

      Ioff = Iuhf*T2ln_aa
      If (Iuhf .Ne. 0) Call Getall(Work(I010+Ioff),T2ln_bb,Irrepx,
     +                             Listr2+1)

      Ioff = T2ln_aa + Iuhf*T2ln_bb
      Call Getall(Work(I010+Ioff),T2ln_ab,Irrepx,Listr2+2)

      Call Dzero(Doo,Lenoo(1)+Iuhf*Lenoo(2))
      Call Dzero(Dvv,LenvV(1)+Iuhf*Lenvv(2))

      Call T2in_densoo(Doo,Work(I010),Work(Iend),Maxcor,T2ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,Irrepx,Iuhf) 
      Call T2in_densvv(Dvv,Work(I010),Work(Iend),Maxcor,T2ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,Irrepx,Iuhf) 

      Call T1in_densoo(Doo,Work(I000),Work(Iend),Maxcor,T1ln,T1ln_aa,
     +                 T1ln_bb,Irrepx,Iuhf) 
      Call T1in_densvv(Dvv,Work(I000),Work(Iend),Maxcor,T1ln,T1ln_aa,
     +                 T1ln_bb,Irrepx,Iuhf) 

      Return
      End 



