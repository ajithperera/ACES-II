#include "flags.h"
      Subroutine Get_2ndord_t1s(T2amp,Work,Maxcor,Tln,T2ln,T1ln,
     +                          T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +                          T1ln_bb,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Integer T2ln,T2ln_aa,T2ln_bb,T2ln_ab
      Integer T1ln,T1ln_aa,T1ln_bb
      Integer Tln
      Integer T1off
      Dimension Work(Maxcor),T2amp(Tln)

      Data Ione /1/

#include "machsp.com"
#include "syminf.com"
#include "sym.com"
#include "sympop.com"

      Ioff_t1 = T2ln + Ione
      Ioff_t2 = Ione

      Call Dzero(T2amp(Ioff_t1),T1ln_aa)
      If (Iuhf .Eq. 0) Call Dzero(T2amp(Ioff_t1+T1ln_aa),T1ln_bb)

      If (Iuhf .Eq. 1) Then 
         Call T2t1_aa1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,1),Vrt(1,1),1,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

         Call T2t1_aa1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,2),Vrt(1,2),2,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

      Endif 

      Call T2t1_ab1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +              Pop(1,1),Pop(1,2),Vrt(1,1),Vrt(1,2),1,Iuhf,
     +              T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +              T1ln_bb)

      If (Iuhf .Eq. 1) Then
          Call T2t1_ab1(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                  Pop(1,2),Pop(1,1),Vrt(1,2),Vrt(1,1),2,Iuhf,
     +                  T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +                  T1ln_bb)
      Endif 

      If (Iuhf .Eq. 1) Then 

         Call T2t1_aa2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,1),Vrt(1,1),1,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

         Call T2t1_aa2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                 Pop(1,2),Vrt(1,2),2,T2ln,T1ln,T2ln_aa,
     +                 T2ln_bb,T2ln_ab,T1ln_aa,T1ln_bb)

      Endif

      Call T2t1_ab2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +              Pop(1,1),Pop(1,2),Vrt(1,1),Vrt(1,2),1,Iuhf,
     +              T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +              T1ln_bb)

      If (Iuhf .Eq. 1) Then
          Call T2t1_ab2(T2amp(Ioff_t2),T2amp(Ioff_t1),Work,Maxcor,
     +                  Pop(1,2),Pop(1,1),Vrt(1,2),Vrt(1,1),2,Iuhf,
     +                  T2ln,T1ln,T2ln_aa,T2ln_bb,T2ln_ab,T1ln_aa,
     +                  T1ln_bb)
      Endif

#ifdef _DEBUG_LVL0
      Call checksum("T1AA :",T2amp(t2ln+Ione),t1ln_aa)
      If (Iuhf .Ne.0) Call checksum("T1BB :",T2amp(T2ln+T1ln_aa+Ione),
     +                               t1ln_bb)
#endif 

      Do Ispin = 1, Iuhf+1
         T1off = Ione + T2ln + (Ispin-1)*T1ln_aa 
         Call Form_t1(T2amp(T1off),Work,Maxcor,Pop(1,Ispin),
     +                 Vrt(1,Ispin),Ispin,T1ln,T1ln_aa)
      Enddo 

#ifdef _DEBUG_LVL0
      Call checksum("T1AA :",T2amp(t2ln+Ione),t1ln_aa)
      If (Iuhf .Ne.0) Call checksum("T1BB :",T2amp(T2ln+T1ln_aa+Ione),
     +                               t1ln_bb)
#endif 
   
      Return
      End 

      
      
