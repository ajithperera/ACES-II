#include "flags.h"
      Subroutine Tdcc_form_mutilde_dot_d(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)

      Double Precision Mone, Mut_dot 
      Dimension Ioff_vo(8,2)
 
      Data Zero, One, Mone, Two /0.0D0, 1.0D0, -1.0D0, 2.0D0/

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "info.com"
#include "symloc.com"

      Call Tdcc_load_mutilde_vo(Work,Memleft,Memleft_modf,Iuhf,
     +                          Ioff_vo,Irrepx)

C NDep paper Eqn. 29 first term(^t indicate transpose)
C
C UHF: Z(AB,IJ) = P(AB) {Mu_t(B,M)^t * Hbar(AM,IJ)^t}
C    : Z(ab,ij) = P(ab) {Mu_t(b,m)^t * Hbar(am,ij)^t}
C    : Z(aB,iJ) = P(aB) {Mu_t(B,M)^t * Hbar(aM,iJ)^t}

      Call Tdcc_form_mutilde_dot_da(Work,Memleft_modf,Irrepx,
     +                              Iuhf,Ioff_vo)
#ifdef _DEBUG_LVL0
      call tdcc_mutilde_dot_da_debug(Work,Memleft,Iuhf,Irrepx)
#endif 

C NDep paper Eqn. 29 second term
C
C UHF: Z(AB,IJ) = P(IJ) {Mu_t(E,J)^t * Hbar(AB,IE)^t}
C    : Z(ab,ij) = P(ab) {Mu_t(e,j)^t * Hbar(AB,ie)^t}
C    : Z(aB,iJ) = P(aB) {Mu_t(E,J)^t * Hbar(aB,iJ)^t}
   
      Call Tdcc_form_mutilde_dot_db(Work,Memleft_modf,Irrepx,
     +                              Iuhf,Ioff_vo)
#ifdef _DEBUG_LVL0
      call tdcc_mutilde_dot_db_debug(Work,Memleft,Iuhf,Irrepx)
#endif 

C NDep paper Eqn. 29 third term.
C Z(AB,IJ) = P(AB)P(IJ) Mu_t(I,A) F(B,J)

      If (Iuhf .Ne. 0) Then
         Call Tdcc_form_mutilde_dot_dc_uhf(Work,Memleft,Irrepx)
     +                                     
      Else
         Call Tdcc_form_mutilde_dot_dc_rhf(Work,Memleft,Irrepx)
      Endif 

C NDep paper Eqn. 29 Forth-term.
C Z(AB,IJ) = P(AB)P(IJ) Mu_t(I,A) F(B,J)

      If (Irrepx .EQ. 1) Then
         Call Tdcc_form_mutilde_dot_dd(Work,Memleft,Irrepx,
     +                                 Iuhf,Mut_dot)
      Else
         Mut_dot = 0.0D0
      Endif 

      Return
      End
     
           
        
