#include "flags.h"
      Subroutine Tdcc_form_mu_dot_d(Work,Memleft,Irrepx,Iuhf)

      Implicit Double Precision (A-H,O-Z)

      Dimension Work(Memleft)

      Double Precision Mone
      Dimension Ioff_vo(8,2)
 
      Data Zero, One, Mone, Two /0.0D0, 1.0D0, -1.0D0, 2.0D0/

      Call Tdcc_load_mu_vo(Work,Memleft,Memleft_modf,Iuhf,Ioff_vo,
     +                      Irrepx) 

C NDep paper Eqn. 26 first term
C
C UHF: Z(AB,IJ) = P(AB) {Mu(B,M) * Hbar(AM,IJ)}
C    : Z(ab,ij) = P(ab) {Mu(b,m) * Hbar(am,ij)}
C    : Z(aB,iJ) = P(aB) {Mu(B,M) * Hbar(aM,iJ)}

      Call Tdcc_form_mu_dot_da(Work,Memleft_modf,Irrepx,Iuhf,
     +                         Ioff_vo)
#ifdef _DEBUG_LVL0
      call tdcc_mu_dot_da_debug(Work,Memleft,Iuhf,Irrepx)
#endif
C NDep paper Eqn. 26 second term
C
C UHF: Z(AB,IJ) = P(IJ) {Mu(E,J) * Hbar(AB,IE)}
C    : Z(ab,ij) = P(ab) {Mu(e,j) * Hbar(AB,ie)}
C    : Z(aB,iJ) = P(aB) {Mu(E,J) * Hbar(aB,iJ)}
   
      Call Tdcc_form_mu_dot_db(Work,Memleft_modf,Irrepx,Iuhf,
     +                         Ioff_vo)
#ifdef _DEBUG_LVL0
      call tdcc_mu_dot_db_debug(Work,Memleft,Iuhf,Irrepx)
#endif

      Return
      End
     
           
        
