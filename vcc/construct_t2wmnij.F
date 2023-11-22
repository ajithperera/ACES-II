#include "flags.h"
      Subroutine construct_t2wmnij(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

      Imode = 0
#ifdef _DEBUG_LVL0
      Write(6,"(a,a)") " Construction of T2 piece of the", 
     &                 " W(mn,ij) = 1/4t2(ij,ef)*<mn||ef>."
      Write(6,*)
#endif 

      Call Inipck(1,14,14,253,Imode,0,1)

      If (Iuhf .ne. 0) Then
         Call Inipck(1,3,3,251,Imode,0,1)
         Call Inipck(1,4,4,252,Imode,0,1)
      Endif
C
C No need to do anything here for methods that do not have
C singles, but needs the 251-253 lists.
C
      If (Iflags(h_IFLAGS_calc) .lt. 9) Return
C
C Construct the W(mn,ij) intermediate as it would be for
C standard Cc calc.
C
      Intpck = 5
      Fact   = -(Pargamma - 1.0D0)

      Call Pdcc_quad1(Work,Length,Intpck,Iuhf,Fact)
      Call pdcc_drlad(Work,Length,Iuhf,1)

      Return 
      End 

