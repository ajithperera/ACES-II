      Subroutine form_modf_hbar_mbej_4ee_corr(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H, O-Z)
      Dimension Work(Maxcor)

      LOGICAL MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      LOGICAL  CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES

      COMMON /REFTYPE/ MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD,CC2
      COMMON /METH/ CIS,RPA,EOMCC,CISD,FULDIAG,INCORE,READGUES

#include "flags.com"
#include "flags2.com"
#include "machsp.com"
#include "sym.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"

      If (Drccd) Then
         CALL Drcl_dwmbej_4ee_corr(Work,Maxcor,Iuhf)
      Elseif (Rccd) Then
         If (Iuhf .Eq.0) Then
            CALL Rcl_dwmbej_r_4ee_corr(Work,Maxcor,Iuhf)
         Else
            CALL Rcl_dwmbej_u_4ee_corr(Work,Maxcor,Iuhf)
         Endif
      Endif

      Return
      End

