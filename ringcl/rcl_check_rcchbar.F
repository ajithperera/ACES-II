      Subroutine Rcl_check_rcchbar(Work,Maxcor,Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Logical MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD,RLE
   
      Dimension Work(Maxcor)

#include "machsp.com"
#include "syminf.com" 
#include "sym.com"
#include "flags.com"
#include "flags2.com"
#include "bwcc.com"
#include "sympop.com"

      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC,
     &            CC2,RCCD,DRCCD
   

C This built drCCD(RPA) Hbar matrix and diagonalize in singles space.
C If we restrict one particle contributions only to Fock diagonals this
C should give the RPA eigenvalues and vectors back. 

      Nmult = 1
      If (RCCD .AND. IUHF .EQ. 0) Nmult = 2

      Do Imult = 1, Nmult 
         Call Rcl_heffxc(Work,Maxcor,Iuhf,Imult)

         Call Rcl_drvrpa(Work,Maxcor,Iuhf)
      Enddo 

      Return 
      End 
