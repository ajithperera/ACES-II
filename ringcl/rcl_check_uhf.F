      Subroutine Rcl_check_uhf(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)

      Dimension Work(Maxcor)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"

      Data One, Onem, Half /1.0D0, -1.0D0, 0.50D0/

      Do Irrep = 1, Nirrep

         If (Iuhf .eq. 0) Then

         Nrow  = Irpdpd(Irrep,13)
         Ncol  = Irpdpd(Irrep,14)
         Length = Nrow * Ncol

         I000 = 1
         I010 = I000 + Length 
         I020 = I010 + Length 
         I030 = I020 + Length
         I040 = I030 + Length
         I040 = Iend 

         Call Getlst(Work(I000),1,Ncol,1,Irrep,144)
         Call Getlst(Work(I010),1,Ncol,1,Irrep,146)
#ifdef _DEBUG_LVL0
         Call checksum("144:",Work(I000),Length)
         Call checksum("146:",Work(I010),Length)
#endif 

C Form 3^T(AB,AB) + 1^T(ABAB) = 2 * T(AAAA) and 1^T(ABAB) - 3^T(AB,AB) = 
C 2 * T*(AB,AB)

         Call Vadd(Work(I020),Work(I000),Work(I010),Length,One)
         Call Dscal(Length,Half,Work(I020),1)

         Call Vadd(Work(I030),Work(I010),Work(I000),Length,Onem)
         Call Dscal(Length,Half,Work(I030),1)

#ifdef _DEBUG_LVL0
         Write(6,"(a)") "@-Rcl_check_uhf; RHF T(AAAA) and T(ABAB)"
         Call checksum("T(AAAA):",Work(I020),Length)
CSSS         call output(Work(i020),1,Nrow,1,Ncol,Nrow,Ncol_aa)
         Call checksum("T(ABAB):",Work(I030),Length)
CSSS         call output(Work(i030),1,Nrow,1,Ncol,Nrow,Ncol_aa)
#endif 
         Else
          
         Nrow_aa = Irpdpd(Irrep,19)
         Nrow_bb = Irpdpd(Irrep,20)
         Nrow_ab = Irpdpd(Irrep,13)

         Ncol_aa = Irpdpd(Irrep,21)
         Ncol_bb = Irpdpd(Irrep,22)
         Ncol_ab = Irpdpd(Irrep,14)

         Length_aa = Nrow_aa * Ncol_aa
         Length_bb = Nrow_bb * Ncol_bb
         Length_ab = Nrow_ab * Ncol_ab

         I000 = 1
         I010 = I000 + Length_aa
         I020 = I010 + Length_bb
         I030 = I020 + Length_ab
         Iend = I030

         Call Getlst(Work(I000),1,Ncol_aa,1,Irrep,144)
         Call Getlst(Work(I010),1,Ncol_bb,1,Irrep,145)
         Call Getlst(Work(I020),1,Ncol_ab,1,Irrep,146)

#ifdef _DEBUG_LVL0
        Write(6,"(a)") "@-Rcl_check_uhf; UHF T(AAAA),T(BBBB),TABAB"
        Call checksum("T(AAAA):",Work(I000),Length_aa)
CSSS        call output(Work(i000),1,Nrow_aa,1,Ncol_aa,Nrow_aa,
CSSS     +              Ncol_aa,1)
        Call checksum("T(BBBB):",Work(I010),Length_bb)
        Call checksum("T(ABAB):",Work(I020),Length_ab)
CSSS        call output(Work(i000),1,Nrow_ab,1,Ncol_ab,Nrow_ab,
CSSS     +              Ncol_ab,1)
#endif 
        Endif 

       Enddo
   
       Return
       End 
