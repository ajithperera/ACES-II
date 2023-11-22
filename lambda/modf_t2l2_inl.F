      Subroutine Modf_T2L2_INL(Work,Length,Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Length)

#include "machsp.com"
#include "sympop.com"
#include "syminf.com"
#include "sym.com"
#include "info.com"
#include "parcc.com"

      Data Onem,One /-1.0D0,1.0D0/

C First T2 lists. COnstruct the T2(AAAA) and T2(BBBB) from 
C T2(ABAB) list. This is valid for RHF (or UHF=RHF) situations.

C Create lists to store the T2 (AAAA,BBBB) and L2(AAAA,BBBB).
C We can not simply overwrite the existing 34,35,134 and 135 
C list.

      Imode  = 0
      Irrepx = 1
      Call Inipck(Irrepx,9,9,234,Imode,0,1)
      Call Inipck(Irrepx,10,10,235,Imode,0,1)
      Call Inipck(Irrepx,9,9,236,Imode,0,1)
      Call Inipck(Irrepx,10,10,237,Imode,0,1)
 
      Do Irrepr = 1, Nirrep
         Irrepl = DIRPRD(Irrepr,Irrepx)

C Read T2(Ij,Ab) list 37 stored as (AI,bj)
C Read T2(Ij,Ab) list 38 stores as (bI,Aj)

         Len_ph_aa_l = Irpdpd(Irrepl,9)
         Len_ph_aa_r = Irpdpd(Irrepr,9)
         Len_ph_bb_r = Irpdpd(Irrepr,10)
         Len_ph_ab_l = Irpdpd(Irrepl,11)
         Len_ph_ba_r = Irpdpd(Irrepr,12)
         Len_pp_ab_l = Irpdpd(Irrepr,15)
         Len_hh_ab_r = Irpdpd(Irrepr,14)

         Ndim_aabb = Len_ph_aa_l * Len_ph_bb_r
         Ndim_abba = Len_ph_ab_l * Len_ph_ba_r
         Ndim_aaaa = Len_ph_aa_l * Len_ph_aa_r
         Ndim_abab = Len_pp_ab_l * Len_hh_ab_r

         If (Ndim_aabb .NE. Ndim_abba) Then 
            Write(6,"(a)") " The two lists have different dimensions"
            Call Errex
         Endif 

         I000 = 1
         I010 = I000 + Ndim_aabb
         I020 = I010 + Ndim_abba
         I030 = I020 + Ndim_aaaa
         Iend = I030
 
         If (Iend .GE. Length) Call Insmem("@-Mody_TL2",Iend,Length)
 
         Call Getlst(Work(I000),1,Len_ph_bb_r,1,Irrepr,37)
         Call Getlst(Work(I010),1,Len_ph_ba_r,1,Irrepr,39)

#ifdef _DEBUG_LVL0
      call checksum("List-37:",Work(I000),Ndim_aabb)
      call checksum("List-39:",Work(I010),Ndim_abba)
#endif
  
C Form the antisymmetrized T2.

         Call Daxpy(Ndim_aabb,Onem,Work(I010),1,Work(I000),1)

C Store these T2 vectors as AAAA and BBBB in to special lists

         Call Putlst(Work(I000),1,Len_ph_bb_r,1,Irrepr,234)
         Call Putlst(Work(I000),1,Len_ph_bb_r,1,Irrepr,235)

#ifdef _DEBUG_LVL0
      call checksum("Check1 :",Work(I000),Ndim_aabb)
C Read the T2(IJ,AB) list 34 stored as (AJ,BI)
CSSS      call output(Work(I000),1,Len_ph_aa_l,1,Len_ph_aa_r,
CSSS     +            Len_ph_aa_l,Len_ph_aa_r,1)
      Call Getlst(Work(I020),1,Len_ph_aa_r,1,Irrepr,34)
      Call Daxpy(Ndim_aaaa,One,Work(I000),1,Work(I020),1)
      call checksum("Check2 :",Work(I020),Ndim_aaaa)
#endif  

C Do exactly the same things for the Lambda ring lists 

         Call Getlst(Work(I000),1,Len_ph_bb_r,1,Irrepr,137)
         Call Getlst(Work(I010),1,Len_ph_ba_r,1,Irrepr,139)
C
         Call Daxpy(Ndim_aabb,Onem,Work(I010),1,Work(I000),1)

C Store these L2 vector as L2 AAAA and BBBB in to special lists

         Call Putlst(Work(I000),1,Len_ph_bb_r,1,Irrepr,236)
         Call Putlst(Work(I000),1,Len_ph_bb_r,1,Irrepr,237)

#ifdef _DEBUG_LVL0
C Another way to built antisymmetric (AJ,BI). Read (Ab,Ij) list
      I040 = I030 + Ndim_abab
      I050 = I040 + Ndim_aaaa
      I060 = I050 + Len_pp_ab_l
      Iend = I060 
      If (Iend .GE. Length) Call Insmem("@-Mody_TL2",Iend,Length)

      Call Getlst(Work(I030),1,Len_hh_ab_r,1,Irrepr,46)
      Call Assym2(Irrepr,Pop(1,1),Len_pp_ab_l,Work(I030))
C This leads to  (Ab,I<j),expand it to (Ab,ij)
      Call Symexp(Irrepr,Pop(1,1),Len_pp_ab_l,Work(I030))
C Reorder from (Ab,Ij) to Aj,bI)
      Call Sstgen(Work(I030),Work(I040),Ndim_abab,Vrt(1,1),
     +            Vrt(1,2),Pop(1,1),Pop(1,2),Work(I050),
     +            Irrepx,"1324")
      Call Getlst(Work(I000),1,Len_ph_aa_r,1,Irrepr,34)
      Call Daxpy(Ndim_aaaa,One,Work(I000),1,Work(I040),1)
      Call Checksum("Check3 :",Work(I040),Ndim_aaaa)
#endif 
      Enddo 

      Return 
      End 
