      Subroutine check_rhf_rcc(Work,Maxcor,Iuhf)

      Implicit Double precision (A-H,O-Z)
      Dimension Work(Maxcor)
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "machsp.com"
#include "syminf.com" /* for nirrep */
#include "sym.com"
#include "jodaflags.com"
#include "bwcc.com"
#include "sympop.com"

C Form spin-adapted integral list and write rhem to lists 117 and 118.

      Imode  = 0
      Irrepx = 1
      CALL INIPCK(1,19,21,197,IMODE,0,1)
      CALL INIPCK(1,19,21,198,IMODE,0,1)
C
      List_aa  = 61
      List_ab  = 63 
      Nsize_aa = Idsymsz(Irrepx,Isytyp(1,LIst_aa),Isytyp(2,List_aa))
      Nsize_ab = Idsymsz(Irrepx,Isytyp(1,LIst_ab),Isytyp(2,List_ab))
      If (Nsize_aa .NE. Nsize_ab) Then
         Write(6,"(a)") "Different AAAA and BBBB blocks"
         call Errex
      Endif 

      I000 = 1
      I010 = I000 + Nsize_aa
      I020 = I010 + Nsize_ab
      I030 = I020 + Nsize_ab

      Call Getall(Work(I000),Nsize_AA,Irrepx,List_AA)
      Call Dcopy(Nsize_aa,Work(I000),1,Work(I010),1)
      Call Getall(Work(I020),Nsize_AA,Irrepx,List_AB)

c Form T2(AAAA) + T2(ABAB)
     
       Call Daxpy(Nsize_aa,1.0D0,Work(I020),1,Work(I000),1)

#ifdef _DEBUG_LVL0
      call checksum("2J-K:",Work(I000),NSIZE_ab)
#endif 
c Form T2(AAAA) - T2(ABAB) 

       Call Daxpy(Nsize_aa,-1.0D0,Work(I020),1,Work(I010),1)
CSSS       Call Dcopy(Nsize_aa,Work(I020),1,Work(I010),1)
#ifdef _DEBUG_LVL0
      call checksum("-K  :",Work(I010),NSIZE_ab)
#endif 
       Nsize_ajbi = Idsymsz(Irrepx,Isytyp(1,21),Isytyp(2,21))
       Nsize_abij = Idsymsz(Irrepx,Isytyp(1,16),Isytyp(2,16))

      I030  = I020 + NSIZE_AJBI 
      I040  = I030 + NSIZE_ABIJ 
      I050  = I040 + NSIZE_ABIJ
      IEND  = I050

      If ((Nsize_aa .NE. Nsize_ab)     .or. 
     +    (NSIZE_AJBI .ne. NSIZE_ABIJ) .or.
     +    (Nsize_aa .NE. NSIZE_AJBI))   Then
         Write(6,"(a)") "Different AAAA and BBBB blocks"
         call Errex
      Endif 
      If (Iend .GE. MAXcor) Call Insmem("check_rhf_rcc",Iend,
     +                                   Maxcor)

C Read <Aj|bi> integrals

      CALL GETALL(Work(I020),NSIZE_AJBI,IRREPX,21)

C <Aj|bi> -> <Ab|ji>

      CALL SSTGEN(Work(I020),Work(I030),NSIZE_ABIJ,VRT(1,1),
     +            POP(1,2),VRT(1,2),POP(1,1),Work(I040),IRREPX,
     +           "1324")

C These are exchange integrals of <Ab|Ij> (RHF only)

C Form the triplet part of the spin-adapted <Ab|ji> by -<Ab|ji>

      CALL DSCAL(NSIZE_ABIJ,-1.0D0,Work(I030),1)
      CALL PUTALL(Work(I030),NSIZE_ABIJ,IRREPX,198)

#ifdef _DEBUG_LVL0
      call checksum("List-198:",Work(I030),NSIZE_ABIJ)
#endif

c Form the singlet part of the spin-adapted <Ab|ij> by 2<Ab|Ij>-<Ab|jI>

      CALL GETALL(Work(I020),NSIZE_ABIJ,IRREPX,16)
      CALL DSCAL(NSIZE_ABIJ,2.0D0,Work(I020),1)
      CALL DAXPY(NSIZE_ABIJ,1.0D0,Work(I030),1,Work(I020),1)
      CALL PUTALL(Work(I020),NSIZE_ABIJ,IRREPX,197)
#ifdef _DEBUG_LVL0
      call checksum("List-197:",Work(I020),NSIZE_ABIJ)
#endif 

c Compute the triplet and singlet spin-adapted energies.

      Call Getall(Work(I020),Nsize_abij,Irrepx,197)
      E_s = Ddot(Nsize_abij,Work(I020),1,Work(I000),1)

      Call Getall(Work(I020),Nsize_abij,Irrepx,198)
      E_t = Ddot(Nsize_abij,Work(I020),1,Work(I010),1)
     
      Write(6,"(a,F15.9)") "Singlet-spin adapted energy: ", E_s
      Write(6,"(a,F15.9)") "Triplet-spin adapted energy: ", E_t

      E_total = 0.75D0*E_t + 0.25D0*E_s
      
      Write(6,"(a,F15.9)") "Total energy: ", E_total

      Return
      End
